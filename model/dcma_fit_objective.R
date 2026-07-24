dcma_objective <- function(
  data,
  objective = c("wgr"),
  epsm_dim = 4, epsy_dim = 4,
  hidden_dim = 64, num_layer = 5,
  hidden_dim_M = NULL, num_layer_M = NULL,
  hidden_dim_Y = NULL, num_layer_Y = NULL,
  m_generator = c("plain", "shared"),
  m_shared_noise_dim = NULL,
  y_generator = c("plain", "locscale", "xgated"),
  epochs_M = 2000, epochs_Y = 2000,
  lr = 5e-4,
  silent = FALSE, standardize = TRUE,
  val_p = 0.2, val_freq = 1L,
  patience = 100L, min_delta = 0.0,
  c_vars = "C",
  critic_steps = 5L,
  gp_lambda = 10,
  critic_hidden_dim = NULL,
  critic_num_layer = 3L,
  lambda_w = 0.9,
  lambda_l = 0.1,
  wgr_J_size = NULL,
  device = "mps",
  y_family = c("continuous", "binary")
) {
  objective <- match.arg(objective)
  y_family <- match.arg(y_family)
  m_generator <- match.arg(m_generator)
  y_generator <- match.arg(y_generator)
  critic_steps <- as.integer(critic_steps)
  if (critic_steps < 1L) stop("critic_steps must be >= 1.")
  stopifnot(is.numeric(lambda_w), length(lambda_w) == 1L, is.finite(lambda_w), lambda_w >= 0)
  stopifnot(is.numeric(lambda_l), length(lambda_l) == 1L, is.finite(lambda_l), lambda_l >= 0)
  if (is.null(wgr_J_size)) wgr_J_size <- 4L
  wgr_J_size <- as.integer(wgr_J_size)
  if (wgr_J_size < 1L) stop("wgr_J_size must be >= 1.")

  if (!identical(y_family, "continuous")) {
    stop("dcma_objective() currently supports only continuous outcomes.")
  }
  if (identical(m_generator, "shared")) {
    if (epsm_dim < 2L) stop("m_generator='shared' requires epsm_dim >= 2.")
    if (!is.null(m_shared_noise_dim)) {
      stopifnot(is.numeric(m_shared_noise_dim), length(m_shared_noise_dim) == 1L, is.finite(m_shared_noise_dim))
      m_shared_noise_dim <- as.integer(m_shared_noise_dim)
      if (m_shared_noise_dim < 1L || m_shared_noise_dim >= epsm_dim) {
        stop("m_shared_noise_dim must be between 1 and epsm_dim - 1.")
      }
    }
  }
  if (identical(y_generator, "xgated")) {
    if (!identical(y_family, "continuous")) {
      stop("dcma_objective() currently supports only continuous outcomes.")
    }
  }

  stopifnot(is.data.frame(data))
  stopifnot(all(c("X", "Y") %in% names(data)))
  stopifnot(all(c_vars %in% names(data)))

  n_all <- nrow(data)
  idx <- sample(seq_len(n_all), size = floor((1 - val_p) * n_all))
  train_raw <- data[idx, , drop = FALSE]
  val_raw <- data[-idx, , drop = FALSE]

  m_vars <- grep("^M", names(train_raw), value = TRUE)
  S <- length(m_vars)
  if (S < 1L) stop("No mediator columns found.")

  vars_core <- c(m_vars, "Y", c_vars)
  stats <- lapply(train_raw[vars_core], function(v) c(mean = mean(v), sd = sd(v)))
  mu <- sapply(stats, `[[`, "mean")
  sdd <- sapply(stats, `[[`, "sd")
  bad_sd <- !is.finite(sdd) | sdd == 0
  sdd[bad_sd] <- 1

  std_df <- function(df) {
    core <- df[vars_core]
    if (standardize) {
      core <- as.data.frame(scale(core, center = mu[vars_core], scale = sdd[vars_core]))
    }
    if (is.factor(df$X)) {
      core$X <- as.numeric(df$X) - 1L
    } else {
      core$X <- as.numeric(df$X)
    }
    core
  }

  train_std <- std_df(train_raw)
  val_std <- std_df(val_raw)

  X <- as.matrix(train_std[, "X", drop = FALSE])
  M_all <- as.matrix(train_std[, m_vars, drop = FALSE])
  Y <- as.matrix(train_std[, "Y", drop = FALSE])
  C <- as.matrix(train_std[, c_vars, drop = FALSE])
  colnames(C) <- c_vars
  C_raw <- as.matrix(train_raw[, c_vars, drop = FALSE])
  colnames(C_raw) <- c_vars

  Xv <- as.matrix(val_std[, "X", drop = FALSE])
  Mv_all <- as.matrix(val_std[, m_vars, drop = FALSE])
  Yv <- as.matrix(val_std[, "Y", drop = FALSE])
  Cv <- as.matrix(val_std[, c_vars, drop = FALSE])
  colnames(Cv) <- c_vars
  Cv_raw <- as.matrix(val_raw[, c_vars, drop = FALSE])
  colnames(Cv_raw) <- c_vars

  p <- ncol(C)
  n <- nrow(X)
  nv <- nrow(Xv)

  X_t <- torch_tensor(X, dtype = torch_float(), device = device)
  M_all_t <- torch_tensor(M_all, dtype = torch_float(), device = device)
  Y_t <- torch_tensor(Y, dtype = torch_float(), device = device)
  C_t <- torch_tensor(C, dtype = torch_float(), device = device)
  XC_t <- torch_cat(list(X_t, C_t), dim = 2)

  if (nv > 0L) {
    Xv_t <- torch_tensor(Xv, dtype = torch_float(), device = device)
    Mv_all_t <- torch_tensor(Mv_all, dtype = torch_float(), device = device)
    Yv_t <- torch_tensor(Yv, dtype = torch_float(), device = device)
    Cv_t <- torch_tensor(Cv, dtype = torch_float(), device = device)
    XCv_t <- torch_cat(list(Xv_t, Cv_t), dim = 2)
  } else {
    Xv_t <- Mv_all_t <- Yv_t <- Cv_t <- XCv_t <- NULL
  }

  in_f <- as.integer(1L + p)
  in_g <- as.integer(1L + p + S)
  hd_M <- if (is.null(hidden_dim_M)) hidden_dim else hidden_dim_M
  nl_M <- if (is.null(num_layer_M)) num_layer else num_layer_M
  hd_Y <- if (is.null(hidden_dim_Y)) hidden_dim else hidden_dim_Y
  nl_Y <- if (is.null(num_layer_Y)) num_layer else num_layer_Y
  critic_hd <- if (is.null(critic_hidden_dim)) hidden_dim else critic_hidden_dim

  if (identical(m_generator, "shared")) {
    gen_f <- nn_model_shared_resid_mediator(
      in_dim = in_f,
      noise_dim = epsm_dim,
      hidden_dim = hd_M,
      out_dim = S,
      shared_noise_dim = m_shared_noise_dim
    )$to(device = device)
  } else {
    gen_f <- nn_model(in_dim = in_f, noise_dim = epsm_dim, hidden_dim = hd_M, out_dim = S, num_layer = nl_M)$to(device = device)
  }
  if (identical(y_generator, "locscale")) {
    gen_g <- nn_model_locscale_outcome(
      in_dim = in_g,
      noise_dim = epsy_dim,
      hidden_dim = hd_Y
    )$to(device = device)
  } else if (identical(y_generator, "xgated")) {
    gen_g <- nn_model_xgated_heads(
      in_dim = in_g,
      noise_dim = epsy_dim,
      hidden_dim = hd_Y,
      out_dim = 1L,
      num_layer = nl_Y
    )$to(device = device)
  } else {
    gen_g <- nn_model(in_dim = in_g, noise_dim = epsy_dim, hidden_dim = hd_Y, out_dim = 1L, num_layer = nl_Y)$to(device = device)
  }

  opt_f <- optim_adam(gen_f$parameters, lr = lr)
  opt_g <- optim_adam(gen_g$parameters, lr = lr)

  sample_eps <- function(nrow, dim) {
    torch_randn(c(nrow, dim), device = device)
  }

  loss_M_hist <- numeric(epochs_M)
  val_loss_M_hist <- rep(NA_real_, epochs_M)
  loss_Y_hist <- numeric(epochs_Y)
  val_loss_Y_hist <- rep(NA_real_, epochs_Y)

  state_list <- list()

  grad_fn <- function() {
    if (exists("autograd_grad", mode = "function", inherits = TRUE)) {
      return(get("autograd_grad", mode = "function", inherits = TRUE))
    }
    if (exists("torch_autograd_grad", mode = "function", inherits = TRUE)) {
      return(get("torch_autograd_grad", mode = "function", inherits = TRUE))
    }
    NULL
  }

  gradient_penalty <- function(critic, real_xy, fake_xy) {
    if (gp_lambda <= 0) return(torch_tensor(0, device = device, dtype = real_xy$dtype))
    grad_fun <- grad_fn()
    if (is.null(grad_fun)) {
      stop("Gradient penalty requested but autograd_grad is unavailable in this torch build.")
    }
    bsz <- real_xy$size(1)
    eps <- torch_rand(bsz, 1, device = device, dtype = real_xy$dtype)
    eps <- eps$expand_as(real_xy)
    xy_hat <- eps * real_xy + (1 - eps) * fake_xy
    xy_hat$requires_grad_(TRUE)
    d_hat <- critic(xy_hat)
    grads <- grad_fun(
      outputs = d_hat,
      inputs = xy_hat,
      grad_outputs = torch_ones_like(d_hat),
      create_graph = TRUE,
      retain_graph = TRUE
    )[[1]]
    grad_norm <- torch_sqrt(torch_sum(grads$pow(2), dim = 2) + 1e-12)
    torch_mean((grad_norm - 1)$pow(2)) * gp_lambda
  }

  mean_l2_regularizer <- function(real_t, cond_t, noise_dim, generator, J_size) {
    fake_list <- lapply(seq_len(J_size), function(j) {
      generator(torch_cat(list(cond_t, sample_eps(real_t$size(1), noise_dim)), dim = 2))
    })
    fake_mean <- torch_stack(fake_list, dim = 3L)$mean(dim = 3L)
    torch_mean((fake_mean - real_t)$pow(2))
  }

  train_stage_wgr <- function(real_t, cond_t, noise_dim, out_dim, generator, optimizer, epochs, stage_name, val_real_t = NULL, val_cond_t = NULL) {
    critic <- nn_critic(
      input_dim = cond_t$size(2) + out_dim,
      hidden_dim = critic_hd,
      num_layer = critic_num_layer
    )$to(device = device)
    opt_c <- optim_adam(critic$parameters, lr = lr, betas = c(0.5, 0.9))

    loss_hist <- numeric(epochs)
    val_hist <- rep(NA_real_, epochs)
    best_val <- Inf
    best_epoch <- 0L
    best_gen_state <- NULL
    best_critic_state <- NULL

    for (epoch in seq_len(epochs)) {
      d_last <- NULL
      for (k in seq_len(critic_steps)) {
        opt_c$zero_grad()
        eps <- sample_eps(n, noise_dim)
        fake_t <- generator(torch_cat(list(cond_t, eps), dim = 2))$detach()
        real_xy <- torch_cat(list(cond_t, real_t), dim = 2)
        fake_xy <- torch_cat(list(cond_t, fake_t), dim = 2)
        d_real <- critic(real_xy)
        d_fake <- critic(fake_xy)
        wasserstein_part <- torch_mean(d_fake) - torch_mean(d_real)
        gp_term <- gradient_penalty(critic, real_xy, fake_xy)
        d_loss <- wasserstein_part + gp_term
        d_loss$backward()
        opt_c$step()
        d_last <- d_loss
      }

      optimizer$zero_grad()
      eps <- sample_eps(n, noise_dim)
      fake_t <- generator(torch_cat(list(cond_t, eps), dim = 2))
      fake_xy <- torch_cat(list(cond_t, fake_t), dim = 2)
      wasserstein_gen <- -torch_mean(critic(fake_xy))
      if (identical(objective, "wgr") && lambda_l > 0) {
        reg_loss <- mean_l2_regularizer(real_t, cond_t, noise_dim, generator, wgr_J_size)
        g_loss <- lambda_w * wasserstein_gen + lambda_l * reg_loss
      } else {
        g_loss <- wasserstein_gen
      }
      g_loss$backward()
      optimizer$step()
      loss_hist[epoch] <- g_loss$item()

      if (!is.null(val_real_t) && epoch %% val_freq == 0L) {
        with_no_grad({
          epsv <- sample_eps(nv, noise_dim)
          fake_v <- generator(torch_cat(list(val_cond_t, epsv), dim = 2))
          fake_xy_v <- torch_cat(list(val_cond_t, fake_v), dim = 2)
          val_w <- -torch_mean(critic(fake_xy_v))
          if (identical(objective, "wgr") && lambda_l > 0) {
            val_reg <- mean_l2_regularizer(val_real_t, val_cond_t, noise_dim, generator, wgr_J_size)
            cur_val <- val_reg$item()
          } else {
            cur_val <- val_w$item()
          }
          val_hist[epoch] <- cur_val
          if (is.infinite(best_val) || cur_val + min_delta < best_val) {
            best_val <- cur_val
            best_epoch <- epoch
            best_gen_state <- generator$state_dict()
            best_critic_state <- critic$state_dict()
          } else if (patience > 0L && (epoch - best_epoch) >= patience) {
            if (!silent) {
              cat(sprintf("[%s] early stopping at epoch %d best epoch %d val obj = %.4f\n", stage_name, epoch, best_epoch, best_val))
            }
            break
          }
        })
      }
    }

    if (!is.null(best_gen_state)) {
      generator$load_state_dict(best_gen_state)
    }
    if (!is.null(best_critic_state)) {
      critic$load_state_dict(best_critic_state)
    }

    list(
      generator = generator,
      critic = critic,
      loss_hist = loss_hist,
      val_hist = val_hist
    )
  }

  res_f <- train_stage_wgr(
    real_t = M_all_t,
    cond_t = XC_t,
    noise_dim = epsm_dim,
    out_dim = S,
    generator = gen_f,
    optimizer = opt_f,
    epochs = epochs_M,
    stage_name = "M",
    val_real_t = Mv_all_t,
    val_cond_t = XCv_t
  )
  gen_f <- res_f$generator
  critic_f <- res_f$critic
  loss_M_hist <- res_f$loss_hist
  val_loss_M_hist <- res_f$val_hist

  res_g <- train_stage_wgr(
    real_t = Y_t,
    cond_t = torch_cat(list(XC_t, M_all_t), dim = 2),
    noise_dim = epsy_dim,
    out_dim = 1L,
    generator = gen_g,
    optimizer = opt_g,
    epochs = epochs_Y,
    stage_name = "Y",
    val_real_t = Yv_t,
    val_cond_t = if (nv > 0L) torch_cat(list(XCv_t, Mv_all_t), dim = 2) else NULL
  )
  gen_g <- res_g$generator
  critic_g <- res_g$critic
  loss_Y_hist <- res_g$loss_hist
  val_loss_Y_hist <- res_g$val_hist

  DCMA <- list(
    gen_f = gen_f,
    gen_g = gen_g,
    critic_f = critic_f,
    critic_g = critic_g,
    loss_M = loss_M_hist,
    loss_Y = loss_Y_hist,
    val_loss_M = val_loss_M_hist,
    val_loss_Y = val_loss_Y_hist,
    state_list = state_list,
    X = X,
    C = C,
    C_raw = C_raw,
    Cv = Cv,
    Cv_raw = Cv_raw,
    M_obs = M_all,
    Y = Y,
    m_vars = m_vars,
    epsm_dim = epsm_dim,
    epsy_dim = epsy_dim,
    mu = mu,
    sdd = sdd,
    vars_core = vars_core,
    standardize = standardize,
    K_es = wgr_J_size,
    y_family = y_family,
    config = list(
      c_vars = c_vars,
      hidden_dim = hidden_dim,
      num_layer = num_layer,
      hidden_dim_M = hd_M,
      num_layer_M = nl_M,
      hidden_dim_Y = hd_Y,
      num_layer_Y = nl_Y,
      m_generator = m_generator,
      m_shared_noise_dim = m_shared_noise_dim,
      y_generator = y_generator,
      epochs_M = epochs_M,
      epochs_Y = epochs_Y,
      lr = lr,
      val_p = val_p,
      val_freq = val_freq,
      patience = patience,
      min_delta = min_delta,
      device = device,
      y_family = y_family,
      objective = objective,
      critic_steps = critic_steps,
      gp_lambda = gp_lambda,
      lambda_w = lambda_w,
      lambda_l = lambda_l,
      wgr_J_size = wgr_J_size
    )
  )
  class(DCMA) <- "DCMA"
  DCMA
}
