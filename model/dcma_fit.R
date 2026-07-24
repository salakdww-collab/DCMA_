dcma <- function(
  data,
  epsm_dim = 4, epsy_dim = 4,
  hidden_dim = 64, num_layer = 5,
  hidden_dim_M = NULL, num_layer_M = NULL,
  hidden_dim_Y = NULL, num_layer_Y = NULL,
  epochs_M = 2000, epochs_Y = 2000,
  lr = 5e-4, beta = 1, beta_M = NULL, beta_Y = NULL,
  silent = FALSE, standardize = TRUE,
  val_p = 0.2, val_freq = 1L,
  patience = 100L, min_delta = 0.0,
  patience_M = NULL, patience_Y = NULL,
  min_delta_M = NULL, min_delta_Y = NULL,
  c_vars = "C",
  save_epochs = c(1, 10, 100, 1000),
  K_es = 2L,
  m_generator = c("plain", "shared", "separate"),
  m_shared_noise_dim = NULL,
  m_reinject_mode = c("none", "obs"),
  y_balance_groups = FALSE,
  grad_clip_norm = NA_real_,
  y_reinject_x = FALSE,
  y_generator = c("plain", "locscale", "xgated", "additive_homosk", "additive_heterosk"),
  y_reinject_mode = c("none", "x", "obs", "obs2", "film2", "ln2", "adaln2"),
  split_seed = NULL,
  init_seed = NULL,
  device = "mps",
  y_family = c("continuous", "binary", "categorical"),
  y_levels = NULL
) {
  ## -----------------------------
  ## 0. Basic checks and setup
  ## -----------------------------
  stopifnot(is.data.frame(data))
  stopifnot(all(c("X", "Y") %in% names(data)))
  stopifnot(all(c_vars %in% names(data)))
  stopifnot(is.numeric(val_p), val_p > 0, val_p < 1)
  stopifnot(is.integer(K_es) || is.numeric(K_es))
  K_es <- as.integer(K_es)
  stopifnot(K_es >= 2L)
  if (is.null(beta_M)) beta_M <- beta
  if (is.null(beta_Y)) beta_Y <- beta
  stopifnot(is.numeric(beta_M), length(beta_M) == 1L, is.finite(beta_M), beta_M > 0)
  stopifnot(is.numeric(beta_Y), length(beta_Y) == 1L, is.finite(beta_Y), beta_Y > 0)
  m_generator <- match.arg(m_generator)
  m_reinject_mode <- match.arg(m_reinject_mode)
  stopifnot(is.logical(y_balance_groups), length(y_balance_groups) == 1L)
  stopifnot(is.logical(y_reinject_x), length(y_reinject_x) == 1L)
  y_generator <- match.arg(y_generator)
  y_reinject_mode <- match.arg(y_reinject_mode)
  if (isTRUE(y_reinject_x) && identical(y_reinject_mode, "none")) y_reinject_mode <- "x"
  if (!is.null(split_seed)) stopifnot(is.numeric(split_seed), length(split_seed) == 1L, is.finite(split_seed))
  if (!is.null(init_seed)) stopifnot(is.numeric(init_seed), length(init_seed) == 1L, is.finite(init_seed))
  if (!is.na(grad_clip_norm)) {
    stopifnot(is.numeric(grad_clip_norm), length(grad_clip_norm) == 1L, is.finite(grad_clip_norm), grad_clip_norm > 0)
  }
  if (is.null(patience_M)) patience_M <- patience
  if (is.null(patience_Y)) patience_Y <- patience
  if (is.null(min_delta_M)) min_delta_M <- min_delta
  if (is.null(min_delta_Y)) min_delta_Y <- min_delta
  patience_M <- as.integer(patience_M)
  patience_Y <- as.integer(patience_Y)
  stopifnot(patience_M >= 0L, patience_Y >= 0L)
  stopifnot(is.numeric(min_delta_M), length(min_delta_M) == 1L, is.finite(min_delta_M))
  stopifnot(is.numeric(min_delta_Y), length(min_delta_Y) == 1L, is.finite(min_delta_Y))
  y_family <- match.arg(y_family)
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
  if (identical(y_generator, "locscale") || grepl("^additive_", y_generator)) {
    if (!identical(y_family, "continuous")) stop("Structured/additive y_generator variants currently support only continuous outcomes.")
    if (!identical(y_reinject_mode, "none") || isTRUE(y_reinject_x)) {
      stop("Structured/additive y_generator variants currently do not combine with y_reinject_mode/y_reinject_x.")
    }
  }
  if (identical(y_generator, "xgated")) {
    if (!identical(y_reinject_mode, "none") || isTRUE(y_reinject_x)) {
      stop("y_generator='xgated' currently does not combine with y_reinject_mode/y_reinject_x.")
    }
  }

  if (!is.null(split_seed)) set.seed(as.integer(split_seed))

  if (identical(y_family, "binary")) {
    y_raw <- as.numeric(data$Y)
    if (any(is.na(y_raw)) || !all(y_raw %in% c(0, 1))) {
      stop("y_family='binary' requires Y in {0,1} without NA.")
    }
  } else if (identical(y_family, "categorical")) {
    y_raw <- as.numeric(data$Y)
    if (any(is.na(y_raw))) stop("y_family='categorical' requires Y without NA.")
    if (is.null(y_levels)) {
      y_levels <- sort(unique(y_raw))
    } else {
      y_levels <- as.numeric(y_levels)
    }
    if (length(y_levels) < 2L || any(is.na(y_levels)) || anyDuplicated(y_levels)) {
      stop("y_family='categorical' requires at least two unique, non-missing y_levels.")
    }
    if (any(is.na(match(y_raw, y_levels)))) {
      stop("Some observed Y values are not present in y_levels.")
    }
  }

  n <- nrow(data)
  idx <- sample(seq_len(n), size = floor((1 - val_p) * n))

  train_raw <- data[idx, , drop = FALSE]
  val_raw   <- data[-idx, , drop = FALSE]

  ## -----------------------------
  ## 1. Detect mediators and standardize (optional)
  ## -----------------------------
  m_vars <- grep("^M", names(train_raw), value = TRUE)
  S <- length(m_vars)
  if (S < 1L) stop("No mediator columns found: expect names starting with 'M' (e.g., M1, M2, ...).")

  vars_core <- c(m_vars, "Y", c_vars)

  stats <- lapply(train_raw[vars_core], function(v) c(mean = mean(v), sd = sd(v)))
  mu  <- sapply(stats, `[[`, "mean")
  sdd <- sapply(stats, `[[`, "sd")
  bad_sd <- !is.finite(sdd) | sdd == 0
  if (any(bad_sd)) {
    sdd[bad_sd] <- 1
    if (!silent) {
      cat(
        sprintf(
          "dcma(): replaced non-finite/zero sd with 1 for: %s\n",
          paste(names(sdd)[bad_sd], collapse = ", ")
        )
      )
    }
  }

  std_df <- function(df) {
    core <- df[vars_core]
    if (standardize) {
      core <- as.data.frame(scale(core, center = mu[vars_core], scale = sdd[vars_core]))
    }
    ## Keep X unstandardized and robustly coerce to 0/1 numeric
    if (is.factor(df$X)) {
      core$X <- as.numeric(df$X) - 1L
    } else {
      core$X <- as.numeric(df$X)
    }
    if (identical(y_family, "binary")) {
      # Keep binary Y on original 0/1 scale (do not standardize).
      core$Y <- as.numeric(df$Y)
    } else if (identical(y_family, "categorical")) {
      # Keep categorical Y as 1-based class ids for torch cross-entropy.
      core$Y <- as.integer(match(as.numeric(df$Y), y_levels))
    }
    core
  }

  train_std <- std_df(train_raw)
  val_std   <- std_df(val_raw)

  ## -----------------------------
  ## 2. Build torch tensors for training and validation
  ## -----------------------------
  X      <- as.matrix(train_std[, "X", drop = FALSE])
  M_all  <- as.matrix(train_std[, m_vars, drop = FALSE])  # n x S
  Y      <- if (identical(y_family, "categorical")) {
    as.integer(train_std$Y)
  } else {
    as.matrix(train_std[, "Y", drop = FALSE])
  }

  ## IMPORTANT: keep covariate matrix colnames
  C <- as.matrix(train_std[, c_vars, drop = FALSE])
  colnames(C) <- c_vars

  ## Also store raw-scale covariates (useful for reconstruction defaults)
  C_raw <- as.matrix(train_raw[, c_vars, drop = FALSE])
  colnames(C_raw) <- c_vars

  p <- ncol(C)
  n <- nrow(X)

  X_t     <- torch_tensor(X,     dtype = torch_float(), device = device)
  M_all_t <- torch_tensor(M_all, dtype = torch_float(), device = device)
  Y_t     <- if (identical(y_family, "categorical")) {
    torch_tensor(Y, dtype = torch_long(), device = device)
  } else {
    torch_tensor(Y, dtype = torch_float(), device = device)
  }
  C_t     <- torch_tensor(C,     dtype = torch_float(), device = device)
  XC_t    <- torch_cat(list(X_t, C_t), dim = 2)

  ## Validation tensors
  Xv     <- as.matrix(val_std$X)
  Mv_all <- as.matrix(val_std[, m_vars, drop = FALSE])
  Yv     <- if (identical(y_family, "categorical")) {
    as.integer(val_std$Y)
  } else {
    as.matrix(val_std$Y)
  }

  Cv <- as.matrix(val_std[, c_vars, drop = FALSE])
  colnames(Cv) <- c_vars

  Cv_raw <- as.matrix(val_raw[, c_vars, drop = FALSE])
  colnames(Cv_raw) <- c_vars

  nv     <- nrow(Xv)

  if (nv > 0L) {
    Xv_t     <- torch_tensor(Xv,     dtype = torch_float(), device = device)
    Mv_all_t <- torch_tensor(Mv_all, dtype = torch_float(), device = device)
    Yv_t     <- if (identical(y_family, "categorical")) {
      torch_tensor(Yv, dtype = torch_long(), device = device)
    } else {
      torch_tensor(Yv, dtype = torch_float(), device = device)
    }
    Cv_t     <- torch_tensor(Cv,     dtype = torch_float(), device = device)
    XCv_t    <- torch_cat(list(Xv_t, Cv_t), dim = 2)
  } else {
    Xv_t <- Mv_all_t <- Yv_t <- Cv_t <- XCv_t <- NULL
  }

  compute_y_stage_loss <- function(target_t, pred_list, x_group_num) {
    if (!isTRUE(y_balance_groups)) {
      if (identical(y_family, "binary")) {
        loss_terms <- lapply(pred_list, function(logits) {
          nnf_binary_cross_entropy_with_logits(logits, target_t)
        })
        return(torch_stack(loss_terms)$mean())
      } else if (identical(y_family, "categorical")) {
        loss_terms <- lapply(pred_list, function(logits) {
          nnf_cross_entropy(logits, target_t)
        })
        return(torch_stack(loss_terms)$mean())
      }
      return(do.call(energyloss_es, c(list(target_t), pred_list, list(beta = beta_Y))))
    }

    idx0 <- which(x_group_num == 0)
    idx1 <- which(x_group_num == 1)
    if (length(idx0) == 0L || length(idx1) == 0L) {
      if (identical(y_family, "binary")) {
        loss_terms <- lapply(pred_list, function(logits) {
          nnf_binary_cross_entropy_with_logits(logits, target_t)
        })
        return(torch_stack(loss_terms)$mean())
      } else if (identical(y_family, "categorical")) {
        loss_terms <- lapply(pred_list, function(logits) {
          nnf_cross_entropy(logits, target_t)
        })
        return(torch_stack(loss_terms)$mean())
      }
      return(do.call(energyloss_es, c(list(target_t), pred_list, list(beta = beta_Y))))
    }

    idx0_t <- torch_tensor(as.integer(idx0), dtype = torch_long(), device = device)
    idx1_t <- torch_tensor(as.integer(idx1), dtype = torch_long(), device = device)
    target0 <- if (identical(y_family, "categorical")) target_t[idx0_t] else target_t[idx0_t, ]
    target1 <- if (identical(y_family, "categorical")) target_t[idx1_t] else target_t[idx1_t, ]
    pred0 <- lapply(pred_list, function(x) x[idx0_t, ])
    pred1 <- lapply(pred_list, function(x) x[idx1_t, ])

    if (identical(y_family, "binary")) {
      loss0 <- torch_stack(lapply(pred0, function(logits) {
        nnf_binary_cross_entropy_with_logits(logits, target0)
      }))$mean()
      loss1 <- torch_stack(lapply(pred1, function(logits) {
        nnf_binary_cross_entropy_with_logits(logits, target1)
      }))$mean()
    } else if (identical(y_family, "categorical")) {
      loss0 <- torch_stack(lapply(pred0, function(logits) {
        nnf_cross_entropy(logits, target0)
      }))$mean()
      loss1 <- torch_stack(lapply(pred1, function(logits) {
        nnf_cross_entropy(logits, target1)
      }))$mean()
    } else {
      loss0 <- do.call(energyloss_es, c(list(target0), pred0, list(beta = beta_Y)))
      loss1 <- do.call(energyloss_es, c(list(target1), pred1, list(beta = beta_Y)))
    }
    0.5 * (loss0 + loss1)
  }

  ## -----------------------------
  ## 3. Define generator architectures
  ## -----------------------------
  in_f <- as.integer(1L + p)      # inputs: X + C
  in_g <- as.integer(1L + p + S)  # inputs: X + C + M
  out_g <- if (identical(y_family, "categorical")) length(y_levels) else 1L

  hd_M <- if (is.null(hidden_dim_M)) hidden_dim else hidden_dim_M
  nl_M <- if (is.null(num_layer_M)) num_layer else num_layer_M
  hd_Y <- if (is.null(hidden_dim_Y)) hidden_dim else hidden_dim_Y
  nl_Y <- if (is.null(num_layer_Y)) num_layer else num_layer_Y
  if (y_reinject_mode %in% c("obs2", "film2", "ln2", "adaln2") && nl_Y < 3L) {
    stop("y_reinject_mode='", y_reinject_mode, "' requires num_layer_Y >= 3 (or num_layer >= 3) so the second-layer modification is not silently skipped.")
  }

  if (!is.null(init_seed)) {
    torch::torch_manual_seed(as.integer(init_seed))
  }

  if (identical(m_generator, "shared")) {
    gen_f <- nn_model_shared_resid_mediator(
      in_dim = in_f,
      noise_dim = epsm_dim,
      hidden_dim = hd_M,
      out_dim = S,
      shared_noise_dim = m_shared_noise_dim
    )$to(device = device)
  } else if (identical(m_generator, "separate")) {
    gen_f <- nn_model_separate_mediator(
      in_dim = in_f,
      noise_dim = epsm_dim,
      hidden_dim = hd_M,
      out_dim = S,
      num_layer = nl_M
    )$to(device = device)
  } else if (identical(m_reinject_mode, "obs")) {
    gen_f <- nn_model_reinject_obs(
      in_dim     = in_f,
      noise_dim  = epsm_dim,
      hidden_dim = hd_M,
      out_dim    = S,
      num_layer  = nl_M,
      obs_dim    = in_f
    )$to(device = device)
  } else {
    gen_f <- nn_model(
      in_dim     = in_f,
      noise_dim  = epsm_dim,
      hidden_dim = hd_M,
      out_dim    = S,
      num_layer  = nl_M
    )$to(device = device)
  }

  if (identical(y_generator, "locscale")) {
    gen_g <- nn_model_locscale_outcome(
      in_dim = in_g,
      noise_dim = epsy_dim,
      hidden_dim = hd_Y
    )$to(device = device)
  } else if (identical(y_generator, "additive_homosk")) {
    gen_g <- nn_model_additive_outcome(
      in_dim = in_g,
      hidden_dim = hd_Y,
      heteroskedastic = FALSE
    )$to(device = device)
  } else if (identical(y_generator, "additive_heterosk")) {
    gen_g <- nn_model_additive_outcome(
      in_dim = in_g,
      hidden_dim = hd_Y,
      heteroskedastic = TRUE
    )$to(device = device)
  } else if (identical(y_generator, "xgated")) {
    gen_g <- nn_model_xgated_heads(
      in_dim = in_g,
      noise_dim = epsy_dim,
      hidden_dim = hd_Y,
      out_dim = out_g,
      num_layer = nl_Y
    )$to(device = device)
  } else if (identical(y_reinject_mode, "x")) {
    gen_g <- nn_model_reinject_x(
      in_dim     = in_g,
      noise_dim  = epsy_dim,
      hidden_dim = hd_Y,
      out_dim    = out_g,
      num_layer  = nl_Y,
      x_dim      = 1L
    )$to(device = device)
  } else if (identical(y_reinject_mode, "obs")) {
    gen_g <- nn_model_reinject_obs(
      in_dim     = in_g,
      noise_dim  = epsy_dim,
      hidden_dim = hd_Y,
      out_dim    = out_g,
      num_layer  = nl_Y,
      obs_dim    = in_g
    )$to(device = device)
  } else if (identical(y_reinject_mode, "obs2")) {
    gen_g <- nn_model_reinject_obs_second(
      in_dim     = in_g,
      noise_dim  = epsy_dim,
      hidden_dim = hd_Y,
      out_dim    = out_g,
      num_layer  = nl_Y,
      obs_dim    = in_g
    )$to(device = device)
  } else if (identical(y_reinject_mode, "film2")) {
    gen_g <- nn_model_film_obs_second(
      in_dim     = in_g,
      noise_dim  = epsy_dim,
      hidden_dim = hd_Y,
      out_dim    = out_g,
      num_layer  = nl_Y,
      obs_dim    = in_g
    )$to(device = device)
  } else if (identical(y_reinject_mode, "ln2")) {
    gen_g <- nn_model_ln_second(
      in_dim     = in_g,
      noise_dim  = epsy_dim,
      hidden_dim = hd_Y,
      out_dim    = out_g,
      num_layer  = nl_Y
    )$to(device = device)
  } else if (identical(y_reinject_mode, "adaln2")) {
    gen_g <- nn_model_adaln_obs_second(
      in_dim     = in_g,
      noise_dim  = epsy_dim,
      hidden_dim = hd_Y,
      out_dim    = out_g,
      num_layer  = nl_Y,
      obs_dim    = in_g
    )$to(device = device)
  } else {
    gen_g <- nn_model(
      in_dim     = in_g,
      noise_dim  = epsy_dim,
      hidden_dim = hd_Y,
      out_dim    = out_g,
      num_layer  = nl_Y
    )$to(device = device)
  }

  opt_f <- optim_adam(gen_f$parameters, lr = lr)
  opt_g <- optim_adam(gen_g$parameters, lr = lr)

  loss_M_hist     <- rep(NA_real_, epochs_M)
  val_loss_M_hist <- rep(NA_real_, epochs_M)
  loss_Y_hist     <- rep(NA_real_, epochs_Y)
  val_loss_Y_hist <- rep(NA_real_, epochs_Y)

  state_list <- list()

  sample_eps <- function(nrow, dim) {
    torch_randn(c(nrow, dim), device = device)
  }
  epsm_input_dim <- if (identical(m_generator, "separate")) epsm_dim * S else epsm_dim

  mediator_stage_loss <- function(target_t, pred_list) {
    if (!identical(m_generator, "separate")) {
      return(do.call(energyloss_es, c(list(target_t), pred_list, list(beta = beta_M))))
    }
    loss_terms <- vector("list", S)
    for (j in seq_len(S)) {
      target_j <- target_t[, j, drop = FALSE]
      pred_j <- lapply(pred_list, function(pred) pred[, j, drop = FALSE])
      loss_terms[[j]] <- do.call(energyloss_es, c(list(target_j), pred_j, list(beta = beta_M)))
    }
    torch_stack(loss_terms)$mean()
  }

  ## -----------------------------
  ## 4. Train mediator generator (M stage) with ES + early stopping
  ## -----------------------------
  best_val_M   <- Inf
  best_epoch_M <- 0L
  best_state_f <- NULL

  for (epoch in seq_len(epochs_M)) {
    opt_f$zero_grad()

    eps_list <- lapply(seq_len(K_es), function(k) sample_eps(n, epsm_input_dim))
    M_hat_list <- vector("list", K_es)
    for (k in seq_len(K_es)) {
      input_fk <- torch_cat(list(XC_t, eps_list[[k]]), dim = 2)
      M_hat_list[[k]] <- gen_f(input_fk)
    }

    if (epoch %in% save_epochs) {
      key <- as.character(epoch)
      if (is.null(state_list[[key]])) state_list[[key]] <- list()
      state_list[[key]][["M_hat"]] <- as.matrix(M_hat_list[[1]]$to(device = "cpu"))
      if (K_es >= 2L) {
        state_list[[key]][["M_hat1"]] <- as.matrix(M_hat_list[[2]]$to(device = "cpu"))
      }
    }

    loss_joint <- mediator_stage_loss(M_all_t, M_hat_list)
    loss_joint$backward()
    if (!is.na(grad_clip_norm)) {
      nn_utils_clip_grad_norm_(gen_f$parameters, max_norm = grad_clip_norm)
    }
    opt_f$step()

    loss_M_hist[epoch] <- loss_joint$item()

    if (nv > 0L && epoch %% val_freq == 0L) {
      with_no_grad({
        epsv_list <- lapply(seq_len(K_es), function(k) sample_eps(nv, epsm_input_dim))
        Mv_hat_list <- vector("list", K_es)
        for (k in seq_len(K_es)) {
          input_fvk <- torch_cat(list(XCv_t, epsv_list[[k]]), dim = 2)
          Mv_hat_list[[k]] <- gen_f(input_fvk)
        }
        val_loss_M <- mediator_stage_loss(Mv_all_t, Mv_hat_list)
        cur_val <- val_loss_M$item()
        val_loss_M_hist[epoch] <- cur_val

        if (is.infinite(best_val_M) || cur_val + min_delta_M < best_val_M) {
          best_val_M   <- cur_val
          best_epoch_M <- epoch
          best_state_f <- gen_f$state_dict()
        } else if (patience_M > 0L && (epoch - best_epoch_M) >= patience_M) {
          if (!silent) {
            cat(sprintf("[M] early stopping at epoch %d best epoch %d val ES = %.4f\n",
                        epoch, best_epoch_M, best_val_M))
          }
          break
        }
      })
    }

    if (!silent && epoch %% 50L == 0L) {
      cat(sprintf("[M] epoch %d train ES = %.4f val ES = %.4f\n",
                  epoch,
                  loss_M_hist[epoch],
                  ifelse(is.na(val_loss_M_hist[epoch]), NA_real_, val_loss_M_hist[epoch])))
    }
  }

  if (!is.null(best_state_f)) {
    gen_f$load_state_dict(best_state_f)
  }

  ## -----------------------------
  ## 5. Train outcome generator (Y stage) with ES + early stopping
  ## -----------------------------
  best_val_Y   <- Inf
  best_epoch_Y <- 0L
  best_state_g <- NULL

  for (epoch in seq_len(epochs_Y)) {
    opt_g$zero_grad()

    eps_list <- lapply(seq_len(K_es), function(k) sample_eps(n, epsy_dim))
    Y_hat_list <- vector("list", K_es)
    for (k in seq_len(K_es)) {
      input_gk <- torch_cat(list(XC_t, M_all_t, eps_list[[k]]), dim = 2)
      Y_hat_list[[k]] <- gen_g(input_gk)
    }

    loss_y <- compute_y_stage_loss(Y_t, Y_hat_list, as.numeric(train_std$X))
    loss_y$backward()
    if (!is.na(grad_clip_norm)) {
      nn_utils_clip_grad_norm_(gen_g$parameters, max_norm = grad_clip_norm)
    }
    opt_g$step()

    loss_Y_hist[epoch] <- loss_y$item()

    if (nv > 0L && epoch %% val_freq == 0L) {
      with_no_grad({
        epsv_list <- lapply(seq_len(K_es), function(k) sample_eps(nv, epsy_dim))
        Yv_hat_list <- vector("list", K_es)
        for (k in seq_len(K_es)) {
          input_gvk <- torch_cat(list(XCv_t, Mv_all_t, epsv_list[[k]]), dim = 2)
          Yv_hat_list[[k]] <- gen_g(input_gvk)
        }
        val_loss_Y <- compute_y_stage_loss(Yv_t, Yv_hat_list, as.numeric(val_std$X))
        cur_val <- val_loss_Y$item()
        val_loss_Y_hist[epoch] <- cur_val

        if (is.infinite(best_val_Y) || cur_val + min_delta_Y < best_val_Y) {
          best_val_Y   <- cur_val
          best_epoch_Y <- epoch
          best_state_g <- gen_g$state_dict()
        } else if (patience_Y > 0L && (epoch - best_epoch_Y) >= patience_Y) {
          if (!silent) {
            cat(sprintf("[Y] early stopping at epoch %d best epoch %d val ES = %.4f\n",
                        epoch, best_epoch_Y, best_val_Y))
          }
          break
        }
      })
    }

    if (!silent && epoch %% 50L == 0L) {
      cat(sprintf("[Y] epoch %d train ES = %.4f val ES = %.4f\n",
                  epoch,
                  loss_Y_hist[epoch],
                  ifelse(is.na(val_loss_Y_hist[epoch]), NA_real_, val_loss_Y_hist[epoch])))
    }
  }

  if (!is.null(best_state_g)) {
    gen_g$load_state_dict(best_state_g)
  }

  ## -----------------------------
  ## 6. Return fitted object
  ## -----------------------------
  DCMA <- list(
    gen_f       = gen_f,
    gen_g       = gen_g,
    loss_M      = loss_M_hist,
    loss_Y      = loss_Y_hist,
    val_loss_M  = val_loss_M_hist,
    val_loss_Y  = val_loss_Y_hist,
    state_list  = state_list,

    X           = X,
    C           = C,       # training-scale covariates (standardized if standardize=TRUE)
    C_raw       = C_raw,   # raw-scale covariates (always original scale, preserves colnames)
    Cv          = Cv,
    Cv_raw      = Cv_raw,

    M_obs       = M_all,
    Y           = Y,
    m_vars      = m_vars,
    epsm_dim    = epsm_dim,
    epsm_input_dim = epsm_input_dim,
    epsy_dim    = epsy_dim,
    mu          = mu,
    sdd         = sdd,
    vars_core   = vars_core,
    standardize = standardize,
    K_es        = K_es,
    y_family    = y_family,
    y_levels    = if (identical(y_family, "categorical")) y_levels else NULL,
    config      = list(
      c_vars     = c_vars,
      hidden_dim = hidden_dim,
      num_layer  = num_layer,
      hidden_dim_M = hd_M,
      num_layer_M = nl_M,
      hidden_dim_Y = hd_Y,
      num_layer_Y = nl_Y,
      m_generator = m_generator,
      m_shared_noise_dim = m_shared_noise_dim,
      m_reinject_mode = m_reinject_mode,
      epochs_M   = epochs_M,
      epochs_Y   = epochs_Y,
      lr         = lr,
      beta       = beta,
      beta_M     = beta_M,
      beta_Y     = beta_Y,
      val_p      = val_p,
      val_freq   = val_freq,
      patience   = patience,
      min_delta  = min_delta,
      patience_M = patience_M,
      patience_Y = patience_Y,
      min_delta_M = min_delta_M,
      min_delta_Y = min_delta_Y,
      y_balance_groups = y_balance_groups,
      grad_clip_norm = grad_clip_norm,
      y_reinject_x = y_reinject_x,
      y_generator = y_generator,
      y_reinject_mode = y_reinject_mode,
      split_seed = split_seed,
      init_seed = init_seed,
      device     = device,
      y_family   = y_family,
      y_levels   = if (identical(y_family, "categorical")) y_levels else NULL
    )
  )
  class(DCMA) <- "DCMA"
  DCMA
}
