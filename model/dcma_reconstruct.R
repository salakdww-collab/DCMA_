# Reconstruct interventional outcome distributions via forward simulation (DCMA)
dcma_reconstruct_interventions <- function(
  fit,
  B      = 50L,
  a0     = 0,
  a1     = 1,
  a_ipse = 1,
  C_mat  = NULL,
  seed   = NULL,
  device = "mps",
  product_mode = c("permute", "redraw")
) {
  ## -----------------------------
  ## 0. Basic checks and setup
  ## -----------------------------
  stopifnot(inherits(fit, "DCMA"))
  B <- as.integer(B)
  product_mode <- match.arg(product_mode)
  stopifnot(B >= 1L)
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1L, is.finite(seed))
    set.seed(as.integer(seed))
    torch::torch_manual_seed(as.integer(seed))
  }

  if (is.null(fit$config) || is.null(fit$config$c_vars)) {
    stop("fit$config$c_vars is missing. Please ensure dcma() stores c_vars in fit$config$c_vars.")
  }
  c_vars <- fit$config$c_vars

  ## -----------------------------
  ## 0.1 Resolve C_mat: default to fit-stored covariates
  ## -----------------------------
  if (is.null(C_mat)) {
    if (!is.null(fit$C_raw)) {
      C_mat <- fit$C_raw
    } else if (!is.null(fit$C)) {
      C_mat <- fit$C
    } else {
      stop("C_mat is NULL and neither fit$C_raw nor fit$C is available. Please provide C_mat.")
    }
  }

  ## Preserve colnames when coercing to matrix
  if (is.data.frame(C_mat)) {
    if (is.null(colnames(C_mat))) colnames(C_mat) <- c_vars
    C_df <- C_mat
    C_mat <- as.matrix(C_mat)
    colnames(C_mat) <- colnames(C_df)
  } else {
    C_mat <- as.matrix(C_mat)
  }

  if (is.null(colnames(C_mat))) {
    colnames(C_mat) <- c_vars
  }

  ## Make sure covariate columns match training covariates
  if (!all(c_vars %in% colnames(C_mat))) {
    stop("C_mat is missing some covariates used in training: ",
         paste(setdiff(c_vars, colnames(C_mat)), collapse = ", "))
  }
  ## Reorder to training order and drop extras if any
  C_mat <- C_mat[, c_vars, drop = FALSE]
  storage.mode(C_mat) <- "numeric"
  if (!is.numeric(C_mat)) stop("C_mat must be numeric after coercion.")

  ## Extract fitted generators and dimensions
  gen_f <- fit$gen_f$to(device = device)
  gen_g <- fit$gen_g$to(device = device)

  n <- nrow(C_mat)
  p <- ncol(C_mat)
  S <- length(fit$m_vars)

  epsm_dim <- fit$epsm_dim
  epsm_input_dim <- if (!is.null(fit$epsm_input_dim)) {
    fit$epsm_input_dim
  } else if (!is.null(fit$config$m_generator) && identical(fit$config$m_generator, "separate")) {
    fit$epsm_dim * S
  } else {
    fit$epsm_dim
  }
  epsy_dim <- fit$epsy_dim
  y_family <- if (!is.null(fit$y_family)) fit$y_family else "continuous"
  y_levels <- if (!is.null(fit$y_levels)) fit$y_levels else fit$config$y_levels
  if (identical(y_family, "categorical")) {
    if (is.null(y_levels) || length(y_levels) < 2L) {
      stop("fit$y_levels is required for y_family='categorical'.")
    }
    y_levels <- as.numeric(y_levels)
  }

  ## -----------------------------
  ## 1. Standardize covariates consistently with training (if needed)
  ## -----------------------------
  ## Important: C_mat at this point is on the original covariate scale unless user passed standardized.
  ## We standardize here whenever fit$standardize=TRUE using fit$mu/sdd.
  if (isTRUE(fit$standardize)) {
    if (is.null(fit$mu) || is.null(fit$sdd)) stop("fit$mu/fit$sdd missing.")

    miss_mu <- setdiff(c_vars, names(fit$mu))
    miss_sd <- setdiff(c_vars, names(fit$sdd))
    if (length(miss_mu) > 0L || length(miss_sd) > 0L) {
      stop("fit$mu/fit$sdd missing covariate names: ",
           paste(unique(c(miss_mu, miss_sd)), collapse = ", "))
    }

    muC <- as.numeric(fit$mu[c_vars])
    sdC <- as.numeric(fit$sdd[c_vars])

    sdC[is.na(sdC) | sdC == 0] <- 1

    C_mat <- sweep(C_mat, 2, muC, "-")
    C_mat <- sweep(C_mat, 2, sdC, "/")
    colnames(C_mat) <- c_vars
  }

  ## -----------------------------
  ## 2. Helper: draw mediators M_a (joint), returns n x B x S
  ## -----------------------------
  draw_M_given_a <- function(a) {
    X_vec <- matrix(a, nrow = n, ncol = 1)
    X_rep <- X_vec[rep(seq_len(n), each = B), , drop = FALSE]     # (nB) x 1
    C_rep <- C_mat[rep(seq_len(n), each = B), , drop = FALSE]     # (nB) x p

    XC_rep <- torch_tensor(
      cbind(X_rep, C_rep),
      dtype  = torch_float(),
      device = device
    )

    epsM  <- torch_randn(c(n * B, epsm_input_dim), device = device)
    input <- torch_cat(list(XC_rep, epsM), dim = 2)

    M_flat <- gen_f(input)  # (nB) x S
    M_arr  <- M_flat$to(device = "cpu")$view(c(n, B, S))
    as.array(M_arr)
  }

  ## -----------------------------
  ## 3. Helper: draw outcomes given mediators and exposure, returns n x B
  ## -----------------------------
  flatten_m_i_major <- function(M_array) {
    stopifnot(length(dim(M_array)) == 3L)
    S_local <- dim(M_array)[3]
    do.call(cbind, lapply(seq_len(S_local), function(j) as.vector(t(M_array[, , j]))))
  }

  draw_Y_given_M_a <- function(M_array, a) {
    stopifnot(length(dim(M_array)) == 3L)
    if (dim(M_array)[1] != n) stop("M_array first dimension must equal n.")

    B_local <- dim(M_array)[2]
    M_flat <- flatten_m_i_major(M_array)

    X_rep <- matrix(a, nrow = n * B_local, ncol = 1)
    C_rep <- C_mat[rep(seq_len(n), each = B_local), , drop = FALSE]

    XC_rep <- torch_tensor(
      cbind(X_rep, C_rep),
      dtype  = torch_float(),
      device = device
    )
    M_t  <- torch_tensor(M_flat, dtype = torch_float(), device = device)
    epsY <- torch_randn(c(n * B_local, epsy_dim), device = device)

    input  <- torch_cat(list(XC_rep, M_t, epsY), dim = 2)
    logits <- gen_g(input)  # (nB) x 1 or (nB) x K
    if (identical(y_family, "binary")) {
      probs <- torch_sigmoid(logits)
      Y_flat <- torch_bernoulli(probs)
    } else if (identical(y_family, "categorical")) {
      probs <- nnf_softmax(logits, dim = 2)
      class_idx <- as.integer(torch_multinomial(probs, num_samples = 1, replacement = TRUE)$to(device = "cpu"))
      Y_flat <- torch_tensor(matrix(y_levels[class_idx], ncol = 1), dtype = torch_float())
    } else {
      Y_flat <- logits
    }

    Y_arr <- Y_flat$to(device = "cpu")$view(c(n, B_local))
    as.array(Y_arr)
  }

  ## -----------------------------
  ## 4. Helper: build pathway-specific mediator pairs for multiple mediators
  ## -----------------------------
  build_pathway_M_pair <- function(M0, M1, s) {
    stopifnot(length(dim(M0)) == 3L, length(dim(M1)) == 3L)
    stopifnot(all(dim(M0) == dim(M1)))
    stopifnot(s >= 1L, s <= dim(M0)[3])
    stopifnot(dim(M0)[3] >= 2L)

    n0 <- dim(M0)[1]
    B0 <- dim(M0)[2]
    S0 <- dim(M0)[3]

    M0s <- array(NA_real_, dim = c(n0, B0, S0))
    M1s <- array(NA_real_, dim = c(n0, B0, S0))

    idx_low  <- if (s > 1L) 1:(s - 1L) else integer(0)
    idx_high <- if (s < S0) (s + 1L):S0 else integer(0)

    for (i in seq_len(n0)) {
      M0_i <- matrix(M0[i, , , drop = FALSE], nrow = B0, ncol = S0)
      M1_i <- matrix(M1[i, , , drop = FALSE], nrow = B0, ncol = S0)

      parts1 <- list()
      if (length(idx_low)) {
        parts1[[length(parts1) + 1L]] <- M0_i[sample(B0), idx_low, drop = FALSE]
      }
      parts1[[length(parts1) + 1L]] <- M1_i[sample(B0), s, drop = FALSE]
      if (length(idx_high)) {
        parts1[[length(parts1) + 1L]] <- M1_i[sample(B0), idx_high, drop = FALSE]
      }
      M1s[i, , ] <- do.call(cbind, parts1)

      parts0 <- list()
      if (length(idx_low)) {
        parts0[[length(parts0) + 1L]] <- M0_i[sample(B0), idx_low, drop = FALSE]
      }
      parts0[[length(parts0) + 1L]] <- M0_i[sample(B0), s, drop = FALSE]
      if (length(idx_high)) {
        parts0[[length(parts0) + 1L]] <- M1_i[sample(B0), idx_high, drop = FALSE]
      }
      M0s[i, , ] <- do.call(cbind, parts0)
    }

    list(M0s = M0s, M1s = M1s)
  }

  build_pathway_M_pair_redraw <- function(s) {
    stopifnot(s >= 1L, s <= S)
    stopifnot(S >= 2L)

    idx_low  <- if (s > 1L) 1:(s - 1L) else integer(0)
    idx_high <- if (s < S) (s + 1L):S else integer(0)

    M1_parts <- list()
    if (length(idx_low)) {
      M1_parts[[length(M1_parts) + 1L]] <- draw_M_given_a(a0)[, , idx_low, drop = FALSE]
    }
    M1_parts[[length(M1_parts) + 1L]] <- draw_M_given_a(a1)[, , s, drop = FALSE]
    if (length(idx_high)) {
      M1_parts[[length(M1_parts) + 1L]] <- draw_M_given_a(a1)[, , idx_high, drop = FALSE]
    }

    M0_parts <- list()
    if (length(idx_low)) {
      M0_parts[[length(M0_parts) + 1L]] <- draw_M_given_a(a0)[, , idx_low, drop = FALSE]
    }
    M0_parts[[length(M0_parts) + 1L]] <- draw_M_given_a(a0)[, , s, drop = FALSE]
    if (length(idx_high)) {
      M0_parts[[length(M0_parts) + 1L]] <- draw_M_given_a(a1)[, , idx_high, drop = FALSE]
    }

    bind_parts <- function(parts) {
      out <- array(NA_real_, dim = c(n, B, S))
      pos <- 1L
      for (part in parts) {
        width <- dim(part)[3]
        out[, , pos:(pos + width - 1L)] <- part
        pos <- pos + width
      }
      out
    }

    list(M0s = bind_parts(M0_parts), M1s = bind_parts(M1_parts))
  }

  ## -----------------------------
  ## 5. Step 1: draw joint mediators under a0 and a1
  ## -----------------------------
  M0 <- draw_M_given_a(a0)  # n x B x S
  M1 <- draw_M_given_a(a1)  # n x B x S

  ## -----------------------------
  ## 6. Step 2: draw joint interventional outcomes
  ## -----------------------------
  Y_1M1 <- draw_Y_given_M_a(M1, a1)
  Y_0M0 <- draw_Y_given_M_a(M0, a0)
  Y_1M0 <- draw_Y_given_M_a(M0, a1)
  Y_0M1 <- draw_Y_given_M_a(M1, a0)

  ## -----------------------------
  ## 7. Step 3: pathway-specific outcomes (per mediator)
  ## -----------------------------
  path_mediators <- vector("list", S)
  path_outcomes  <- vector("list", S)

  if (S == 1L) {
    Y_a_M0s <- draw_Y_given_M_a(M0, a_ipse)
    Y_a_M1s <- draw_Y_given_M_a(M1, a_ipse)

    path_mediators[[1L]] <- list(M0_s = M0, M1_s = M1)
    path_outcomes[[1L]]  <- list(Y_a_M0s = Y_a_M0s, Y_a_M1s = Y_a_M1s)
  } else {
    for (s in seq_len(S)) {
      pair_s <- if (identical(product_mode, "redraw")) {
        build_pathway_M_pair_redraw(s)
      } else {
        build_pathway_M_pair(M0, M1, s)
      }
      M0_s   <- pair_s$M0s
      M1_s   <- pair_s$M1s

      Y_a_M0s <- draw_Y_given_M_a(M0_s, a_ipse)
      Y_a_M1s <- draw_Y_given_M_a(M1_s, a_ipse)

      path_mediators[[s]] <- list(M0_s = M0_s, M1_s = M1_s)
      path_outcomes[[s]]  <- list(Y_a_M0s = Y_a_M0s, Y_a_M1s = Y_a_M1s)
    }
  }

  names(path_mediators) <- fit$m_vars
  names(path_outcomes)  <- fit$m_vars

  ## -----------------------------
  ## 8. Assemble output object
  ## -----------------------------
  out <- list(
    B      = B,
    a0     = a0,
    a1     = a1,
    a_ipse = a_ipse,
    m_vars = fit$m_vars,
    mediators = list(
      M_a0 = M0,
      M_a1 = M1,
      path = path_mediators
    ),
    outcomes = list(
      Y_1M1 = Y_1M1,
      Y_0M0 = Y_0M0,
      Y_1M0 = Y_1M0,
      Y_0M1 = Y_0M1,
      path  = path_outcomes
    ),
    scale_info = list(
      covariates_standardized = isTRUE(fit$standardize),
      mediators_on_original_scale = !isTRUE(fit$standardize),
      outcomes_on_original_scale = !isTRUE(fit$standardize) || !identical(y_family, "continuous"),
      product_mode = product_mode
    )
  )

  ## -----------------------------
  ## 9. Unstandardize outcomes back to original scale (if trained on standardized Y)
  ## -----------------------------
  if (isTRUE(fit$standardize) && identical(y_family, "continuous")) {
    if (!("Y" %in% names(fit$mu)) || !("Y" %in% names(fit$sdd))) {
      stop("fit$mu/fit$sdd missing Y for unstandardization.")
    }
    sdY <- as.numeric(fit$sdd["Y"])
    muY <- as.numeric(fit$mu["Y"])
    if (is.na(sdY) || sdY == 0) sdY <- 1

    unstd_Y <- function(mat) mat * sdY + muY

    out$outcomes$Y_1M1 <- unstd_Y(out$outcomes$Y_1M1)
    out$outcomes$Y_0M0 <- unstd_Y(out$outcomes$Y_0M0)
    out$outcomes$Y_1M0 <- unstd_Y(out$outcomes$Y_1M0)
    out$outcomes$Y_0M1 <- unstd_Y(out$outcomes$Y_0M1)

    out$outcomes$path <- lapply(out$outcomes$path, function(obj) {
      list(
        Y_a_M0s = unstd_Y(obj$Y_a_M0s),
        Y_a_M1s = unstd_Y(obj$Y_a_M1s)
      )
    })
    out$scale_info$outcomes_on_original_scale <- TRUE
  }

  ## -----------------------------
  ## 10. Unstandardize mediators back to original scale (if trained on standardized M)
  ## -----------------------------
  if (isTRUE(fit$standardize)) {
    if (!all(fit$m_vars %in% names(fit$mu)) || !all(fit$m_vars %in% names(fit$sdd))) {
      stop("fit$mu/fit$sdd missing mediator names for unstandardization.")
    }

    muM <- as.numeric(fit$mu[fit$m_vars])
    sdM <- as.numeric(fit$sdd[fit$m_vars])
    sdM[is.na(sdM) | sdM == 0] <- 1

    unstd_M <- function(arr3) {
      stopifnot(length(dim(arr3)) == 3L, dim(arr3)[3] == length(muM))
      out_arr <- arr3
      for (j in seq_along(muM)) {
        out_arr[, , j] <- out_arr[, , j] * sdM[j] + muM[j]
      }
      out_arr
    }

    out$mediators$M_a0 <- unstd_M(out$mediators$M_a0)
    out$mediators$M_a1 <- unstd_M(out$mediators$M_a1)
    out$mediators$path <- lapply(out$mediators$path, function(obj) {
      list(
        M0_s = unstd_M(obj$M0_s),
        M1_s = unstd_M(obj$M1_s)
      )
    })
    out$scale_info$mediators_on_original_scale <- TRUE
  }

  class(out) <- "DCMA_interventions"
  out
}
