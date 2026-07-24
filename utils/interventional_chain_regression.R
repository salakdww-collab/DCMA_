# Interventional mediation via parametric chain regression + Monte Carlo
#
# This implements a practical regression-based g-computation estimator for
# joint interventional effects with multiple mediators.

is_binary01 <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  length(x) > 0L && all(x %in% c(0, 1))
}

fit_interventional_chain <- function(
  data,
  c_vars = "Z",
  x_var = "X",
  y_var = "Y",
  m_vars = NULL,
  mediator_formulas = NULL,
  outcome_formula = NULL
) {
  stopifnot(is.data.frame(data))
  stopifnot(all(c(x_var, y_var, c_vars) %in% names(data)))

  if (is.null(m_vars)) {
    m_vars <- grep("^M", names(data), value = TRUE)
  }
  if (length(m_vars) < 1L) {
    stop("No mediators found. Expect columns named M1, M2, ...")
  }

  if (is.null(mediator_formulas)) {
    mediator_formulas <- vector("list", length(m_vars))
    names(mediator_formulas) <- m_vars
    for (j in seq_along(m_vars)) {
      rhs <- c(x_var, c_vars, m_vars[seq_len(j - 1L)])
      mediator_formulas[[j]] <- stats::as.formula(
        paste(m_vars[j], "~", paste(rhs, collapse = " + "))
      )
    }
  }

  mediator_models <- vector("list", length(m_vars))
  names(mediator_models) <- m_vars
  mediator_family <- setNames(rep("gaussian", length(m_vars)), m_vars)
  mediator_sigma <- setNames(rep(NA_real_, length(m_vars)), m_vars)

  for (j in seq_along(m_vars)) {
    mj <- m_vars[j]
    fm <- mediator_formulas[[j]]
    if (is_binary01(data[[mj]])) {
      mediator_family[[mj]] <- "binomial"
      mediator_models[[mj]] <- stats::glm(fm, data = data, family = stats::binomial(link = "logit"))
      mediator_sigma[[mj]] <- NA_real_
    } else {
      mediator_family[[mj]] <- "gaussian"
      fit_j <- stats::lm(fm, data = data)
      mediator_models[[mj]] <- fit_j
      s <- summary(fit_j)$sigma
      if (!is.finite(s) || s <= 0) s <- stats::sd(stats::residuals(fit_j))
      if (!is.finite(s) || s <= 0) s <- 1e-6
      mediator_sigma[[mj]] <- s
    }
  }

  if (is.null(outcome_formula)) {
    rhs <- c(x_var, c_vars, m_vars)
    outcome_formula <- stats::as.formula(
      paste(y_var, "~", paste(rhs, collapse = " + "))
    )
  }
  if (is_binary01(data[[y_var]])) {
    outcome_family <- "binomial"
    outcome_model <- stats::glm(outcome_formula, data = data, family = stats::binomial(link = "logit"))
    outcome_sigma <- NA_real_
  } else {
    outcome_family <- "gaussian"
    outcome_model <- stats::lm(outcome_formula, data = data)
    outcome_sigma <- summary(outcome_model)$sigma
    if (!is.finite(outcome_sigma) || outcome_sigma <= 0) {
      outcome_sigma <- stats::sd(stats::residuals(outcome_model))
    }
    if (!is.finite(outcome_sigma) || outcome_sigma <= 0) outcome_sigma <- 1e-6
  }

  list(
    x_var = x_var,
    y_var = y_var,
    c_vars = c_vars,
    m_vars = m_vars,
    mediator_models = mediator_models,
    mediator_family = mediator_family,
    mediator_sigma = mediator_sigma,
    outcome_model = outcome_model,
    outcome_family = outcome_family,
    outcome_sigma = outcome_sigma
  )
}

simulate_mediators_chain <- function(fit_obj, C_df, a, B = 50L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(all(fit_obj$c_vars %in% names(C_df)))

  n <- nrow(C_df)
  S <- length(fit_obj$m_vars)
  B <- as.integer(B)
  out <- array(
    NA_real_,
    dim = c(n, B, S),
    dimnames = list(NULL, NULL, fit_obj$m_vars)
  )

  for (b in seq_len(B)) {
    nd <- as.data.frame(C_df[, fit_obj$c_vars, drop = FALSE])
    nd[[fit_obj$x_var]] <- as.numeric(a)
    for (j in seq_along(fit_obj$m_vars)) {
      mj <- fit_obj$m_vars[j]
      fam <- fit_obj$mediator_family[[mj]]
      if (identical(fam, "binomial")) {
        prob <- as.numeric(stats::predict(fit_obj$mediator_models[[mj]], newdata = nd, type = "response"))
        prob <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
        draw <- stats::rbinom(n, size = 1L, prob = prob)
      } else {
        pred <- as.numeric(stats::predict(fit_obj$mediator_models[[mj]], newdata = nd))
        draw <- pred + stats::rnorm(n, mean = 0, sd = fit_obj$mediator_sigma[[mj]])
      }
      nd[[mj]] <- draw
      out[, b, j] <- draw
    }
  }
  out
}

simulate_outcomes_chain <- function(fit_obj, C_df, M_array, a, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(dim(M_array)) == 3L)

  n <- nrow(C_df)
  B <- dim(M_array)[2]
  out <- matrix(NA_real_, nrow = n, ncol = B)

  for (b in seq_len(B)) {
    nd <- as.data.frame(C_df[, fit_obj$c_vars, drop = FALSE])
    nd[[fit_obj$x_var]] <- as.numeric(a)
    for (j in seq_along(fit_obj$m_vars)) {
      nd[[fit_obj$m_vars[j]]] <- M_array[, b, j]
    }
    if (identical(fit_obj$outcome_family, "binomial")) {
      prob <- as.numeric(stats::predict(fit_obj$outcome_model, newdata = nd, type = "response"))
      prob <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
      out[, b] <- stats::rbinom(n, size = 1L, prob = prob)
    } else {
      pred <- as.numeric(stats::predict(fit_obj$outcome_model, newdata = nd))
      out[, b] <- pred + stats::rnorm(n, mean = 0, sd = fit_obj$outcome_sigma)
    }
  }
  out
}

build_pathway_M_pair_chain <- function(M0, M1, s) {
  stopifnot(length(dim(M0)) == 3L, length(dim(M1)) == 3L)
  stopifnot(all(dim(M0) == dim(M1)))
  stopifnot(s >= 1L, s <= dim(M0)[3])

  n0 <- dim(M0)[1]
  B0 <- dim(M0)[2]
  S0 <- dim(M0)[3]

  M0s <- array(NA_real_, dim = c(n0, B0, S0))
  M1s <- array(NA_real_, dim = c(n0, B0, S0))

  idx_low <- if (s > 1L) 1:(s - 1L) else integer(0)
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

reconstruct_interventions_chain <- function(
  fit_obj,
  data,
  B = 50L,
  seed = NULL,
  a0 = 0,
  a1 = 1,
  a_ipse = 1
) {
  stopifnot(all(fit_obj$c_vars %in% names(data)))
  C_df <- as.data.frame(data[, fit_obj$c_vars, drop = FALSE])

  base_seed <- if (is.null(seed)) sample.int(1e7, 1) else as.integer(seed)

  M0 <- simulate_mediators_chain(fit_obj, C_df, a = a0, B = B, seed = base_seed + 11L)
  M1 <- simulate_mediators_chain(fit_obj, C_df, a = a1, B = B, seed = base_seed + 13L)

  Y_0M0 <- simulate_outcomes_chain(fit_obj, C_df, M0, a = a0, seed = base_seed + 21L)
  Y_1M1 <- simulate_outcomes_chain(fit_obj, C_df, M1, a = a1, seed = base_seed + 23L)
  Y_1M0 <- simulate_outcomes_chain(fit_obj, C_df, M0, a = a1, seed = base_seed + 25L)
  Y_0M1 <- simulate_outcomes_chain(fit_obj, C_df, M1, a = a0, seed = base_seed + 27L)

  S <- length(fit_obj$m_vars)
  path_outcomes <- vector("list", S)
  names(path_outcomes) <- fit_obj$m_vars

  if (S == 1L) {
    path_outcomes[[1L]] <- list(
      Y_a_M0s = simulate_outcomes_chain(fit_obj, C_df, M0, a = a_ipse, seed = base_seed + 31L),
      Y_a_M1s = simulate_outcomes_chain(fit_obj, C_df, M1, a = a_ipse, seed = base_seed + 33L)
    )
  } else {
    for (s in seq_len(S)) {
      pair <- build_pathway_M_pair_chain(M0, M1, s)
      path_outcomes[[s]] <- list(
        Y_a_M0s = simulate_outcomes_chain(fit_obj, C_df, pair$M0s, a = a_ipse, seed = base_seed + 100L + 2L * s),
        Y_a_M1s = simulate_outcomes_chain(fit_obj, C_df, pair$M1s, a = a_ipse, seed = base_seed + 101L + 2L * s)
      )
    }
  }

  list(
    mediators = list(
      M_a0 = M0,
      M_a1 = M1
    ),
    outcomes = list(
      Y_1M1 = Y_1M1,
      Y_0M0 = Y_0M0,
      Y_1M0 = Y_1M0,
      Y_0M1 = Y_0M1,
      path = path_outcomes
    )
  )
}

compute_effect_pair_summary <- function(
  v1_true,
  v0_true,
  v1_est,
  v0_est,
  probs = c(0.1, 0.5, 0.9),
  n_ed = 2000L
) {
  v1_true <- as.numeric(v1_true)
  v0_true <- as.numeric(v0_true)
  v1_est <- as.numeric(v1_est)
  v0_est <- as.numeric(v0_est)

  q_true <- stats::quantile(v1_true, probs = probs) - stats::quantile(v0_true, probs = probs)
  q_est <- stats::quantile(v1_est, probs = probs) - stats::quantile(v0_est, probs = probs)

  idx1t <- sample.int(length(v1_true), min(n_ed, length(v1_true)))
  idx0t <- sample.int(length(v0_true), min(n_ed, length(v0_true)))
  idx1e <- sample.int(length(v1_est), min(n_ed, length(v1_est)))
  idx0e <- sample.int(length(v0_est), min(n_ed, length(v0_est)))

  tibble::tibble(
    mean_true = mean(v1_true) - mean(v0_true),
    mean_est = mean(v1_est) - mean(v0_est),
    q_true = list(q_true),
    q_est = list(q_est),
    ED_true = energy_distance(v1_true[idx1t], v0_true[idx0t]),
    ED_est = energy_distance(v1_est[idx1e], v0_est[idx0e])
  )
}

compute_interventional_basic_effects <- function(
  truth_outcomes,
  est_outcomes,
  probs = c(0.1, 0.5, 0.9),
  n_ed = 2000L,
  ed_seed = NULL
) {
  stopifnot(all(c("Y_1M1", "Y_0M0", "Y_1M0") %in% names(truth_outcomes)))
  stopifnot(all(c("Y_1M1", "Y_0M0", "Y_1M0") %in% names(est_outcomes)))
  if (!is.null(ed_seed)) {
    stopifnot(is.numeric(ed_seed), length(ed_seed) == 1L, is.finite(ed_seed))
    set.seed(as.integer(ed_seed))
  }

  dplyr::bind_rows(
    compute_effect_pair_summary(
      truth_outcomes$Y_1M1, truth_outcomes$Y_0M0,
      est_outcomes$Y_1M1, est_outcomes$Y_0M0,
      probs = probs, n_ed = n_ed
    ) |>
      dplyr::mutate(effect = "ITE", .before = 1),
    compute_effect_pair_summary(
      truth_outcomes$Y_1M0, truth_outcomes$Y_0M0,
      est_outcomes$Y_1M0, est_outcomes$Y_0M0,
      probs = probs, n_ed = n_ed
    ) |>
      dplyr::mutate(effect = "IDE", .before = 1),
    compute_effect_pair_summary(
      truth_outcomes$Y_1M1, truth_outcomes$Y_1M0,
      est_outcomes$Y_1M1, est_outcomes$Y_1M0,
      probs = probs, n_ed = n_ed
    ) |>
      dplyr::mutate(effect = "IIE", .before = 1)
  )
}

compute_path_specific_effects <- function(
  truth_outcomes,
  est_outcomes,
  probs = c(0.1, 0.5, 0.9),
  n_ed = 2000L,
  ed_seed = NULL,
  effect_prefix = "IPSE"
) {
  stopifnot("path" %in% names(truth_outcomes), "path" %in% names(est_outcomes))
  stopifnot(length(truth_outcomes$path) == length(est_outcomes$path))
  if (!is.null(ed_seed)) {
    stopifnot(is.numeric(ed_seed), length(ed_seed) == 1L, is.finite(ed_seed))
    set.seed(as.integer(ed_seed))
  }

  pieces <- vector("list", length(truth_outcomes$path))
  for (s in seq_along(truth_outcomes$path)) {
    t_obj <- truth_outcomes$path[[s]]
    e_obj <- est_outcomes$path[[s]]
    stopifnot(all(c("Y_a_M0s", "Y_a_M1s") %in% names(t_obj)))
    stopifnot(all(c("Y_a_M0s", "Y_a_M1s") %in% names(e_obj)))
    pieces[[s]] <- compute_effect_pair_summary(
      t_obj$Y_a_M1s, t_obj$Y_a_M0s,
      e_obj$Y_a_M1s, e_obj$Y_a_M0s,
      probs = probs, n_ed = n_ed
    ) |>
      dplyr::mutate(effect = paste0(effect_prefix, s), .before = 1)
  }
  dplyr::bind_rows(pieces)
}
