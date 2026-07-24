# Rich distributional effect summaries for interventional outcome pairs
#
# This file provides effect-level summaries used in the IHDP semi-synthetic
# benchmark. For each pair of distributions (v1, v0), we compute:
# - mean effect
# - quantile effects
# - exceedance effects at user-specified thresholds
# - divergence-style summaries: energy distance, Wasserstein-1, KL

metric_label_from_prob <- function(p) {
  sprintf("q%02d", as.integer(round(100 * p)))
}

wasserstein1_distance_1d <- function(x, y, n_prob = 512L) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) < 1L || length(y) < 1L) return(NA_real_)
  probs <- seq(0, 1, length.out = as.integer(n_prob))
  qx <- as.numeric(stats::quantile(x, probs = probs, names = FALSE, type = 8))
  qy <- as.numeric(stats::quantile(y, probs = probs, names = FALSE, type = 8))
  mean(abs(qx - qy))
}

kl_divergence_kde_1d <- function(x, y, n_grid = 512L, eps = 1e-8) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)

  xr <- range(c(x, y), finite = TRUE)
  if (!all(is.finite(xr)) || diff(xr) <= 0) return(0)

  dx <- stats::density(x, from = xr[1], to = xr[2], n = as.integer(n_grid))
  dy <- stats::density(y, from = xr[1], to = xr[2], n = as.integer(n_grid))

  grid_x <- dx$x
  p_raw <- pmax(dx$y, eps)
  q_raw <- pmax(dy$y, eps)
  step <- mean(diff(grid_x))

  p <- p_raw / sum(p_raw * step)
  q <- q_raw / sum(q_raw * step)

  sum(p * log(p / q)) * step
}

derive_tail_thresholds <- function(
  truth_outcomes,
  baseline_key = "Y_0M0",
  probs = c(0.5, 0.75)
) {
  stopifnot(baseline_key %in% names(truth_outcomes))
  vals <- as.numeric(truth_outcomes[[baseline_key]])
  thr <- as.numeric(stats::quantile(vals, probs = probs, names = FALSE, type = 8))
  names(thr) <- paste0("tail_q", sprintf("%02d", as.integer(round(100 * probs))))
  thr
}

effect_metric_rows <- function(
  effect,
  v1_true,
  v0_true,
  v1_est,
  v0_est,
  probs = c(0.1, 0.5, 0.9),
  thresholds = NULL,
  n_div = 5000L,
  seed = NULL
) {
  v1_true <- as.numeric(v1_true)
  v0_true <- as.numeric(v0_true)
  v1_est <- as.numeric(v1_est)
  v0_est <- as.numeric(v0_est)
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1L, is.finite(seed))
    set.seed(as.integer(seed))
  }

  out <- list()
  idx <- 1L

  out[[idx]] <- tibble::tibble(
    effect = effect,
    metric = "mean",
    truth = mean(v1_true) - mean(v0_true),
    est = mean(v1_est) - mean(v0_est)
  )
  idx <- idx + 1L

  for (p in probs) {
    lab <- metric_label_from_prob(p)
    out[[idx]] <- tibble::tibble(
      effect = effect,
      metric = lab,
      truth = as.numeric(stats::quantile(v1_true, probs = p, names = FALSE, type = 8) -
                           stats::quantile(v0_true, probs = p, names = FALSE, type = 8)),
      est = as.numeric(stats::quantile(v1_est, probs = p, names = FALSE, type = 8) -
                         stats::quantile(v0_est, probs = p, names = FALSE, type = 8))
    )
    idx <- idx + 1L
  }

  if (!is.null(thresholds) && length(thresholds) > 0L) {
    if (is.null(names(thresholds))) {
      names(thresholds) <- paste0("tail_", seq_along(thresholds))
    }
    for (nm in names(thresholds)) {
      thr <- as.numeric(thresholds[[nm]])
      out[[idx]] <- tibble::tibble(
        effect = effect,
        metric = nm,
        truth = mean(v1_true > thr) - mean(v0_true > thr),
        est = mean(v1_est > thr) - mean(v0_est > thr)
      )
      idx <- idx + 1L
    }
  }

  take_n <- function(x, n_take) {
  if (length(x) <= n_take) return(x)
  x[sample.int(length(x), size = n_take)]
}
  n_div <- as.integer(n_div)
  n_div <- max(1L, min(n_div, min(length(v1_true), length(v0_true), length(v1_est), length(v0_est))))
  x1t <- take_n(v1_true, n_div)
  x0t <- take_n(v0_true, n_div)
  x1e <- take_n(v1_est, n_div)
  x0e <- take_n(v0_est, n_div)

  out[[idx]] <- tibble::tibble(
    effect = effect,
    metric = "ED",
    truth = energy_distance(x1t, x0t),
    est = energy_distance(x1e, x0e)
  )
  idx <- idx + 1L

  out[[idx]] <- tibble::tibble(
    effect = effect,
    metric = "W1",
    truth = wasserstein1_distance_1d(x1t, x0t),
    est = wasserstein1_distance_1d(x1e, x0e)
  )
  idx <- idx + 1L

  out[[idx]] <- tibble::tibble(
    effect = effect,
    metric = "KL",
    truth = kl_divergence_kde_1d(x1t, x0t),
    est = kl_divergence_kde_1d(x1e, x0e)
  )

  dplyr::bind_rows(out)
}

compute_rich_interventional_metrics <- function(
  truth_outcomes,
  est_outcomes,
  probs = c(0.1, 0.5, 0.9),
  thresholds = NULL,
  n_div = 5000L,
  effect_prefix = "IPSE",
  seed = NULL
) {
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1L, is.finite(seed))
    base_seed <- as.integer(seed)
  } else {
    base_seed <- NULL
  }

  stopifnot(all(c("Y_1M1", "Y_0M0", "Y_1M0") %in% names(truth_outcomes)))
  stopifnot(all(c("Y_1M1", "Y_0M0", "Y_1M0") %in% names(est_outcomes)))

  pieces <- list(
    effect_metric_rows(
      effect = "ITE",
      v1_true = truth_outcomes$Y_1M1,
      v0_true = truth_outcomes$Y_0M0,
      v1_est = est_outcomes$Y_1M1,
      v0_est = est_outcomes$Y_0M0,
      probs = probs,
      thresholds = thresholds,
      n_div = n_div,
      seed = if (is.null(base_seed)) NULL else (base_seed + 11L)
    ),
    effect_metric_rows(
      effect = "IDE",
      v1_true = truth_outcomes$Y_1M0,
      v0_true = truth_outcomes$Y_0M0,
      v1_est = est_outcomes$Y_1M0,
      v0_est = est_outcomes$Y_0M0,
      probs = probs,
      thresholds = thresholds,
      n_div = n_div,
      seed = if (is.null(base_seed)) NULL else (base_seed + 23L)
    ),
    effect_metric_rows(
      effect = "IIE",
      v1_true = truth_outcomes$Y_1M1,
      v0_true = truth_outcomes$Y_1M0,
      v1_est = est_outcomes$Y_1M1,
      v0_est = est_outcomes$Y_1M0,
      probs = probs,
      thresholds = thresholds,
      n_div = n_div,
      seed = if (is.null(base_seed)) NULL else (base_seed + 37L)
    )
  )

  if ("path" %in% names(truth_outcomes) && "path" %in% names(est_outcomes)) {
    stopifnot(length(truth_outcomes$path) == length(est_outcomes$path))
    for (s in seq_along(truth_outcomes$path)) {
      path_seed <- if (is.null(base_seed)) NULL else (base_seed + 500L + 100L * s)
      pieces[[length(pieces) + 1L]] <- effect_metric_rows(
        effect = paste0(effect_prefix, s),
        v1_true = truth_outcomes$path[[s]]$Y_a_M1s,
        v0_true = truth_outcomes$path[[s]]$Y_a_M0s,
        v1_est = est_outcomes$path[[s]]$Y_a_M1s,
        v0_est = est_outcomes$path[[s]]$Y_a_M0s,
        probs = probs,
        thresholds = thresholds,
        n_div = n_div,
        seed = path_seed
      )
    }
  }

  dplyr::bind_rows(pieces)
}
