# Compute ITE, IDE, and IPSE_s for one scenario (truth vs. estimated interventions)
#
# Compares truth and estimated interventional outcome samples for a given scenario and
# computes summary functionals for:
#   ITE:  Y_1M1 vs Y_0M0
#   IDE:  Y_1M0 vs Y_0M0
#   IPSEs: Y_a_M1s vs Y_a_M0s for each mediator pathway s
# across all replications.
#
# Args:
#   truth_outcomes_list: list containing truth interventional outcomes by scenario
#   est_outcomes_rep_list: list of replication outputs, each containing estimated outcomes by scenario
#   scenario_label: scenario key used to index the above lists
#   probs: quantile levels for quantile effects
#   tail_threshold: if not NULL, compute tail probability effects
#   n_ed: subsample size used for ED computation
#
# Returns:
#   A tibble with one row per (replication × effect), including scenario and replication id.

compute_effects <- function(
  truth_outcomes_list,
  est_outcomes_rep_list,
  scenario_label,
  probs          = c(0.1, 0.5, 0.9),
  tail_threshold = NULL,
  n_ed           = 2000L
) {
  # Compute distributional effect functionals for two distributions
  #
  # Given samples from two distributions under "1" and "0" (truth and estimate),
  # computes mean effect, quantile effects, optional tail probability effect,
  # and energy distance (ED) between the two distributions (computed on subsamples).
  #
  # Args:
  #   v1_true, v0_true: numeric vectors of true samples
  #   v1_est,  v0_est:  numeric vectors of estimated samples
  #   probs: quantile levels
  #   tail_threshold: if not NULL, compute P(Y > threshold) difference
  #   n_ed: subsample size for ED computation
  #
  # Returns:
  #   A tibble with mean effects, quantile effects, optional tail effects, and ED.
  compute_effect_functionals <- function(
    v1_true, v0_true,
    v1_est,  v0_est,
    probs  = c(0.1, 0.5, 0.9),
    tail_threshold = NULL,
    n_ed = 2000L
  ) {
    mean_true <- mean(v1_true) - mean(v0_true)
    mean_est  <- mean(v1_est)  - mean(v0_est)

    q_true <- stats::quantile(v1_true, probs) - stats::quantile(v0_true, probs)
    q_est  <- stats::quantile(v1_est,  probs) - stats::quantile(v0_est,  probs)

    if (!is.null(tail_threshold)) {
      tail_true <- mean(v1_true > tail_threshold) - mean(v0_true > tail_threshold)
      tail_est  <- mean(v1_est  > tail_threshold) - mean(v0_est  > tail_threshold)
    } else {
      tail_true <- NA_real_
      tail_est  <- NA_real_
    }

    idx1t <- sample.int(length(v1_true), min(n_ed, length(v1_true)))
    idx0t <- sample.int(length(v0_true), min(n_ed, length(v0_true)))
    idx1e <- sample.int(length(v1_est),  min(n_ed, length(v1_est)))
    idx0e <- sample.int(length(v0_est),  min(n_ed, length(v0_est)))

    ed_true <- energy_distance(v1_true[idx1t], v0_true[idx0t])
    ed_est  <- energy_distance(v1_est[idx1e],  v0_est[idx0e])

    tibble::tibble(
      mean_true = mean_true,
      mean_est  = mean_est,
      q_true    = list(q_true),
      q_est     = list(q_est),
      tail_true = tail_true,
      tail_est  = tail_est,
      ED_true   = ed_true,
      ED_est    = ed_est
    )
  }

  R_rep <- length(est_outcomes_rep_list)

  truth_obj <- truth_outcomes_list[[scenario_label]]$outcomes
  S <- length(truth_obj$path)

  truth_basic <- list(
    Y_1M1 = truth_obj$Y_1M1,
    Y_1M0 = truth_obj$Y_1M0,
    Y_0M0 = truth_obj$Y_0M0
  )

  truth_paths <- lapply(seq_len(S), function(s) {
    list(
      M0s = truth_obj$path[[s]]$Y_a_M0s,
      M1s = truth_obj$path[[s]]$Y_a_M1s
    )
  })

  purrr::map_dfr(seq_len(R_rep), function(r) {
    est_obj <- est_outcomes_rep_list[[r]][[scenario_label]]$outcomes

    est_basic <- list(
      Y_1M1 = est_obj$Y_1M1,
      Y_1M0 = est_obj$Y_1M0,
      Y_0M0 = est_obj$Y_0M0
    )

    est_paths <- lapply(seq_len(S), function(s) {
      list(
        M0s = est_obj$path[[s]]$Y_a_M0s,
        M1s = est_obj$path[[s]]$Y_a_M1s
      )
    })

    out_list_rep <- list()

    ite_tbl <- compute_effect_functionals(
      v1_true = truth_basic$Y_1M1,
      v0_true = truth_basic$Y_0M0,
      v1_est  = est_basic$Y_1M1,
      v0_est  = est_basic$Y_0M0,
      probs          = probs,
      tail_threshold = tail_threshold,
      n_ed           = n_ed
    ) |>
      dplyr::mutate(effect = "ITE")
    out_list_rep[[length(out_list_rep) + 1L]] <- ite_tbl

    ide_tbl <- compute_effect_functionals(
      v1_true = truth_basic$Y_1M0,
      v0_true = truth_basic$Y_0M0,
      v1_est  = est_basic$Y_1M0,
      v0_est  = est_basic$Y_0M0,
      probs          = probs,
      tail_threshold = tail_threshold,
      n_ed           = n_ed
    ) |>
      dplyr::mutate(effect = "IDE")
    out_list_rep[[length(out_list_rep) + 1L]] <- ide_tbl

    for (s in seq_len(S)) {
      ipse_tbl <- compute_effect_functionals(
        v1_true = truth_paths[[s]]$M1s,
        v0_true = truth_paths[[s]]$M0s,
        v1_est  = est_paths[[s]]$M1s,
        v0_est  = est_paths[[s]]$M0s,
        probs          = probs,
        tail_threshold = tail_threshold,
        n_ed           = n_ed
      ) |>
        dplyr::mutate(effect = paste0("IPSE", s))

      out_list_rep[[length(out_list_rep) + 1L]] <- ipse_tbl
    }

    dplyr::bind_rows(out_list_rep) |>
      dplyr::mutate(
        scenario = scenario_label,
        rep      = r,
        .before  = 1
      )
  })
}
