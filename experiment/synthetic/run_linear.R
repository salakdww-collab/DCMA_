suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

root <- Sys.getenv("DCMA_CODE_ROOT", ".")
bench_dir <- file.path(root, "experiment", "synthetic")
manifest_path <- Sys.getenv("BENCH_MANIFEST_PATH", file.path(bench_dir, "data", "manifest.csv"))
out_dir <- Sys.getenv("BENCH_OUT_DIR", file.path(bench_dir, "results", "linear"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(root, "utils", "energy_distance.R"))
source(file.path(root, "utils", "interventional_chain_regression.R"))
source(file.path(root, "experiment", "synthetic", "dcma_truth.R"))

parse_chr <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])

B_est <- as.integer(Sys.getenv("BENCH_B_EST", "200"))
n_ed <- as.integer(Sys.getenv("BENCH_N_ED", "10000"))
N_truth <- as.integer(Sys.getenv("BENCH_N_TRUTH", "30000"))
B_truth <- as.integer(Sys.getenv("BENCH_B_TRUTH", "100"))
seed_truth <- as.integer(Sys.getenv("BENCH_SEED_TRUTH", "123"))

raw_path <- file.path(out_dir, "raw.csv")
sum_path <- file.path(out_dir, "summary.csv")
runtime_path <- file.path(out_dir, "runtime.csv")

manifest <- read_csv(manifest_path, show_col_types = FALSE)
scenario_filter <- Sys.getenv("BENCH_SCENARIOS", "")
if (nzchar(scenario_filter)) {
  manifest <- manifest |> filter(scenario %in% parse_chr(scenario_filter))
}
max_reps <- as.integer(Sys.getenv("BENCH_MAX_REPS", "0"))
if (is.finite(max_reps) && max_reps > 0L) {
  manifest <- manifest |> group_by(scenario) |> slice_head(n = max_reps) |> ungroup()
}
raw <- if (file.exists(raw_path)) read_csv(raw_path, show_col_types = FALSE) else tibble()
runtime <- if (file.exists(runtime_path)) read_csv(runtime_path, show_col_types = FALSE) else tibble()
done_keys <- if (nrow(runtime)) paste(runtime$scenario, runtime$rep, sep = "_") else character()

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  key <- paste(row$scenario, row$rep, sep = "_")
  if (key %in% done_keys) next

  dat <- read_csv(row$data_csv, show_col_types = FALSE)
  truth <- dcma_truth_interventions(
    scenario = row$scenario,
    N = N_truth,
    B = B_truth,
    S = row$S,
    seed = seed_truth
  )
  c_vars <- if ("Z" %in% names(dat)) "Z" else paste0("C", 1:5)
  m_vars <- grep("^M", names(dat), value = TRUE)

  fit_t <- system.time({
    fit <- fit_interventional_chain(
      data = as.data.frame(dat),
      c_vars = c_vars,
      x_var = "X",
      y_var = "Y",
      m_vars = m_vars
    )
  })

  recon_t <- system.time({
    est <- reconstruct_interventions_chain(
      fit_obj = fit,
      data = as.data.frame(dat),
      B = B_est,
      seed = row$seed_data + 900000L
    )
  })

  one_rep <- compute_interventional_basic_effects(
    truth_outcomes = truth$outcomes,
    est_outcomes = est$outcomes,
    n_ed = n_ed
  )

  if ("path" %in% names(truth$outcomes) && "path" %in% names(est$outcomes)) {
    one_rep <- bind_rows(
      one_rep,
      compute_path_specific_effects(
        truth_outcomes = truth$outcomes,
        est_outcomes = est$outcomes,
        n_ed = n_ed
      )
    )
  }

  one_rep <- one_rep |>
    mutate(
      method = "Linear",
      scenario = row$scenario,
      rep = row$rep,
      n_est = row$n_est,
      fit_sec = unname(fit_t[["elapsed"]]),
      recon_sec = unname(recon_t[["elapsed"]]),
      total_sec = unname(fit_t[["elapsed"]]) + unname(recon_t[["elapsed"]]),
      .before = 1
    ) |>
    select(method, scenario, rep, n_est, effect, mean_true, mean_est, ED_true, ED_est, fit_sec, recon_sec, total_sec)

  raw <- bind_rows(raw, one_rep) |> arrange(scenario, rep, effect)
  runtime <- bind_rows(runtime, one_rep |> distinct(method, scenario, rep, n_est, fit_sec, recon_sec, total_sec)) |>
    arrange(scenario, rep)

  write_csv(raw, raw_path)
  write_csv(runtime, runtime_path)
}

summary_mean <- raw |>
  group_by(method, scenario, effect) |>
  summarise(
    metric_family = "mean",
    n_rep = n(),
    truth = mean(mean_true),
    bias = mean(mean_est - mean_true),
    rmse = sqrt(mean((mean_est - mean_true)^2)),
    sd = sd(mean_est),
    .groups = "drop"
  )

summary_ed <- raw |>
  group_by(method, scenario, effect) |>
  summarise(
    metric_family = "ED",
    n_rep = n(),
    truth = mean(ED_true),
    bias = mean(ED_est - ED_true),
    rmse = sqrt(mean((ED_est - ED_true)^2)),
    sd = sd(ED_est),
    .groups = "drop"
  )

write_csv(bind_rows(summary_mean, summary_ed), sum_path)
cat("Wrote:\n")
cat(raw_path, "\n")
cat(sum_path, "\n")
