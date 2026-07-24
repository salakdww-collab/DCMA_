suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(torch)
})

root <- Sys.getenv("DCMA_CODE_ROOT", ".")
bench_dir <- file.path(root, "experiment", "synthetic")

source(file.path(root, "experiment", "synthetic", "simulate_mediation_data.R"))
source(file.path(root, "experiment", "synthetic", "dcma_truth.R"))
source(file.path(root, "utils", "nn_model.R"))
source(file.path(root, "utils", "energyloss_es.R"))
source(file.path(root, "utils", "energy_distance.R"))
source(file.path(root, "utils", "interventional_chain_regression.R"))
source(file.path(root, "model", "dcma_fit.R"))
source(file.path(root, "model", "dcma_reconstruct.R"))

summarise_raw <- function(raw) {
  summary_mean <- raw |>
    group_by(ablation, variant, scenario, effect, metric_family = "mean") |>
    summarise(
      n_data_rep = n_distinct(rep),
      n_row = n(),
      truth = mean(mean_true),
      bias = mean(mean_est - mean_true),
      rmse = sqrt(mean((mean_est - mean_true)^2)),
      sd_est = sd(mean_est),
      var_est = var(mean_est),
      .groups = "drop"
    )

  summary_ed <- raw |>
    group_by(ablation, variant, scenario, effect, metric_family = "ED") |>
    summarise(
      n_data_rep = n_distinct(rep),
      n_row = n(),
      truth = mean(ED_true),
      bias = mean(ED_est - ED_true),
      rmse = sqrt(mean((ED_est - ED_true)^2)),
      sd_est = sd(ED_est),
      var_est = var(ED_est),
      .groups = "drop"
    )

  bind_rows(summary_mean, summary_ed) |>
    arrange(ablation, variant, metric_family, effect)
}

ablation <- "noise"
R_rep <- as.integer(Sys.getenv("S1_NOISE_R_REP", "100"))
n_est <- as.integer(Sys.getenv("S1_NOISE_N_EST", "5000"))
N_truth <- as.integer(Sys.getenv("S1_NOISE_N_TRUTH", "30000"))
B_truth <- as.integer(Sys.getenv("S1_NOISE_B_TRUTH", "100"))
B_est <- as.integer(Sys.getenv("S1_NOISE_B_EST", "200"))
n_ed <- as.integer(Sys.getenv("S1_NOISE_N_ED", "10000"))
seed_truth <- as.integer(Sys.getenv("S1_NOISE_SEED_TRUTH", "123"))
seed_data_base <- as.integer(Sys.getenv("S1_NOISE_SEED_DATA_BASE", "500000"))
s1_mode_shift <- as.numeric(Sys.getenv("S1_NOISE_MODE_SHIFT", "2"))

epsm_dim <- as.integer(Sys.getenv("S1_NOISE_EPSM_DIM", "4"))
epsy_dim <- as.integer(Sys.getenv("S1_NOISE_EPSY_DIM", "8"))
fit_device <- Sys.getenv("DCMA_DEVICE", "cpu")
recon_device <- Sys.getenv("DCMA_RECON_DEVICE", fit_device)

variants <- tibble(
  variant = c("Internal noise", "Additive homosk.", "Additive heterosk."),
  y_generator = c("plain", "additive_homosk", "additive_heterosk")
)

out_suffix <- Sys.getenv("S1_NOISE_OUT_SUFFIX", "noise_s1_shift2_100rep")
out_dir <- file.path(bench_dir, "results", paste0("s1_ablation_", out_suffix))
data_dir <- file.path(out_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

raw_path <- file.path(out_dir, "raw.csv")
runtime_path <- file.path(out_dir, "runtime.csv")
summary_path <- file.path(out_dir, "summary.csv")
raw <- if (file.exists(raw_path)) read_csv(raw_path, show_col_types = FALSE) else tibble()
runtime <- if (file.exists(runtime_path)) read_csv(runtime_path, show_col_types = FALSE) else tibble()
has_raw_rows <- function(...) {
  if (!nrow(raw)) return(FALSE)
  nrow(filter(raw, ...)) > 0L
}

truth <- dcma_truth_interventions(
  scenario = "S1",
  N = N_truth,
  B = B_truth,
  s1_mode_shift = s1_mode_shift,
  seed = seed_truth
)

cat("S1 outcome-noise ablation runner\n")
cat("R_rep =", R_rep, "| n =", n_est, "| s1_mode_shift =", s1_mode_shift, "\n")
cat("out_dir =", out_dir, "\n")
print(variants)

for (r in seq_len(R_rep)) {
  seed_data <- seed_data_base + r
  data_path <- file.path(data_dir, sprintf("s1_rep_%03d.csv", r))
  if (file.exists(data_path)) {
    dat <- read_csv(data_path, show_col_types = FALSE)
  } else {
    dat <- simulate_mediation_data(
      n = n_est,
      scenario = "S1",
      S = 1L,
      s1_mode_shift = s1_mode_shift,
      seed = seed_data
    )
    write_csv(dat, data_path)
  }

  for (i in seq_len(nrow(variants))) {
    v <- variants[i, ]
    if (has_raw_rows(rep == r, variant == v$variant)) next
    seed_split <- seed_data + 1000L * i + 101L
    seed_init <- seed_data + 1000L * i + 102L
    seed_recon <- seed_data + 1000L * i + 103L
    seed_ed <- seed_data + 1000L * i + 104L

    cat(sprintf("[run] S1 noise rep=%03d variant=%s\n", r, v$variant))
    fit_t <- system.time({
      fit <- dcma(
        data = as.data.frame(dat),
        epsm_dim = epsm_dim,
        epsy_dim = epsy_dim,
        c_vars = "Z",
        standardize = TRUE,
        silent = TRUE,
        m_generator = "plain",
        y_generator = v$y_generator,
        split_seed = seed_split,
        init_seed = seed_init,
        device = fit_device
      )
    })
    recon_t <- system.time({
      est <- dcma_reconstruct_interventions(
        fit = fit,
        B = B_est,
        seed = seed_recon,
        device = recon_device
      )
    })

    one_rep <- compute_interventional_basic_effects(
      truth_outcomes = truth$outcomes,
      est_outcomes = est$outcomes,
      n_ed = n_ed,
      ed_seed = seed_ed
    ) |>
      mutate(
        ablation = ablation,
        variant = v$variant,
        scenario = "S1",
        rep = r,
        recon_rep = 1L,
        n_est = n_est,
        m_generator = "plain",
        y_generator = v$y_generator,
        fit_sec = unname(fit_t[["elapsed"]]),
        recon_sec = unname(recon_t[["elapsed"]]),
        .before = 1
      )

    raw <- bind_rows(raw, one_rep) |> arrange(rep, variant, effect)
    runtime <- bind_rows(
      runtime,
      tibble(
        ablation = ablation,
        variant = v$variant,
        rep = r,
        recon_rep = 1L,
        fit_sec = unname(fit_t[["elapsed"]]),
        recon_sec = unname(recon_t[["elapsed"]])
      )
    )
    write_csv(raw, raw_path)
    write_csv(runtime, runtime_path)
    write_csv(summarise_raw(raw), summary_path)
    cat(sprintf("[done] S1 noise rep=%03d variant=%s\n", r, v$variant))
  }
}

write_csv(summarise_raw(raw), summary_path)
cat("Wrote S1 outcome-noise ablation results to:\n")
cat(out_dir, "\n")
