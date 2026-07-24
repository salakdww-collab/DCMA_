## =========================================================
## NHANES liver: full-data bootstrap (row resampling + refit)
## Methods:
##   - DCMA
##
## Outputs:
##   - bootstrap raw estimates (per replicate)
##   - runtime per replicate/method
##   - summary table with percentile CI across bootstrap replicates
## =========================================================

suppressPackageStartupMessages({
  library(torch)
  library(dplyr)
  library(tibble)
})

source("utils/nn_model.R")
source("utils/energyloss_es.R")
source("utils/energy_distance.R")
source("model/dcma_fit.R")
source("model/dcma_reconstruct.R")
source(file.path("experiment", "nhanes_liver", "nhanes_liver_covariates.R"))

cfg <- list(
  in_csv = Sys.getenv(
    "NHANES_LIVER_IN_CSV",
    Sys.getenv("NHANES_LIVER_IN_CSV", file.path("data", "nhanes_2017_2018_liver", "nhanes_2017_2018_liver_analysis_fasting.csv"))
  ),
  out_dir = Sys.getenv("NHANES_LIVER_FULLBOOT_OUT", "results/nhanes_liver_threeway_fullbootstrap_completecase_threshold"),
  y_var = Sys.getenv("NHANES_LIVER_Y", "LUXSMED"),
  y_log1p = as.logical(Sys.getenv("NHANES_LIVER_Y_LOG1P", "FALSE")),
  threshold = as.numeric(Sys.getenv("NHANES_LIVER_THRESHOLD", "8")),
  seed_base = as.integer(Sys.getenv("NHANES_LIVER_FULLBOOT_SEED", "20260329")),
  n_boot = as.integer(Sys.getenv("NHANES_LIVER_FULLBOOT_N", "100")),
  conf_level = as.numeric(Sys.getenv("NHANES_LIVER_CONF", "0.95")),
  B_est = as.integer(Sys.getenv("NHANES_LIVER_B_EST", "1000")),
  ed_n = as.integer(Sys.getenv("NHANES_LIVER_ED_N", "5000")),
  device = Sys.getenv("NHANES_LIVER_DCMA_DEVICE", "mps")
)

if (!file.exists(cfg$in_csv)) stop("Input file not found: ", cfg$in_csv)
if (!cfg$y_var %in% c("LUXCAPM", "LUXSMED")) stop("NHANES_LIVER_Y must be LUXCAPM or LUXSMED.")
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

if (cfg$device == "mps" && !torch::backends_mps_is_available()) {
  message("MPS unavailable, fallback to CPU.")
  cfg$device <- "cpu"
}

raw_file <- file.path(cfg$out_dir, "nhanes_liver_threeway_fullbootstrap_raw.csv")
rt_file <- file.path(cfg$out_dir, "nhanes_liver_threeway_fullbootstrap_runtime.csv")
sum_file <- file.path(cfg$out_dir, "nhanes_liver_threeway_fullbootstrap_summary.csv")
calc_ed <- function(v1, v0, n_ed, seed) {
  v1 <- as.numeric(v1)
  v0 <- as.numeric(v0)
  n1 <- min(length(v1), n_ed)
  n0 <- min(length(v0), n_ed)
  set.seed(seed)
  i1 <- sample.int(length(v1), n1)
  i0 <- sample.int(length(v0), n0)
  energy_distance(v1[i1], v0[i0])
}

calc_exceedance_draws <- function(v1, v0, threshold) {
  stopifnot(length(dim(v1)) == 2L, length(dim(v0)) == 2L)
  colMeans(v1 >= threshold) - colMeans(v0 >= threshold)
}

## -----------------------------
## Base data
## -----------------------------
dat0 <- read.csv(cfg$in_csv, stringsAsFactors = FALSE)
dat_base <- nhanes_liver_make_analysis_dat(dat0, cfg$y_var)

if (cfg$y_log1p) dat_base$Y <- log1p(dat_base$Y)

n <- nrow(dat_base)
cov_vars <- nhanes_liver_covariate_vars()
effects <- c("ITE", "IDE", "IIE", "IPSE1", "IPSE2")

existing_raw <- if (file.exists(raw_file)) read.csv(raw_file, stringsAsFactors = FALSE) else tibble()
existing_rt <- if (file.exists(rt_file)) read.csv(rt_file, stringsAsFactors = FALSE) else tibble()
if ("method" %in% names(existing_raw)) {
  existing_raw <- existing_raw |> filter(method == "DCMA")
}
if ("method" %in% names(existing_rt)) {
  existing_rt <- existing_rt |> filter(method == "DCMA")
}
done_reps <- if (nrow(existing_raw) > 0L) sort(unique(existing_raw$boot_rep)) else integer(0)
todo <- setdiff(seq_len(cfg$n_boot), done_reps)

if (length(todo) == 0L) {
  message("All bootstrap replicates already completed. Rebuilding summary.")
} else {
  message(sprintf("Running full bootstrap reps: %d to %d (n=%d)", min(todo), max(todo), length(todo)))
}

for (b in todo) {
  seed_b <- cfg$seed_base + b * 1000L
  set.seed(seed_b)
  idx <- sample.int(n, n, replace = TRUE)
  db <- dat_base[idx, , drop = FALSE]

  rep_rows <- list()
  rep_rt <- list()

  ## -------------------------
  ## DCMA (mean + ED)
  ## -------------------------
  ok_dcma <- TRUE
  t_dcma_fit <- NA_real_
  t_dcma_rec <- NA_real_
  try({
    ddc <- nhanes_liver_make_dcma_dat(db)
    cov_names <- setdiff(names(ddc), c("X", "Y", "M1", "M2"))

    set.seed(seed_b + 101L)
    torch::torch_manual_seed(seed_b + 101L)
    t1 <- system.time({
      fit_dcma <- dcma(
        data = ddc,
        c_vars = cov_names,
        standardize = TRUE,
        silent = TRUE,
        device = cfg$device,
        y_family = "continuous"
      )
    })
    t_dcma_fit <- unname(t1[["elapsed"]])

    set.seed(seed_b + 151L)
    torch::torch_manual_seed(seed_b + 151L)
    t2 <- system.time({
      out_dcma <- dcma_reconstruct_interventions(
        fit = fit_dcma,
        B = cfg$B_est,
        C_mat = ddc[, cov_names, drop = FALSE],
        device = cfg$device
      )
    })
    t_dcma_rec <- unname(t2[["elapsed"]])

    draw_dc <- list(
      ITE = colMeans(out_dcma$outcomes$Y_1M1 - out_dcma$outcomes$Y_0M0),
      IDE = colMeans(out_dcma$outcomes$Y_1M0 - out_dcma$outcomes$Y_0M0),
      IIE = colMeans(out_dcma$outcomes$Y_1M1 - out_dcma$outcomes$Y_1M0),
      IPSE1 = colMeans(out_dcma$outcomes$path[[1]]$Y_a_M1s - out_dcma$outcomes$path[[1]]$Y_a_M0s),
      IPSE2 = colMeans(out_dcma$outcomes$path[[2]]$Y_a_M1s - out_dcma$outcomes$path[[2]]$Y_a_M0s)
    )

    pair_dc <- list(
      ITE = list(v1 = out_dcma$outcomes$Y_1M1, v0 = out_dcma$outcomes$Y_0M0),
      IDE = list(v1 = out_dcma$outcomes$Y_1M0, v0 = out_dcma$outcomes$Y_0M0),
      IIE = list(v1 = out_dcma$outcomes$Y_1M1, v0 = out_dcma$outcomes$Y_1M0),
      IPSE1 = list(v1 = out_dcma$outcomes$path[[1]]$Y_a_M1s, v0 = out_dcma$outcomes$path[[1]]$Y_a_M0s),
      IPSE2 = list(v1 = out_dcma$outcomes$path[[2]]$Y_a_M1s, v0 = out_dcma$outcomes$path[[2]]$Y_a_M0s)
    )

    rep_rows[[length(rep_rows) + 1L]] <- bind_rows(
      bind_rows(lapply(effects, function(eff) {
        tibble(
          boot_rep = b,
          method = "DCMA",
          metric_family = "mean",
          effect = eff,
          estimate = mean(draw_dc[[eff]]),
          n_draw = length(draw_dc[[eff]])
        )
      })),
      bind_rows(lapply(effects, function(eff) {
        tibble(
          boot_rep = b,
          method = "DCMA",
          metric_family = "exceedance",
          effect = eff,
          estimate = mean(calc_exceedance_draws(pair_dc[[eff]]$v1, pair_dc[[eff]]$v0, cfg$threshold)),
          n_draw = ncol(pair_dc[[eff]]$v1)
        )
      })),
      bind_rows(lapply(effects, function(eff) {
        tibble(
          boot_rep = b,
          method = "DCMA",
          metric_family = "ED",
          effect = eff,
          estimate = calc_ed(pair_dc[[eff]]$v1, pair_dc[[eff]]$v0, n_ed = cfg$ed_n, seed = seed_b + 211L + match(eff, effects)),
          n_draw = min(cfg$ed_n, length(as.numeric(pair_dc[[eff]]$v1)), length(as.numeric(pair_dc[[eff]]$v0)))
        )
      }))
    )
  }, silent = TRUE)
  if (is.na(t_dcma_fit) || is.na(t_dcma_rec)) ok_dcma <- FALSE

  rep_rt[[length(rep_rt) + 1L]] <- tibble(
    boot_rep = b,
    method = "DCMA",
    ok = ok_dcma,
    fit_sec = ifelse(ok_dcma, t_dcma_fit, NA_real_),
    recon_sec = ifelse(ok_dcma, t_dcma_rec, NA_real_),
    total_sec = ifelse(ok_dcma, t_dcma_fit + t_dcma_rec, NA_real_)
  )

  add_raw <- if (length(rep_rows) > 0L) bind_rows(rep_rows) else tibble()
  add_rt <- bind_rows(rep_rt)

  existing_raw <- bind_rows(existing_raw, add_raw)
  existing_rt <- bind_rows(existing_rt, add_rt)
  write.csv(existing_raw, raw_file, row.names = FALSE)
  write.csv(existing_rt, rt_file, row.names = FALSE)

  message(sprintf(
    "[boot %03d/%03d] done | ok: DCMA=%s",
    b, cfg$n_boot, ok_dcma
  ))
}

## -----------------------------
## Summary with percentile CI
## -----------------------------
summary_tbl <- existing_raw |>
  group_by(method, metric_family, effect) |>
  summarise(
    n_boot_eff = n(),
    mean_est = mean(estimate, na.rm = TRUE),
    sd_est = sd(estimate, na.rm = TRUE),
    lcl = as.numeric(quantile(estimate, probs = 0.025, na.rm = TRUE)),
    ucl = as.numeric(quantile(estimate, probs = 0.975, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  mutate(
    metric_family = factor(metric_family, levels = c("mean", "exceedance", "ED")),
    effect = factor(effect, levels = effects),
    method = factor(method, levels = c("DCMA"))
  ) |>
  arrange(metric_family, effect, method) |>
  mutate(
    metric_family = as.character(metric_family),
    effect = as.character(effect),
    method = as.character(method),
    mean_est = ifelse(metric_family == "ED", pmax(mean_est, 0), mean_est),
    lcl = ifelse(metric_family == "ED", pmax(lcl, 0), lcl),
    ucl = ifelse(metric_family == "ED", pmax(ucl, 0), ucl),
    interval_type = "full_refit_bootstrap_CI"
  )

write.csv(summary_tbl, sum_file, row.names = FALSE)

rt_summary <- existing_rt |>
  group_by(method) |>
  summarise(
    n_rep = n(),
    n_ok = sum(ok, na.rm = TRUE),
    mean_fit_sec = mean(fit_sec, na.rm = TRUE),
    mean_recon_sec = mean(recon_sec, na.rm = TRUE),
    mean_total_sec = mean(total_sec, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(rt_summary, file.path(cfg$out_dir, "nhanes_liver_threeway_fullbootstrap_runtime_summary.csv"), row.names = FALSE)

cat("Wrote:\n")
cat(" -", raw_file, "\n")
cat(" -", rt_file, "\n")
cat(" -", sum_file, "\n")
cat(" -", file.path(cfg$out_dir, "nhanes_liver_threeway_fullbootstrap_runtime_summary.csv"), "\n")
cat("\nSummary:\n")
print(summary_tbl)
cat("\nRuntime summary:\n")
print(rt_summary)
