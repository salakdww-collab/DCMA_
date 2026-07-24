suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

root <- Sys.getenv("DCMA_CODE_ROOT", ".")
source(file.path(root, "experiment", "ihdp_semisynthetic", "ihdp_semisynth_3m.R"))

base_csv <- Sys.getenv("IHDP_BASE_CSV", file.path(root, "data", "ihdp", "ihdp_npci_1.csv"))
out_dir <- Sys.getenv("IHDP_SEMISYNTH_OUT_DIR", file.path(root, "data", "ihdp_semisynth_3m"))
R_rep <- as.integer(Sys.getenv("IHDP_R_REP", "100"))
seed_base <- as.integer(Sys.getenv("IHDP_SEED_BASE", "820000"))
base_n <- as.integer(Sys.getenv("IHDP_BASE_N", "0"))
scenario_env <- Sys.getenv("IHDP_SCENARIO_LIST", "IHDPMean,IHDPDist")
scenario_list <- trimws(strsplit(scenario_env, ",", fixed = TRUE)[[1]])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base_df <- load_ihdp_base(csv_path = base_csv)
if (!is.na(base_n) && base_n > 0L && base_n != nrow(base_df)) {
  set.seed(seed_base + 17L)
  idx_base <- sample(seq_len(nrow(base_df)), size = base_n, replace = TRUE)
  base_df <- base_df[idx_base, , drop = FALSE]
  base_df$row_id <- seq_len(nrow(base_df))
}
write_csv(base_df, file.path(out_dir, "ihdp_base_covariates.csv"))

manifest_rows <- list()
row_id <- 1L

for (scenario in scenario_list) {
  sc_dir <- file.path(out_dir, scenario)
  dir.create(sc_dir, recursive = TRUE, showWarnings = FALSE)
  for (r in seq_len(R_rep)) {
    seed_offset <- match(scenario, scenario_list) * 100000L
    seed <- seed_base + seed_offset + r
    dat <- generate_ihdp_semisynth_3m(
      base_df = base_df,
      scenario = scenario,
      seed = seed
    )
    out_csv <- file.path(sc_dir, sprintf("%s_rep_%03d.csv", tolower(scenario), r))
    write_csv(dat, out_csv)
    m_vars <- grep("^M", names(dat), value = TRUE)
    manifest_rows[[row_id]] <- tibble(
      scenario = scenario,
      S = length(m_vars),
      rep = r,
      n_est = nrow(dat),
      seed_data = seed,
      data_csv = out_csv
    )
    row_id <- row_id + 1L
  }
}

manifest <- bind_rows(manifest_rows) |>
  arrange(scenario, rep)
write_csv(manifest, file.path(out_dir, "manifest.csv"))

dict_df <- tribble(
  ~variable, ~type, ~description,
  "X", "binary", "Treatment copied from IHDP treatment column",
  "C1-C10", "mixed", "Selected IHDP covariates from x1-x10",
  "M1", "continuous", "Mediator 1, primarily mean-driving",
  "M2", "continuous", "Mediator 2, primarily distribution-shape-driving in the two-mediator display scenario",
  "M3", "continuous", "Mediator 3, used only in three-mediator scenarios",
  "S_latent", "binary", "Latent mixture indicator used only for generation",
  "Y", "continuous", "Semi-synthetic outcome"
)
write_csv(dict_df, file.path(out_dir, "dictionary.csv"))

summary_df <- bind_rows(
  lapply(scenario_list, function(sc) {
    dat <- read_csv(file.path(out_dir, sc, sprintf("%s_rep_%03d.csv", tolower(sc), 1L)), show_col_types = FALSE)
    tibble(
      scenario = sc,
      n = nrow(dat),
      rate_X1 = mean(dat$X),
      mean_Y = mean(dat$Y),
      sd_Y = sd(dat$Y),
      mean_M1 = mean(dat$M1),
      mean_M2 = mean(dat$M2),
      mean_M3 = if ("M3" %in% names(dat)) mean(dat$M3) else NA_real_,
      mix_rate = mean(dat$S_latent)
    )
  })
)
write_csv(summary_df, file.path(out_dir, "summary_preview_rep1.csv"))

cat("Wrote IHDP semi-synthetic data to:\n")
cat(out_dir, "\n")
