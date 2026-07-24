suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

root <- Sys.getenv("DCMA_CODE_ROOT", ".")
bench_dir <- file.path(root, "experiment", "synthetic")
data_dir <- file.path(bench_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(root, "experiment", "synthetic", "simulate_mediation_data.R"))
source(file.path(root, "experiment", "synthetic", "dcma_truth.R"))

parse_chr <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])
scenarios <- parse_chr(Sys.getenv("BENCH_SCENARIOS", "S1,S2"))
n_est <- as.integer(Sys.getenv("BENCH_N_EST", "5000"))
R_rep <- as.integer(Sys.getenv("BENCH_R_REP", "100"))
N_truth <- as.integer(Sys.getenv("BENCH_N_TRUTH", "30000"))
B_truth <- as.integer(Sys.getenv("BENCH_B_TRUTH", "100"))
seed_truth <- as.integer(Sys.getenv("BENCH_SEED_TRUTH", "123"))
seed_data_base <- as.integer(Sys.getenv("BENCH_SEED_DATA_BASE", "500000"))

manifest_rows <- vector("list", length(scenarios) * R_rep)
row_id <- 1L

for (sc in scenarios) {
  S <- if (identical(sc, "S1")) 1L else 5L

  truth_obj <- dcma_truth_interventions(
    scenario = sc,
    N = N_truth,
    B = B_truth,
    S = S,
    seed = seed_truth
  )
  sc_dir <- file.path(data_dir, sc)
  dir.create(sc_dir, recursive = TRUE, showWarnings = FALSE)

  for (r in seq_len(R_rep)) {
    seed_data <- seed_data_base + if (identical(sc, "S1")) r else 100000L + r
    dat <- simulate_mediation_data(
      n = n_est,
      scenario = sc,
      S = S,
      seed = seed_data
    )
    out_csv <- file.path(sc_dir, sprintf("%s_rep_%03d.csv", tolower(sc), r))
    write_csv(dat, out_csv)

    manifest_rows[[row_id]] <- tibble(
      scenario = sc,
      S = S,
      rep = r,
      n_est = n_est,
      seed_data = seed_data,
      data_csv = out_csv
    )
    row_id <- row_id + 1L
  }
}

manifest <- bind_rows(manifest_rows)
write_csv(manifest, file.path(data_dir, "manifest.csv"))
cat("Wrote shared data manifest:\n")
cat(file.path(data_dir, "manifest.csv"), "\n")
