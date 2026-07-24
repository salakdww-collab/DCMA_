suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

root <- Sys.getenv("DCMA_CODE_ROOT", ".")
bench_dir <- file.path(root, "experiment", "synthetic")
results_dir <- file.path(bench_dir, "results")
out_path <- file.path(results_dir, "combined_summary.csv")

summary_files <- c(
  file.path(results_dir, "dcma_es", "summary.csv"),
  file.path(results_dir, "dcma_wgr", "summary.csv"),
  file.path(results_dir, "linear", "summary.csv")
)

existing <- summary_files[file.exists(summary_files)]
if (!length(existing)) stop("No benchmark summary files found.")

combined <- bind_rows(lapply(existing, function(p) read_csv(p, show_col_types = FALSE))) |>
  arrange(scenario, metric_family, effect, method)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(combined, out_path)
cat("Wrote combined benchmark summary:\n")
cat(out_path, "\n")
