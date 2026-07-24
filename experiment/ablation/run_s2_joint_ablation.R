suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
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

parse_bool <- function(x, default = FALSE) {
  x <- tolower(trimws(x))
  if (x == "") return(default)
  if (x %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (x %in% c("0", "false", "f", "no", "n")) return(FALSE)
  warning(sprintf("Unrecognized boolean env value '%s'; fallback to default=%s", x, default))
  default
}

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
      .groups = "drop"
    )

  bind_rows(summary_mean, summary_ed) |>
    arrange(ablation, variant, metric_family, effect)
}

write_product_comparison <- function(raw, path) {
  if (!all(c("permute", "redraw") %in% unique(raw$variant))) return(invisible(NULL))
  cmp <- raw |>
    filter(effect %in% c("IDE", "IIE", "ITE", paste0("IPSE", 1:5))) |>
    select(rep, recon_rep, effect, variant, mean_est, ED_est) |>
    pivot_wider(names_from = variant, values_from = c(mean_est, ED_est)) |>
    mutate(
      mean_est_diff = mean_est_permute - mean_est_redraw,
      ED_est_diff = ED_est_permute - ED_est_redraw
    ) |>
    group_by(effect) |>
    summarise(
      n = n(),
      mean_abs_mean_diff = mean(abs(mean_est_diff), na.rm = TRUE),
      rmse_mean_diff = sqrt(mean(mean_est_diff^2, na.rm = TRUE)),
      mean_abs_ED_diff = mean(abs(ED_est_diff), na.rm = TRUE),
      rmse_ED_diff = sqrt(mean(ED_est_diff^2, na.rm = TRUE)),
      sd_ED_diff = sd(ED_est_diff, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(match(effect, c("IDE", "IIE", "ITE", paste0("IPSE", 1:5))))
  write_csv(cmp, path)
}

ablation <- Sys.getenv("S2_ABLATION_KIND", "joint_separate")
valid_ablation <- c("joint_separate")
if (!ablation %in% valid_ablation) {
  stop("S2_ABLATION_KIND must be one of: ", paste(valid_ablation, collapse = ", "))
}

R_rep <- as.integer(Sys.getenv("S2_ABLAT_R_REP", "100"))
n_est <- as.integer(Sys.getenv("S2_ABLAT_N_EST", "5000"))
N_truth <- as.integer(Sys.getenv("S2_ABLAT_N_TRUTH", "30000"))
B_truth <- as.integer(Sys.getenv("S2_ABLAT_B_TRUTH", "100"))
B_est <- as.integer(Sys.getenv("S2_ABLAT_B_EST", "200"))
n_ed <- as.integer(Sys.getenv("S2_ABLAT_N_ED", "10000"))
seed_truth <- as.integer(Sys.getenv("S2_ABLAT_SEED_TRUTH", "123"))
seed_data_base <- as.integer(Sys.getenv("S2_ABLAT_SEED_DATA_BASE", "980000"))
product_recon_rep <- as.integer(Sys.getenv("S2_ABLAT_PRODUCT_RECON_REP", "3"))

epsm_dim <- as.integer(Sys.getenv("S2_ABLAT_EPSM_DIM", "4"))
epsy_dim <- as.integer(Sys.getenv("S2_ABLAT_EPSY_DIM", "8"))
fit_device <- Sys.getenv("DCMA_DEVICE", "cpu")
recon_device <- Sys.getenv("DCMA_RECON_DEVICE", fit_device)

s2_noise <- Sys.getenv("S2_ABLAT_S2_NOISE", "heavyMY")
s2_df_y <- as.numeric(Sys.getenv("S2_ABLAT_S2_DF_Y", "3"))
s2_df_m <- as.numeric(Sys.getenv("S2_ABLAT_S2_DF_M", "5"))
s2_scale_y <- as.numeric(Sys.getenv("S2_ABLAT_S2_SCALE_Y", "1"))
s2_scale_m <- as.numeric(Sys.getenv("S2_ABLAT_S2_SCALE_M", "1"))
s2_mediator_scale <- as.numeric(Sys.getenv("S2_ABLAT_S2_MEDIATOR_SCALE", "0.5"))
s2_outcome_scale <- as.numeric(Sys.getenv("S2_ABLAT_S2_OUTCOME_SCALE", "0.5"))
s2_standardize_y <- parse_bool(Sys.getenv("S2_ABLAT_S2_STANDARDIZE_Y", "true"), default = TRUE)
s2_standardize_m <- parse_bool(Sys.getenv("S2_ABLAT_S2_STANDARDIZE_M", "true"), default = TRUE)

variants <- switch(
  ablation,
  joint_separate = tibble(
    variant = c("Joint-M", "Separate-M"),
    m_generator = c("plain", "separate"),
    y_generator = c("plain", "plain"),
    product_mode = c("permute", "permute")
  ),
  product = tibble(
    variant = c("permute", "redraw"),
    m_generator = c("plain", "plain"),
    y_generator = c("plain", "plain"),
    product_mode = c("permute", "redraw")
  ),
  noise = tibble(
    variant = c("Internal noise", "Additive homosk.", "Additive heterosk."),
    m_generator = c("plain", "plain", "plain"),
    y_generator = c("plain", "additive_homosk", "additive_heterosk"),
    product_mode = c("permute", "permute", "permute")
  )
)

out_suffix <- Sys.getenv("S2_ABLAT_OUT_SUFFIX", paste0(ablation, "_s2_heavyMY_100rep"))
out_dir <- file.path(bench_dir, "results", paste0("s2_ablation_", out_suffix))
data_dir <- file.path(out_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

raw_path <- file.path(out_dir, "raw.csv")
runtime_path <- file.path(out_dir, "runtime.csv")
summary_path <- file.path(out_dir, "summary.csv")
product_cmp_path <- file.path(out_dir, "product_mode_difference.csv")
raw <- if (file.exists(raw_path)) read_csv(raw_path, show_col_types = FALSE) else tibble()
runtime <- if (file.exists(runtime_path)) read_csv(runtime_path, show_col_types = FALSE) else tibble()
has_raw_rows <- function(...) {
  if (!nrow(raw)) return(FALSE)
  nrow(filter(raw, ...)) > 0L
}

truth <- dcma_truth_interventions(
  scenario = "S2",
  N = N_truth,
  B = B_truth,
  S = 5L,
  seed = seed_truth,
  s2_noise = s2_noise,
  s2_df_y = s2_df_y,
  s2_df_m = s2_df_m,
  s2_scale_y = s2_scale_y,
  s2_scale_m = s2_scale_m,
  s2_mediator_scale = s2_mediator_scale,
  s2_outcome_scale = s2_outcome_scale,
  s2_standardize_y = s2_standardize_y,
  s2_standardize_m = s2_standardize_m
)

cat("S2 ablation runner\n")
cat("ablation =", ablation, "| R_rep =", R_rep, "| n =", n_est, "| noise =", s2_noise, "\n")
cat("out_dir =", out_dir, "\n")
print(variants)

for (r in seq_len(R_rep)) {
  seed_data <- seed_data_base + r
  data_path <- file.path(data_dir, sprintf("s2_rep_%03d.csv", r))
  if (file.exists(data_path)) {
    dat <- read_csv(data_path, show_col_types = FALSE)
  } else {
    dat <- simulate_mediation_data(
      n = n_est,
      scenario = "S2",
      S = 5L,
      seed = seed_data,
      s2_noise = s2_noise,
      s2_df_y = s2_df_y,
      s2_df_m = s2_df_m,
      s2_scale_y = s2_scale_y,
      s2_scale_m = s2_scale_m,
      s2_mediator_scale = s2_mediator_scale,
      s2_outcome_scale = s2_outcome_scale,
      s2_standardize_y = s2_standardize_y,
      s2_standardize_m = s2_standardize_m
    )
    write_csv(dat, data_path)
  }

  if (identical(ablation, "product")) {
    rep_rows <- if (nrow(raw)) raw |> filter(rep == r) else tibble()
    done_variants <- unique(rep_rows$variant)
    done_recon <- if (nrow(rep_rows)) unique(rep_rows$recon_rep) else integer()
    if (all(variants$variant %in% done_variants) && length(done_recon) >= product_recon_rep) next

    seed_split <- seed_data + 101L
    seed_init <- seed_data + 102L
    fit_t <- system.time({
      fit <- dcma(
        data = as.data.frame(dat),
        epsm_dim = epsm_dim,
        epsy_dim = epsy_dim,
        c_vars = "Z",
        standardize = TRUE,
        silent = TRUE,
        m_generator = "plain",
        y_generator = "plain",
        split_seed = seed_split,
        init_seed = seed_init,
        device = fit_device
      )
    })

    for (q in seq_len(product_recon_rep)) {
      for (i in seq_len(nrow(variants))) {
        v <- variants[i, ]
        if (has_raw_rows(rep == r, recon_rep == q, variant == v$variant)) next
        seed_recon <- seed_data + 1000L * q + ifelse(v$variant == "permute", 201L, 301L)
        seed_ed <- seed_data + 1000L * q + ifelse(v$variant == "permute", 202L, 302L)
        recon_t <- system.time({
          est <- dcma_reconstruct_interventions(
            fit = fit,
            B = B_est,
            seed = seed_recon,
            device = recon_device,
            product_mode = v$product_mode
          )
        })
        one_rep <- compute_interventional_basic_effects(
          truth_outcomes = truth$outcomes,
          est_outcomes = est$outcomes,
          n_ed = n_ed,
          ed_seed = seed_ed
        )
        one_rep <- bind_rows(
          one_rep,
          compute_path_specific_effects(
            truth_outcomes = truth$outcomes,
            est_outcomes = est$outcomes,
            n_ed = n_ed,
            ed_seed = seed_ed + 1000L
          )
        ) |>
          mutate(
            ablation = ablation,
            variant = v$variant,
            scenario = "S2",
            rep = r,
            recon_rep = q,
            n_est = n_est,
            m_generator = v$m_generator,
            y_generator = v$y_generator,
            product_mode = v$product_mode,
            fit_sec = unname(fit_t[["elapsed"]]),
            recon_sec = unname(recon_t[["elapsed"]]),
            .before = 1
          )
        raw <- bind_rows(raw, one_rep) |> arrange(rep, recon_rep, variant, effect)
        runtime <- bind_rows(
          runtime,
          tibble(
            ablation = ablation,
            variant = v$variant,
            rep = r,
            recon_rep = q,
            fit_sec = unname(fit_t[["elapsed"]]),
            recon_sec = unname(recon_t[["elapsed"]])
          )
        )
        write_csv(raw, raw_path)
        write_csv(runtime, runtime_path)
        write_csv(summarise_raw(raw), summary_path)
        write_product_comparison(raw, product_cmp_path)
      }
    }
    cat(sprintf("[done] product rep=%03d\n", r))
    next
  }

  for (i in seq_len(nrow(variants))) {
    v <- variants[i, ]
    if (has_raw_rows(rep == r, variant == v$variant)) next
    seed_split <- seed_data + 1000L * i + 101L
    seed_init <- seed_data + 1000L * i + 102L
    seed_recon <- seed_data + 1000L * i + 103L
    seed_ed <- seed_data + 1000L * i + 104L

    cat(sprintf("[run] %s rep=%03d variant=%s\n", ablation, r, v$variant))
    fit_t <- system.time({
      fit <- dcma(
        data = as.data.frame(dat),
        epsm_dim = epsm_dim,
        epsy_dim = epsy_dim,
        c_vars = "Z",
        standardize = TRUE,
        silent = TRUE,
        m_generator = v$m_generator,
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
        device = recon_device,
        product_mode = v$product_mode
      )
    })
    one_rep <- compute_interventional_basic_effects(
      truth_outcomes = truth$outcomes,
      est_outcomes = est$outcomes,
      n_ed = n_ed,
      ed_seed = seed_ed
    )
    one_rep <- bind_rows(
      one_rep,
      compute_path_specific_effects(
        truth_outcomes = truth$outcomes,
        est_outcomes = est$outcomes,
        n_ed = n_ed,
        ed_seed = seed_ed + 1000L
      )
    ) |>
      mutate(
        ablation = ablation,
        variant = v$variant,
        scenario = "S2",
        rep = r,
        recon_rep = 1L,
        n_est = n_est,
        m_generator = v$m_generator,
        y_generator = v$y_generator,
        product_mode = v$product_mode,
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
    cat(sprintf("[done] %s rep=%03d variant=%s\n", ablation, r, v$variant))
  }
}

write_csv(summarise_raw(raw), summary_path)
if (identical(ablation, "product")) write_product_comparison(raw, product_cmp_path)
cat("Wrote ablation results to:\n")
cat(out_dir, "\n")
