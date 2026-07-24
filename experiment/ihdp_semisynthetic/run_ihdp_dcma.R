suppressPackageStartupMessages({
  library(torch)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

root <- Sys.getenv("DCMA_CODE_ROOT", ".")
objective <- tolower(Sys.getenv("TNNLS_COVTAIL_OBJECTIVE", "es"))
if (!objective %in% c("es", "wgr")) {
  stop("TNNLS_COVTAIL_OBJECTIVE must be one of: es, wgr")
}
manifest_path <- Sys.getenv(
  "TNNLS_COVTAIL_REP100_MANIFEST",
  file.path(root, "data", "ihdp_covtail_skew_tnnls_v7d_n5000_rep100", "manifest.csv")
)
out_dir <- Sys.getenv(
  "TNNLS_COVTAIL_REP100_OUT_DIR",
  if (identical(objective, "es")) {
    file.path(root, "results", "covtail_skew_tnnls_v7d_estimated_density_rep100_smooth")
  } else {
    file.path(root, "results", paste0("covtail_skew_tnnls_v7d_", objective, "_rep100_smooth"))
  }
)
ref_dir <- Sys.getenv(
  "TNNLS_COVTAIL_REF_DIR",
  file.path(root, "results", "covtail_skew_tnnls_v7d_estimated_density_rep5_smooth")
)
scenario <- Sys.getenv("TNNLS_COVTAIL_SCENARIO", "IHDPMechanismThreeChannelSmoothTail")
family_label <- Sys.getenv("IHDP_S2_SMOOTHTAIL_FAMILY", "covtail")
reps <- as.integer(Sys.getenv("TNNLS_COVTAIL_REPS", "100"))
B_est <- as.integer(Sys.getenv("TNNLS_COVTAIL_B_EST", "160"))
n_sample <- as.integer(Sys.getenv("TNNLS_COVTAIL_SAMPLE_N", "70000"))
n_grid <- as.integer(Sys.getenv("TNNLS_COVTAIL_GRID_N", "512"))
density_adjust <- as.numeric(Sys.getenv("TNNLS_COVTAIL_DENSITY_ADJUST", "0.95"))
epsm_dim <- as.integer(Sys.getenv("TNNLS_COVTAIL_EPSM", "4"))
epsy_dim <- as.integer(Sys.getenv("TNNLS_COVTAIL_EPSY", "32"))
beta_Y <- as.numeric(Sys.getenv("TNNLS_COVTAIL_BETA_Y", "1"))
device <- Sys.getenv("TNNLS_COVTAIL_DEVICE", "cpu")
lambda_w <- as.numeric(Sys.getenv("TNNLS_COVTAIL_WGR_LAMBDA_W", "0.9"))
lambda_l <- as.numeric(Sys.getenv("TNNLS_COVTAIL_WGR_LAMBDA_L", "0.1"))
wgr_J_size <- as.integer(Sys.getenv("TNNLS_COVTAIL_WGR_J_SIZE", "4"))
model_source_label <- if (identical(objective, "es")) "Full DCMA" else paste0("Full DCMA-", toupper(objective))
oracle_source_label <- if (identical(objective, "es")) "Oracle M + Model Y" else paste0("Oracle M + Model Y-", toupper(objective))

Sys.setenv(
  IHDP_S2_SMOOTHTAIL_FAMILY = Sys.getenv("IHDP_S2_SMOOTHTAIL_FAMILY", "covtail"),
  IHDP_S2_SMOOTHTAIL_M1_X = Sys.getenv("IHDP_S2_SMOOTHTAIL_M1_X", "0.55"),
  IHDP_S2_SMOOTHTAIL_M2_X = Sys.getenv("IHDP_S2_SMOOTHTAIL_M2_X", "0.60"),
  IHDP_S2_SMOOTHTAIL_M3_X = Sys.getenv("IHDP_S2_SMOOTHTAIL_M3_X", "0.65"),
  IHDP_S2_SMOOTHTAIL_M_NOISE = Sys.getenv("IHDP_S2_SMOOTHTAIL_M_NOISE", "1.10"),
  IHDP_S2_SMOOTHTAIL_M2_SIGMA_COEF = Sys.getenv("IHDP_S2_SMOOTHTAIL_M2_SIGMA_COEF", "0.55"),
  IHDP_S2_COVTAIL_Q = Sys.getenv("IHDP_S2_COVTAIL_Q", "0.88"),
  IHDP_S2_COVTAIL_C_SLOPE = Sys.getenv("IHDP_S2_COVTAIL_C_SLOPE", "5.00"),
  IHDP_S2_COVTAIL_AMP_MAX = Sys.getenv("IHDP_S2_COVTAIL_AMP_MAX", "10.00"),
  IHDP_S2_COVTAIL_SLOPE = Sys.getenv("IHDP_S2_COVTAIL_SLOPE", "5.00"),
  IHDP_S2_COVTAIL_CENTER = Sys.getenv("IHDP_S2_COVTAIL_CENTER", "-0.05"),
  IHDP_S2_COVTAIL_SCALE_GROUP = Sys.getenv("IHDP_S2_COVTAIL_SCALE_GROUP", "false"),
  IHDP_S2_COVTAIL_NORMALIZE = Sys.getenv("IHDP_S2_COVTAIL_NORMALIZE", "false")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
part_dir <- file.path(out_dir, "rep_parts")
dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(root, "utils", "nn_model.R"))
source(file.path(root, "utils", "nn_critic.R"))
source(file.path(root, "utils", "energyloss_es.R"))
source(file.path(root, "model", "dcma_fit.R"))
source(file.path(root, "model", "dcma_fit_objective.R"))
source(file.path(root, "model", "dcma_reconstruct.R"))
source(file.path(root, "experiment", "ihdp_semisynthetic", "ihdp_semisynth_3m.R"))

path_labs <- c(
  IPSE1 = "M1 path: location",
  IPSE2 = "M2 path: scale",
  IPSE3 = "M3 path: right-tail shape"
)
path_levels <- unname(path_labs[c("IPSE1", "IPSE2", "IPSE3")])

sample_vec <- function(x, n, seed) {
  x <- as.numeric(x)
  if (length(x) <= n) return(x)
  set.seed(as.integer(seed))
  x[sample.int(length(x), n)]
}

standardize_M_array <- function(fit, M_array) {
  if (!isTRUE(fit$standardize)) return(M_array)
  out <- M_array
  muM <- as.numeric(fit$mu[fit$m_vars])
  sdM <- as.numeric(fit$sdd[fit$m_vars])
  sdM[!is.finite(sdM) | sdM == 0] <- 1
  for (j in seq_along(muM)) out[, , j] <- (out[, , j] - muM[[j]]) / sdM[[j]]
  out
}

flatten_m_i_major <- function(M_array) {
  S <- dim(M_array)[3]
  do.call(cbind, lapply(seq_len(S), function(j) as.vector(t(M_array[, , j]))))
}

draw_model_Y_given_M_a <- function(fit, M_array_raw, a, C_raw, seed) {
  set.seed(as.integer(seed))
  torch::torch_manual_seed(as.integer(seed))
  c_vars <- fit$config$c_vars
  C_mat <- as.matrix(C_raw[, c_vars, drop = FALSE])
  if (isTRUE(fit$standardize)) {
    muC <- as.numeric(fit$mu[c_vars])
    sdC <- as.numeric(fit$sdd[c_vars])
    sdC[!is.finite(sdC) | sdC == 0] <- 1
    C_mat <- sweep(C_mat, 2, muC, "-")
    C_mat <- sweep(C_mat, 2, sdC, "/")
  }
  M_array <- standardize_M_array(fit, M_array_raw)
  n <- dim(M_array)[1]
  B_local <- dim(M_array)[2]
  M_flat <- flatten_m_i_major(M_array)
  C_rep <- C_mat[rep(seq_len(n), each = B_local), , drop = FALSE]
  X_rep <- matrix(a, nrow = n * B_local, ncol = 1)
  input <- torch_cat(list(
    torch_tensor(cbind(X_rep, C_rep), dtype = torch_float(), device = device),
    torch_tensor(M_flat, dtype = torch_float(), device = device),
    torch_randn(c(n * B_local, fit$epsy_dim), device = device)
  ), dim = 2)
  Y <- matrix(
    as.numeric(as.array(fit$gen_g$to(device = device)(input)$to(device = "cpu"))),
    nrow = n,
    ncol = B_local,
    byrow = TRUE
  )
  if (isTRUE(fit$standardize)) {
    sdY <- as.numeric(fit$sdd["Y"])
    muY <- as.numeric(fit$mu["Y"])
    if (!is.finite(sdY) || sdY == 0) sdY <- 1
    Y <- Y * sdY + muY
  }
  Y
}

path_outcomes <- function(M0, M1, draw_Y, seed_offset) {
  S <- dim(M0)[3]
  path <- vector("list", S)
  for (s in seq_len(S)) {
    pair <- build_ihdp_pathway_M_pair(M0, M1, s)
    path[[s]] <- list(
      Y_a_M0s = draw_Y(pair$M0s, 1, seed_offset + 2L * s),
      Y_a_M1s = draw_Y(pair$M1s, 1, seed_offset + 2L * s + 1L)
    )
  }
  names(path) <- paste0("M", seq_len(S))
  list(path = path)
}

collect_density_samples <- function(outcomes, source_label, rep_id, source_order) {
  bind_rows(lapply(seq_along(outcomes$path), function(s) {
    bind_rows(
      tibble(
        rep = rep_id,
        source = source_label,
        effect = paste0("IPSE", s),
        distribution = "Before switch",
        value = sample_vec(outcomes$path[[s]]$Y_a_M0s, n_sample, 100000L * rep_id + 1000L * source_order + 10L + s)
      ),
      tibble(
        rep = rep_id,
        source = source_label,
        effect = paste0("IPSE", s),
        distribution = "After switch",
        value = sample_vec(outcomes$path[[s]]$Y_a_M1s, n_sample, 100000L * rep_id + 1000L * source_order + 20L + s)
      )
    )
  }))
}

stats_vec <- function(y) {
  y <- as.numeric(y)
  s <- stats::sd(y)
  if (!is.finite(s) || s <= 0) s <- 1
  z <- (y - mean(y)) / s
  q <- stats::quantile(y, c(0.10, 0.90), type = 8, names = FALSE)
  c(
    mean = mean(y),
    sd = s,
    iqr90 = q[[2]] - q[[1]],
    upper25 = mean(z > 2.5),
    lower25 = mean(z < -2.5),
    asym25 = mean(z > 2.5) - mean(z < -2.5)
  )
}

signature_from_pair <- function(before, after) {
  b <- stats_vec(before)
  a <- stats_vec(after)
  tibble(
    metric = c("Mean shift", "Log SD ratio", "Rel. Q90-Q10", "Std P(Z>2.5)", "Std P(Z<-2.5)", "Std tail asym."),
    value = c(
      a[["mean"]] - b[["mean"]],
      log(a[["sd"]] / b[["sd"]]),
      (a[["iqr90"]] - b[["iqr90"]]) / b[["iqr90"]],
      a[["upper25"]] - b[["upper25"]],
      a[["lower25"]] - b[["lower25"]],
      a[["asym25"]] - b[["asym25"]]
    )
  )
}

make_range_df <- function() {
  ref_grid <- file.path(ref_dir, "estimated_density_grid_average.csv")
  if (file.exists(ref_grid)) {
    read_csv(ref_grid, show_col_types = FALSE) |>
      mutate(effect_label = as.character(effect_label)) |>
      group_by(effect_label) |>
      summarise(
        xmin = min(x, na.rm = TRUE),
        xmax = max(x, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    tibble(
      effect_label = path_levels,
      xmin = c(-12, -20, -12),
      xmax = c(28, 35, 28)
    )
  }
}

range_df <- make_range_df()

make_density <- function(df) {
  effect_name <- unique(as.character(df$effect_label))
  rg <- range_df |> filter(.data$effect_label == effect_name)
  if (nrow(rg) != 1L) stop("Missing density range for ", effect_name)
  den <- density(df$value, from = rg$xmin, to = rg$xmax, n = n_grid, adjust = density_adjust)
  tibble(x = den$x, density = den$y)
}

process_samples <- function(samples) {
  samples <- samples |>
    mutate(
      source = factor(source, levels = c("Truth", model_source_label, oracle_source_label)),
      distribution = factor(distribution, levels = c("Before switch", "After switch")),
      effect_label = recode(effect, !!!path_labs),
      effect_label = factor(effect_label, levels = path_levels)
    )

  density_grid <- samples |>
    group_by(rep, source, effect, effect_label, distribution) |>
    group_modify(~make_density(bind_cols(.y, .x))) |>
    ungroup()

  signature <- samples |>
    group_by(rep, source, effect) |>
    group_modify(~{
      before <- .x |> filter(distribution == "Before switch") |> pull(value)
      after <- .x |> filter(distribution == "After switch") |> pull(value)
      signature_from_pair(before, after)
    }) |>
    ungroup()

  before_stats <- samples |>
    filter(effect == "IPSE3", distribution == "Before switch") |>
    group_by(rep, source) |>
    summarise(mu0 = mean(value), sd0 = sd(value), .groups = "drop")
  thresholds <- seq(1.25, 4.00, by = 0.25)
  tail <- samples |>
    filter(effect == "IPSE3") |>
    left_join(before_stats, by = c("rep", "source")) |>
    mutate(z = (value - mu0) / sd0) |>
    group_by(rep, source, distribution) |>
    summarise(
      tail_tbl = list(tibble(threshold = thresholds, prob = vapply(thresholds, function(t) mean(z > t), numeric(1)))),
      .groups = "drop"
    ) |>
    unnest(tail_tbl)

  list(density_grid = density_grid, signature = signature, tail = tail)
}

write_final_outputs <- function() {
  density_files <- sort(list.files(part_dir, pattern = "_density_grid\\.csv$", full.names = TRUE))
  signature_files <- sort(list.files(part_dir, pattern = "_signature\\.csv$", full.names = TRUE))
  tail_files <- sort(list.files(part_dir, pattern = "_tail\\.csv$", full.names = TRUE))
  runtime_files <- sort(list.files(part_dir, pattern = "_runtime\\.csv$", full.names = TRUE))
  loss_files <- sort(list.files(part_dir, pattern = "_loss\\.csv$", full.names = TRUE))

  density_grid <- bind_rows(lapply(density_files, read_csv, show_col_types = FALSE))
  density_avg <- density_grid |>
    group_by(source, effect, effect_label, distribution, x) |>
    summarise(
      lo = quantile(density, 0.10, type = 8, names = FALSE),
      hi = quantile(density, 0.90, type = 8, names = FALSE),
      density = mean(density),
      n_rep = n_distinct(rep),
      .groups = "drop"
    )
  write_csv(density_grid, file.path(out_dir, "estimated_density_grid_by_rep.csv"))
  write_csv(density_avg, file.path(out_dir, "estimated_density_grid_average.csv"))

  raw_sig <- bind_rows(lapply(signature_files, read_csv, show_col_types = FALSE))
  avg_sig <- raw_sig |>
    group_by(source, effect, metric) |>
    summarise(
      lo = stats::quantile(value, 0.10, type = 8, names = FALSE),
      hi = stats::quantile(value, 0.90, type = 8, names = FALSE),
      value = mean(value),
      n_rep = n_distinct(rep),
      .groups = "drop"
    ) |>
    mutate(
      effect = recode(effect,
        IPSE1 = "M1 / IPSE1: location",
        IPSE2 = "M2 / IPSE2: scale",
        IPSE3 = "M3 / IPSE3: right-tail shape"
      ),
      effect = factor(effect, levels = c("M1 / IPSE1: location", "M2 / IPSE2: scale", "M3 / IPSE3: right-tail shape")),
      metric = factor(metric, levels = c("Mean shift", "Log SD ratio", "Rel. Q90-Q10", "Std P(Z>2.5)", "Std P(Z<-2.5)", "Std tail asym.")),
      label = if_else(grepl("P\\(|tail", as.character(metric)), sprintf("%+.1f pp", 100 * value), sprintf("%+.2f", value))
    ) |>
    group_by(metric) |>
    mutate(
      scale_floor = case_when(
        metric == "Mean shift" ~ 0.45,
        metric == "Log SD ratio" ~ 0.18,
        metric == "Rel. Q90-Q10" ~ 0.16,
        metric == "Std tail asym." ~ 0.004,
        TRUE ~ 0.004
      ),
      scale = pmax(max(abs(value), na.rm = TRUE), scale_floor),
      value_scaled = pmax(pmin(value / scale, 1), -1)
    ) |>
    ungroup()
  write_csv(raw_sig, file.path(out_dir, "estimated_signature_by_rep.csv"))
  write_csv(avg_sig, file.path(out_dir, "estimated_signature_heatmap_data.csv"))

  tail <- bind_rows(lapply(tail_files, read_csv, show_col_types = FALSE))
  tail_avg <- tail |>
    group_by(source, distribution, threshold) |>
    summarise(
      lo = quantile(prob, 0.10, type = 8, names = FALSE),
      hi = quantile(prob, 0.90, type = 8, names = FALSE),
      prob = mean(prob),
      n_rep = n_distinct(rep),
      .groups = "drop"
    )
  write_csv(tail, file.path(out_dir, "m3_tail_probability_by_rep.csv"))
  write_csv(tail_avg, file.path(out_dir, "m3_tail_probability_average.csv"))

  runtime <- bind_rows(lapply(runtime_files, read_csv, show_col_types = FALSE))
  write_csv(runtime, file.path(out_dir, "estimated_density_runtime.csv"))
  if (length(loss_files)) {
    loss <- bind_rows(lapply(loss_files, read_csv, show_col_types = FALSE))
    write_csv(loss, file.path(out_dir, "estimated_density_loss_history.csv"))
  }

  cat(sprintf("Finalized %d reps in %s\n", n_distinct(runtime$rep), out_dir))
}

manifest <- read_csv(manifest_path, show_col_types = FALSE) |>
  filter(.data$scenario == .env$scenario) |>
  slice_head(n = reps)
if (nrow(manifest) == 0L) stop("No manifest rows found for scenario: ", scenario)
manifest <- manifest |>
  mutate(
    data_csv = if_else(
      file.exists(.data$data_csv),
      .data$data_csv,
      file.path(
        root,
        sub("^.*?/DCMA_code/", "", .data$data_csv)
      )
    )
  )

cat(sprintf(
  "TNNLS covtail v7d light runner | objective=%s | reps=%d | B=%d | sample=%d | device=%s | out=%s\n",
  objective, nrow(manifest), B_est, n_sample, device, out_dir
))

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  rep_id <- as.integer(row$rep)
  prefix <- file.path(part_dir, sprintf("rep_%03d", rep_id))
  done_path <- paste0(prefix, "_runtime.csv")
  if (file.exists(done_path)) {
    cat(sprintf("[rep %03d/%03d] skip existing\n", rep_id, nrow(manifest)))
    next
  }

  seed_base <- as.integer(row$seed_data) + round(1000 * beta_Y)
  dat <- read_csv(row$data_csv, show_col_types = FALSE)
  c_vars <- grep("^C", names(dat), value = TRUE)
  cat(sprintf("[rep %03d/%03d] fit and reconstruct\n", rep_id, nrow(manifest)))

  set.seed(seed_base)
  torch::torch_manual_seed(seed_base)
  t_fit <- system.time({
    if (identical(objective, "es")) {
      fit <- dcma(
        data = as.data.frame(dat),
        c_vars = c_vars,
        epsm_dim = epsm_dim,
        epsy_dim = epsy_dim,
        beta_M = 1,
        beta_Y = beta_Y,
        y_generator = "plain",
        standardize = TRUE,
        silent = TRUE,
        device = device
      )
    } else {
      fit <- dcma_objective(
        data = as.data.frame(dat),
        objective = objective,
        c_vars = c_vars,
        epsm_dim = epsm_dim,
        epsy_dim = epsy_dim,
        y_generator = "plain",
        lambda_w = lambda_w,
        lambda_l = lambda_l,
        wgr_J_size = wgr_J_size,
        standardize = TRUE,
        silent = TRUE,
        device = device
      )
    }
  })

  c_df <- as_tibble(fit$C_raw)
  colnames(c_df) <- fit$config$c_vars
  t_eval <- system.time({
    full <- dcma_reconstruct_interventions(
      fit,
      B = B_est,
      seed = as.integer(row$seed_data) + 2000000L + round(1000 * beta_Y),
      device = device
    )
    M0_true <- draw_ihdp_semisynth_mediators_given_a(c_df, a = 0, B = B_est, scenario = scenario, seed = as.integer(row$seed_data) + 2100000L)
    M1_true <- draw_ihdp_semisynth_mediators_given_a(c_df, a = 1, B = B_est, scenario = scenario, seed = as.integer(row$seed_data) + 2200000L)
    truth <- path_outcomes(
      M0_true,
      M1_true,
      draw_Y = function(M, a, seed) {
        draw_ihdp_semisynth_outcomes_given_M(c_df, M, a = a, scenario = scenario, seed = as.integer(row$seed_data) + 2300000L + seed)
      },
      seed_offset = 100L
    )
    oracleM_modelY <- path_outcomes(
      M0_true,
      M1_true,
      draw_Y = function(M, a, seed) {
        draw_model_Y_given_M_a(fit, M, a = a, C_raw = c_df, seed = as.integer(row$seed_data) + 2400000L + seed)
      },
      seed_offset = 200L
    )
  })

  samples <- bind_rows(
    collect_density_samples(truth, "Truth", rep_id, 1L),
    collect_density_samples(full$outcomes, model_source_label, rep_id, 2L),
    collect_density_samples(oracleM_modelY, oracle_source_label, rep_id, 3L)
  )
  processed <- process_samples(samples)

  write_csv(processed$density_grid, paste0(prefix, "_density_grid.csv"))
  write_csv(processed$signature, paste0(prefix, "_signature.csv"))
  write_csv(processed$tail, paste0(prefix, "_tail.csv"))
  write_csv(bind_rows(
    tibble(rep = rep_id, stage = "M", epoch = seq_along(fit$loss_M), train_loss = as.numeric(fit$loss_M), val_loss = as.numeric(fit$val_loss_M)),
    tibble(rep = rep_id, stage = "Y", epoch = seq_along(fit$loss_Y), train_loss = as.numeric(fit$loss_Y), val_loss = as.numeric(fit$val_loss_Y))
  ), paste0(prefix, "_loss.csv"))
  write_csv(tibble(rep = rep_id, fit_sec = unname(t_fit[["elapsed"]]), eval_sec = unname(t_eval[["elapsed"]])), done_path)
  cat(sprintf("[rep %03d/%03d] done fit=%.1fs eval=%.1fs\n", rep_id, nrow(manifest), unname(t_fit[["elapsed"]]), unname(t_eval[["elapsed"]])))

  rm(fit, full, M0_true, M1_true, truth, oracleM_modelY, dat, c_df, samples, processed)
  gc()
}

write_final_outputs()
