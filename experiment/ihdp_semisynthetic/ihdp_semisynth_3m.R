suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Load one IHDP realization as a covariate/treatment base
#
# Returns a data frame with:
# - X: binary treatment
# - C1..Cp: selected covariates
load_ihdp_base <- function(
  csv_path = Sys.getenv("IHDP_NPCI_CSV", file.path("data", "ihdp", "ihdp_npci_1.csv")),
  covariate_idx = c(1:6, 7:10)
) {
  raw <- read_csv(csv_path, col_names = FALSE, show_col_types = FALSE)
  stopifnot(ncol(raw) == 30L)

  covariate_idx <- as.integer(covariate_idx)
  stopifnot(length(covariate_idx) >= 3L, all(covariate_idx >= 1L), all(covariate_idx <= 25L))

  x_cols <- 5L + covariate_idx
  cov_df <- raw[, x_cols, drop = FALSE]
  names(cov_df) <- paste0("C", seq_along(covariate_idx))

  tibble(
    row_id = seq_len(nrow(raw)),
    X = as.integer(raw[[1]])
  ) |>
    bind_cols(cov_df)
}

softplus_num <- function(x) log1p(exp(pmin(x, 30)))

ihdp_env_num <- function(name, default) {
  value <- Sys.getenv(name, "")
  if (nzchar(value)) as.numeric(value) else default
}

capped_exp_centered <- function(n, cap = 1.25) {
  cap <- as.numeric(cap)
  pmin(stats::rexp(n, rate = 1) - 1, cap) + exp(-(cap + 1))
}

bounded_quadratic_skew_noise <- function(n, cap = 3.0) {
  # Centered, variance-standardized right-skewed noise with bounded extreme tail.
  # Constants are Monte Carlo approximations for cap = 3.0.
  q <- (stats::rnorm(n)^2 - 1) / sqrt(2)
  q <- pmin(q, cap)
  if (abs(cap - 3.0) < 1e-8) {
    return((q + 0.02749) / 0.87271)
  }
  as.numeric(scale(q))
}

right_hinge_skew_noise <- function(n) {
  # Centered, variance-standardized positive-half residual. This creates a
  # smooth right-skewness channel without relying on rare extreme tail events.
  z <- stats::rnorm(n)
  mu <- 1 / sqrt(2 * pi)
  sd <- sqrt(0.5 - mu^2)
  (pmax(z, 0) - mu) / sd
}

standardized_bernoulli_skew_noise <- function(p) {
  # Two-point residual with conditional mean zero and variance one. Lower event
  # probabilities create stronger right-skewness without changing the variance.
  p <- pmin(pmax(as.numeric(p), 0.03), 0.97)
  b <- stats::rbinom(length(p), size = 1L, prob = p)
  (b - p) / sqrt(p * (1 - p))
}

standardized_skew_mixture <- function(p, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- pmin(pmax(as.numeric(p), 1e-4), 1 - 1e-4)
  n <- length(p)
  hi_mean <- 1.10
  hi_sd <- 0.25
  lo_mean <- -0.45
  lo_sd <- 0.70
  z <- stats::rbinom(n, size = 1L, prob = p)
  u <- ifelse(
    z == 1L,
    stats::rnorm(n, mean = hi_mean, sd = hi_sd),
    stats::rnorm(n, mean = lo_mean, sd = lo_sd)
  )
  cond_mean <- p * hi_mean + (1 - p) * lo_mean
  cond_var <- p * (hi_sd^2 + (hi_mean - cond_mean)^2) +
    (1 - p) * (lo_sd^2 + (lo_mean - cond_mean)^2)
  (u - cond_mean) / sqrt(pmax(cond_var, 1e-8))
}

standardized_skew_normal <- function(alpha, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  alpha <- as.numeric(alpha)
  delta <- alpha / sqrt(1 + alpha^2)
  n <- length(alpha)
  z0 <- stats::rnorm(n)
  z1 <- stats::rnorm(n)
  u <- delta * abs(z0) + sqrt(pmax(1 - delta^2, 1e-8)) * z1
  cond_mean <- sqrt(2 / pi) * delta
  cond_var <- 1 - 2 * delta^2 / pi
  (u - cond_mean) / sqrt(pmax(cond_var, 1e-8))
}

# Generate a 3-mediator IHDP-based semi-synthetic dataset
#
# Scenarios:
# - IHDPMean: primarily mean-shift, mild heteroskedasticity / shape change
# - IHDPDist: stronger distributional complexity via heteroskedasticity + mixture
# - IHDPMechanismStrong: pathway-mechanism stress test with distinct
#   location, scale/tail, and mixture/shape mediator channels
# - IHDPMechanismSafe: reviewer-safe version using synthetic balanced
#   treatment assignment and correlated mediator errors, without latent
#   variables that could be mistaken for unmeasured confounding
# - IHDPMechanismSkew/Tail/Threshold/Display: variants of the reviewer-safe
#   mechanism that replace the M3 mixture channel with skewness, tail-risk,
#   threshold, or stronger display-oriented channels.
# - IHDPMechanismTwoChannelSimple: two-mediator display mechanism where M1
#   shifts location and M2 changes scale/right-tail shape.
# - IHDPMechanismThreeChannelDisplay: three-mediator display mechanism where
#   M1 shifts location, M2 changes scale, and M3 changes right-tail risk.
# - IHDPMechanismThreeChannelSimple: three-mediator display mechanism with a
#   low-order outcome model: M1 shifts location, M2 changes scale, and M3
#   changes right-skewness through a centered quadratic noise term.
# - IHDPMechanismThreeChannelBernoulliSkew: same mediator design, but M3
#   changes the probability of a centered upper-tail component.
# - IHDPMechanismThreeChannelFinal: reviewer-facing three-channel mechanism
#   with location, symmetric spread, and smooth upper-tail pathways.
# - IHDPMechanismThreeChannelOrthogonal: cleaner three-channel mechanism with
#   weakly correlated mediators and more separated location/spread/upper-tail
#   pathway signatures.
# - IHDPMechanismThreeChannelTailShift: clean three-channel mechanism where
#   M1 shifts location, M2 changes spread, and M3 shifts standardized tail
#   risk from the left tail to the right tail while preserving mean/variance.
# - IHDPMechanismThreeChannelSmoothTail: learnable three-channel mechanism
#   where M3 controls a smooth standardized upper-tail deformation instead of
#   a rare discrete tail jump.
generate_ihdp_semisynth_3m <- function(
  base_df,
  scenario = c(
    "IHDPMean", "IHDPDist", "IHDPDistStrong",
    "IHDPMechanismStrong", "IHDPMechanismSafe",
    "IHDPMechanismSkew", "IHDPMechanismTail", "IHDPMechanismThreshold",
    "IHDPMechanismDisplay", "IHDPMechanismDisplayLearnable",
    "IHDPMechanismThreeChannel", "IHDPMechanismTwoChannelSimple",
    "IHDPMechanismThreeChannelDisplay", "IHDPMechanismThreeChannelSimple",
    "IHDPMechanismThreeChannelBernoulliSkew",
    "IHDPMechanismThreeChannelFinal",
    "IHDPMechanismThreeChannelOrthogonal",
    "IHDPMechanismThreeChannelTailShift",
    "IHDPMechanismThreeChannelSmoothTail"
  ),
  seed = NULL
) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)

  stopifnot(is.data.frame(base_df), "X" %in% names(base_df))
  c_vars <- grep("^C", names(base_df), value = TRUE)
  stopifnot(length(c_vars) >= 8L)

  C <- as.matrix(base_df[, c_vars, drop = FALSE])
  n <- nrow(C)
  X <- as.numeric(base_df$X)

  mechanism_safe_family <- c(
    "IHDPMechanismSafe",
    "IHDPMechanismSkew",
    "IHDPMechanismTail",
    "IHDPMechanismThreshold",
    "IHDPMechanismDisplay",
    "IHDPMechanismDisplayLearnable",
    "IHDPMechanismThreeChannel",
    "IHDPMechanismTwoChannelSimple",
    "IHDPMechanismThreeChannelDisplay",
    "IHDPMechanismThreeChannelSimple",
    "IHDPMechanismThreeChannelBernoulliSkew",
    "IHDPMechanismThreeChannelFinal",
    "IHDPMechanismThreeChannelOrthogonal",
    "IHDPMechanismThreeChannelTailShift",
    "IHDPMechanismThreeChannelSmoothTail"
  )

  if (scenario %in% mechanism_safe_family) {
    px <- plogis(0.10 + 0.45 * C[, 1] - 0.35 * C[, 2] +
                   0.25 * C[, 3] - 0.20 * C[, 4] + 0.15 * C[, 5])
    px <- pmin(pmax(px, 0.12), 0.88)
    X <- rbinom(n, size = 1L, prob = px)
  }

  # Shared latent factors induce mediator dependence beyond observed C.
  H1 <- rnorm(n)
  H2 <- rnorm(n)

  if (identical(scenario, "IHDPMechanismTwoChannelSimple")) {
    Sigma_M2 <- matrix(c(1.00, 0.12, 0.12, 1.00), nrow = 2L, byrow = TRUE)
    E2ch <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L) %*% chol(Sigma_M2)
    M1 <- 0.20 + 1.00 * X +
      0.30 * C[, 1] - 0.20 * C[, 2] + 0.45 * E2ch[, 1]
    M2 <- -0.75 + 1.30 * X +
      0.25 * C[, 5] + 0.45 * E2ch[, 2]

    mu <- 1.00 + 0.08 * X + 0.95 * M1 + 0.15 * C[, 1] - 0.10 * C[, 2]
    sigma <- exp(0.22 * M2)
    tail_strength <- 0.15 + 1.05 * plogis(2.20 * M2)
    z_skew <- rnorm(n)
    q_skew <- (z_skew^2 - 1) / sqrt(2)
    Y <- mu + sigma * rnorm(n) + tail_strength * q_skew
    S_latent <- as.numeric(M2 > stats::median(M2))

    return(
      base_df |>
        transmute(
          row_id = row_id,
          X = as.integer(.env$X),
          !!!as.data.frame(C),
          M1 = M1,
          M2 = M2,
          S_latent = S_latent,
          Y = Y
        )
    )
  }

  if (scenario %in% c(
    "IHDPMechanismThreeChannelDisplay",
    "IHDPMechanismThreeChannelSimple",
    "IHDPMechanismThreeChannelBernoulliSkew",
    "IHDPMechanismThreeChannelFinal",
    "IHDPMechanismThreeChannelOrthogonal",
    "IHDPMechanismThreeChannelTailShift",
    "IHDPMechanismThreeChannelSmoothTail"
  )) {
    Sigma_M3 <- if (identical(scenario, "IHDPMechanismThreeChannelOrthogonal") ||
                    identical(scenario, "IHDPMechanismThreeChannelTailShift") ||
                    identical(scenario, "IHDPMechanismThreeChannelSmoothTail")) {
      diag(3L)
    } else {
      matrix(c(
        1.00, 0.12, 0.08,
        0.12, 1.00, 0.10,
        0.08, 0.10, 1.00
      ), nrow = 3L, byrow = TRUE)
    }
    E3ch <- matrix(rnorm(n * 3L), nrow = n, ncol = 3L) %*% chol(Sigma_M3)

    if (identical(scenario, "IHDPMechanismThreeChannelTailShift") ||
        identical(scenario, "IHDPMechanismThreeChannelSmoothTail")) {
      smooth_m1_x <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_M1_X", 0.90)
      smooth_m2_x <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_M2_X", 1.00)
      smooth_m3_x <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_M3_X", 2.40)
      smooth_m_noise <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_M_NOISE", 0.34)
      M1 <- 0.10 + smooth_m1_x * X +
        0.20 * C[, 1] - 0.15 * C[, 2] + smooth_m_noise * E3ch[, 1]
      M2 <- -0.15 + smooth_m2_x * X +
        0.10 * C[, 4] + smooth_m_noise * E3ch[, 2]
      m3_x <- if (identical(scenario, "IHDPMechanismThreeChannelSmoothTail")) smooth_m3_x else 2.60
      M3 <- -0.70 + m3_x * X +
        0.08 * C[, 5] + smooth_m_noise * E3ch[, 3]
    } else if (identical(scenario, "IHDPMechanismThreeChannelOrthogonal")) {
      M1 <- 0.15 + 0.95 * X +
        0.22 * C[, 1] - 0.16 * C[, 2] + 0.36 * E3ch[, 1]
      M2 <- -0.30 + 1.05 * X +
        0.12 * C[, 4] + 0.34 * E3ch[, 2]
      M3 <- -0.85 + 1.35 * X +
        0.12 * C[, 5] + 0.34 * E3ch[, 3]
    } else if (identical(scenario, "IHDPMechanismThreeChannelFinal")) {
      M1 <- 0.20 + 0.90 * X +
        0.25 * C[, 1] - 0.18 * C[, 2] + 0.42 * E3ch[, 1]
      M2 <- -0.45 + 1.05 * X +
        0.16 * C[, 4] + 0.45 * E3ch[, 2]
      M3 <- -0.70 + 1.25 * X +
        0.18 * C[, 5] + 0.45 * E3ch[, 3]
    } else {
      M1 <- 0.20 + 0.95 * X +
        0.30 * C[, 1] - 0.20 * C[, 2] + 0.45 * E3ch[, 1]
      M2 <- -0.60 + 1.25 * X +
        0.20 * C[, 4] + 0.45 * E3ch[, 2]
      M3 <- -1.00 + 1.65 * X +
        0.25 * C[, 5] + 0.45 * E3ch[, 3]
    }

    if (identical(scenario, "IHDPMechanismThreeChannelSmoothTail")) {
      mu <- 1.00 + 0.98 * M1 + 0.08 * C[, 1] - 0.06 * C[, 2]
      sigma <- exp(ihdp_env_num("IHDP_S2_SMOOTHTAIL_M2_SIGMA_COEF", 0.22) * M2)
      smooth_family <- Sys.getenv("IHDP_S2_SMOOTHTAIL_FAMILY", "hermite")
      gamma_max <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_GAMMA_MAX", 0.80)
      shape_slope <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_SHAPE_SLOPE", 3.00)
      tail_center <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_CENTER", 0.50)
      gamma <- gamma_max * (2 * plogis(shape_slope * (M3 - tail_center)) - 1)
      if (identical(smooth_family, "additive")) {
        tail_lambda <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_TAIL_LAMBDA", 0.70)
        z_scale <- rnorm(n)
        z_tail <- rnorm(n)
        q_tail <- (z_tail^2 - 1) / sqrt(2)
        Y <- mu + sigma * z_scale + tail_lambda * gamma * q_tail
      } else if (identical(smooth_family, "covskew")) {
        shape_lambda <- ihdp_env_num("IHDP_S2_COVSKEW_LAMBDA", 0.55)
        h_shape <- as.numeric(scale(C[, 5]^2))
        z_tail <- rnorm(n)
        Y <- mu + sigma * z_tail + shape_lambda * gamma * h_shape
      } else if (identical(smooth_family, "bimodal")) {
        noise_sd <- ihdp_env_num("IHDP_S2_BIMODAL_NOISE_SD", 0.42)
        sep <- exp(gamma)
        signs <- ifelse(runif(n) < 0.5, -1, 1)
        raw_tail <- signs * sep + rnorm(n, sd = noise_sd)
        z_tail <- raw_tail / sqrt(sep^2 + noise_sd^2)
      } else if (identical(smooth_family, "s1bimodal")) {
        noise_sd <- ihdp_env_num("IHDP_S2_S1BIMODAL_NOISE_SD", 0.35)
        sep_min <- ihdp_env_num("IHDP_S2_S1BIMODAL_SEP_MIN", 0.10)
        sep_amp <- ihdp_env_num("IHDP_S2_S1BIMODAL_SEP_AMP", 1.35)
        sep_slope <- ihdp_env_num("IHDP_S2_S1BIMODAL_SLOPE", 2.80)
        sep_center <- ihdp_env_num("IHDP_S2_S1BIMODAL_CENTER", 0.50)
        sep <- sep_min + sep_amp * plogis(sep_slope * (M3 - sep_center))
        signs <- ifelse(runif(n) < 0.5, -1, 1)
        raw_tail <- signs * sep + rnorm(n, sd = noise_sd)
        z_tail <- raw_tail / sqrt(mean(sep^2 + noise_sd^2))
      } else if (identical(smooth_family, "covbimodal")) {
        split_var <- C[, 5]
        split_group <- ifelse(split_var >= median(split_var), 1, -1)
        split_group <- split_group - mean(split_group)
        split_group <- split_group / sd(split_group)
        amp_max <- ihdp_env_num("IHDP_S2_COVBIMODAL_AMP_MAX", 1.35)
        amp_slope <- ihdp_env_num("IHDP_S2_COVBIMODAL_SLOPE", 3.00)
        amp_center <- ihdp_env_num("IHDP_S2_COVBIMODAL_CENTER", 0.40)
        normalize_covbimodal <- tolower(Sys.getenv("IHDP_S2_COVBIMODAL_NORMALIZE", "true")) %in% c("1", "true", "yes")
        amp <- amp_max * plogis(amp_slope * (M3 - amp_center))
        raw_tail <- rnorm(n) + amp * split_group
        z_tail <- if (normalize_covbimodal) raw_tail / sqrt(pmax(1 + amp^2, 1e-8)) else raw_tail
      } else if (identical(smooth_family, "covtail")) {
        split_var <- C[, 5]
        tail_q <- ihdp_env_num("IHDP_S2_COVTAIL_Q", 0.75)
        tail_cut <- as.numeric(quantile(split_var, tail_q, type = 8, names = FALSE))
        tail_slope <- ihdp_env_num("IHDP_S2_COVTAIL_C_SLOPE", 5.00)
        tail_gate <- plogis(tail_slope * (split_var - tail_cut))
        tail_group <- tail_gate - mean(tail_gate)
        scale_tail_group <- tolower(Sys.getenv("IHDP_S2_COVTAIL_SCALE_GROUP", "false")) %in% c("1", "true", "yes")
        if (scale_tail_group) tail_group <- tail_group / sd(tail_group)
        amp_max <- ihdp_env_num("IHDP_S2_COVTAIL_AMP_MAX", 2.80)
        amp_slope <- ihdp_env_num("IHDP_S2_COVTAIL_SLOPE", 5.00)
        amp_center <- ihdp_env_num("IHDP_S2_COVTAIL_CENTER", -0.50)
        normalize_covtail <- tolower(Sys.getenv("IHDP_S2_COVTAIL_NORMALIZE", "false")) %in% c("1", "true", "yes")
        amp <- amp_max * plogis(amp_slope * (M3 - amp_center))
        raw_tail <- rnorm(n) + amp * tail_group
        z_tail <- if (normalize_covtail) raw_tail / sqrt(pmax(1 + amp^2 * stats::var(tail_group), 1e-8)) else raw_tail
      } else if (identical(smooth_family, "coretail")) {
        p0 <- ihdp_env_num("IHDP_S2_CORETAIL_P0", 0.15)
        core_sd <- ihdp_env_num("IHDP_S2_CORETAIL_CORE_SD", 0.55)
        tail_sd <- ihdp_env_num("IHDP_S2_CORETAIL_TAIL_SD", 2.30)
        p_tail <- plogis(qlogis(p0) + gamma)
        is_tail <- runif(n) < p_tail
        raw_tail <- ifelse(is_tail, rnorm(n, sd = tail_sd), rnorm(n, sd = core_sd))
        tail_var <- p_tail * tail_sd^2 + (1 - p_tail) * core_sd^2
        z_tail <- raw_tail / sqrt(pmax(tail_var, 1e-8))
      } else if (identical(smooth_family, "twopiece")) {
        z <- rnorm(n)
        right_scale <- exp(gamma)
        left_scale <- exp(-gamma)
        raw_tail <- ifelse(z >= 0, right_scale * z, left_scale * z)
        tail_mean <- (right_scale - left_scale) / sqrt(2 * pi)
        tail_second <- 0.5 * (right_scale^2 + left_scale^2)
        z_tail <- (raw_tail - tail_mean) / sqrt(pmax(tail_second - tail_mean^2, 1e-8))
      } else if (identical(smooth_family, "skewnormal")) {
        z_tail <- standardized_skew_normal(gamma)
      } else if (identical(smooth_family, "hinge")) {
        z <- rnorm(n)
        hinge_mu <- 1 / sqrt(2 * pi)
        hinge_sd <- sqrt(0.5 - hinge_mu^2)
        hinge_cov <- 0.5 / hinge_sd
        h <- (pmax(z, 0) - hinge_mu) / hinge_sd
        z_tail <- (z + gamma * h) / sqrt(pmax(1 + gamma^2 + 2 * gamma * hinge_cov, 1e-8))
      } else if (identical(smooth_family, "residualshape")) {
        z_scale <- rnorm(n)
        z_sym <- rnorm(n)
        z_shape <- right_hinge_skew_noise(n)
        shape_amp <- ihdp_env_num("IHDP_S2_RESIDUALSHAPE_AMP", 0.65)
        residual_slope <- ihdp_env_num("IHDP_S2_RESIDUALSHAPE_SLOPE", 3.40)
        p_shape <- plogis(residual_slope * M3)
        shape_noise <- ((1 - p_shape) * z_sym + p_shape * z_shape) /
          sqrt(pmax((1 - p_shape)^2 + p_shape^2, 1e-8))
        Y <- mu + sigma * z_scale + shape_amp * shape_noise
      } else {
        z <- rnorm(n)
        z_tail <- (z + gamma * ((z^2 - 1) / sqrt(2))) / sqrt(1 + gamma^2)
      }
      if (!identical(smooth_family, "additive") && !identical(smooth_family, "residualshape")) {
        Y <- mu + sigma * z_tail
      }
    } else if (identical(scenario, "IHDPMechanismThreeChannelTailShift")) {
      mu <- 1.00 + 0.98 * M1 + 0.08 * C[, 1] - 0.06 * C[, 2]
      sigma <- exp(ihdp_env_num("IHDP_S2_TAILSHIFT_M2_SIGMA_COEF", 0.22) * M2)
      q_tail <- ihdp_env_num("IHDP_S2_TAILSHIFT_Q_TAIL", 0.06)
      shape_slope <- ihdp_env_num("IHDP_S2_TAILSHIFT_SHAPE_SLOPE", 9.00)
      core_sd <- ihdp_env_num("IHDP_S2_TAILSHIFT_CORE_SD", 0.45)
      tail_loc <- ihdp_env_num("IHDP_S2_TAILSHIFT_TAIL_LOC", 6.20)
      tail_sd <- ihdp_env_num("IHDP_S2_TAILSHIFT_TAIL_SD", 0.25)
      tail_center <- ihdp_env_num("IHDP_S2_TAILSHIFT_CENTER", 0.60)
      right_prob <- plogis(shape_slope * (M3 - tail_center))
      p_right <- q_tail * right_prob
      p_left <- q_tail * (1 - right_prob)
      u_tail <- runif(n)
      raw_tail <- rnorm(n, mean = 0, sd = core_sd)
      is_right <- u_tail < p_right
      is_left <- u_tail >= p_right & u_tail < p_right + p_left
      raw_tail[is_right] <- rnorm(sum(is_right), mean = tail_loc, sd = tail_sd)
      raw_tail[is_left] <- rnorm(sum(is_left), mean = -tail_loc, sd = tail_sd)
      mean_tail <- (p_right - p_left) * tail_loc
      second_tail <- (1 - q_tail) * core_sd^2 +
        p_right * (tail_sd^2 + tail_loc^2) +
        p_left * (tail_sd^2 + tail_loc^2)
      z_tail <- (raw_tail - mean_tail) / sqrt(pmax(second_tail - mean_tail^2, 1e-8))
      Y <- mu + sigma * z_tail
    } else if (identical(scenario, "IHDPMechanismThreeChannelOrthogonal")) {
      mu <- 1.00 + 0.05 * X + 0.90 * M1 + 0.10 * C[, 1] - 0.08 * C[, 2]
      sigma <- exp(ihdp_env_num("IHDP_S2_ORTHO_M2_SIGMA_COEF", 0.18) * M2)
      tail_base <- ihdp_env_num("IHDP_S2_ORTHO_M3_TAIL_BASE", 0.04)
      tail_amp <- ihdp_env_num("IHDP_S2_ORTHO_M3_TAIL_AMP", 1.35)
      tail_slope <- ihdp_env_num("IHDP_S2_ORTHO_M3_TAIL_SLOPE", 2.35)
      tail_strength <- tail_base + tail_amp * plogis(tail_slope * M3)
      tail_noise <- right_hinge_skew_noise(n)
      tail_noise <- tail_noise - mean(tail_noise)
      Y <- mu + sigma * rnorm(n) + tail_strength * tail_noise
    } else if (identical(scenario, "IHDPMechanismThreeChannelFinal")) {
      mu <- 1.00 + 0.06 * X + 0.85 * M1 + 0.12 * C[, 1] - 0.08 * C[, 2]
      sigma <- exp(ihdp_env_num("IHDP_S2_FINAL_M2_SIGMA_COEF", 0.22) * M2)
      tail_base <- ihdp_env_num("IHDP_S2_FINAL_M3_TAIL_BASE", 0.12)
      tail_amp <- ihdp_env_num("IHDP_S2_FINAL_M3_TAIL_AMP", 1.00)
      tail_slope <- ihdp_env_num("IHDP_S2_FINAL_M3_TAIL_SLOPE", 1.80)
      tail_strength <- tail_base + tail_amp * plogis(tail_slope * M3)
      Y <- mu + sigma * rnorm(n) + tail_strength * right_hinge_skew_noise(n)
    } else {
      mu <- 1.00 + 0.08 * X + 0.85 * M1 + 0.15 * C[, 1] - 0.10 * C[, 2]
      sigma <- exp(ihdp_env_num("IHDP_S2_M2_SIGMA_COEF", 0.34) * M2)
    }
    if (identical(scenario, "IHDPMechanismThreeChannelFinal") ||
        identical(scenario, "IHDPMechanismThreeChannelTailShift") ||
        identical(scenario, "IHDPMechanismThreeChannelSmoothTail") ||
        identical(scenario, "IHDPMechanismThreeChannelOrthogonal")) {
      # Y has already been generated by the final reviewer-facing mechanism.
    } else if (identical(scenario, "IHDPMechanismThreeChannelSimple")) {
      z_scale <- rnorm(n)
      skew_base <- ihdp_env_num("IHDP_S2_M3_HINGE_BASE", 0.10)
      skew_amp <- ihdp_env_num("IHDP_S2_M3_HINGE_AMP", 0.90)
      skew_slope <- ihdp_env_num("IHDP_S2_M3_HINGE_SLOPE", 2.25)
      skew_strength <- skew_base + skew_amp * plogis(skew_slope * M3)
      skew_component <- skew_strength * right_hinge_skew_noise(n)
      Y <- mu + sigma * z_scale + skew_component
    } else if (identical(scenario, "IHDPMechanismThreeChannelBernoulliSkew")) {
      z_scale <- rnorm(n)
      p_base <- ihdp_env_num("IHDP_S2_M3_TAIL_BASE", 0.05)
      p_amp <- ihdp_env_num("IHDP_S2_M3_TAIL_AMP", 0.45)
      p_slope <- ihdp_env_num("IHDP_S2_M3_TAIL_SLOPE", 2.00)
      tail_size <- ihdp_env_num("IHDP_S2_M3_TAIL_SIZE", 2.00)
      p_tail <- p_base + p_amp * plogis(p_slope * M3)
      tail_event <- stats::rbinom(n, size = 1L, prob = p_tail)
      skew_component <- tail_size * (tail_event - p_tail)
      Y <- mu + sigma * z_scale + skew_component
    } else {
      tail_strength <- 0.10 + 1.25 * plogis(2.80 * M3)
      tail_noise <- capped_exp_centered(n, cap = 2.25)
      Y <- mu + sigma * rnorm(n) + tail_strength * tail_noise
    }
    S_latent <- as.numeric(M3 > stats::median(M3))

    return(
      base_df |>
        transmute(
          row_id = row_id,
          X = as.integer(.env$X),
          !!!as.data.frame(C),
          M1 = M1,
          M2 = M2,
          M3 = M3,
          S_latent = S_latent,
          Y = Y
        )
    )
  }

  if (scenario %in% mechanism_safe_family) {
    Sigma_M <- matrix(c(
      1.00, 0.12, 0.08,
      0.12, 1.00, 0.10,
      0.08, 0.10, 1.00
    ), nrow = 3L, byrow = TRUE)
    E <- matrix(rnorm(n * 3L), nrow = n, ncol = 3L) %*% chol(Sigma_M)

    m1_x <- if (identical(scenario, "IHDPMechanismDisplay")) {
      1.25
    } else if (identical(scenario, "IHDPMechanismDisplayLearnable")) {
      1.10
    } else if (identical(scenario, "IHDPMechanismThreeChannel")) {
      1.05
    } else {
      0.95
    }
    M1 <- 0.20 + m1_x * X +
      0.30 * C[, 1] - 0.20 * C[, 2] + 0.45 * E[, 1]

    if (identical(scenario, "IHDPMechanismDisplayLearnable")) {
      M2 <- -0.80 + 1.15 * X + 0.25 * C[, 5] + 0.45 * E[, 2]
      mu <- 1.00 + 0.12 * X + 1.05 * M1 + 0.15 * C[, 1] - 0.10 * C[, 2]
      sigma <- 0.90
      skew_strength <- pmax(0.20, 1.05 + 0.75 * M2)
      skew_component <- skew_strength * capped_exp_centered(n, cap = 2.00)
      S_latent <- as.numeric(skew_strength > stats::median(skew_strength))
      Y <- mu + sigma * rnorm(n) + skew_component

      return(
        base_df |>
          transmute(
            row_id = row_id,
            X = as.integer(.env$X),
            !!!as.data.frame(C),
            M1 = M1,
            M2 = M2,
            S_latent = S_latent,
            Y = Y
          )
      )
    }

    if (identical(scenario, "IHDPMechanismThreeChannel")) {
      M2 <- 0.60 * X + 0.20 * C[, 4] + 0.45 * E[, 2]
      M3 <- 0.60 * X + 0.20 * C[, 5] + 0.45 * E[, 3]
      mu <- 1.00 + 0.10 * X + 1.00 * M1 + 0.15 * C[, 1] - 0.10 * C[, 2]
      sigma <- exp(0.35 * M2)
      gamma_skew <- 0.85 * M3
      z_skew <- rnorm(n)
      q_skew <- (z_skew^2 - 1) / sqrt(2)
      skew_component <- gamma_skew * q_skew
      S_latent <- as.numeric(gamma_skew > stats::median(gamma_skew))
      Y <- mu + sigma * rnorm(n) + skew_component

      return(
        base_df |>
          transmute(
            row_id = row_id,
            X = as.integer(.env$X),
            !!!as.data.frame(C),
            M1 = M1,
            M2 = M2,
            M3 = M3,
            S_latent = S_latent,
            Y = Y
          )
      )
    }

    scale2 <- if (identical(scenario, "IHDPMechanismDisplay")) {
      1.45 - 0.95 * X + 0.05 * softplus_num(C[, 4])
    } else {
      1.16 - 0.54 * X + 0.06 * softplus_num(C[, 4])
    }
    M2 <- -0.05 + 0.04 * X - 0.18 * C[, 3] + 0.18 * C[, 4] + scale2 * E[, 2]

    m3_x <- if (identical(scenario, "IHDPMechanismDisplay")) {
      1.65
    } else {
      1.25
    }
    M3 <- -1.20 + m3_x * X +
      0.12 * C[, 5] - 0.10 * C[, 8] +
      0.12 * C[, 6] + 0.34 * E[, 3]
  } else if (identical(scenario, "IHDPMechanismStrong")) {
    M1 <- 0.20 + 1.35 * X + 0.30 * C[, 1] - 0.20 * C[, 2] +
      0.25 * H1 + rnorm(n, sd = 0.45)

    scale2 <- 0.55 + 0.85 * X + 0.12 * softplus_num(C[, 4])
    M2 <- -0.05 + 0.08 * X - 0.18 * C[, 3] + 0.18 * C[, 4] +
      0.18 * H1 + scale2 * rnorm(n)

    p3 <- plogis(-1.35 + 2.65 * X + 0.35 * C[, 5] - 0.20 * C[, 8])
    S3 <- rbinom(n, size = 1L, prob = p3)
    M3 <- ifelse(
      S3 == 1L,
      rnorm(n, mean = 1.20, sd = 0.35),
      rnorm(n, mean = -1.05, sd = 0.40)
    ) + 0.12 * C[, 6] + 0.10 * H1
  } else {
    e1 <- rnorm(n, sd = 0.45)
    e2 <- rnorm(n, sd = 0.50)
    e3 <- rnorm(n, sd = 0.45)

    # Three mediators with different roles:
    # M1: mainly shifts location
    # M2: mainly controls scale / spread
    # M3: mainly controls shape / mixture tendency
    M1 <- 0.20 + 0.90 * X + 0.30 * C[, 1] - 0.20 * C[, 2] + 0.15 * C[, 7] + 0.55 * H1 + 0.15 * H2 + e1
    M2 <- -0.10 + 0.55 * X - 0.15 * C[, 3] + 0.20 * C[, 4] + 0.20 * C[, 8] + 0.20 * H1 + 0.60 * H2 + e2
    M3 <- 0.15 + 0.35 * X + 0.20 * C[, 5] - 0.20 * C[, 6] + 0.10 * C[, 9] + 0.40 * H1 - 0.30 * H2 + e3
  }

  if (scenario %in% mechanism_safe_family) {
    mu <- 1.00 + 0.12 * X +
      ifelse(identical(scenario, "IHDPMechanismDisplay"), 0.95, 0.75) * M1 +
      0.15 * C[, 1] - 0.10 * C[, 2]
    sigma <- if (identical(scenario, "IHDPMechanismDisplay")) {
      0.18 + 0.95 * softplus_num(-0.55 + 1.80 * abs(M2) + 0.04 * C[, 4])
    } else {
      0.22 + 0.70 * softplus_num(-0.70 + 1.35 * abs(M2) + 0.05 * C[, 4])
    }

    if (identical(scenario, "IHDPMechanismSafe")) {
      tail_noise <- rnorm(n)
      p_mix <- plogis(0.35 - 2.60 * abs(M3) + 0.04 * C[, 5])
      S_latent <- rbinom(n, size = 1L, prob = p_mix)
      shape_component <- ifelse(
        S_latent == 1L,
        rnorm(n, mean = 1.55, sd = 0.18),
        rnorm(n, mean = -1.10, sd = 0.22)
      )
      Y <- mu + shape_component + sigma * tail_noise
    } else if (identical(scenario, "IHDPMechanismSkew")) {
      p_skew <- plogis(-0.45 + 1.35 * M3 + 0.05 * C[, 5])
      S_latent <- rbinom(n, size = 1L, prob = p_skew)
      skew_noise <- ifelse(
        S_latent == 1L,
        stats::rexp(n, rate = 1) - 1,
        -(stats::rexp(n, rate = 1) - 1)
      )
      Y <- mu + sigma * rnorm(n) + 0.55 * skew_noise
    } else if (identical(scenario, "IHDPMechanismTail")) {
      p_tail <- plogis(-0.85 + 1.45 * M3 + 0.05 * C[, 5])
      S_latent <- rbinom(n, size = 1L, prob = p_tail)
      tail_component <- ifelse(
        S_latent == 1L,
        stats::rt(n, df = 3) / sqrt(3),
        rnorm(n, sd = 0.55)
      )
      Y <- mu + sigma * rnorm(n) + 0.85 * tail_component
    } else {
      p_jump <- if (identical(scenario, "IHDPMechanismDisplay")) {
        plogis(-0.30 + 2.05 * M3 + 0.05 * C[, 5])
      } else {
        plogis(-0.25 + 1.55 * M3 + 0.05 * C[, 5])
      }
      S_latent <- rbinom(n, size = 1L, prob = p_jump)
      threshold_component <- if (identical(scenario, "IHDPMechanismDisplay")) {
        ifelse(
          S_latent == 1L,
          rnorm(n, mean = 1.20, sd = 0.20),
          rnorm(n, mean = -0.40, sd = 0.20)
        )
      } else {
        ifelse(
          S_latent == 1L,
          rnorm(n, mean = 0.85, sd = 0.22),
          rnorm(n, mean = -0.25, sd = 0.22)
        )
      }
      Y <- mu + sigma * rnorm(n) + threshold_component
    }

    return(
      base_df |>
        transmute(
          row_id = row_id,
          X = as.integer(.env$X),
          !!!as.data.frame(C),
          M1 = M1,
          M2 = M2,
          M3 = M3,
          S_latent = S_latent,
          Y = Y
        )
    )
  }

  if (identical(scenario, "IHDPMechanismStrong")) {
    mu <- 1.20 + 0.20 * X + 0.95 * M1 +
      0.18 * C[, 1] - 0.12 * C[, 2] + 0.08 * C[, 7]
    sigma <- 0.22 + 0.88 * softplus_num(-1.10 + 1.85 * abs(M2) + 0.10 * C[, 4])
    tail_noise <- rnorm(n)
    p_mix <- plogis(-1.25 + 3.35 * (M3 > 0) + 0.45 * M3 + 0.10 * C[, 5])
    S_latent <- rbinom(n, size = 1L, prob = p_mix)
    shape_component <- ifelse(
      S_latent == 1L,
      rnorm(n, mean = 1.45, sd = 0.22),
      rnorm(n, mean = -0.95, sd = 0.28)
    )
    Y <- mu + shape_component + sigma * tail_noise

    return(
      base_df |>
        transmute(
          row_id = row_id,
          X = as.integer(.env$X),
          !!!as.data.frame(C),
          M1 = M1,
          M2 = M2,
          M3 = M3,
          S_latent = S_latent,
          Y = Y
        )
    )
  }

  if (identical(scenario, "IHDPMean")) {
    mu <- 1.00 + 0.55 * X + 0.45 * M1 + 0.12 * M2 + 0.10 * M3 +
      0.18 * C[, 1] - 0.10 * C[, 2] + 0.08 * C[, 7]
    sigma <- 0.70 + 0.12 * softplus_num(-0.25 + 0.15 * X + 0.35 * M2)
    p_mix <- plogis(-0.90 + 0.35 * X + 0.35 * M3 + 0.10 * C[, 3])
    shift <- 0.55
    sd_pos <- 0.80
    sd_neg <- 0.80
  } else {
    if (identical(scenario, "IHDPDist") || identical(scenario, "IHDPDistStrong")) {
      dist_boost <- if (identical(scenario, "IHDPDistStrong")) 1.25 else 1.00
      mu <- 1.00 + 0.40 * X + 0.32 * M1 + 0.08 * M2 + 0.12 * M3 +
        0.16 * C[, 1] - 0.08 * C[, 2] + 0.06 * C[, 7]
      sigma <- 0.55 + 0.35 * softplus_num(-0.35 + 0.35 * X + 0.70 * M2 - 0.20 * C[, 4])
      p_mix <- plogis(-0.60 + dist_boost * (0.70 * X + 0.95 * M3) + 0.18 * C[, 5] - 0.12 * C[, 8])
      shift <- dist_boost * 1.20
      sd_pos <- dist_boost * 0.80
      sd_neg <- dist_boost * 0.80
    } else {
      m3_pos <- pmax(M3, 0)
      m3_neg <- pmax(-M3, 0)
      mu <- 0.90 + 0.28 * X + 0.28 * M1 + 0.05 * M2 + 0.06 * M3 +
        0.16 * C[, 1] - 0.08 * C[, 2] + 0.05 * C[, 7] +
        0.18 * X * M1 + 0.22 * m3_pos - 0.10 * m3_neg
      sigma <- 0.35 + 0.45 * softplus_num(
        -0.40 + 0.95 * abs(M2) + 0.55 * X * (M2 > 0) + 0.20 * C[, 4]
      )
      p_mix <- plogis(
        -1.10 + 1.20 * X + 1.45 * (M3 > 0) + 0.90 * X * (M3 > 0) +
          0.20 * C[, 5] - 0.15 * C[, 8]
      )
      shift <- 1.65
      sd_pos <- 0.55
      sd_neg <- 1.05
    }
  }

  S_latent <- rbinom(n, size = 1L, prob = p_mix)
  mix_noise <- ifelse(
    S_latent == 1L,
    rnorm(n, mean = shift, sd = sd_pos),
    rnorm(n, mean = -0.75 * shift, sd = sd_neg)
  )

  Y <- mu + sigma * mix_noise

  out <- base_df |>
    transmute(
      row_id = row_id,
      X = as.integer(.env$X),
      !!!as.data.frame(C),
      M1 = M1,
      M2 = M2,
      M3 = M3,
      S_latent = S_latent,
      Y = Y
    )

  out
}

draw_ihdp_semisynth_mediators_given_a <- function(
  base_df,
  a,
  B = 200L,
  scenario = c(
    "IHDPMean", "IHDPDist", "IHDPDistStrong",
    "IHDPMechanismStrong", "IHDPMechanismSafe",
    "IHDPMechanismSkew", "IHDPMechanismTail", "IHDPMechanismThreshold",
    "IHDPMechanismDisplay", "IHDPMechanismDisplayLearnable",
    "IHDPMechanismThreeChannel", "IHDPMechanismTwoChannelSimple",
    "IHDPMechanismThreeChannelDisplay", "IHDPMechanismThreeChannelSimple",
    "IHDPMechanismThreeChannelBernoulliSkew",
    "IHDPMechanismThreeChannelFinal",
    "IHDPMechanismThreeChannelOrthogonal",
    "IHDPMechanismThreeChannelTailShift",
    "IHDPMechanismThreeChannelSmoothTail"
  ),
  seed = NULL
) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)

  c_vars <- grep("^C", names(base_df), value = TRUE)
  stopifnot(length(c_vars) >= 9L)

  C <- as.matrix(base_df[, c_vars, drop = FALSE])
  n <- nrow(C)
  B <- as.integer(B)

  X_rep <- matrix(as.numeric(a), nrow = n, ncol = B)
  H1 <- matrix(rnorm(n * B), nrow = n, ncol = B)
  H2 <- matrix(rnorm(n * B), nrow = n, ncol = B)

  mechanism_safe_family <- c(
    "IHDPMechanismSafe",
    "IHDPMechanismSkew",
    "IHDPMechanismTail",
    "IHDPMechanismThreshold",
    "IHDPMechanismDisplay",
    "IHDPMechanismDisplayLearnable",
    "IHDPMechanismThreeChannel",
    "IHDPMechanismTwoChannelSimple",
    "IHDPMechanismThreeChannelDisplay",
    "IHDPMechanismThreeChannelSimple",
    "IHDPMechanismThreeChannelBernoulliSkew",
    "IHDPMechanismThreeChannelFinal",
    "IHDPMechanismThreeChannelOrthogonal",
    "IHDPMechanismThreeChannelTailShift",
    "IHDPMechanismThreeChannelSmoothTail"
  )

  if (scenario %in% c(
    "IHDPMechanismThreeChannelDisplay",
    "IHDPMechanismThreeChannelSimple",
    "IHDPMechanismThreeChannelBernoulliSkew",
    "IHDPMechanismThreeChannelFinal",
    "IHDPMechanismThreeChannelOrthogonal",
    "IHDPMechanismThreeChannelTailShift",
    "IHDPMechanismThreeChannelSmoothTail"
  )) {
    Sigma_M3 <- if (identical(scenario, "IHDPMechanismThreeChannelOrthogonal") ||
                    identical(scenario, "IHDPMechanismThreeChannelTailShift") ||
                    identical(scenario, "IHDPMechanismThreeChannelSmoothTail")) {
      diag(3L)
    } else {
      matrix(c(
        1.00, 0.12, 0.08,
        0.12, 1.00, 0.10,
        0.08, 0.10, 1.00
      ), nrow = 3L, byrow = TRUE)
    }
    E <- matrix(rnorm(n * B * 3L), nrow = n * B, ncol = 3L) %*% chol(Sigma_M3)
    E1 <- matrix(E[, 1], nrow = n, ncol = B)
    E2 <- matrix(E[, 2], nrow = n, ncol = B)
    E3 <- matrix(E[, 3], nrow = n, ncol = B)
    if (identical(scenario, "IHDPMechanismThreeChannelTailShift") ||
        identical(scenario, "IHDPMechanismThreeChannelSmoothTail")) {
      smooth_m1_x <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_M1_X", 0.90)
      smooth_m2_x <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_M2_X", 1.00)
      smooth_m3_x <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_M3_X", 2.40)
      smooth_m_noise <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_M_NOISE", 0.34)
      M1 <- 0.10 + smooth_m1_x * X_rep +
        0.20 * C[, 1] - 0.15 * C[, 2] + smooth_m_noise * E1
      M2 <- -0.15 + smooth_m2_x * X_rep +
        0.10 * C[, 4] + smooth_m_noise * E2
      m3_x <- if (identical(scenario, "IHDPMechanismThreeChannelSmoothTail")) smooth_m3_x else 2.60
      M3 <- -0.70 + m3_x * X_rep +
        0.08 * C[, 5] + smooth_m_noise * E3
    } else if (identical(scenario, "IHDPMechanismThreeChannelOrthogonal")) {
      M1 <- 0.15 + 0.95 * X_rep +
        0.22 * C[, 1] - 0.16 * C[, 2] + 0.36 * E1
      M2 <- -0.30 + 1.05 * X_rep +
        0.12 * C[, 4] + 0.34 * E2
      M3 <- -0.85 + 1.35 * X_rep +
        0.12 * C[, 5] + 0.34 * E3
    } else if (identical(scenario, "IHDPMechanismThreeChannelFinal")) {
      M1 <- 0.20 + 0.90 * X_rep +
        0.25 * C[, 1] - 0.18 * C[, 2] + 0.42 * E1
      M2 <- -0.45 + 1.05 * X_rep +
        0.16 * C[, 4] + 0.45 * E2
      M3 <- -0.70 + 1.25 * X_rep +
        0.18 * C[, 5] + 0.45 * E3
    } else {
      M1 <- 0.20 + 0.95 * X_rep +
        0.30 * C[, 1] - 0.20 * C[, 2] + 0.45 * E1
      M2 <- -0.60 + 1.25 * X_rep +
        0.20 * C[, 4] + 0.45 * E2
      M3 <- -1.00 + 1.65 * X_rep +
        0.25 * C[, 5] + 0.45 * E3
    }
    out <- array(NA_real_, dim = c(n, B, 3L))
    out[, , 1] <- M1
    out[, , 2] <- M2
    out[, , 3] <- M3
    dimnames(out) <- list(NULL, NULL, c("M1", "M2", "M3"))
    return(out)
  }

  if (identical(scenario, "IHDPMechanismTwoChannelSimple")) {
    Sigma_M2 <- matrix(c(1.00, 0.12, 0.12, 1.00), nrow = 2L, byrow = TRUE)
    E <- matrix(rnorm(n * B * 2L), nrow = n * B, ncol = 2L) %*% chol(Sigma_M2)
    E1 <- matrix(E[, 1], nrow = n, ncol = B)
    E2 <- matrix(E[, 2], nrow = n, ncol = B)
    M1 <- 0.20 + 1.00 * X_rep +
      0.30 * C[, 1] - 0.20 * C[, 2] + 0.45 * E1
    M2 <- -0.75 + 1.30 * X_rep +
      0.25 * C[, 5] + 0.45 * E2
    out <- array(NA_real_, dim = c(n, B, 2L))
    out[, , 1] <- M1
    out[, , 2] <- M2
    dimnames(out) <- list(NULL, NULL, c("M1", "M2"))
    return(out)
  }

  if (scenario %in% mechanism_safe_family) {
    Sigma_M <- matrix(c(
      1.00, 0.25, 0.18,
      0.25, 1.00, 0.22,
      0.18, 0.22, 1.00
    ), nrow = 3L, byrow = TRUE)
    E <- matrix(rnorm(n * B * 3L), nrow = n * B, ncol = 3L) %*% chol(Sigma_M)
    E1 <- matrix(E[, 1], nrow = n, ncol = B)
    E2 <- matrix(E[, 2], nrow = n, ncol = B)
    E3 <- matrix(E[, 3], nrow = n, ncol = B)

    m1_x <- if (identical(scenario, "IHDPMechanismDisplay")) {
      1.25
    } else if (identical(scenario, "IHDPMechanismDisplayLearnable")) {
      1.10
    } else if (identical(scenario, "IHDPMechanismThreeChannel")) {
      1.05
    } else {
      0.95
    }
    M1 <- 0.20 + m1_x * X_rep +
      0.30 * C[, 1] - 0.20 * C[, 2] + 0.45 * E1

    if (identical(scenario, "IHDPMechanismDisplayLearnable")) {
      M2 <- -0.80 + 1.15 * X_rep + 0.25 * C[, 5] + 0.45 * E2
      out <- array(NA_real_, dim = c(n, B, 2L))
      out[, , 1] <- M1
      out[, , 2] <- M2
      dimnames(out) <- list(NULL, NULL, c("M1", "M2"))
      return(out)
    }

    if (identical(scenario, "IHDPMechanismThreeChannel")) {
      M2 <- -0.45 + 1.05 * X_rep + 0.20 * C[, 4] + 0.45 * E2
      M3 <- -0.65 + 1.25 * X_rep + 0.25 * C[, 5] + 0.45 * E3
      out <- array(NA_real_, dim = c(n, B, 3L))
      out[, , 1] <- M1
      out[, , 2] <- M2
      out[, , 3] <- M3
      dimnames(out) <- list(NULL, NULL, c("M1", "M2", "M3"))
      return(out)
    }

    scale2 <- if (identical(scenario, "IHDPMechanismDisplay")) {
      1.45 - 0.95 * X_rep + 0.05 * softplus_num(C[, 4])
    } else {
      1.16 - 0.54 * X_rep + 0.06 * softplus_num(C[, 4])
    }
    M2 <- -0.05 + 0.04 * X_rep - 0.18 * C[, 3] + 0.18 * C[, 4] + scale2 * E2

    m3_x <- if (identical(scenario, "IHDPMechanismDisplay")) {
      1.65
    } else {
      1.25
    }
    M3 <- -1.20 + m3_x * X_rep +
        0.12 * C[, 5] - 0.10 * C[, 8] +
        0.12 * C[, 6] + 0.34 * E3
  } else if (identical(scenario, "IHDPMechanismStrong")) {
    M1 <- 0.20 + 1.35 * X_rep + 0.30 * C[, 1] - 0.20 * C[, 2] +
      0.25 * H1 + matrix(rnorm(n * B, sd = 0.45), nrow = n, ncol = B)

    scale2 <- 0.55 + 0.85 * X_rep + 0.12 * softplus_num(C[, 4])
    M2 <- -0.05 + 0.08 * X_rep - 0.18 * C[, 3] + 0.18 * C[, 4] +
      0.18 * H1 + scale2 * matrix(rnorm(n * B), nrow = n, ncol = B)

    p3 <- plogis(-1.35 + 2.65 * X_rep + 0.35 * C[, 5] - 0.20 * C[, 8])
    S3 <- matrix(rbinom(n * B, size = 1L, prob = as.vector(p3)), nrow = n, ncol = B)
    M3 <- ifelse(
      S3 == 1L,
      matrix(rnorm(n * B, mean = 1.20, sd = 0.35), nrow = n, ncol = B),
      matrix(rnorm(n * B, mean = -1.05, sd = 0.40), nrow = n, ncol = B)
    ) + 0.12 * C[, 6] + 0.10 * H1
  } else {
    e1 <- matrix(rnorm(n * B, sd = 0.45), nrow = n, ncol = B)
    e2 <- matrix(rnorm(n * B, sd = 0.50), nrow = n, ncol = B)
    e3 <- matrix(rnorm(n * B, sd = 0.45), nrow = n, ncol = B)

    M1 <- 0.20 + 0.90 * X_rep + 0.30 * C[, 1] - 0.20 * C[, 2] + 0.15 * C[, 7] + 0.55 * H1 + 0.15 * H2 + e1
    M2 <- -0.10 + 0.55 * X_rep - 0.15 * C[, 3] + 0.20 * C[, 4] + 0.20 * C[, 8] + 0.20 * H1 + 0.60 * H2 + e2
    M3 <- 0.15 + 0.35 * X_rep + 0.20 * C[, 5] - 0.20 * C[, 6] + 0.10 * C[, 9] + 0.40 * H1 - 0.30 * H2 + e3
  }

  out <- array(NA_real_, dim = c(n, B, 3L))
  out[, , 1] <- M1
  out[, , 2] <- M2
  out[, , 3] <- M3
  dimnames(out) <- list(NULL, NULL, c("M1", "M2", "M3"))
  out
}

draw_ihdp_semisynth_outcomes_given_M <- function(
  base_df,
  M_array,
  a,
  scenario = c(
    "IHDPMean", "IHDPDist", "IHDPDistStrong",
    "IHDPMechanismStrong", "IHDPMechanismSafe",
    "IHDPMechanismSkew", "IHDPMechanismTail", "IHDPMechanismThreshold",
    "IHDPMechanismDisplay", "IHDPMechanismDisplayLearnable",
    "IHDPMechanismThreeChannel", "IHDPMechanismTwoChannelSimple",
    "IHDPMechanismThreeChannelDisplay", "IHDPMechanismThreeChannelSimple",
    "IHDPMechanismThreeChannelBernoulliSkew",
    "IHDPMechanismThreeChannelFinal",
    "IHDPMechanismThreeChannelOrthogonal",
    "IHDPMechanismThreeChannelTailShift",
    "IHDPMechanismThreeChannelSmoothTail"
  ),
  seed = NULL
) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)

  c_vars <- grep("^C", names(base_df), value = TRUE)
  C <- as.matrix(base_df[, c_vars, drop = FALSE])
  n <- nrow(C)
  stopifnot(length(dim(M_array)) == 3L, dim(M_array)[1] == n)
  if (identical(scenario, "IHDPMechanismDisplayLearnable") ||
      identical(scenario, "IHDPMechanismTwoChannelSimple")) {
    stopifnot(dim(M_array)[3] == 2L)
  } else {
    stopifnot(dim(M_array)[3] == 3L)
  }
  B <- dim(M_array)[2]

  X_rep <- matrix(as.numeric(a), nrow = n, ncol = B)
  M1 <- M_array[, , 1]
  M2 <- M_array[, , 2]
  M3 <- if (dim(M_array)[3] >= 3L) M_array[, , 3] else NULL

  mechanism_safe_family <- c(
    "IHDPMechanismSafe",
    "IHDPMechanismSkew",
    "IHDPMechanismTail",
    "IHDPMechanismThreshold",
    "IHDPMechanismDisplay",
    "IHDPMechanismDisplayLearnable",
    "IHDPMechanismThreeChannel",
    "IHDPMechanismTwoChannelSimple",
    "IHDPMechanismThreeChannelDisplay",
    "IHDPMechanismThreeChannelSimple",
    "IHDPMechanismThreeChannelBernoulliSkew",
    "IHDPMechanismThreeChannelFinal",
    "IHDPMechanismThreeChannelOrthogonal",
    "IHDPMechanismThreeChannelTailShift",
    "IHDPMechanismThreeChannelSmoothTail"
  )

  if (scenario %in% c(
    "IHDPMechanismThreeChannelDisplay",
    "IHDPMechanismThreeChannelSimple",
    "IHDPMechanismThreeChannelBernoulliSkew",
    "IHDPMechanismThreeChannelFinal",
    "IHDPMechanismThreeChannelOrthogonal",
    "IHDPMechanismThreeChannelTailShift",
    "IHDPMechanismThreeChannelSmoothTail"
  )) {
    if (identical(scenario, "IHDPMechanismThreeChannelSmoothTail")) {
      mu <- 1.00 + 0.98 * M1 + 0.08 * C[, 1] - 0.06 * C[, 2]
      sigma <- exp(ihdp_env_num("IHDP_S2_SMOOTHTAIL_M2_SIGMA_COEF", 0.22) * M2)
      smooth_family <- Sys.getenv("IHDP_S2_SMOOTHTAIL_FAMILY", "hermite")
      gamma_max <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_GAMMA_MAX", 0.80)
      shape_slope <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_SHAPE_SLOPE", 3.00)
      tail_center <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_CENTER", 0.50)
      gamma <- gamma_max * (2 * plogis(shape_slope * (M3 - tail_center)) - 1)
      if (identical(smooth_family, "additive")) {
        tail_lambda <- ihdp_env_num("IHDP_S2_SMOOTHTAIL_TAIL_LAMBDA", 0.70)
        z_scale <- matrix(rnorm(n * B), nrow = n, ncol = B)
        z_tail <- matrix(rnorm(n * B), nrow = n, ncol = B)
        q_tail <- (z_tail^2 - 1) / sqrt(2)
        return(mu + sigma * z_scale + tail_lambda * gamma * q_tail)
      }
      if (identical(smooth_family, "covskew")) {
        shape_lambda <- ihdp_env_num("IHDP_S2_COVSKEW_LAMBDA", 0.55)
        h_shape <- as.numeric(scale(C[, 5]^2))
        z_tail <- matrix(rnorm(n * B), nrow = n, ncol = B)
        return(mu + sigma * z_tail + shape_lambda * gamma * h_shape)
      }
      if (identical(smooth_family, "bimodal")) {
        noise_sd <- ihdp_env_num("IHDP_S2_BIMODAL_NOISE_SD", 0.42)
        sep <- exp(gamma)
        signs <- ifelse(matrix(runif(n * B), nrow = n, ncol = B) < 0.5, -1, 1)
        raw_tail <- signs * sep + matrix(rnorm(n * B, sd = noise_sd), nrow = n, ncol = B)
        z_tail <- raw_tail / sqrt(sep^2 + noise_sd^2)
        return(mu + sigma * z_tail)
      }
      if (identical(smooth_family, "s1bimodal")) {
        noise_sd <- ihdp_env_num("IHDP_S2_S1BIMODAL_NOISE_SD", 0.35)
        sep_min <- ihdp_env_num("IHDP_S2_S1BIMODAL_SEP_MIN", 0.10)
        sep_amp <- ihdp_env_num("IHDP_S2_S1BIMODAL_SEP_AMP", 1.35)
        sep_slope <- ihdp_env_num("IHDP_S2_S1BIMODAL_SLOPE", 2.80)
        sep_center <- ihdp_env_num("IHDP_S2_S1BIMODAL_CENTER", 0.50)
        sep <- sep_min + sep_amp * plogis(sep_slope * (M3 - sep_center))
        signs <- ifelse(matrix(runif(n * B), nrow = n, ncol = B) < 0.5, -1, 1)
        raw_tail <- signs * sep + matrix(rnorm(n * B, sd = noise_sd), nrow = n, ncol = B)
        z_tail <- raw_tail / sqrt(mean(sep^2 + noise_sd^2))
        return(mu + sigma * z_tail)
      }
      if (identical(smooth_family, "covbimodal")) {
        split_var <- C[, 5]
        split_group <- ifelse(split_var >= median(split_var), 1, -1)
        split_group <- split_group - mean(split_group)
        split_group <- split_group / sd(split_group)
        amp_max <- ihdp_env_num("IHDP_S2_COVBIMODAL_AMP_MAX", 1.35)
        amp_slope <- ihdp_env_num("IHDP_S2_COVBIMODAL_SLOPE", 3.00)
        amp_center <- ihdp_env_num("IHDP_S2_COVBIMODAL_CENTER", 0.40)
        normalize_covbimodal <- tolower(Sys.getenv("IHDP_S2_COVBIMODAL_NORMALIZE", "true")) %in% c("1", "true", "yes")
        amp <- amp_max * plogis(amp_slope * (M3 - amp_center))
        raw_tail <- matrix(rnorm(n * B), nrow = n, ncol = B) + amp * split_group
        z_tail <- if (normalize_covbimodal) raw_tail / sqrt(pmax(1 + amp^2, 1e-8)) else raw_tail
        return(mu + sigma * z_tail)
      }
      if (identical(smooth_family, "covtail")) {
        split_var <- C[, 5]
        tail_q <- ihdp_env_num("IHDP_S2_COVTAIL_Q", 0.75)
        tail_cut <- as.numeric(quantile(split_var, tail_q, type = 8, names = FALSE))
        tail_slope <- ihdp_env_num("IHDP_S2_COVTAIL_C_SLOPE", 5.00)
        tail_gate <- plogis(tail_slope * (split_var - tail_cut))
        tail_group <- tail_gate - mean(tail_gate)
        scale_tail_group <- tolower(Sys.getenv("IHDP_S2_COVTAIL_SCALE_GROUP", "false")) %in% c("1", "true", "yes")
        if (scale_tail_group) tail_group <- tail_group / sd(tail_group)
        amp_max <- ihdp_env_num("IHDP_S2_COVTAIL_AMP_MAX", 2.80)
        amp_slope <- ihdp_env_num("IHDP_S2_COVTAIL_SLOPE", 5.00)
        amp_center <- ihdp_env_num("IHDP_S2_COVTAIL_CENTER", -0.50)
        normalize_covtail <- tolower(Sys.getenv("IHDP_S2_COVTAIL_NORMALIZE", "false")) %in% c("1", "true", "yes")
        amp <- amp_max * plogis(amp_slope * (M3 - amp_center))
        raw_tail <- matrix(rnorm(n * B), nrow = n, ncol = B) + amp * tail_group
        z_tail <- if (normalize_covtail) raw_tail / sqrt(pmax(1 + amp^2 * stats::var(tail_group), 1e-8)) else raw_tail
        return(mu + sigma * z_tail)
      }
      if (identical(smooth_family, "coretail")) {
        p0 <- ihdp_env_num("IHDP_S2_CORETAIL_P0", 0.15)
        core_sd <- ihdp_env_num("IHDP_S2_CORETAIL_CORE_SD", 0.55)
        tail_sd <- ihdp_env_num("IHDP_S2_CORETAIL_TAIL_SD", 2.30)
        p_tail <- plogis(qlogis(p0) + gamma)
        is_tail <- matrix(runif(n * B), nrow = n, ncol = B) < p_tail
        raw_tail <- ifelse(
          is_tail,
          matrix(rnorm(n * B, sd = tail_sd), nrow = n, ncol = B),
          matrix(rnorm(n * B, sd = core_sd), nrow = n, ncol = B)
        )
        tail_var <- p_tail * tail_sd^2 + (1 - p_tail) * core_sd^2
        z_tail <- raw_tail / sqrt(pmax(tail_var, 1e-8))
        return(mu + sigma * z_tail)
      }
      if (identical(smooth_family, "twopiece")) {
        z <- matrix(rnorm(n * B), nrow = n, ncol = B)
        right_scale <- exp(gamma)
        left_scale <- exp(-gamma)
        raw_tail <- ifelse(z >= 0, right_scale * z, left_scale * z)
        tail_mean <- (right_scale - left_scale) / sqrt(2 * pi)
        tail_second <- 0.5 * (right_scale^2 + left_scale^2)
        z_tail <- (raw_tail - tail_mean) / sqrt(pmax(tail_second - tail_mean^2, 1e-8))
        return(mu + sigma * z_tail)
      }
      if (identical(smooth_family, "skewnormal")) {
        gamma_vec <- as.numeric(gamma)
        if (length(gamma_vec) == n) {
          gamma_vec <- rep(gamma_vec, times = B)
        }
        z_tail <- matrix(standardized_skew_normal(gamma_vec), nrow = n, ncol = B)
      } else if (identical(smooth_family, "hinge")) {
        z <- matrix(rnorm(n * B), nrow = n, ncol = B)
        hinge_mu <- 1 / sqrt(2 * pi)
        hinge_sd <- sqrt(0.5 - hinge_mu^2)
        hinge_cov <- 0.5 / hinge_sd
        h <- (pmax(z, 0) - hinge_mu) / hinge_sd
        z_tail <- (z + gamma * h) / sqrt(pmax(1 + gamma^2 + 2 * gamma * hinge_cov, 1e-8))
      } else {
        z <- matrix(rnorm(n * B), nrow = n, ncol = B)
        z_tail <- (z + gamma * ((z^2 - 1) / sqrt(2))) / sqrt(1 + gamma^2)
      }
      return(mu + sigma * z_tail)
    }
    if (identical(scenario, "IHDPMechanismThreeChannelTailShift")) {
      mu <- 1.00 + 0.98 * M1 + 0.08 * C[, 1] - 0.06 * C[, 2]
      sigma <- exp(ihdp_env_num("IHDP_S2_TAILSHIFT_M2_SIGMA_COEF", 0.22) * M2)
      q_tail <- ihdp_env_num("IHDP_S2_TAILSHIFT_Q_TAIL", 0.06)
      shape_slope <- ihdp_env_num("IHDP_S2_TAILSHIFT_SHAPE_SLOPE", 9.00)
      core_sd <- ihdp_env_num("IHDP_S2_TAILSHIFT_CORE_SD", 0.45)
      tail_loc <- ihdp_env_num("IHDP_S2_TAILSHIFT_TAIL_LOC", 6.20)
      tail_sd <- ihdp_env_num("IHDP_S2_TAILSHIFT_TAIL_SD", 0.25)
      tail_center <- ihdp_env_num("IHDP_S2_TAILSHIFT_CENTER", 0.60)
      right_prob <- plogis(shape_slope * (M3 - tail_center))
      p_right <- q_tail * right_prob
      p_left <- q_tail * (1 - right_prob)
      u_tail <- matrix(runif(n * B), nrow = n, ncol = B)
      raw_tail <- matrix(rnorm(n * B, mean = 0, sd = core_sd), nrow = n, ncol = B)
      is_right <- u_tail < p_right
      is_left <- u_tail >= p_right & u_tail < p_right + p_left
      raw_tail[is_right] <- rnorm(sum(is_right), mean = tail_loc, sd = tail_sd)
      raw_tail[is_left] <- rnorm(sum(is_left), mean = -tail_loc, sd = tail_sd)
      mean_tail <- (p_right - p_left) * tail_loc
      second_tail <- (1 - q_tail) * core_sd^2 +
        p_right * (tail_sd^2 + tail_loc^2) +
        p_left * (tail_sd^2 + tail_loc^2)
      z_tail <- (raw_tail - mean_tail) / sqrt(pmax(second_tail - mean_tail^2, 1e-8))
      return(mu + sigma * z_tail)
    }
    if (identical(scenario, "IHDPMechanismThreeChannelOrthogonal")) {
      mu <- 1.00 + 0.05 * X_rep + 0.90 * M1 + 0.10 * C[, 1] - 0.08 * C[, 2]
      sigma <- exp(ihdp_env_num("IHDP_S2_ORTHO_M2_SIGMA_COEF", 0.18) * M2)
      tail_base <- ihdp_env_num("IHDP_S2_ORTHO_M3_TAIL_BASE", 0.04)
      tail_amp <- ihdp_env_num("IHDP_S2_ORTHO_M3_TAIL_AMP", 1.35)
      tail_slope <- ihdp_env_num("IHDP_S2_ORTHO_M3_TAIL_SLOPE", 2.35)
      tail_strength <- tail_base + tail_amp * plogis(tail_slope * M3)
      tail_noise <- matrix(right_hinge_skew_noise(n * B), nrow = n, ncol = B)
      tail_noise <- tail_noise - mean(tail_noise)
      return(mu + sigma * matrix(rnorm(n * B), nrow = n, ncol = B) +
               tail_strength * tail_noise)
    }
    if (identical(scenario, "IHDPMechanismThreeChannelFinal")) {
      mu <- 1.00 + 0.06 * X_rep + 0.85 * M1 + 0.12 * C[, 1] - 0.08 * C[, 2]
      sigma <- exp(ihdp_env_num("IHDP_S2_FINAL_M2_SIGMA_COEF", 0.22) * M2)
      tail_base <- ihdp_env_num("IHDP_S2_FINAL_M3_TAIL_BASE", 0.12)
      tail_amp <- ihdp_env_num("IHDP_S2_FINAL_M3_TAIL_AMP", 1.00)
      tail_slope <- ihdp_env_num("IHDP_S2_FINAL_M3_TAIL_SLOPE", 1.80)
      tail_strength <- tail_base + tail_amp * plogis(tail_slope * M3)
      return(mu + sigma * matrix(rnorm(n * B), nrow = n, ncol = B) +
               tail_strength * matrix(right_hinge_skew_noise(n * B), nrow = n, ncol = B))
    }
    mu <- 1.00 + 0.08 * X_rep + 0.85 * M1 + 0.15 * C[, 1] - 0.10 * C[, 2]
    sigma <- exp(ihdp_env_num("IHDP_S2_M2_SIGMA_COEF", 0.34) * M2)
    if (identical(scenario, "IHDPMechanismThreeChannelSimple")) {
      z_scale <- matrix(rnorm(n * B), nrow = n, ncol = B)
      skew_base <- ihdp_env_num("IHDP_S2_M3_HINGE_BASE", 0.10)
      skew_amp <- ihdp_env_num("IHDP_S2_M3_HINGE_AMP", 0.90)
      skew_slope <- ihdp_env_num("IHDP_S2_M3_HINGE_SLOPE", 2.25)
      skew_strength <- skew_base + skew_amp * plogis(skew_slope * M3)
      skew_component <- skew_strength *
        matrix(right_hinge_skew_noise(n * B), nrow = n, ncol = B)
      return(mu + sigma * z_scale + skew_component)
    } else if (identical(scenario, "IHDPMechanismThreeChannelBernoulliSkew")) {
      z_scale <- matrix(rnorm(n * B), nrow = n, ncol = B)
      p_base <- ihdp_env_num("IHDP_S2_M3_TAIL_BASE", 0.05)
      p_amp <- ihdp_env_num("IHDP_S2_M3_TAIL_AMP", 0.45)
      p_slope <- ihdp_env_num("IHDP_S2_M3_TAIL_SLOPE", 2.00)
      tail_size <- ihdp_env_num("IHDP_S2_M3_TAIL_SIZE", 2.00)
      p_tail <- p_base + p_amp * plogis(p_slope * M3)
      tail_event <- matrix(
        stats::rbinom(n * B, size = 1L, prob = as.vector(p_tail)),
        nrow = n,
        ncol = B
      )
      skew_component <- tail_size * (tail_event - p_tail)
      return(mu + sigma * z_scale + skew_component)
    }
    tail_strength <- 0.10 + 1.25 * plogis(2.80 * M3)
    tail_noise <- matrix(capped_exp_centered(n * B, cap = 2.25), nrow = n, ncol = B)
    return(mu + sigma * matrix(rnorm(n * B), nrow = n, ncol = B) + tail_strength * tail_noise)
  }

  if (identical(scenario, "IHDPMechanismTwoChannelSimple")) {
    mu <- 1.00 + 0.08 * X_rep + 0.95 * M1 + 0.15 * C[, 1] - 0.10 * C[, 2]
    sigma <- exp(0.22 * M2)
    tail_strength <- 0.15 + 1.05 * plogis(2.20 * M2)
    z_skew <- matrix(rnorm(n * B), nrow = n, ncol = B)
    q_skew <- (z_skew^2 - 1) / sqrt(2)
    return(mu + sigma * matrix(rnorm(n * B), nrow = n, ncol = B) + tail_strength * q_skew)
  }

  if (identical(scenario, "IHDPMechanismDisplayLearnable")) {
    mu <- 1.00 + 0.12 * X_rep + 1.05 * M1 + 0.15 * C[, 1] - 0.10 * C[, 2]
    sigma <- 0.90
    skew_strength <- pmax(0.20, 1.05 + 0.75 * M2)
    skew_component <- skew_strength * matrix(capped_exp_centered(n * B, cap = 2.00), nrow = n, ncol = B)
    return(mu + sigma * matrix(rnorm(n * B), nrow = n, ncol = B) + skew_component)
  }

  if (identical(scenario, "IHDPMechanismThreeChannel")) {
    mu <- 1.00 + 0.10 * X_rep + 1.00 * M1 + 0.15 * C[, 1] - 0.10 * C[, 2]
    sigma <- exp(0.35 * M2)
    gamma_skew <- 0.85 * M3
    z_skew <- matrix(rnorm(n * B), nrow = n, ncol = B)
    q_skew <- (z_skew^2 - 1) / sqrt(2)
    skew_component <- gamma_skew * q_skew
    return(mu + sigma * matrix(rnorm(n * B), nrow = n, ncol = B) + skew_component)
  }

  if (scenario %in% mechanism_safe_family) {
    mu <- 1.00 + 0.12 * X_rep +
      ifelse(identical(scenario, "IHDPMechanismDisplay"), 0.95, 0.75) * M1 +
      0.15 * C[, 1] - 0.10 * C[, 2]
    sigma <- if (identical(scenario, "IHDPMechanismDisplay")) {
      0.18 + 0.95 * softplus_num(-0.55 + 1.80 * abs(M2) + 0.04 * C[, 4])
    } else {
      0.22 + 0.70 * softplus_num(-0.70 + 1.35 * abs(M2) + 0.05 * C[, 4])
    }

    if (identical(scenario, "IHDPMechanismSafe")) {
      tail_noise <- matrix(rnorm(n * B), nrow = n, ncol = B)
      p_mix <- plogis(0.35 - 2.60 * abs(M3) + 0.04 * C[, 5])
      S_latent <- matrix(stats::rbinom(n * B, size = 1L, prob = as.vector(p_mix)), nrow = n, ncol = B)
      shape_component <- ifelse(
        S_latent == 1L,
        matrix(rnorm(n * B, mean = 1.55, sd = 0.18), nrow = n, ncol = B),
        matrix(rnorm(n * B, mean = -1.10, sd = 0.22), nrow = n, ncol = B)
      )
      return(mu + shape_component + sigma * tail_noise)
    }

    if (identical(scenario, "IHDPMechanismSkew")) {
      p_skew <- plogis(-0.45 + 1.35 * M3 + 0.05 * C[, 5])
      S_latent <- matrix(stats::rbinom(n * B, size = 1L, prob = as.vector(p_skew)), nrow = n, ncol = B)
      skew_draw <- stats::rexp(n * B, rate = 1) - 1
      skew_noise <- matrix(ifelse(as.vector(S_latent) == 1L, skew_draw, -skew_draw), nrow = n, ncol = B)
      return(mu + sigma * matrix(rnorm(n * B), nrow = n, ncol = B) + 0.55 * skew_noise)
    }

    if (identical(scenario, "IHDPMechanismTail")) {
      p_tail <- plogis(-0.85 + 1.45 * M3 + 0.05 * C[, 5])
      S_latent <- matrix(stats::rbinom(n * B, size = 1L, prob = as.vector(p_tail)), nrow = n, ncol = B)
      tail_component <- matrix(
        ifelse(
          as.vector(S_latent) == 1L,
          stats::rt(n * B, df = 3) / sqrt(3),
          rnorm(n * B, sd = 0.55)
        ),
        nrow = n,
        ncol = B
      )
      return(mu + sigma * matrix(rnorm(n * B), nrow = n, ncol = B) + 0.85 * tail_component)
    }

    p_jump <- if (identical(scenario, "IHDPMechanismDisplay")) {
      plogis(-0.30 + 2.05 * M3 + 0.05 * C[, 5])
    } else {
      plogis(-0.25 + 1.55 * M3 + 0.05 * C[, 5])
    }
    S_latent <- matrix(stats::rbinom(n * B, size = 1L, prob = as.vector(p_jump)), nrow = n, ncol = B)
    threshold_component <- if (identical(scenario, "IHDPMechanismDisplay")) {
      matrix(
        ifelse(
          as.vector(S_latent) == 1L,
          rnorm(n * B, mean = 1.20, sd = 0.20),
          rnorm(n * B, mean = -0.40, sd = 0.20)
        ),
        nrow = n,
        ncol = B
      )
    } else {
      matrix(
        ifelse(
          as.vector(S_latent) == 1L,
          rnorm(n * B, mean = 0.85, sd = 0.22),
          rnorm(n * B, mean = -0.25, sd = 0.22)
        ),
        nrow = n,
        ncol = B
      )
    }
    return(mu + sigma * matrix(rnorm(n * B), nrow = n, ncol = B) + threshold_component)
  }

  if (identical(scenario, "IHDPMechanismStrong")) {
    mu <- 1.20 + 0.20 * X_rep + 0.95 * M1 +
      0.18 * C[, 1] - 0.12 * C[, 2] + 0.08 * C[, 7]
    sigma <- 0.22 + 0.88 * softplus_num(-1.10 + 1.85 * abs(M2) + 0.10 * C[, 4])
    tail_noise <- matrix(rnorm(n * B), nrow = n, ncol = B)
    p_mix <- plogis(-1.25 + 3.35 * (M3 > 0) + 0.45 * M3 + 0.10 * C[, 5])
    S_latent <- matrix(stats::rbinom(n * B, size = 1L, prob = as.vector(p_mix)), nrow = n, ncol = B)
    shape_component <- ifelse(
      S_latent == 1L,
      matrix(rnorm(n * B, mean = 1.45, sd = 0.22), nrow = n, ncol = B),
      matrix(rnorm(n * B, mean = -0.95, sd = 0.28), nrow = n, ncol = B)
    )
    return(mu + shape_component + sigma * tail_noise)
  }

  if (identical(scenario, "IHDPMean")) {
    mu <- 1.00 + 0.55 * X_rep + 0.45 * M1 + 0.12 * M2 + 0.10 * M3 +
      0.18 * C[, 1] - 0.10 * C[, 2] + 0.08 * C[, 7]
    sigma <- 0.70 + 0.12 * softplus_num(-0.25 + 0.15 * X_rep + 0.35 * M2)
    p_mix <- plogis(-0.90 + 0.35 * X_rep + 0.35 * M3 + 0.10 * C[, 3])
    shift <- 0.55
    sd_pos <- 0.80
    sd_neg <- 0.80
  } else {
    if (identical(scenario, "IHDPDist") || identical(scenario, "IHDPDistStrong")) {
      dist_boost <- if (identical(scenario, "IHDPDistStrong")) 1.25 else 1.00
      mu <- 1.00 + 0.40 * X_rep + 0.32 * M1 + 0.08 * M2 + 0.12 * M3 +
        0.16 * C[, 1] - 0.08 * C[, 2] + 0.06 * C[, 7]
      sigma <- 0.55 + 0.35 * softplus_num(-0.35 + 0.35 * X_rep + 0.70 * M2 - 0.20 * C[, 4])
      p_mix <- plogis(-0.60 + dist_boost * (0.70 * X_rep + 0.95 * M3) + 0.18 * C[, 5] - 0.12 * C[, 8])
      shift <- dist_boost * 1.20
      sd_pos <- dist_boost * 0.80
      sd_neg <- dist_boost * 0.80
    } else {
      m3_pos <- pmax(M3, 0)
      m3_neg <- pmax(-M3, 0)
      mu <- 0.90 + 0.28 * X_rep + 0.28 * M1 + 0.05 * M2 + 0.06 * M3 +
        0.16 * C[, 1] - 0.08 * C[, 2] + 0.05 * C[, 7] +
        0.18 * X_rep * M1 + 0.22 * m3_pos - 0.10 * m3_neg
      sigma <- 0.35 + 0.45 * softplus_num(
        -0.40 + 0.95 * abs(M2) + 0.55 * X_rep * (M2 > 0) + 0.20 * C[, 4]
      )
      p_mix <- plogis(
        -1.10 + 1.20 * X_rep + 1.45 * (M3 > 0) + 0.90 * X_rep * (M3 > 0) +
          0.20 * C[, 5] - 0.15 * C[, 8]
      )
      shift <- 1.65
      sd_pos <- 0.55
      sd_neg <- 1.05
    }
  }

  S_latent <- matrix(stats::rbinom(n * B, size = 1L, prob = as.vector(p_mix)), nrow = n, ncol = B)
  mix_noise <- matrix(
    ifelse(
      as.vector(S_latent) == 1L,
      rnorm(n * B, mean = shift, sd = sd_pos),
      rnorm(n * B, mean = -0.75 * shift, sd = sd_neg)
    ),
    nrow = n,
    ncol = B
  )

  mu + sigma * mix_noise
}

build_ihdp_pathway_M_pair <- function(M0, M1, s) {
  stopifnot(length(dim(M0)) == 3L, length(dim(M1)) == 3L, all(dim(M0) == dim(M1)))
  n0 <- dim(M0)[1]
  B0 <- dim(M0)[2]
  S0 <- dim(M0)[3]
  stopifnot(s >= 1L, s <= S0)

  M0s <- array(NA_real_, dim = c(n0, B0, S0))
  M1s <- array(NA_real_, dim = c(n0, B0, S0))
  idx_low <- if (s > 1L) 1:(s - 1L) else integer(0)
  idx_high <- if (s < S0) (s + 1L):S0 else integer(0)

  for (i in seq_len(n0)) {
    M0_i <- matrix(M0[i, , , drop = FALSE], nrow = B0, ncol = S0)
    M1_i <- matrix(M1[i, , , drop = FALSE], nrow = B0, ncol = S0)

    parts1 <- list()
    if (length(idx_low)) parts1[[length(parts1) + 1L]] <- M0_i[sample(B0), idx_low, drop = FALSE]
    parts1[[length(parts1) + 1L]] <- M1_i[sample(B0), s, drop = FALSE]
    if (length(idx_high)) parts1[[length(parts1) + 1L]] <- M1_i[sample(B0), idx_high, drop = FALSE]
    M1s[i, , ] <- do.call(cbind, parts1)

    parts0 <- list()
    if (length(idx_low)) parts0[[length(parts0) + 1L]] <- M0_i[sample(B0), idx_low, drop = FALSE]
    parts0[[length(parts0) + 1L]] <- M0_i[sample(B0), s, drop = FALSE]
    if (length(idx_high)) parts0[[length(parts0) + 1L]] <- M1_i[sample(B0), idx_high, drop = FALSE]
    M0s[i, , ] <- do.call(cbind, parts0)
  }

  list(M0s = M0s, M1s = M1s)
}

ihdp_semisynth_truth_interventions <- function(
  base_df,
  scenario = c(
    "IHDPMean", "IHDPDist", "IHDPDistStrong",
    "IHDPMechanismStrong", "IHDPMechanismSafe",
    "IHDPMechanismSkew", "IHDPMechanismTail", "IHDPMechanismThreshold",
    "IHDPMechanismDisplay", "IHDPMechanismDisplayLearnable",
    "IHDPMechanismThreeChannel", "IHDPMechanismTwoChannelSimple",
    "IHDPMechanismThreeChannelDisplay", "IHDPMechanismThreeChannelSimple",
    "IHDPMechanismThreeChannelBernoulliSkew",
    "IHDPMechanismThreeChannelFinal",
    "IHDPMechanismThreeChannelOrthogonal",
    "IHDPMechanismThreeChannelTailShift",
    "IHDPMechanismThreeChannelSmoothTail"
  ),
  B = 2000L,
  seed = 123L,
  a0 = 0,
  a1 = 1,
  a_ipse = 1
) {
  scenario <- match.arg(scenario)
  B <- as.integer(B)
  seed <- as.integer(seed)

  M0 <- draw_ihdp_semisynth_mediators_given_a(base_df, a = a0, B = B, scenario = scenario, seed = seed + 11L)
  M1 <- draw_ihdp_semisynth_mediators_given_a(base_df, a = a1, B = B, scenario = scenario, seed = seed + 13L)

  Y_0M0 <- draw_ihdp_semisynth_outcomes_given_M(base_df, M0, a = a0, scenario = scenario, seed = seed + 21L)
  Y_1M1 <- draw_ihdp_semisynth_outcomes_given_M(base_df, M1, a = a1, scenario = scenario, seed = seed + 23L)
  Y_1M0 <- draw_ihdp_semisynth_outcomes_given_M(base_df, M0, a = a1, scenario = scenario, seed = seed + 25L)
  Y_0M1 <- draw_ihdp_semisynth_outcomes_given_M(base_df, M1, a = a0, scenario = scenario, seed = seed + 27L)

  S <- dim(M0)[3]
  m_vars <- dimnames(M0)[[3]]
  if (is.null(m_vars)) m_vars <- paste0("M", seq_len(S))

  path_outcomes <- vector("list", S)
  for (s in seq_len(S)) {
    pair_s <- build_ihdp_pathway_M_pair(M0, M1, s)
    path_outcomes[[s]] <- list(
      Y_a_M0s = draw_ihdp_semisynth_outcomes_given_M(base_df, pair_s$M0s, a = a_ipse, scenario = scenario, seed = seed + 100L + 2L * s),
      Y_a_M1s = draw_ihdp_semisynth_outcomes_given_M(base_df, pair_s$M1s, a = a_ipse, scenario = scenario, seed = seed + 101L + 2L * s)
    )
  }
  names(path_outcomes) <- m_vars

  list(
    B = B,
    a0 = a0,
    a1 = a1,
    a_ipse = a_ipse,
    m_vars = m_vars,
    outcomes = list(
      Y_1M1 = Y_1M1,
      Y_0M0 = Y_0M0,
      Y_1M0 = Y_1M0,
      Y_0M1 = Y_0M1,
      path = path_outcomes
    )
  )
}
