# Generate simulation data for DCMA experiments
#
# Generates synthetic datasets under two scenarios used in DCMA experiments.
#
# Scenario S1 (single mediator, mixture outcome when X=1):
#   - One mediator M1 depends on treatment X and covariate Z.
#   - When X=1, Y follows a two-component mixture with component chosen by sign(Z).
#   - When X=0, Y is unimodal.
#
# Scenario S2 (multiple correlated mediators, nonlinear outcome):
#   - S mediators are generated with correlated Gaussian noise (AR(1)-like correlation).
#   - Y depends on X, Z, sum of mediators, and a nonlinear interaction sin(M1 * M2).
#
# Args:
#   n: sample size
#   scenario: "S1" or "S2"
#   S: number of mediators for S2 (must be >= 2)
#   s1_mode_shift: half-distance between the two treated modes in S1
#   seed: optional random seed
#
# Returns:
#   A data.frame with columns:
#     - X: binary treatment
#     - Z: continuous covariate
#     - M1...MS: mediator(s)
#     - Y: outcome
#     - S_latent (S1 only): indicator for whether the observation is in the mixture regime (X=1)


simulate_mediation_data <- function(
  n,
  scenario = c("S1", "S2", "S3"),
  S = 5L,
  s1_mode_shift = 2,
  s2_noise = c("gaussian", "heavyY", "heavyMY"),
  s2_df_y = 3,
  s2_df_m = 5,
  s2_scale_y = 1,
  s2_scale_m = 1,
  s2_mediator_scale = 0.5,
  s2_outcome_scale = 0.5,
  s2_standardize_y = TRUE,
  s2_standardize_m = TRUE,
  seed = NULL
) {
  scenario <- match.arg(scenario)
  s2_noise <- match.arg(s2_noise)
  if (!is.numeric(s1_mode_shift) || length(s1_mode_shift) != 1L || !is.finite(s1_mode_shift) || s1_mode_shift < 0) {
    stop("s1_mode_shift must be one finite non-negative number.")
  }
  if (!is.null(seed)) set.seed(seed)

  n <- as.integer(n)
  if (n <= 0L) stop("n must be a positive integer")

  Z <- rnorm(n, mean = 0, sd = 1)
  X <- rbinom(n, size = 1, prob = 0.5)

  if (scenario == "S1") {
    ## -----------------------------
    ## S1: single mediator + mixture outcome when X=1
    ## -----------------------------
    eM <- rnorm(n, mean = 0, sd = 0.5)
    M1 <- 0.5 + 1.0 * X + 0.3 * Z + eM

    ## Baseline mean: ensures intercept 4 at M1=0 when X=0, Z=0
    mu_base <- 4 + 0.3 * X + 0.2 * Z + 0.5 * M1

    ## Latent component indicators
    S_x1 <- as.integer(X == 1)          # whether in "mixture regime"
    S_z  <- as.integer(Z < 0)           # selects left/right mode when X=1

    ## Noise scales
    sd0 <- 1
    sd1 <- 1
    sd2 <- 1

    Y <- numeric(n)

    idx0 <- which(S_x1 == 1L & S_z == 0L)
    if (length(idx0) > 0L) {
      Y[idx0] <- (mu_base[idx0] - s1_mode_shift) + rnorm(length(idx0), 0, sd0)
    }

    idx1 <- which(S_x1 == 1L & S_z == 1L)
    if (length(idx1) > 0L) {
      Y[idx1] <- (mu_base[idx1] + s1_mode_shift) + rnorm(length(idx1), 0, sd1)
    }

    idx2 <- which(S_x1 == 0L)
    if (length(idx2) > 0L) {
      Y[idx2] <- (mu_base[idx2] + 0.3) + rnorm(length(idx2), 0, sd2)
    }

    dat <- data.frame(
      X = X,
      Z = Z,
      M1 = M1,
      S_latent = S_x1,   # keep your original output name, but now clearly defined
      Y = Y
    )

  } else if (scenario == "S2") {
    ## -----------------------------
    ## S2: multiple mediators with correlated noise + nonlinear outcome
    ## -----------------------------
    S <- as.integer(S)
    if (S < 2L) stop("For scenario S2, S must be >= 2 (needs M1 and M2 for sin(M1*M2)).")

    ## Coefficients (truncate/recycle safely)
    b0 <- rep(0.5, S)

    if (!is.numeric(s2_mediator_scale) || length(s2_mediator_scale) != 1L || !is.finite(s2_mediator_scale) || s2_mediator_scale <= 0) {
      stop("s2_mediator_scale must be one positive finite number.")
    }
    if (!is.numeric(s2_outcome_scale) || length(s2_outcome_scale) != 1L || !is.finite(s2_outcome_scale) || s2_outcome_scale <= 0) {
      stop("s2_outcome_scale must be one positive finite number.")
    }

    bA_full <- c(1.0, 0.8, 0.6, 0.4, 0.2) * s2_mediator_scale
    bZ_full <- c(0.3, 0.3, 0.2, 0.2, 0.1) * s2_mediator_scale
    bA <- rep(bA_full, length.out = S)
    bZ <- rep(bZ_full, length.out = S)

    ## Correlation Σ(i,j) = rho^{|i-j|}
    rho <- 0.6
    Sigma <- outer(seq_len(S), seq_len(S), function(i, j) rho^abs(i - j))
    L <- chol(Sigma)

    ## Correlated noise for mediators
    if (identical(s2_noise, "heavyMY")) {
      if (!is.finite(s2_df_m) || s2_df_m <= 2) stop("s2_df_m must be > 2 for heavyMY.")
      z_raw <- matrix(rnorm(n * S), n, S) %*% L
      w <- rchisq(n, df = s2_df_m) / s2_df_m
      epsM <- z_raw / sqrt(w)
      if (isTRUE(s2_standardize_m)) {
        ## standardize each margin to unit variance
        epsM <- epsM / sqrt(s2_df_m / (s2_df_m - 2))
      }
      epsM <- s2_scale_m * epsM
    } else {
      epsM <- matrix(rnorm(n * S), n, S) %*% L
    }

    ## Mediators
    M_mat <- sweep(epsM, 2, b0, "+") +
      sweep(matrix(X, n, S), 2, bA, "*") +
      sweep(matrix(Z, n, S), 2, bZ, "*")
    colnames(M_mat) <- paste0("M", seq_len(S))

    ## Outcome
    if (identical(s2_noise, "gaussian")) {
      eY <- rnorm(n, mean = 0, sd = 1)
    } else {
      if (!is.finite(s2_df_y) || s2_df_y <= 2) stop("s2_df_y must be > 2 for heavy-tailed Y.")
      eY <- rt(n, df = s2_df_y)
      if (isTRUE(s2_standardize_y)) {
        eY <- eY / sqrt(s2_df_y / (s2_df_y - 2))
      }
      eY <- s2_scale_y * eY
    }
    Y <- 1 +
      0.6 * s2_outcome_scale * X +
      0.2 * Z +
      0.2 * s2_outcome_scale * rowSums(M_mat) +
      s2_outcome_scale * sin(M_mat[, 1] * M_mat[, 2]) +
      eY

    dat <- data.frame(
      X = X,
      Z = Z,
      M_mat,
      Y = Y
    )
  } else if (scenario == "S3") {
    ## -----------------------------
    ## S3: simple-but-challenging setting
    ## - 5 non-Gaussian confounders
    ## - binary treatment X ~ Bernoulli(0.5)
    ## - multiple/high-dimensional mediators with shared latent factor
    ## - nonlinear outcome via M1*M2 interaction
    ## -----------------------------
    S <- as.integer(S)
    if (S < 2L) stop("For scenario S3, S must be >= 2.")

    ## Non-Gaussian confounders
    C1 <- rbinom(n, size = 1, prob = 0.5)
    C2 <- rt(n, df = 4)
    C3 <- rgamma(n, shape = 2, rate = 1) - 2
    C4 <- runif(n, min = -2, max = 2)
    C5 <- ifelse(
      rbinom(n, 1, 0.5) == 1,
      rnorm(n, mean = -1.0, sd = 1.0),
      rnorm(n, mean = 1.0, sd = 1.0)
    )

    ## Binary treatment with fixed propensity 0.5
    X <- rbinom(n, size = 1, prob = 0.5)

    ## Mediator model: simple linear structure + shared latent factor
    a <- ifelse(seq_len(S) <= 5L, 0.6, 0.3)
    u <- rnorm(n, mean = 0, sd = 1)  # shared latent factor induces dependence
    eta <- matrix(rt(n * S, df = 5) / sqrt(5 / 3), nrow = n, ncol = S)
    epsM <- 0.5 * matrix(u, nrow = n, ncol = S) + eta

    C5_pos <- as.numeric(C5 > 0)
    M_mat <- matrix(NA_real_, nrow = n, ncol = S)
    for (j in seq_len(S)) {
      M_mat[, j] <-
        0.3 +
        a[j] * X +
        0.2 * C1 +
        0.15 * C2 +
        0.1 * C3 -
        0.1 * C4 +
        0.1 * C5_pos +
        epsM[, j]
    }
    colnames(M_mat) <- paste0("M", seq_len(S))

    ## Outcome model
    M_bar <- rowMeans(M_mat)
    eY <- rt(n, df = 5) / sqrt(5 / 3)
    Y <- 1.5 +
      0.5 * X +
      0.4 * C1 +
      0.2 * C2 +
      0.1 * C3 -
      0.10 * C4 +
      0.2 * as.numeric(C5 > 0) +
      0.15 * M_bar +
      0.2 * sin(M_mat[, 1] * M_mat[, 2]) +
      eY

    dat <- data.frame(
      X = X,
      C1 = C1,
      C2 = C2,
      C3 = C3,
      C4 = C4,
      C5 = C5,
      M_mat,
      Y = Y
    )
  }

  dat
}
