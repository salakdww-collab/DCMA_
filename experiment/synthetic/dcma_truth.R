# Generate truth interventional distributions for DCMA simulations
#
# Scenario S1: single mediator with bimodal outcome under X = 1.
#   - M1 depends on (X, Z) with Gaussian noise.
#   - When X = 1, Y is bimodal and the mixture component is selected by Z.
#   - When X = 0, Y is unimodal.
#
# Scenario S2: multiple correlated mediators with nonlinear outcome.
#   - S mediators are generated with correlated Gaussian noise.
#   - Y depends on rowSums(M) and a nonlinear interaction sin(M1 * M2).
#
# Args:
#   scenario: simulation scenario ("S1" or "S2")
#   N: number of covariate profiles (draws of Z)
#   B: Monte Carlo draws per profile
#   S: number of mediators for scenario S2
#   s1_mode_shift: half-distance between the two treated modes in S1
#   a0, a1: baseline and treated exposure levels
#   a_ipse: exposure level used for pathway specific outcomes
#   seed: optional random seed
#
# Returns:
#   A list with joint interventional outcomes and pathway specific outcomes
#   in the same structure used by downstream effect computation.

dcma_truth_interventions <- function(
  scenario = c("S1", "S2", "S3"),
  N = 20000L,
  B = 200L,
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
  a0 = 0,
  a1 = 1,
  a_ipse = 1,
  seed = NULL
) {
  scenario <- match.arg(scenario)
  s2_noise <- match.arg(s2_noise)
  if (!is.numeric(s1_mode_shift) || length(s1_mode_shift) != 1L || !is.finite(s1_mode_shift) || s1_mode_shift < 0) {
    stop("s1_mode_shift must be one finite non-negative number.")
  }
  if (!is.null(seed)) set.seed(seed)

  ## Number of mediators: S for S2/S3, single mediator otherwise
  S_eff <- if (scenario %in% c("S2", "S3")) S else 1L
  if (scenario == "S3" && S_eff < 2L) stop("For scenario S3, S must be >= 2.")
  m_vars <- paste0("M", seq_len(S_eff))

  ## Covariate profiles
  if (scenario == "S3") {
    C1 <- rbinom(N, size = 1, prob = 0.5)
    C2 <- rt(N, df = 4)
    C3 <- rgamma(N, shape = 2, rate = 1) - 2
    C4 <- runif(N, min = -2, max = 2)
    C5 <- ifelse(
      rbinom(N, 1, 0.5) == 1,
      rnorm(N, mean = -1.0, sd = 1.0),
      rnorm(N, mean = 1.0, sd = 1.0)
    )
    C_mat <- cbind(C1, C2, C3, C4, C5)
    colnames(C_mat) <- paste0("C", 1:5)
  } else {
    Z <- rnorm(N, 0, 1)
  }

  ## ========================================================
  ## 1. True mediator model: draw M | A, Z
  ## ========================================================
  draw_M_given_a_truth <- function(a) {
    n <- N
    B_local <- B

    X_vec <- rep(a, times = n * B_local)
    if (scenario != "S3") {
      if (!is.numeric(s2_mediator_scale) || length(s2_mediator_scale) != 1L || !is.finite(s2_mediator_scale) || s2_mediator_scale <= 0) {
        stop("s2_mediator_scale must be one positive finite number.")
      }
      if (!is.numeric(s2_outcome_scale) || length(s2_outcome_scale) != 1L || !is.finite(s2_outcome_scale) || s2_outcome_scale <= 0) {
        stop("s2_outcome_scale must be one positive finite number.")
      }
      Z_rep <- rep(Z, each = B_local)
    }

    if (scenario %in% c("S1")) {
      eM <- rnorm(n * B_local, mean = 0, sd = 0.5)
      M1_flat <- 0.5 + 1 * X_vec + 0.3 * Z_rep + eM
      M_mat <- cbind(M1_flat)

    } else if (scenario == "S2") {
      S_local <- S_eff

      b0 <- rep(0.5, S_local)
      bA <- c(1, 0.8, 0.6, 0.4, rep(0.2, max(0, S_local - 4))) * s2_mediator_scale
      bZ <- c(0.3, 0.3, 0.2, 0.2, rep(0.1, max(0, S_local - 4))) * s2_mediator_scale

      ## Correlation structure: Sigma(i,j) = rho^{|i-j|}
      rho <- 0.6
      Sigma <- outer(1:S_local, 1:S_local, function(i, j) rho^abs(i - j))
      L <- chol(Sigma)

      ## Correlated noise for mediators
      if (identical(s2_noise, "heavyMY")) {
        if (!is.finite(s2_df_m) || s2_df_m <= 2) stop("s2_df_m must be > 2 for heavyMY.")
        eps_raw <- matrix(rnorm(n * B_local * S_local), n * B_local, S_local) %*% L
        w <- rchisq(n * B_local, df = s2_df_m) / s2_df_m
        epsM_flat <- eps_raw / sqrt(w)
        if (isTRUE(s2_standardize_m)) {
          epsM_flat <- epsM_flat / sqrt(s2_df_m / (s2_df_m - 2))
        }
        epsM_flat <- s2_scale_m * epsM_flat
      } else {
        eps_raw <- matrix(rnorm(n * B_local * S_local), n * B_local, S_local)
        epsM_flat <- eps_raw %*% L
      }

      M_mat <- matrix(NA_real_, n * B_local, S_local)
      for (j in seq_len(S_local)) {
        M_mat[, j] <- b0[j] + bA[j] * X_vec + bZ[j] * Z_rep + epsM_flat[, j]
      }
    } else if (scenario == "S3") {
      S_local <- S_eff
      C_rep <- C_mat[rep(seq_len(N), each = B_local), , drop = FALSE]
      a_j <- ifelse(seq_len(S_local) <= 5L, 0.6, 0.3)
      u <- rnorm(n * B_local, mean = 0, sd = 1)
      eta <- matrix(rt(n * B_local * S_local, df = 5) / sqrt(5 / 3),
                    nrow = n * B_local, ncol = S_local)
      epsM <- 0.5 * matrix(u, nrow = n * B_local, ncol = S_local) + eta
      C5_pos <- as.numeric(C_rep[, "C5"] > 0)

      M_mat <- matrix(NA_real_, n * B_local, S_local)
      for (j in seq_len(S_local)) {
        M_mat[, j] <-
          0.3 +
          a_j[j] * X_vec +
          0.2 * C_rep[, "C1"] +
          0.15 * C_rep[, "C2"] +
          0.1 * C_rep[, "C3"] -
          0.1 * C_rep[, "C4"] +
          0.1 * C5_pos +
          epsM[, j]
      }
    }

    array(
      M_mat,
      dim = c(n, B_local, S_eff),
      dimnames = list(NULL, NULL, m_vars)
    )
  }

  ## ========================================================
  ## 2. True outcome model: draw Y | A, M, Z
  ## ========================================================
  draw_Y_given_M_a_truth <- function(M_array, a) {
    stopifnot(length(dim(M_array)) == 3L)

    n <- dim(M_array)[1]
    B_local <- dim(M_array)[2]
    S_local <- dim(M_array)[3]

    M_mat <- matrix(M_array, nrow = n * B_local, ncol = S_local)
    X_vec <- rep(a, times = n * B_local)
    if (scenario != "S3") {
      Z_rep <- rep(Z, each = B_local)
    } else {
      C_rep <- C_mat[rep(seq_len(N), each = B_local), , drop = FALSE]
    }

    if (scenario == "S1") {
      ## Bimodal outcome when X = 1, unimodal when X = 0; component selected by Z
      M1 <- M_mat[, 1]

      ## Baseline mean (centered so X=0, Z=0, M=0 has mean 4)
      mu_base <- 4 + 0.3 * X_vec + 0.2 * Z_rep + 0.5 * M1

      ## Regime indicators
      Scls  <- ifelse(X_vec == 1, 1, 0)
      Scls1 <- ifelse(Z_rep > 0, 0, 1)

      Y_flat <- numeric(n * B_local)

      sd0 <- 1
      sd1 <- 1
      sd2 <- 1

      ## Left mode (X=1, Z>0)
      idx0 <- which(Scls == 1 & Scls1 == 0)
      if (length(idx0) > 0) {
        eY0 <- rnorm(length(idx0), 0, sd0)
        Y_flat[idx0] <- (mu_base[idx0] - s1_mode_shift) + eY0
      }

      ## Right mode (X=1, Z<=0)
      idx1 <- which(Scls == 1 & Scls1 == 1)
      if (length(idx1) > 0) {
        eY1 <- rnorm(length(idx1), 0, sd1)
        Y_flat[idx1] <- (mu_base[idx1] + s1_mode_shift) + eY1
      }

      ## Unimodal regime (X=0)
      idx2 <- which(Scls == 0)
      if (length(idx2) > 0) {
        eY2 <- rnorm(length(idx2), 0, sd2)
        Y_flat[idx2] <- (mu_base[idx2] + 0.3) + eY2
      }

    } else if (scenario == "S2") {
      ## Multiple mediators with nonlinear outcome: rowSums(M) and sin(M1 * M2)
      S_local <- S_eff
      if (identical(s2_noise, "gaussian")) {
        eY <- rnorm(n * B_local, 0, 1)
      } else {
        if (!is.finite(s2_df_y) || s2_df_y <= 2) stop("s2_df_y must be > 2 for heavy-tailed Y.")
        eY <- rt(n * B_local, df = s2_df_y)
        if (isTRUE(s2_standardize_y)) {
          eY <- eY / sqrt(s2_df_y / (s2_df_y - 2))
        }
        eY <- s2_scale_y * eY
      }

      Y_flat <-
        1 +
        0.6 * s2_outcome_scale * X_vec +
        0.2 * Z_rep +
        0.2 * s2_outcome_scale * rowSums(M_mat) +
        s2_outcome_scale * sin((M_mat[, 1] * M_mat[, 2])) +
        eY
    } else if (scenario == "S3") {
      m_bar <- rowMeans(M_mat)
      eY <- rt(n * B_local, df = 5) / sqrt(5 / 3)

      Y_flat <-
        1.5 +
        0.5 * X_vec +
        0.4 * C_rep[, "C1"] +
        0.2 * C_rep[, "C2"] +
        0.1 * C_rep[, "C3"] -
        0.10 * C_rep[, "C4"] +
        0.2 * as.numeric(C_rep[, "C5"] > 0) +
        0.15 * m_bar +
        0.2 * sin((M_mat[, 1] * M_mat[, 2])) +
        eY
    }

    array(Y_flat, dim = c(n, B_local))
  }

  ## ========================================================
  ## 3. Block independent recombination for pathway specific mediators
  ## ========================================================
  build_pathway_M_pair_truth <- function(M0, M1, s) {
    stopifnot(length(dim(M0)) == 3L, length(dim(M1)) == 3L)
    stopifnot(all(dim(M0) == dim(M1)))

    n <- dim(M0)[1]
    B_local <- dim(M0)[2]
    S_local <- dim(M0)[3]

    stopifnot(s >= 1L, s <= S_local)

    if (S_local == 1L) {
      return(list(M0s = M0, M1s = M1))
    }

    M0s <- array(NA_real_, dim = c(n, B_local, S_local))
    M1s <- array(NA_real_, dim = c(n, B_local, S_local))

    idx_low  <- if (s > 1L) 1:(s - 1L) else integer(0)
    idx_high <- if (s < S_local) (s + 1L):S_local else integer(0)

    for (i in seq_len(n)) {
      M0_i <- matrix(M0[i, , , drop = FALSE], nrow = B_local)
      M1_i <- matrix(M1[i, , , drop = FALSE], nrow = B_local)

      ## M1s: take upstream from M0, mediator s from M1, downstream from M1
      parts1 <- list()
      if (length(idx_low))
        parts1[[length(parts1) + 1L]] <- M0_i[sample(B_local), idx_low, drop = FALSE]
      parts1[[length(parts1) + 1L]] <- M1_i[sample(B_local), s, drop = FALSE]
      if (length(idx_high))
        parts1[[length(parts1) + 1L]] <- M1_i[sample(B_local), idx_high, drop = FALSE]
      M1s[i, , ] <- do.call(cbind, parts1)

      ## M0s: take upstream from M0, mediator s from M0, downstream from M1
      parts0 <- list()
      if (length(idx_low))
        parts0[[length(parts0) + 1L]] <- M0_i[sample(B_local), idx_low, drop = FALSE]
      parts0[[length(parts0) + 1L]] <- M0_i[sample(B_local), s, drop = FALSE]
      if (length(idx_high))
        parts0[[length(parts0) + 1L]] <- M1_i[sample(B_local), idx_high, drop = FALSE]
      M0s[i, , ] <- do.call(cbind, parts0)
    }

    list(M0s = M0s, M1s = M1s)
  }

  ## ========================================================
  ## 4. Joint interventional quantities
  ## ========================================================
  M0 <- draw_M_given_a_truth(a0)
  M1 <- draw_M_given_a_truth(a1)

  Y_1M1 <- draw_Y_given_M_a_truth(M1, a1)
  Y_0M0 <- draw_Y_given_M_a_truth(M0, a0)
  Y_1M0 <- draw_Y_given_M_a_truth(M0, a1)
  Y_0M1 <- draw_Y_given_M_a_truth(M1, a0)

  ## ========================================================
  ## 5. Pathway specific (block independent) quantities
  ## ========================================================
  path_mediators <- vector("list", S_eff)
  path_outcomes  <- vector("list", S_eff)

  for (s in seq_len(S_eff)) {
    pair_s <- build_pathway_M_pair_truth(M0, M1, s)

    M0_s <- pair_s$M0s
    M1_s <- pair_s$M1s

    Y_a_M0s <- draw_Y_given_M_a_truth(M0_s, a_ipse)
    Y_a_M1s <- draw_Y_given_M_a_truth(M1_s, a_ipse)

    path_mediators[[s]] <- list(M0_s = M0_s, M1_s = M1_s)
    path_outcomes[[s]]  <- list(Y_a_M0s = Y_a_M0s, Y_a_M1s = Y_a_M1s)
  }

  names(path_mediators) <- m_vars
  names(path_outcomes)  <- m_vars

  ## ========================================================
  ## 6. Output object
  ## ========================================================
  out <- list(
    scenario = scenario,
    N = N,
    B = B,
    s1_mode_shift = if (scenario == "S1") s1_mode_shift else NA_real_,
    a0 = a0,
    a1 = a1,
    a_ipse = a_ipse,
    s2_noise = if (scenario == "S2") s2_noise else NA_character_,
    s2_df_y = if (scenario == "S2") s2_df_y else NA_real_,
    s2_df_m = if (scenario == "S2") s2_df_m else NA_real_,
    s2_scale_y = if (scenario == "S2") s2_scale_y else NA_real_,
    s2_scale_m = if (scenario == "S2") s2_scale_m else NA_real_,
    s2_mediator_scale = if (scenario == "S2") s2_mediator_scale else NA_real_,
    s2_outcome_scale = if (scenario == "S2") s2_outcome_scale else NA_real_,
    s2_standardize_y = if (scenario == "S2") s2_standardize_y else NA,
    s2_standardize_m = if (scenario == "S2") s2_standardize_m else NA,
    m_vars = m_vars,
    Z = if (scenario == "S3") NULL else Z,
    C = if (scenario == "S3") C_mat else NULL,
    mediators = list(
      M_a0 = M0,
      M_a1 = M1,
      path = path_mediators
    ),
    outcomes = list(
      Y_1M1 = Y_1M1,
      Y_0M0 = Y_0M0,
      Y_1M0 = Y_1M0,
      Y_0M1 = Y_0M1,
      path  = path_outcomes
    )
  )

  class(out) <- "DCMA_truth_interventions"
  out
}
