# Energy Score loss for training conditional generators
#
# Computes the empirical Energy Score (ES) using K Monte Carlo samples
# per observation. This loss encourages the generator to match the full
# conditional distribution rather than only the mean.
#
# Args:
#   x0: observed samples (n × d torch tensor)
#   ...: at least two generated samples of the same shape
#   beta: exponent for Euclidean distance in ES (beta = 1 recovers the standard ES)
#   verbose: whether to return diagnostic terms
#
# Returns:
#   A scalar torch loss to be minimized (negative ES), or a list if verbose = TRUE.

energyloss_es <- function(x0, ..., beta = 1, verbose = FALSE) {
  if (!is.numeric(beta) || length(beta) != 1L || !is.finite(beta) || beta <= 0) {
    stop("beta must be a positive finite scalar")
  }
  
  gens <- list(...)
  K <- length(gens)
  if (K < 2L) {
    stop("energyloss_es requires at least 2 generated samples")
  }
  # Stack generated samples into n × K × d
  x_stack <- torch_stack(gens, dim = 2)
  
  n <- x_stack$size(1)
  d <- x_stack$size(3)
  
  # Term 2: average distance between generated samples and observations
  x0_exp <- x0$unsqueeze(2)$expand(c(n, K, d))
  diff2  <- x_stack - x0_exp
  norm2  <- torch_norm(diff2, p = 2, dim = 3)
  if (abs(beta - 1) < 1e-12) {
    dist2 <- norm2
  } else {
    dist2 <- torch_pow(torch_clamp(norm2, min = 1e-12), beta)
  }
  term2_i <- torch_mean(dist2, dim = 2)
  
  # Term 1: average pairwise distance among generated samples
  u1 <- x_stack$unsqueeze(3)
  u2 <- x_stack$unsqueeze(2)
  diff1 <- u1 - u2
  dist1 <- torch_norm(diff1, p = 2, dim = 4)
  if (abs(beta - 1) >= 1e-12) {
    dist1 <- torch_pow(torch_clamp(dist1, min = 1e-12), beta)
  }
  
  term1_i <- torch_sum(torch_sum(dist1, dim = 3), dim = 2) / (2 * K * (K - 1))
  
  ES_i <- term1_i - term2_i
  ES   <- torch_mean(ES_i)
  loss <- -ES
  
  if (verbose) {
    return(list(
      total_loss = loss,
      ES         = ES,
      term1_mean = torch_mean(term1_i),
      term2_mean = torch_mean(term2_i)
    ))
  } else {
    return(loss)
  }
}
