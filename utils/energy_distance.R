# 1D energy distance 
#
# Computes the 1D energy distance using an O(n log n) sorting-based implementation.
#
# Args:
#   x: numeric vector
#   y: numeric vector
#
# Returns:
#   Numeric scalar energy distance (NA if x or y is empty).

energy_distance <- function(x, y) {
  sum_abs_pairs_1d <- function(v) {
    v <- sort(as.numeric(v))
    n <- length(v)
    if (n <= 1L) return(0)
    k <- seq_len(n)
    sum((2 * k - n - 1) * v)
  }

  sum_abs_cross_1d <- function(x, y) {
    x <- as.numeric(x)
    y <- sort(as.numeric(y))
    m <- length(y)
    if (m == 0L) stop("y is empty")
    cy <- c(0, cumsum(y))
    Sy <- cy[m + 1L]
    idx <- findInterval(x, y)
    left_sum  <- cy[idx + 1L]
    right_sum <- Sy - left_sum
    left_n  <- idx
    right_n <- m - idx
    sum(x * left_n - left_sum + right_sum - x * right_n)
  }

  x <- as.numeric(x); y <- as.numeric(y)
  n <- length(x); m <- length(y)
  if (n < 1L || m < 1L) return(NA_real_)

  cross <- (2 / (n * m)) * sum_abs_cross_1d(x, y)
  xx <- (2 * sum_abs_pairs_1d(x)) / (n * (n - 1))
  yy <- (2 * sum_abs_pairs_1d(y)) / (m * (m - 1))

  cross - xx - yy
}


