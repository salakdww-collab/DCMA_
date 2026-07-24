# Feedforward generator g(x, ε)
#
# Used as the conditional generator in DCMA. Takes observed inputs x and
# noise ε, and outputs mediators or outcomes.
#
# Args:
#   in_dim: dimension of observed inputs x
#   noise_dim: dimension of noise ε
#   hidden_dim: hidden layer width
#   out_dim: output dimension (S for mediators, 1 for outcome)
#   num_layer: total number of linear layers
#
# Returns:
#   A torch nn_sequential model mapping (x, ε) to outputs

nn_model <- function(in_dim, noise_dim, hidden_dim = 100, out_dim, num_layer = 3) {
  layers <- list()

  layers <- append(layers, nn_linear(in_dim + noise_dim, hidden_dim))
  layers <- append(layers, nn_relu())

  if (num_layer > 2) {
    for (i in 1:(num_layer - 2)) {
      layers <- append(layers, nn_linear(hidden_dim, hidden_dim))
      layers <- append(layers, nn_relu())
    }
  }

  layers <- append(layers, nn_linear(hidden_dim, out_dim))

  return(do.call(nn_sequential, layers))
}

# Feedforward generator with layer-wise reinjection of X
#
# Assumes the first x_dim observed input columns correspond to exposure X and
# concatenates them back to every hidden layer.
nn_model_reinject_x <- function(in_dim, noise_dim, hidden_dim = 100, out_dim, num_layer = 3, x_dim = 1L) {
  x_dim <- as.integer(x_dim)
  stopifnot(x_dim >= 1L, x_dim <= in_dim)

  nn_module(
    initialize = function() {
      self$x_dim <- x_dim
      self$input <- nn_linear(in_dim + noise_dim, hidden_dim)
      self$hidden <- nn_module_list()
      if (num_layer > 2) {
        for (i in seq_len(num_layer - 2L)) {
          self$hidden$append(nn_linear(hidden_dim + x_dim, hidden_dim))
        }
      }
      self$output <- nn_linear(hidden_dim + x_dim, out_dim)
    },
    forward = function(z) {
      x_skip <- z[, seq_len(self$x_dim), drop = FALSE]
      h <- nnf_relu(self$input(z))
      if (length(self$hidden) > 0) {
        for (i in seq_len(length(self$hidden))) {
          h <- nnf_relu(self$hidden[[i]](torch_cat(list(h, x_skip), dim = 2)))
        }
      }
      self$output(torch_cat(list(h, x_skip), dim = 2))
    }
  )()
}

# Feedforward generator with layer-wise reinjection of observed inputs
#
# Concatenates the first obs_dim observed input columns back to every hidden
# layer and the output layer. This keeps deterministic inputs visible even
# when stochastic noise is present.
nn_model_reinject_obs <- function(in_dim, noise_dim, hidden_dim = 100, out_dim, num_layer = 3, obs_dim = in_dim) {
  obs_dim <- as.integer(obs_dim)
  stopifnot(obs_dim >= 1L, obs_dim <= in_dim)

  nn_module(
    initialize = function() {
      self$obs_dim <- obs_dim
      self$input <- nn_linear(in_dim + noise_dim, hidden_dim)
      self$hidden <- nn_module_list()
      if (num_layer > 2) {
        for (i in seq_len(num_layer - 2L)) {
          self$hidden$append(nn_linear(hidden_dim + obs_dim, hidden_dim))
        }
      }
      self$output <- nn_linear(hidden_dim + obs_dim, out_dim)
    },
    forward = function(z) {
      obs_skip <- z[, seq_len(self$obs_dim), drop = FALSE]
      h <- nnf_relu(self$input(z))
      if (length(self$hidden) > 0) {
        for (i in seq_len(length(self$hidden))) {
          h <- nnf_relu(self$hidden[[i]](torch_cat(list(h, obs_skip), dim = 2)))
        }
      }
      self$output(torch_cat(list(h, obs_skip), dim = 2))
    }
  )()
}

# Feedforward generator with a single reinjection at the second hidden layer
#
# This variant concatenates the observed inputs back only once, before the
# second hidden transformation (if it exists), instead of reinjecting at
# every hidden layer.
nn_model_reinject_obs_second <- function(in_dim, noise_dim, hidden_dim = 100, out_dim, num_layer = 3, obs_dim = in_dim) {
  obs_dim <- as.integer(obs_dim)
  stopifnot(obs_dim >= 1L, obs_dim <= in_dim)
  if (num_layer < 3L) {
    stop("nn_model_reinject_obs_second requires num_layer >= 3 so the second-layer reinjection is actually applied.")
  }

  nn_module(
    initialize = function() {
      self$obs_dim <- obs_dim
      self$input <- nn_linear(in_dim + noise_dim, hidden_dim)
      self$hidden <- nn_module_list()
      if (num_layer > 2) {
        for (i in seq_len(num_layer - 2L)) {
          if (i == 1L) {
            self$hidden$append(nn_linear(hidden_dim + obs_dim, hidden_dim))
          } else {
            self$hidden$append(nn_linear(hidden_dim, hidden_dim))
          }
        }
      }
      self$output <- nn_linear(hidden_dim, out_dim)
    },
    forward = function(z) {
      obs_skip <- z[, seq_len(self$obs_dim), drop = FALSE]
      h <- nnf_relu(self$input(z))
      if (length(self$hidden) > 0) {
        for (i in seq_len(length(self$hidden))) {
          if (i == 1L) {
            h <- nnf_relu(self$hidden[[i]](torch_cat(list(h, obs_skip), dim = 2)))
          } else {
            h <- nnf_relu(self$hidden[[i]](h))
          }
        }
      }
      self$output(h)
    }
  )()
}

# Feedforward generator with a single FiLM modulation at the second hidden layer
#
# The observed inputs generate per-feature scale/shift parameters that modulate
# the hidden state once, after the first hidden transformation.
nn_model_film_obs_second <- function(in_dim, noise_dim, hidden_dim = 100, out_dim, num_layer = 3, obs_dim = in_dim) {
  obs_dim <- as.integer(obs_dim)
  stopifnot(obs_dim >= 1L, obs_dim <= in_dim)
  if (num_layer < 3L) {
    stop("nn_model_film_obs_second requires num_layer >= 3 so the second-layer FiLM modulation is actually applied.")
  }

  nn_module(
    initialize = function() {
      self$obs_dim <- obs_dim
      self$input <- nn_linear(in_dim + noise_dim, hidden_dim)
      self$hidden <- nn_module_list()
      if (num_layer > 2) {
        for (i in seq_len(num_layer - 2L)) {
          self$hidden$append(nn_linear(hidden_dim, hidden_dim))
        }
      }
      self$film_scale <- nn_linear(obs_dim, hidden_dim)
      self$film_shift <- nn_linear(obs_dim, hidden_dim)
      self$output <- nn_linear(hidden_dim, out_dim)
    },
    forward = function(z) {
      obs_skip <- z[, seq_len(self$obs_dim), drop = FALSE]
      h <- nnf_relu(self$input(z))
      if (length(self$hidden) > 0) {
        for (i in seq_len(length(self$hidden))) {
          h <- self$hidden[[i]](h)
          if (i == 1L) {
            gamma <- self$film_scale(obs_skip)
            beta <- self$film_shift(obs_skip)
            h <- h * (1 + gamma) + beta
          }
          h <- nnf_relu(h)
        }
      }
      self$output(h)
    }
  )()
}

# Feedforward generator with a single LayerNorm at the second hidden layer
nn_model_ln_second <- function(in_dim, noise_dim, hidden_dim = 100, out_dim, num_layer = 3) {
  if (num_layer < 3L) {
    stop("nn_model_ln_second requires num_layer >= 3 so the second-layer normalization is actually applied.")
  }
  nn_module(
    initialize = function() {
      self$input <- nn_linear(in_dim + noise_dim, hidden_dim)
      self$hidden <- nn_module_list()
      if (num_layer > 2) {
        for (i in seq_len(num_layer - 2L)) {
          self$hidden$append(nn_linear(hidden_dim, hidden_dim))
        }
      }
      self$ln2 <- nn_layer_norm(hidden_dim)
      self$output <- nn_linear(hidden_dim, out_dim)
    },
    forward = function(z) {
      h <- nnf_relu(self$input(z))
      if (length(self$hidden) > 0) {
        for (i in seq_len(length(self$hidden))) {
          h <- self$hidden[[i]](h)
          if (i == 1L) {
            h <- self$ln2(h)
          }
          h <- nnf_relu(h)
        }
      }
      self$output(h)
    }
  )()
}

# Feedforward generator with a single AdaLN modulation at the second hidden layer
nn_model_adaln_obs_second <- function(in_dim, noise_dim, hidden_dim = 100, out_dim, num_layer = 3, obs_dim = in_dim) {
  obs_dim <- as.integer(obs_dim)
  stopifnot(obs_dim >= 1L, obs_dim <= in_dim)
  if (num_layer < 3L) {
    stop("nn_model_adaln_obs_second requires num_layer >= 3 so the second-layer AdaLN modulation is actually applied.")
  }

  nn_module(
    initialize = function() {
      self$obs_dim <- obs_dim
      self$input <- nn_linear(in_dim + noise_dim, hidden_dim)
      self$hidden <- nn_module_list()
      if (num_layer > 2) {
        for (i in seq_len(num_layer - 2L)) {
          self$hidden$append(nn_linear(hidden_dim, hidden_dim))
        }
      }
      self$ln2 <- nn_layer_norm(hidden_dim)
      self$film_scale <- nn_linear(obs_dim, hidden_dim)
      self$film_shift <- nn_linear(obs_dim, hidden_dim)
      self$output <- nn_linear(hidden_dim, out_dim)
    },
    forward = function(z) {
      obs_skip <- z[, seq_len(self$obs_dim), drop = FALSE]
      h <- nnf_relu(self$input(z))
      if (length(self$hidden) > 0) {
        for (i in seq_len(length(self$hidden))) {
          h <- self$hidden[[i]](h)
          if (i == 1L) {
            h <- self$ln2(h)
            gamma <- self$film_scale(obs_skip)
            beta <- self$film_shift(obs_skip)
            h <- h * (1 + gamma) + beta
          }
          h <- nnf_relu(h)
        }
      }
      self$output(h)
    }
  )()
}

# Structured mediator generator with shared and idiosyncratic branches
#
# The deterministic mean branch captures systematic effects of observed inputs,
# the shared branch captures common latent variation across mediators, and the
# idiosyncratic branch absorbs mediator-specific residual randomness.
nn_model_shared_resid_mediator <- function(
  in_dim,
  noise_dim,
  hidden_dim = 100,
  out_dim,
  shared_noise_dim = NULL
) {
  if (noise_dim < 2L) {
    stop("nn_model_shared_resid_mediator requires noise_dim >= 2.")
  }
  if (is.null(shared_noise_dim)) {
    shared_noise_dim <- max(1L, floor(noise_dim / 2))
  }
  shared_noise_dim <- as.integer(shared_noise_dim)
  stopifnot(shared_noise_dim >= 1L, shared_noise_dim < noise_dim)
  idio_noise_dim <- as.integer(noise_dim - shared_noise_dim)

  nn_module(
    initialize = function() {
      self$in_dim <- as.integer(in_dim)
      self$shared_noise_dim <- shared_noise_dim
      self$idio_noise_dim <- idio_noise_dim

      self$mean_net <- nn_sequential(
        nn_linear(in_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, out_dim)
      )

      self$shared_net <- nn_sequential(
        nn_linear(in_dim + shared_noise_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, out_dim)
      )

      self$idio_net <- nn_sequential(
        nn_linear(in_dim + idio_noise_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, out_dim)
      )
    },
    forward = function(z) {
      x <- z[, seq_len(self$in_dim), drop = FALSE]
      eps <- z[, (self$in_dim + 1L):z$size(2), drop = FALSE]
      eps_shared <- eps[, seq_len(self$shared_noise_dim), drop = FALSE]
      eps_idio <- eps[, (self$shared_noise_dim + 1L):eps$size(2), drop = FALSE]

      mu <- self$mean_net(x)
      shared_part <- self$shared_net(torch_cat(list(x, eps_shared), dim = 2))
      idio_part <- self$idio_net(torch_cat(list(x, eps_idio), dim = 2))
      mu + shared_part + idio_part
    }
  )()
}

# Separate mediator generator.
#
# This module keeps one independent generator per mediator. The input is
# x concatenated with S independent noise blocks, each of dimension noise_dim.
# It returns a mediator vector, but the networks do not share parameters or
# noise, so it implements the ablation
#   prod_s P(M_s | A, Z)
# instead of a joint P(M_1, ..., M_S | A, Z).
nn_model_separate_mediator <- function(
  in_dim,
  noise_dim,
  hidden_dim = 100,
  out_dim,
  num_layer = 3
) {
  stopifnot(in_dim >= 1L, noise_dim >= 1L, hidden_dim >= 1L, out_dim >= 1L)

  nn_module(
    initialize = function() {
      self$in_dim <- as.integer(in_dim)
      self$noise_dim <- as.integer(noise_dim)
      self$out_dim <- as.integer(out_dim)
      self$gens <- nn_module_list()
      for (j in seq_len(self$out_dim)) {
        self$gens$append(
          nn_model(
            in_dim = self$in_dim,
            noise_dim = self$noise_dim,
            hidden_dim = hidden_dim,
            out_dim = 1L,
            num_layer = num_layer
          )
        )
      }
    },
    forward = function(z) {
      x <- z[, seq_len(self$in_dim), drop = FALSE]
      eps <- z[, (self$in_dim + 1L):z$size(2), drop = FALSE]
      outs <- vector("list", self$out_dim)
      for (j in seq_len(self$out_dim)) {
        idx <- ((j - 1L) * self$noise_dim + 1L):(j * self$noise_dim)
        eps_j <- eps[, idx, drop = FALSE]
        outs[[j]] <- self$gens[[j]](torch_cat(list(x, eps_j), dim = 2))
      }
      torch_cat(outs, dim = 2)
    }
  )()
}

# Structured outcome generator with location-scale residual decomposition
#
# The generator factorizes output variation into a deterministic location
# component, an input-dependent scale, and a noise-driven residual shape term.
nn_model_locscale_outcome <- function(
  in_dim,
  noise_dim,
  hidden_dim = 100,
  scale_floor = 1e-3
) {
  stopifnot(noise_dim >= 1L)
  stopifnot(is.numeric(scale_floor), length(scale_floor) == 1L, is.finite(scale_floor), scale_floor > 0)

  nn_module(
    initialize = function() {
      self$in_dim <- as.integer(in_dim)
      self$scale_floor <- scale_floor

      self$mu_net <- nn_sequential(
        nn_linear(in_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, 1L)
      )

      self$scale_net <- nn_sequential(
        nn_linear(in_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, 1L)
      )

      self$resid_net <- nn_sequential(
        nn_linear(in_dim + noise_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, 1L)
      )
    },
    forward = function(z) {
      x <- z[, seq_len(self$in_dim), drop = FALSE]
      mu <- self$mu_net(x)
      scale_raw <- self$scale_net(x)
      scale <- nnf_softplus(scale_raw) + self$scale_floor
      resid <- self$resid_net(z)
      mu + scale * resid
    }
  )()
}

# Additive external-noise outcome generator.
#
# The noise is restricted to enter additively through the first noise coordinate:
#   homoskedastic:    Y = mu(x) + eps
#   heteroskedastic:  Y = mu(x) + sigma(x) eps
# This is intentionally less expressive than the main internal-noise generator
# g(x, eps), and is used only for ablation.
nn_model_additive_outcome <- function(
  in_dim,
  hidden_dim = 100,
  heteroskedastic = TRUE,
  scale_floor = 1e-3
) {
  stopifnot(in_dim >= 1L, hidden_dim >= 1L)
  stopifnot(is.logical(heteroskedastic), length(heteroskedastic) == 1L)

  nn_module(
    initialize = function() {
      self$in_dim <- as.integer(in_dim)
      self$heteroskedastic <- heteroskedastic
      self$scale_floor <- scale_floor

      self$mu_net <- nn_sequential(
        nn_linear(in_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, 1L)
      )

      if (isTRUE(self$heteroskedastic)) {
        self$scale_net <- nn_sequential(
          nn_linear(in_dim, hidden_dim),
          nn_relu(),
          nn_linear(hidden_dim, hidden_dim),
          nn_relu(),
          nn_linear(hidden_dim, 1L)
        )
      }
    },
    forward = function(z) {
      x <- z[, seq_len(self$in_dim), drop = FALSE]
      eps <- z[, self$in_dim + 1L, drop = FALSE]
      mu <- self$mu_net(x)
      if (isTRUE(self$heteroskedastic)) {
        scale <- nnf_softplus(self$scale_net(x)) + self$scale_floor
      } else {
        scale <- torch_ones_like(mu)
      }
      mu + scale * eps
    }
  )()
}

# Outcome generator with shared trunk and X-gated expert heads
#
# This is a generic specialization mechanism: one shared representation is
# learned for all conditions, while two lightweight heads are mixed using a
# gate driven by the exposure X. For binary X, this behaves like treatment-
# specific heads; for continuous X, it still provides a smooth interpolation.
nn_model_xgated_heads <- function(
  in_dim,
  noise_dim,
  hidden_dim = 100,
  out_dim = 1L,
  num_layer = 3
) {
  stopifnot(in_dim >= 1L, noise_dim >= 1L, hidden_dim >= 1L, out_dim >= 1L, num_layer >= 2L)

  nn_module(
    initialize = function() {
      self$in_dim <- as.integer(in_dim)
      self$trunk_in <- nn_linear(in_dim + noise_dim, hidden_dim)
      self$trunk_hidden <- nn_module_list()
      if (num_layer > 2L) {
        for (i in seq_len(num_layer - 2L)) {
          self$trunk_hidden$append(nn_linear(hidden_dim, hidden_dim))
        }
      }
      self$head0 <- nn_sequential(
        nn_linear(hidden_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, out_dim)
      )
      self$head1 <- nn_sequential(
        nn_linear(hidden_dim, hidden_dim),
        nn_relu(),
        nn_linear(hidden_dim, out_dim)
      )
      self$gate <- nn_sequential(
        nn_linear(1L, hidden_dim %/% 2L + 1L),
        nn_relu(),
        nn_linear(hidden_dim %/% 2L + 1L, out_dim)
      )
    },
    forward = function(z) {
      x <- z[, 1, drop = FALSE]
      h <- nnf_relu(self$trunk_in(z))
      if (length(self$trunk_hidden) > 0) {
        for (i in seq_len(length(self$trunk_hidden))) {
          h <- nnf_relu(self$trunk_hidden[[i]](h))
        }
      }
      y0 <- self$head0(h)
      y1 <- self$head1(h)
      gate <- torch_sigmoid(self$gate(x))
      (1 - gate) * y0 + gate * y1
    }
  )()
}
