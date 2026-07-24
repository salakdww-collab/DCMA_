nn_critic <- function(input_dim, hidden_dim = 64, num_layer = 3) {
  layers <- list()

  layers <- append(layers, nn_linear(input_dim, hidden_dim))
  layers <- append(layers, nn_relu())

  if (num_layer > 2) {
    for (i in seq_len(num_layer - 2)) {
      layers <- append(layers, nn_linear(hidden_dim, hidden_dim))
      layers <- append(layers, nn_relu())
    }
  }

  layers <- append(layers, nn_linear(hidden_dim, 1L))
  do.call(nn_sequential, layers)
}
