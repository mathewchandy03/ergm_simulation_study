source('gibbs_sampling.R')
sgld <- function(y, epsilon_D, nabla = 0, m = 10, iterations = 10000) {
  theta <- as.matrix(NA, nrow = 2, ncol = iterations + 1)
  theta[0, ] <- c(runif(1, -2.5, -1.5), runif(1, 0, 0.01))
  E_k <- c(kstar(y, 1), kstar(y, 2))
  for(t in 1:iterations) {
    y_tilde <- inner_markov_chain(y, m, theta) # List of adjacency matrices?
    my_sum <- c(0, 0)
    for(i in 1:m) {
      my_sum <- my_sum + c(kstar(y_tilde[[i]], 1), kstar(y_tilde[[i]], 2))
    }
    estimate <- my_sum / m
    nabla <- E_k - estimate + nabla
    theta[1, ] <- as.vector(as.matrix(theta) + as.matrix(epsilon_D / 2) %*% as.matrix(nabla) + as.matrix(rnorm(1, 0, sqrt(epsilon_D))))
  }
  return(theta)
}

inner_markov_chain <- function(y, m, theta) {
  N <- nrow(y)
  degrees <- rowSums(y)
  y_tilde <- vector(mode = "list", length = 10)
  for(k in 1:m) {
    for (i in 1:N) {
      for (j in 1:i) {
        if (j != i) {
          current_value = y[i, j]
          delta = theta[1] + theta[2] * (degrees[i] + degrees[j])
          p = 1/(1+exp(-delta))
          new_value = 1*(runif(1) <= p)
          y[i,j] = new_value
          y[j,i] = new_value
          
          change = new_value - current_value
          degrees[i] = degrees[i] + change
          degrees[j] = degrees[j] + change 
        }
      }
    }
    y_tilde[[k]] <- y
  }
  return(y_tilde)
}
