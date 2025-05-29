source('gibbs_sampling.R')

theta_verify <- function(theta)
{
  condition1 = (theta[1] >= -4) & (theta[1] <= 2)
  condition2 = (theta[2] >= -0.05) & (theta[2] <= 1)
  condition3 = theta[1]^2 + theta[2]^2 <= 16
  condition1 & condition2 & condition3
}

sgld <- function(y, epsilon_D, nabla = 0, m = 10, iterations = 10000,
                 space=1) {
  theta <- matrix(ncol = 2, nrow = iterations + 1)
  theta[1, ] <- c(runif(1, -2.5, -1.5), runif(1, 0, 0.01))
  E_k <- c(kstar(y, 1), kstar(y, 2))
  for(t in 1:iterations) {
    y_tilde <- inner_markov_chain(y, m, theta[t,], space) # List of adjacency matrices?
    my_sum <- c(0, 0)
    for(i in 1:m) {
      my_sum <- my_sum + c(kstar(y_tilde[[i]], 1), kstar(y_tilde[[i]], 2))
    }
    estimate <- my_sum / m
    nabla_cond <- E_k - estimate + nabla
    theta_prime <- theta[t,] + 0.5*epsilon_D %*% nabla_cond + rnorm(1, 0, epsilon_D)
    if (theta_verify(theta_prime))
    {
      theta[t+1, ] = theta_prime
    }
    else
    {
      theta[t+1, ] = theta[t, ]
    }
  }
  return(theta)
}

inner_markov_chain <- function(y, m, theta, space=1)
{
  # from gibbs_sampling.R
  sample_networks(y, m, theta, nrow(y), space)
}
