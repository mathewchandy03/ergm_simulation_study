# Simulating networks with given theta
library(ergm)
library(network)

sweep <- function(Y, theta)
{
  # goes through each possible edge and updates with some probability
  adj_matrix = Y
  N = nrow(adj_matrix)
  degrees = rowSums(adj_matrix)
  
  for (i in 1:N)
  {
    for (j in 1:i)
    {
      if (j != i)
      {
        current_value = adj_matrix[i, j]
        
        # appendix B of supplemental of paper we are recreating
        delta = 2*theta[1] + theta[2] * (degrees[i] + degrees[j])
        # equation 2 of paper 
        # https://csss.uw.edu/Papers/wp39.pdf
        p = 1/(1+exp(-delta))
        
        # i and j share an edge with probability p
        new_value = 1*(runif(1) <= p)
        adj_matrix[i,j] = new_value
        adj_matrix[j,i] = new_value
        
        change = new_value - current_value
        degrees[i] = degrees[i] + change
        degrees[j] = degrees[j] + change 
      }
    }
  }
  adj_matrix
}

sample_network <- function(Y, theta, iters)
{
  # runs a iters sweeps to sample a network, from a starting network Y
  # Y: adjacency matrix
  current_sweep = Y
  for (i in 1:(iters))
  {
    current_sweep = sweep(current_sweep, theta)
  }
  current_sweep
}

sample_networks <- function(Y=NULL, m, theta, nodes, iters)
{
  # net <- network.initialize(nodes, directed=F)
  # model_formula <- net ~ kstar(1:2)
  # coefs <- theta
  # 
  # sim_nets <- simulate(model_formula, coef=coefs, nsim=m, output="network")
  # nets = vector(mode = 'list', length = m)
  # for (i in 1:m)
  # {
  #   nets[[i]] = as.matrix.network(sim_nets[[i]], matrix.type = 'adjacency')
  # }
  
  # n: number of networks to sample
  # iters: iterations between each sample
  # returns a vector with all the networks
  
  
  nets = vector(mode = "list", length = m)

  if (is.null(Y))
  {
    current_net = matrix(0, nrow=nodes, ncol=nodes)
  }
  else
  {
    current_net = Y
  }

  for (i in 1:m)
  {
    current_net = sample_network(current_net, theta, iters)
    nets[[i]] = current_net
  }
  nets
}

kstar <- function(adj_matrix, k)
{
  degrees = rowSums(adj_matrix)
  sum(choose(degrees, k))
}
