# Simulating networks with given parameters
library(ergm)
library(network)

# number of nodes
N <- 20

net <- network.initialize(N, directed=F)
model_formula <- net ~ kstar(1:2)
coefs <- c(-2, 0.0042)

n = 2
sim_nets <- simulate(model_formula, coef=coefs, nsim=n, output="network")

statistics = summary(sim_nets ~ kstar(1:2))
m = as.matrix.network(sim_nets[[1]], matrix.type="adjacency")

sweep <- function(Y, parameters)
{
  # goes through each possible edge and updates with some probability
  # currently I think it ends up with too many edges
  adj_matrix = as.matrix.network(Y, matrix.type="adjacency")
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
        delta = parameters[1] + parameters[2] * (degrees[i] + degrees[j])
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
  network(adj_matrix, directed=F)
}

# for (i in 1:100)
# {
#   sweep_test = sweep(sweep_test, c(-2, 0.0042))
# }

sweep_test = sweep(sim_nets[[1]], c(-2, 0.0042))
plot(sim_nets[[1]])
plot(sweep_test)