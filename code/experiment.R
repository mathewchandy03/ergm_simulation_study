# Running Experiment

source("./gibbs_sampling.R")
source("./sgld.R")

run_experiment <- function(n, theta, nodes, epsilon_D, iterations=10000)
{
  # simulate 50 networks
  networks = sample_networks(m = n,
                             theta = theta, 
                             nodes = nodes, 
                             iters = 10^4)
  
  # estimate theta with SGLD
  estimates = matrix(ncol=2, nrow=n)
  for (i in 1:n)
  {
    estimate_chain = sgld(networks[[i]], epsilon_D, iterations=iterations)
    estimates[i, ] = estimate_chain[iterations+1, ]
  }
  estimates
}

# running an experiment
results = run_experiment(1, c(-2, 0.0042), 20, diag(c(0.004, 0.00012)))

# getting SGLD chain
Y = sample_networks(Y=NULL, 1, c(-2, 0.0042), 20, 10^4)[[1]]
theta_chain = sgld(Y, diag(c(0.004, 0.00012)), iterations=1000)

plot(theta_chain[, 1])

