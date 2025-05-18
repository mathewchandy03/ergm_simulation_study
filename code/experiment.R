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
  times = rep(NA, n)
  for (i in 1:n)
  {
    start.time <- Sys.time()
    estimate_chain = sgld(networks[[i]], epsilon_D, iterations=iterations)
    end.time <- Sys.time()
    estimates[i, ] = estimate_chain[iterations+1, ]
    times[i] = end.time - start.time
  }
 list(estimates, times)
}

# running an experiment
results = run_experiment(1, c(-2, 0.0042), 20, diag(c(0.004, 0.00012)))

# getting SGLD chain
Y = sample_networks(Y=NULL, 1, c(-2, 0.0042), 20, 10^4)[[1]]
theta_chain = sgld(Y, diag(c(0.004, 0.00012)), iterations=1000)

plot(theta_chain[, 1])

data <- data.frame(n = c(), theta_1_bias_mean = c(), theta_1_bias_sd = c(),
                   theta_2_bias_mean = c(), theta_2_bias_msd = c(), 
                   average_cpu_time = c())
# Section 4.2
for (n in c(20, 50, 75, 100)) {
  results = run_experiment(50, c(-2, 0.0042), n, diag(c(0.004, 0.00012)))
  theta_bias <- base::sweep(results[[1]], 2, c(-2, 0.0042))
  theta_1_bias_mean <- mean(theta_bias[,1])
  theta_1_bias_sd <- sd(theta_bias[,1])
  theta_2_bias_mean <- mean(theta_bias[,2])
  theta_2_bias_sd <- sd(theta_bias[,2])
  average_cpu_time <- mean(results[[2]])
  new_row <- c(n, theta_1_bias_mean, theta_1_bias_sd, theta_2_bias_mean, 
               theta_2_bias_sd)
  data <- rbind(data, new_row)
}

