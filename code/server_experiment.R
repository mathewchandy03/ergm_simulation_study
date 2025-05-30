# Running Experiment
library(doRNG)
library(foreach)
doFuture::registerDoFuture()
future::plan("multisession", workers=8)
doRNG::registerDoRNG(123)

source("./gibbs_sampling.R")
source("./sgld.R")

run_experiment <- function(n, theta, nodes, epsilon_D, iterations=10000, 
                           space=1)
{
  # simulate n networks
  my_net = matrix(0, nrow=nodes, ncol=nodes)
  result = foreach(i = 1:n, .packages = c('rpm'), .combine = 'rbind') %dorng% {
    new_net = sample_network(my_net, theta, iters = 1e4)
    start.time <- Sys.time()
    estimate_chain = sgld(new_net, epsilon_D, iterations=iterations,
                          space=space)
    end.time <- Sys.time()
    new_row = c(estimate_chain[iterations+1, ], end.time - start.time)
    new_row
  }
  result
  
  #  # estimate theta with SGLD
  #  estimates = matrix(ncol=2, nrow=n)
  #  times = rep(NA, n)
  #  for (i in 1:n)
  #  {
  #    start.time <- Sys.time()
  #    estimate_chain = sgld(networks[[i]], epsilon_D, iterations=iterations,
  #                          space=space)
  #    end.time <- Sys.time()
  #    estimates[i, ] = estimate_chain[iterations+1, ]
  #    times[i] = end.time - start.time
  #  }
  # list(estimates, times)
}

# # running an experiment
# results = run_experiment(1, c(-2, 0.0042), 20, diag(c(0.004, 0.00012)))
# 
# # getting SGLD chain
# Y = sample_networks(Y=NULL, 1, c(-2, 0.0042), 20, 10^4)[[1]]
# theta_chain = sgld(Y, diag(c(0.004, 0.00012)), iterations=10000)
# 
# plot(theta_chain[, 1])
# plot(theta_chain[, 2])
# 
# # Trace Plots Sweep Test
# net <- network.initialize(50, directed=F)
# model_formula <- net ~ kstar(1:2)
# coefs <- c(-2, 0.0042)
# 
# sim_nets <- simulate(model_formula, coef=coefs, nsim=50, output="stats")
# 
# Y = matrix(0, nrow=50, ncol=50)
# current_sweep = Y
# theta_1_chain = c()
# theta_2_chain = c()
# for (i in 1:(10000))
# {
#   current_sweep = sweep(current_sweep, c(-2, 0.0042))
#   theta_1_chain = c(theta_1_chain, kstar(current_sweep, 1))
#   theta_2_chain = c(theta_2_chain, kstar(current_sweep, 2))
# }
# current_sweep


# Additional parameter choices
# epsilon_D <- diag(0.01 / E_k)
# epsilon_D <- diag(0.01 / estimate)

# tie-no-tie

# dirichlet prior
# normal

# Section 4.2
# data <- data.frame(n = numeric(), theta_1_bias_mean = numeric(),
#                    theta_1_bias_sd = numeric(),
#                    theta_2_bias_mean = numeric(), theta_2_bias_msd = numeric(),
#                    average_cpu_time = numeric(), m = numeric(),
#                    epsilon_D = character())
set.seed(123, kind = "L'Ecuyer-CMRG")
start.time <- Sys.time()
n <- 20
epsilon_D <- diag(c(0.004, 0.00012))
iterations = 1
space = c(1, 5)
data <- foreach(s = space, .packages=c('rpm'), .combine = 'rbind') %dorng% {
  foreach(alt = 0:1, .packages = c('rpm'), .combine = 'rbind') %dorng% {
    if (alt == 0) {
      eD = epsilon_D
    } else {
      eD = NULL
    }
    results = run_experiment(50, c(-2, 0.0042), n, eD, 
                             iterations=iterations, 
                             space=s)
    theta_bias <- base::sweep(results[,1:2], 2, c(-2, 0.0042))
    theta_1_bias_mean <- mean(theta_bias[,1])
    theta_1_bias_sd <- sd(theta_bias[,1])
    theta_2_bias_mean <- mean(theta_bias[,2])
    theta_2_bias_sd <- sd(theta_bias[,2])
    average_cpu_time <- mean(results[,3])
    new_row <- data.frame(n = n, theta_1_bias_mean = theta_1_bias_mean, 
                          theta_1_bias_sd = theta_1_bias_sd, 
                          theta_2_bias_mean = theta_2_bias_mean,
                          theta_2_bias_sd = theta_2_bias_sd, 
                          average_cpu_time = average_cpu_time, s = s, alt = alt,
                          stringsAsFactors = FALSE)
    new_row
  }
}
# colnames(data) = c("n", "theta_1_bias_mean", "theta_1_bias_sd", 
#                    "theta_2_bias_mean", "theta_2_bias_msd", 
#                    "average_cpu_time", "space", "epsilon_D")
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
