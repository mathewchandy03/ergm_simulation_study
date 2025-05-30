# Running Experiment

source("./gibbs_sampling.R")
source("./sgld.R")

library(iterators)
library(parallel)
library(foreach)
library(doParallel)

run_experiment <- function(n, theta, nodes, epsilon_D, iterations=10000, 
                           space=1)
{
  # simulate n networks
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
    estimate_chain = sgld(networks[[i]], epsilon_D, iterations=iterations,
                          space=space)
    end.time <- Sys.time()
    # estimates[i, ] = estimate_chain[iterations+1, ]
    estimates[i, ] = colMeans(estimate_chain[iterations/2:iterations+1,])
    times[i] = end.time - start.time
  }
 list(estimates, times)
}

# running an experiment
results = run_experiment(1, c(-2, 0.0042), 20, diag(c(0.004, 0.00012)))

# getting SGLD chain
Y = sample_networks(Y=NULL, 1, c(-2, 0.0042), 20, 10^4)[[1]]
theta_chain = sgld(Y, diag(c(0.004, 0.00012)), iterations=10000)

plot(theta_chain[, 1])
plot(theta_chain[, 2])

# Trace Plots Sweep Test
net <- network.initialize(50, directed=F)
model_formula <- net ~ kstar(1:2)
coefs <- c(-2, 0.0042)

sim_nets <- simulate(model_formula, coef=coefs, nsim=50, output="stats")

Y = matrix(0, nrow=50, ncol=50)
current_sweep = Y
theta_1_chain = c()
theta_2_chain = c()
for (i in 1:(10000))
{
  current_sweep = sweep(current_sweep, c(-2, 0.0042))
  theta_1_chain = c(theta_1_chain, kstar(current_sweep, 1))
  theta_2_chain = c(theta_2_chain, kstar(current_sweep, 2))
}
current_sweep


# Additional parameter choices
# epsilon_D <- diag(0.01 / E_k)
# epsilon_D <- diag(0.01 / estimate)

# tie-no-tie

# dirichlet prior
# normal

# Section 4.2
data <- data.frame(n = numeric(), theta_1_bias_mean = numeric(),
                   theta_1_bias_sd = numeric(),
                   theta_2_bias_mean = numeric(), theta_2_bias_msd = numeric(),
                   average_cpu_time = numeric(), space = numeric(),
                   epsilon_D = character())

n <- c(20, 50, 75, 100)
epsilon_D <- list(diag(c(0.004, 0.00012)), diag(c(5e-5, 5e-6)),
                  diag(c(5.64e-6, 5.64e-7)), diag(c(3.08e-6, 2.38e-7)))
iterations = 10000
space = c(1, 5)
for (s in space) {
  print(sprintf("space: %d", s))
  for (i in 1:4) {
    results = run_experiment(50, c(-2, 0.0042), n[i], epsilon_D[[i]], 
                             iterations=iterations, 
                             space=s)
    theta_bias <- base::sweep(results[[1]], 2, c(-2, 0.0042))
    theta_1_bias_mean <- mean(theta_bias[,1])
    theta_1_bias_sd <- sd(theta_bias[,1])
    theta_2_bias_mean <- mean(theta_bias[,2])
    theta_2_bias_sd <- sd(theta_bias[,2])
    average_cpu_time <- mean(results[[2]])
    new_row <- c(n[i], theta_1_bias_mean, theta_1_bias_sd, theta_2_bias_mean,
                 theta_2_bias_sd, average_cpu_time, s, "original")
    data <- rbind(data, new_row)
    print(sprintf("original epsilon_D, %d nodes, time: %g" , n[i], 
                  average_cpu_time))
  }
  for (i in 1:4) {
    results = run_experiment(50, c(-2, 0.0042), n[i], epsilon_D=NULL, 
                             iterations=iterations, space=s)
    theta_bias <- base::sweep(results[[1]], 2, c(-2, 0.0042))
    theta_1_bias_mean <- mean(theta_bias[,1])
    theta_1_bias_sd <- sd(theta_bias[,1])
    theta_2_bias_mean <- mean(theta_bias[,2])
    theta_2_bias_sd <- sd(theta_bias[,2])
    average_cpu_time <- mean(results[[2]])
    new_row <- c(n[i], theta_1_bias_mean, theta_1_bias_sd, theta_2_bias_mean,
                 theta_2_bias_sd, average_cpu_time, s, "null")
    data <- rbind(data, new_row)
    print(sprintf("null epsilon_D, %d nodes, time: %g" , n[i], 
                  average_cpu_time))
  }
}
colnames(data) = c("n", "theta_1_bias_mean", "theta_1_bias_sd", 
                   "theta_2_bias_mean", "theta_2_bias_msd", 
                   "average_cpu_time", "space", "epsilon_D")


# parallel
arguments = list()
for (s in space) {
  for (i in 1:4) {
    curr_args = list(n=n[i], epsilon_D=epsilon_D[[i]], space=s, 
                     iterations=iterations)
    arguments = list.append(curr_args)
  }
  for (i in 1:4) {
    curr_args = list(n=n[i], epsilon_D=NULL, space=s, 
                     iterations=iterations)
    arguments = list.append(curr_args)
  }
}

n_experiments = length(arguments)
data_par = matrix(nrow=n_experiments, ncol=8)

data_par <- foreach(1:n_experiments) %dopar% function(x) {
  curr_args = arguments[[x]]
  results = run_experiment(50, c(-2, 0.0042), n=curr_args$n,
                           epsilon_D = curr_args$epsilon_D, 
                           iterations=curr_args$iterations, 
                           space=curr_args$space)
  theta_bias <- base::sweep(results[[1]], 2, c(-2, 0.0042))
  theta_1_bias_mean <- mean(theta_bias[,1])
  theta_1_bias_sd <- sd(theta_bias[,1])
  theta_2_bias_mean <- mean(theta_bias[,2])
  theta_2_bias_sd <- sd(theta_bias[,2])
  average_cpu_time <- mean(results[[2]])
  new_row <- c(curr_args$n, theta_1_bias_mean, theta_1_bias_sd, 
               theta_2_bias_mean,theta_2_bias_sd, average_cpu_time, 
               curr_args$space, is.null(curr_args$epsilon_D))
  new_row
}

registerDoParallel(cores=8)
# single node size, parallel
results_for_node <- function(n, epsilon_D, iterations) {
  results = foreach(s=c(1,5)) %dopar% {
    # loop over number of sweeps between samples (space)
    foreach(x=0:1) %dopar% {
      # loop over which epsilon D value to use (default or NULL)
      if (x == 0) {
        epD = epsilon_D
      }
      else {
        epD = NULL
      }
      results = run_experiment(50, c(-2, 0.0042), n, epsilon_D=epD, 
                               iterations=iterations, 
                               space=s)
      theta_bias <- base::sweep(results[[1]], 2, c(-2, 0.0042))
      theta_1_bias_mean <- mean(theta_bias[,1])
      theta_1_bias_sd <- sd(theta_bias[,1])
      theta_2_bias_mean <- mean(theta_bias[,2])
      theta_2_bias_sd <- sd(theta_bias[,2])
      average_cpu_time <- mean(results[[2]])
      new_row <- c(n, theta_1_bias_mean, theta_1_bias_sd, theta_2_bias_mean,
                   theta_2_bias_sd, average_cpu_time, s, x)
    }
  }
  results = matrix(unlist(unlist(results, recursive=F)), nrow=4, byrow=T)
  colnames(results) = c("n", "theta_1_bias_mean", "theta_1_bias_sd", 
                     "theta_2_bias_mean", "theta_2_bias_msd", 
                     "average_cpu_time", "space", "epsilon_D_null")
  results
}

data20 = results_for_node(20, diag(c(0.004, 0.00012)), 10000)

# Creating Plots
n= 2
iterations=1000
networks = sample_networks(m = n,
                           theta = c(-2, 0.0042),
                           nodes = 50,
                           iters = 10^4)

estimates = matrix(ncol=2, nrow=n)
times = rep(NA, n)
for (i in 1:n)
{
  start.time <- Sys.time()
  estimate_chain = sgld(networks[[i]], diag(c(0.004, 0.00012)), iterations=iterations)
  end.time <- Sys.time()
  estimates[i, ] = estimate_chain[iterations+1, ]
  times[i] = end.time - start.time
  plot(estimate_chain, main=i, xlab="theta1", ylab="theta2")
  points(c(estimate_chain[1, 1], estimate_chain[iterations+1, 1]), 
         c(estimate_chain[1, 2], estimate_chain[iterations+1, 2]),
         col="red")
}
