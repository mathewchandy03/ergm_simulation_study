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
                             iters = 10)
                             # iters = 10^4)
  
  
  # estimate theta with SGLD
  estimates = matrix(ncol=2, nrow=n)
  times = rep(NA, n)
  for (i in 1:n)
  {
    start.time <- Sys.time()
    estimate_chain = sgld(networks[[i]], epsilon_D, iterations=iterations,
                          space=space)
    end.time <- Sys.time()
    estimates[i, ] = colMeans(estimate_chain[iterations/2:iterations+1,])
    times[i] = end.time - start.time
  }
 list(estimates, times)
}

my_local_experiment <- function() {
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
                     "theta_2_bias_mean", "theta_2_bias_sd", 
                     "average_cpu_time", "space", "epsilon_D")
  data
}

# run experiment for single node size in parallel
results_for_node <- function(n, epsilon_D, iterations) {
  results = foreach(s=c(1,5), .packages=c('rpm')) %dorng% {
    # loop over number of sweeps between samples (space)
    foreach(x=0:1, .packages=c('rpm')) %dorng% {
      # loop over which epsilon D value to use (provided or NULL)
      if (x == 0) {
        epD = epsilon_D
      }
      else {
        epD = NULL
      }
      results = run_experiment(5, c(-2, 0.0042), n, epsilon_D=epD, 
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

registerDoParallel(cores=8)
doRNG::registerDoRNG(123)
set.seed(20, kind = "L'Ecuyer-CMRG")

data20 = results_for_node(20, diag(c(0.004, 0.00012)), 10)
# data20 = results_for_node(20, diag(c(0.004, 0.00012)), 10000) 

saveRDS(data20, file="../data/data20.Rds")
# readRDS("../data/data20.Rds")

set.seed(50, kind = "L'Ecuyer-CMRG")

data50 = results_for_node(50, diag(c(0.004, 0.00012)), 10)
# data50 = results_for_node(50, diag(c(0.004, 0.00012)), 10000)

saveRDS(data50, file="../data/data50.Rds")
# readRDS("../data/data50.Rds")

set.seed(75, kind = "L'Ecuyer-CMRG")

data75 = results_for_node(75, diag(c(0.004, 0.00012)), 10)
# data75 = results_for_node(75, diag(c(0.004, 0.00012)), 10000)

saveRDS(data75, file="../data/data75.Rds")
# readRDS("../data/data75.Rds")

set.seed(100, kind = "L'Ecuyer-CMRG")

data75 = results_for_node(75, diag(c(0.004, 0.00012)), 10)
# data100 = results_for_node(100, diag(c(0.004, 0.00012)), 10000)

saveRDS(data100, file="../data/data100.Rds")
# readRDS("../data/data100.Rds")

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
  estimate_chain = sgld(networks[[i]], diag(c(0.004, 0.00012)), 
                        iterations=iterations)
  end.time <- Sys.time()
  estimates[i, ] = estimate_chain[iterations+1, ]
  times[i] = end.time - start.time
  plot(estimate_chain, main=i, xlab="theta1", ylab="theta2")
  points(c(estimate_chain[1, 1], estimate_chain[iterations+1, 1]), 
         c(estimate_chain[1, 2], estimate_chain[iterations+1, 2]),
         col="red")
}


