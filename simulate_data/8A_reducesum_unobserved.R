##### Run models with and without unobserved colonies #####
# Script Initiated: June 17, 2025
# By: J. Melanson
# Goals:
### Test some methods for parallelizing stan code on the server with cmdstan

##### Load packages #####
library(cmdstanr)
library(matrixStats)
library(sp)
library(gstat)
library(ggplot2)
library(reshape2)
library(raster)
library(rasterVis)
library(parallel)
library(future)
library(furrr)
library(dplyr)
library(tidyr)
library(gridExtra)
library(tibble)
library(future.apply)

##### Set Environment #####
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
set_cmdstan_path("/home/melanson/projects/def-ckremen/melanson/cmdstan")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))


##### Load in some simulation data!! #####

## NOTE: this is simulation 74!! this data is for landscape = 4, rho = 50, exponential decay
yobs = readRDS("simulate_data/methods_comparison/data/sim_result_074/yobs.RDS")
colony_data = readRDS("simulate_data/methods_comparison/data/sim_result_074/colonydata.RDS")
trap_data = readRDS("simulate_data/methods_comparison/data/sim_result_074/trapdata.RDS")

ydub = yobs[rowSums(yobs) > 1,]
# Prep data list for Stan
data = list()
data$y = ydub
data$C = nrow(data$y)
data$K = ncol(data$y)
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$lowerbound = 400
data$upperbound = 1100
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$rho_center = 4.5
data$rho_sd = 0.5

# Table of conditions
model = c("reduce_sum_exponential.stan", "reduce_sum_exponentialGQ.stan")
threads = c(4, 8)
grid = expand.grid(model = model,
                 threads = threads)
grid$task = 1:nrow(grid)

current = grid[grid$task == task_id,]

print(current$model)
print(current$threads)

#select stan model to fit
mod_file <- paste("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging/models/", current$model, sep = "")
mod <- cmdstan_model(mod_file, cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)


# Run sampling
threads_per_chain = current$threads

# add grainsize to data list
grainsize <- max(floor(data$C / (threads_per_chain * 5)), 1)
data$grainsize = grainsize


#fit and save model
print("Starting sampling.")
fit <- mod$sample(
  data = data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = threads_per_chain,
  refresh = 100,
  iter_warmup = 1000,
  iter_sampling = 5000
)

saveRDS(fit, paste("methods_comparison/observed_vs_unobserved/reducesum/", current$thread, current$model, "Fit.rds", sep = ""))



###### Post hoc calculations of colony_dist

# function for computing colony dist!
compute_colony_dist_summary <- function(draw) {
  delta <- draw$delta
  rho <- draw$rho
  theta <- draw$theta
  mu <- draw$mu
  zeta <- draw$zeta
  eps <- draw$eps
  tau <- draw$tau
  sigma <- draw$sigma
  
  alpha <- 1e-12
  C <- nrow(delta)
  K <- nrow(trap)
  
  dis <- matrix(NA, nrow = C, ncol = K)
  lambda <- matrix(NA, nrow = C, ncol = K)
  V <- numeric(C)
  colony_dist <- numeric(C)
  
  for(k in 1:K) {
    for(i in 1:C) {
      dis[i, k] <- sqrt((delta[i,1] - trap[k,1])^2 + (delta[i,2] - trap[k,2])^2)
      lambda[i, k] <- dis[i,k]/(-rho * exp(theta * floral[k])) + mu + zeta[i]*sqrt(tau) + eps[k]*sqrt(sigma)
    }
  }
  
  V <- rowSums(exp(lambda))
  
  for(k in 1:K) {
    colony_dist <- colony_dist + (dis[,k] * exp(lambda[,k]) / (V + alpha))
  }
  
  # return mean and sd of colony_dist (averaged across colonies, for one iteration)
  return(c(mean = mean(colony_dist), sd = sd(colony_dist)))
}


# run in parallel over many draws!
plan(multisession, workers = 4)
posterior_draws <- fit$draws(variables = c("delta", "rho", "theta", "mu", "zeta", "eps", "tau", "sigma"),
                             format = "draws_list")


# set some parameters from R
C <- data$C
K <- data$K

summary_stats_list <- future_lapply(posterior_draws, compute_colony_dist_summary)
summary_stats_mat <- do.call(rbind, summary_stats_list)


print(paste("posterior mean = ", apply(summary_stats_mat, 2, mean), sep = ""))
print(paste("posterior sd = ", apply(summary_stats_mat, 2, sd), sep = ""))
print(paste("posterior CIs = ", apply(summary_stats_mat, 2, quantile, probs = c(0.025, 0.975)), sep = ""))


