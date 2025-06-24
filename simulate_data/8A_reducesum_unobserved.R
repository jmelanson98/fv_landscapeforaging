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

saveRDS(fit, paste("simulate_data/methods_comparison/observed_vs_unobserved/reducesum/", current$thread, current$model, "Fit.rds", sep = ""))



###### Post hoc calculations of colony_dist

# function for computing colony dist!
compute_colony_dist_summary <- function(draw,
                                        trap = data$trap,
                                        floral = data$floral,
                                        C = data$C,
                                        K = data$K) {
  # reconstruct delta, zeta
  delta = matrix(NA, nrow = C, ncol = 2)
  zeta = matrix(NA, nrow = C)
  for (i in 1:C) {
    delta[i, 1] <- draw[[paste0("delta[", i, ",1]")]]
    delta[i, 2] <- draw[[paste0("delta[", i, ",2]")]]
    zeta[i] = draw[[paste0("zeta[", i, "]")]]
  }
  
  # reconstruct eps
  eps = matrix(NA, nrow = K)
  for (i in 1:K) {
    eps[k] = draw[[paste0("eps[", i, "]")]]
  }
  
  # pull out other vars
  rho <- draw$rho
  theta <- draw$theta
  mu <- draw$mu
  tau <- draw$tau
  sigma <- draw$sigma
  alpha <- 1e-12

  # temporary declarations
  dis <- matrix(NA, nrow = C, ncol = K)
  lambda <- matrix(NA, nrow = C, ncol = K)
  V <- numeric(C)
  colony_dist <- numeric(C)
  
  # calculate distance and lambda
  for(k in 1:K) {
    for(i in 1:C) {
      dis[i, k] <- sqrt((delta[i,1] - trap[k,1])^2 + (delta[i,2] - trap[k,2])^2)
      lambda[i, k] <- dis[i,k]/(-rho * exp(theta * floral[k])) + mu + zeta[i]*sqrt(tau) + eps[k]*sqrt(sigma)
    }
  }
  
  # calculate total visitation per colony
  V <- rowSums(exp(lambda))
  
  # calculate colony_dist
  for(k in 1:K) {
    colony_dist <- colony_dist + (dis[,k] * exp(lambda[,k]) / (V + alpha))
  }
  
  # return mean and sd of colony_dist (averaged across colonies, for one iteration)
  return(c(mean = mean(colony_dist), sd = sd(colony_dist)))
}


# run in parallel over many draws!
# set up
plan(multisession, workers = 16)

# get posterior draws
posterior_draws <- fit$draws(format = "draws_list")
print(object.size(posterior_draws), units = "auto")

# apply function and summarize (loop over chunks)
print('lapplying')

chunk_size <- 100
n_draws <- length(posterior_draws)
n_chunks <- ceiling(n_draws / chunk_size)
print(paste0("ndraws = ", n_draws))


summary_stats_list <- list()

for (i in seq_len(n_chunks)) {
  start <- (i - 1) * chunk_size + 1
  end <- min(i * chunk_size, n_draws)

  message("Processing draws ", start, " to ", end)

  chunk_draws <- posterior_draws[start:end]
  print(length(chunk_draws))
  
  chunk_summaries <- future_lapply(chunk_draws, compute_colony_dist_summary, trap=data$trap,
  floral=data$floral,C=data$C, K=data$K, future.globals = FALSE)

  summary_stats_list[[i]] <- do.call(rbind, chunk_summaries)
}

summary_stats_mat <- do.call(rbind, summary_stats_list)

print('done lapplying')

print(paste("posterior mean = ", apply(summary_stats_mat, 2, mean), sep = ""))
print(paste("posterior sd = ", apply(summary_stats_mat, 2, sd), sep = ""))
print(paste("posterior CIs = ", apply(summary_stats_mat, 2, quantile, probs = c(0.025, 0.975)), sep = ""))


