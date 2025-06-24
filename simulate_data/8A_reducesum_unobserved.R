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
library(posterior)

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
posterior_draws_matrix <- as_draws_matrix(fit$draws())

# make a function to compute colony_dist for each draw
compute_colony_dist_summary <- function(draw_row,
                                        trap = data$trap,
                                        floral = data$floral,
                                        C = data$C,
                                        K = data$K) {
  # draw_row is a vector (one row of posterior_draws_matrix)
  
  # put delta and zeta in proper format
  delta <- matrix(NA, nrow = C, ncol = 2)
  zeta <- numeric(C)
  for (i in 1:C) {
    delta[i, 1] <- draw_row[paste0("delta[", i, ",1]")]
    delta[i, 2] <- draw_row[paste0("delta[", i, ",2]")]
    zeta[i] <- draw_row[paste0("zeta[", i, "]")]
  }
  
  # put epsilon in proper format
  eps <- numeric(K)
  for (k in 1:K) {
    eps[k] <- draw_row[paste0("eps[", k, "]")]
  }
  
  # get scalars
  rho <- draw_row["rho"]
  theta <- draw_row["theta"]
  mu <- draw_row["mu"]
  tau <- draw_row["tau"]
  sigma <- draw_row["sigma"]
  alpha <- 1e-12
  
  # hold space for temporary declarations
  dis <- matrix(NA, nrow = C, ncol = K)
  lambda <- matrix(NA, nrow = C, ncol = K)
  colony_dist <- numeric(C)
  
  # calculate dis and lambda
  for (k in 1:K) {
    for (i in 1:C) {
      dis[i, k] <- sqrt((delta[i, 1] - trap[k, 1])^2 + (delta[i, 2] - trap[k, 2])^2)
      lambda[i, k] <- dis[i, k] / (-rho * exp(theta * floral[k])) + mu +
        zeta[i] * sqrt(tau) + eps[k] * sqrt(sigma)
    }
  }
  
  # calculate per colony visitation
  V <- rowSums(exp(lambda))
  
  # calculate per colony foraging distance
  for (k in 1:K) {
    colony_dist <- colony_dist + (dis[, k] * exp(lambda[, k]) / (V + alpha))
  }
  
  # return mean foraging distances across colonies, for a single iteration
  return(c(mean = mean(colony_dist), sd = sd(colony_dist)))
}



### Apply the function in parallel across draws!
draws_per_chunk <- 100
total_draws <- nrow(posterior_draws_matrix)
chunk_starts <- seq(1, total_draws, by = draws_per_chunk)

# set up future backend
plan(multisession, workers = 8, gc = TRUE)

# loop over chunks
summary_stats_list <- list()
for (start in chunk_starts) {
  end <- min(start + draws_per_chunk - 1, total_draws)
  print(paste0("Processing draws ", start, " to ", end))
  
  # subset the matrix for this chunk
  chunk_draws <- posterior_draws_matrix[start:end, , drop = FALSE]
  
  # apply function in parallel to each row (draw)
  print("start lapplying")
  chunk_results <- future_lapply(1:nrow(chunk_draws), function(i) {
    compute_colony_dist_summary(
      draw_row = chunk_draws[i, ],
      trap = data$trap,
      floral = data$floral,
      C = data$C,
      K = data$K
    )
    print("done one draw")
  }
  )
  
  # combine results
  summary_stats_list[[length(summary_stats_list) + 1]] <- do.call(rbind, chunk_results)
}
print('done lapplying')

# combine all into one matrix or data frame
summary_stats_mat <- do.call(rbind, summary_stats_list)

print(paste("posterior mean = ", apply(summary_stats_mat, 2, mean), sep = ""))
print(paste("posterior sd = ", apply(summary_stats_mat, 2, sd), sep = ""))
print(paste("posterior CIs = ", apply(summary_stats_mat, 2, quantile, probs = c(0.025, 0.975)), sep = ""))


