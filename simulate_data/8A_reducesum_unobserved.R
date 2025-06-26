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
library(filelock)

##### Set Environment #####
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
set_cmdstan_path("/home/melanson/projects/def-ckremen/melanson/cmdstan")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))


##### Load in some simulation data!! #####
param_grid = readRDS("simulate_data/methods_comparison/param_grid.rds")
inputfilepath <- sprintf("simulate_data/methods_comparison/data/sim_result_%03d", task_id)
yobs = readRDS(paste(inputfilepath, "/yobs.RDS", sep = ""))
colony_data = readRDS(paste(inputfilepath, "/colonydata.RDS", sep = ""))
trap_data = readRDS(paste(inputfilepath, "/trapdata.RDS", sep = ""))
current_params = param_grid[param_grid$task_id == task_id,]


# Prep data list for Stan
data = list()

# how much of simulation data to keep?
if (current_params$model_approach == "all"){
  data$y = yobs
} else if (current_params$model_approach == "singletons"){
  data$y = yobs[rowSums(yobs) >0,]
} else if (current_params$model_approach == "doubletons"){
  data$y = yobs[rowSums(yobs) >1,]
}

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

#select stan model to fit
if (current_params$distance_decay == "exponential"){
  mod_file <- "/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging/models/reduce_sum_exponential.stan"
} else if (current_params$distance_decay == "exponentiated_quadratic"){
  mod_file <- "/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging/models/reduce_sum_expquad.stan"
}
mod <- cmdstan_model(mod_file, cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)


# add grainsize to data list
threads_per_chain = 4
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
  iter_sampling = 5000,
  init = 1
)


###### Post hoc calculations of colony_dist
posterior_draws_matrix <- as_draws_matrix(fit$draws())

# make a function to compute colony_dist for each draw
compute_colony_dist_summary <- function(draw_row,
                                        trap,
                                        floral,
                                        C,
                                        K) {
  # draw_row is a draws matrix with one row
  
  # put delta and zeta in proper format
  delta <- matrix(NA, nrow = C, ncol = 2)
  zeta <- rep(NA,C)
  for (i in 1:C) {
    delta1 <- paste0("delta[", i, ",1]")
    delta2 <- paste0("delta[", i, ",2]")
    delta[i,1] = draw_row[,delta1]
    delta[i,2] = draw_row[,delta2]
    zeta[i] <- draw_row[,paste0("zeta[", i, "]")]
  }
  
  # put epsilon in proper format
  eps <- rep(NA,K)
  for (k in 1:K) {
    eps[k] <- draw_row[,paste0("eps[", k, "]")]
  }
  
  # get scalars
  rho <- draw_row[,"rho"]
  theta <- draw_row[,"theta"]
  mu <- draw_row[,"mu"]
  tau <- draw_row[,"tau"]
  sigma <- draw_row[,"sigma"]
  alpha <- 1e-12

  # hold space for temporary declarations
  dis <- matrix(NA, nrow = C, ncol = K)
  lambda <- matrix(NA, nrow = C, ncol = K)
  colony_dist <- rep(0,C)
  
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
  chunk_draws <- posterior_draws_matrix[start:end,]
  
  # apply function in parallel to each row (draw)
  chunk_start_time = Sys.time()
  chunk_results <- future_lapply(1:nrow(chunk_draws), function(j) {
    compute_colony_dist_summary(
      draw_row = chunk_draws[j, ,drop = FALSE ],
      trap = data$trap,
      floral = data$floral,
      C = data$C,
      K = data$K
    )
  }
  )
  chunk_end_time = Sys.time()
  print(paste0("Chunk time = ", chunk_end_time-chunk_start_time))
  
  # combine results
  summary_stats_list[[length(summary_stats_list) + 1]] <- do.call(rbind, chunk_results)
}
print('done lapplying')

# combine all into one matrix or data frame
summary_stats_mat <- as.data.frame(do.call(rbind, summary_stats_list))

#####  Save colony_dist estimates ######
# use lock to make sure multiple tasks don't write to output at the same time
lockfile <- "output.lock"
rds_file <- "simulate_data/methods_comparison/output.rds"

  
# save the result of *this* simulation only
print("Saving current params")
current_params$model_average_foraging = mean(summary_stats_mat$mean)
current_params$model_sd_foraging = sqrt(mean(summary_stats_mat$sd^2))
current_params$model_mu = fit$summary(variables = "mu")$mean
current_params$model_rho = fit$summary(variables = "rho")$mean
saveRDS(current_params, paste(inputfilepath, "/iteration_output.rds", sep = ""))
  
# save output
# try to acquire the lock (waits up to 60 seconds)
lock <- lock(lockfile, timeout = 60000)

if (!is.null(lock)) {
    df <- readRDS(rds_file)
    
    print("Saving params to output.rds")
    df$model_average_foraging[df$task_id == task_id] = mean(summary_stats_mat$mean)
    df$model_sd_foraging[df$task_id == task_id] = sqrt(mean(summary_stats_mat$sd^2))
    df$model_mu[df$task_id == task_id] = fit$summary(variables = "mu")$mean
    df$model_mu_sd[df$task_id == task_id] = fit$summary(variables = "mu")$sd
    df$model_rho[df$task_id == task_id] = fit$summary(variables = "rho")$mean
    df$model_rho_sd[df$task_id == task_id] = fit$summary(variables = "rho")$sd
    
    # save updated dataframe
    saveRDS(df, rds_file)
    
    # Release lock
    unlock(lock)
  } else {
    stop("Could not acquire lock on output file.")
  }



