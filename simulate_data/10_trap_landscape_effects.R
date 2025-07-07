##### Test accuracy of landscape effects #####
# Script Initiated: July 3, 2025
# By: Jenna Melanson
# Goals:
### Check whether models can accurately return estimates of landscape effects on foraging distance;
### assuming that landscape values are spatially autocorrelated and that we give the model a value for each trap


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
library(posterior)
library(terra)
#library(ggpubr)

##### Set Environment #####
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
#setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
set_cmdstan_path("/home/melanson/projects/def-ckremen/melanson/cmdstan/cmdstan-2.36.0")

##### Source functions #####
source("simulate_data/src/GeneralizedSimFunctions.R")

##### Simulate landscape characteristics "landscape" #####
landscape_char = simulateLandscapeRaster(landscape_size = 1500, resource_range = 200)
landscape_char = terra::rast(landscape_char)

##### Load in floral quality landscape #####
fq <- readRDS("simulate_data/landscapes/landscapes/random_field_range10/landscape_001.rds")
fq = terra::rast(fq)

##### Simulate some data with landscape effects on foraging distance #####
result <- draw_bees_colony_restricted(
  sample_size     = 1000,
  landscape_size  = 1500,
  colonygrid_size = 700,
  trapgrid_size   = 300,
  number_traps    = 25,
  number_colonies = 2000,
  colony_sizes    = rep(100, 2000),
  rho            = 50,
  theta           = 0.2,
  alpha = 0.2,
  resource_landscape = fq,
  configuration = landscape_char,
  nesting_landscape = NULL,
  distance_decay = "exponential"
)

# Save results
saveRDS(result, "simulate_data/methods_comparison/landscape_effects/traplevel_simdata.rds")
#result = readRDS("simulate_data/methods_comparison/landscape_effects/traplevel_simdata.rds")

# Write outputs to variables
yobs = result[[1]]
colony_data = result[[2]]
trap_data = result[[3]]

# # First try for only detected colonies
# yobs_detected = yobs[rowSums(yobs) >0,]
# colony_data_detected = colony_data[rowSums(yobs) > 0,]

# Prep data list for Stan
data = list()
data$y = yobs
data$C = nrow(data$y)
data$K = ncol(data$y)
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$lowerbound = 400
data$upperbound = 1100
data$landscape = trap_data$landscape_metric
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$rho_center = 3.5
data$rho_sd = 0.5

# add grainsize to data list
threads_per_chain = 4
grainsize <- max(floor(data$C / (threads_per_chain * 5)), 1)
data$grainsize = grainsize

# compile model
mod <- cmdstan_model(
  "/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging/models/rs_landscape_traplevel.stan",
  force_recompile = TRUE, cpp_options = list(stan_threads = TRUE)
)

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


saveRDS(fit, "simulate_data/methods_comparison/landscape_effects/lanscape_all_traplevel.rds")
posterior <- fit$draws(format = "df")
write.csv(posterior, "simulate_data/methods_comparison/landscape_effects/traplevel_landscape_all_draws_.csv", row.names = FALSE)

###### Post hoc calculations of colony_dist
# posterior_draws_matrix <- as_draws_matrix(fit$draws())
# 
# # make a function to compute colony_dist for each draw
# compute_colony_dist_summary <- function(draw_row,
#                                         trap,
#                                         floral,
#                                         C,
#                                         K) {
#   # draw_row is a draws matrix with one row
#   
#   # put delta and zeta in proper format
#   delta <- matrix(NA, nrow = C, ncol = 2)
#   zeta <- rep(NA,C)
#   for (i in 1:C) {
#     delta1 <- paste0("delta[", i, ",1]")
#     delta2 <- paste0("delta[", i, ",2]")
#     delta[i,1] = draw_row[,delta1]
#     delta[i,2] = draw_row[,delta2]
#     zeta[i] <- draw_row[,paste0("zeta[", i, "]")]
#   }
#   
#   # put epsilon in proper format
#   eps <- rep(NA,K)
#   for (k in 1:K) {
#     eps[k] <- draw_row[,paste0("eps[", k, "]")]
#   }
#   
#   # get scalars
#   rho <- draw_row[,"rho"]
#   theta <- draw_row[,"theta"]
#   tau <- draw_row[,"tau"]
#   sigma <- draw_row[,"sigma"]
#   alpha <- 1e-12
#   
#   # hold space for temporary declarations
#   dis <- matrix(NA, nrow = C, ncol = K)
#   lambda <- matrix(NA, nrow = C, ncol = K)
#   colony_dist <- rep(0,C)
#   
#   # calculate dis and lambda
#   for (k in 1:K) {
#     for (i in 1:C) {
#       dis[i, k] <- sqrt((delta[i, 1] - trap[k, 1])^2 + (delta[i, 2] - trap[k, 2])^2)
#       lambda[i, k] <- dis[i, k] / (-rho * exp(theta * floral[k])) +
#         zeta[i] * sqrt(tau) + eps[k] * sqrt(sigma)
#     }
#   }
#   
#   # calculate per colony visitation
#   V <- rowSums(exp(lambda))
#   
#   # calculate per colony foraging distance
#   for (k in 1:K) {
#     colony_dist <- colony_dist + (dis[, k] * exp(lambda[, k]) / (V + alpha))
#   }
#   
#   # return mean foraging distances across colonies, for a single iteration
#   return(c(mean = mean(colony_dist), sd = sd(colony_dist)))
# }
# 
# 
# 
# ### Apply the function in parallel across draws!
# draws_per_chunk <- 100
# total_draws <- nrow(posterior_draws_matrix)
# chunk_starts <- seq(1, total_draws, by = draws_per_chunk)
# 
# # set up future backend
# plan(multisession, workers = 8, gc = TRUE)
# 
# # loop over chunks
# summary_stats_list <- list()
# for (start in chunk_starts) {
#   end <- min(start + draws_per_chunk - 1, total_draws)
#   print(paste0("Processing draws ", start, " to ", end))
#   
#   # subset the matrix for this chunk
#   chunk_draws <- posterior_draws_matrix[start:end,]
#   
#   # apply function in parallel to each row (draw)
#   chunk_start_time = Sys.time()
#   chunk_results <- future_lapply(1:nrow(chunk_draws), function(j) {
#     compute_colony_dist_summary(
#       draw_row = chunk_draws[j, ,drop = FALSE ],
#       trap = data$trap,
#       floral = data$floral,
#       C = data$C,
#       K = data$K
#     )
#   }
#   )
#   chunk_end_time = Sys.time()
#   print(paste0("Chunk time = ", chunk_end_time-chunk_start_time))
#   
#   # combine results
#   summary_stats_list[[length(summary_stats_list) + 1]] <- do.call(rbind, chunk_results)
# }
# print('done lapplying')
# 
# # combine all into one matrix or data frame
# summary_stats_mat <- as.data.frame(do.call(rbind, summary_stats_list))
# write.csv(summary_stats_mat, "simulate_data/methods_comparison/landscape_effects/all_summary_stats.csv", row.names = FALSE)
# 
# 
