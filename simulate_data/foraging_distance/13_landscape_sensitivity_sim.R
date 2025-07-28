##### How well can we estimate landscape effects? #####
##### THIS SCRIPT IS FOR THE DATA SIMULATION PHASE ######
# Script Initiated: July 23, 2025
# By: Jenna Melanson
# Goals:
### test the ability of a trap-centric Bayesian model and regression-based centroid approach to
### detect effects of landscapes at (1) different buffer sizes (e.g., scales of landscape 
### spatial autocorrelation), (2) different effect sizes (alpha), and (3) across foraging distances


##### Load packages #####
library(rstan)
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
library(terra)
library(filelock)
library(tibble)

##### Set Environment #####
#setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging") # server
options(mc.cores = parallel::detectCores())

##### Source functions #####
source("simulate_data/src/GeneralizedSimFunctions.R")


##### Create a grid of parameters to test #####
# landscape_ids = 1:10
# autocorrelation = c(200, 300, 400)
# rho = c(20, 50, 100)
# alpha = c(0.2, 0.3, 0.4, 0.5)
# colony_density = c(7000)
# sample_sizes = c(2000)
# distance_decay = c("exponential")
# 
# param_grid <- expand.grid(
#     landscape_id = landscape_ids,
#     autocorrelation = autocorrelation,
#     rho = rho,
#     alpha = alpha,
#     colony_density = colony_density,
#     sample_size = sample_sizes,
#     distance_decay = distance_decay,
#     stringsAsFactors = FALSE
#   )
# param_grid$task_id = 1:nrow(param_grid)
# param_grid$true_average_foraging = NA
# param_grid$true_sd_foraging = NA
# param_grid$model_average_foraging = NA
# param_grid$model_sd_foraging = NA
# param_grid$model_mu = NA
# param_grid$model_rho = NA
# param_grid$model_alpha = NA
# param_grid$model_alpha_min = NA
# param_grid$model_alpha_max = NA
# param_grid$mean_center_dist = NA
# param_grid$mean_center_alpha = NA
# param_grid$mean_center_alpha_sd = NA
# 
# saveRDS(param_grid, "simulate_data/foraging_distance/methods_comparison/landscape_effects/param_grid.rds")
# 

print("Reading param grid.")
param_grid = readRDS("simulate_data/foraging_distance/methods_comparison/landscape_effects/param_grid.rds")

##### Simulate data ######
# Get task ID from SLURM environment variable
print("Setting task id and params.")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
params <- param_grid[param_grid$task_id == task_id, ]

# Get floral resource landscape from saved file
print("Loading floral resource raster.")
fq = readRDS(paste0(sprintf("simulate_data/landscapes/random_field_range10/landscape_%03d", params$landscape_id), ".rds"))

# Get landscape metrics from saved file
print("Loading landscape metric raster.")
repo_name = paste0("simulate_data/landscapes/random_field_range", params$autocorrelation)
landscape_name = sprintf("/landscape_%03d", params$landscape_id)
lmq = readRDS(paste0(repo_name, landscape_name, ".rds"))


# Run simulation
print("Starting simulation.")
result <- draw_bees_colony_restricted(
  sample_size     = params$sample_size,
  landscape_size  = 1500,
  colonygrid_size = 700,
  trapgrid_size   = 300,
  number_traps    = 25,
  number_colonies = params$colony_density,
  colony_sizes    = rep(100, params$colony_density),
  rho            = params$rho,
  alpha          = params$alpha,
  theta           = 0.5,
  resource_landscape = fq,
  configuration = lmq,
  nesting_landscape = NULL,
  distance_decay = params$distance_decay
)


# Write outputs to variables
yobs = result[[1]]
colony_data = result[[2]]
trap_data = result[[3]]

# save main output
outfilepath <- sprintf("simulate_data/foraging_distance/methods_comparison/landscape_effects/data/sim_result_%03d", task_id)
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)

saveRDS(yobs, paste(outfilepath, "/yobs.RDS", sep = ""))
saveRDS(colony_data, paste(outfilepath, "/colonydata.RDS", sep = ""))
saveRDS(trap_data, paste(outfilepath, "/trapdata.RDS", sep = ""))

#save colony metrics / summary statistics to single output file
params$true_average_foraging = mean(colony_data$foraging_range)
params$true_sd_foraging = sd(colony_data$foraging_range)

result = tibble(params)

# use lock to make sure multiple tasks don't write to output at the same time
lockfile <- "output.lock"
rds_file <- "simulate_data/foraging_distance/methods_comparison/landscape_effects/output_sim.rds"

# try to acquire the lock (waits up to 60 seconds)
lock <- lock(lockfile, timeout = 60000)

if (!is.null(lock)) {
  # If the RDS file exists, read it in; else create new
  if (file.exists(rds_file)) {
    df <- readRDS(rds_file)
  } else {
    df <- tibble()  # initialize empty tibble
  }
  
  # append the new row
  df <- dplyr::bind_rows(df, result)
  
  # Save the updated dataframe
  tmp <- tempfile()
  saveRDS(df, tmp)
  file.copy(tmp, rds_file, overwrite = TRUE)
  file.remove(tmp)
  
  # Release lock
  unlock(lock)
} else {
  stop("Could not acquire lock on output file.")
}

