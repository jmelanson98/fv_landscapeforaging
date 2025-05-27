##### Accuracy comparison of different simulation/modeling scenarios #####
##### THIS SCRIPT IS FOR THE SIMULATION PHASE ######
# Script Initiated: May 26, 2025
# By: Jenna Melanson
# Goals:
### use (1) higher sample size, (2) accurate background colony density, (3) larger resource landscape to compare conditions:
##### Exponentiated quadratic (all colonies)
##### Exponential (all colonies)
##### exponentiated quadratic (all nonzero colonies)
##### exponential (all nonzero colonies)
##### exponentiated quadratic (doubleton+ colonies)
##### exponential (doubleton+ colonies)
##### centroid approach on exponential data (doubleton+)
##### centroid approach on exp. quadratic data (doubleton+)


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

##### Set Environment #####
#setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging") # server
options(mc.cores = parallel::detectCores())




##### Source functions #####
source("simulate_data/src/GeneralizedSimFunctions.R")


##### Create a grid of parameters to test #####
if (file.exists("simulate_data/methods_comparison/param_grid.rds")){
  param_grid = readRDS("simulate_data/methods_comparison/param_grid.rds")
} else{
  landscape_ids = 1:10
  rho <- c(50, 75, 100, 125, 150)
  colony_density <- c(8000) # using this because it's 4 * sample_size
  sample_sizes <- c(2000)
  distance_decay = c("exponentiated_quadratic", "exponential")
  model_approach = c("all", "singletons", "doubletons", "centroid")
  
  param_grid <- expand.grid(
    landscape_id = landscape_ids,
    rho = rho,
    sample_size = sample_sizes,
    distance_decay = distance_decay,
    model_approach = model_approach,
    stringsAsFactors = FALSE
  )
  param_grid$true_average_foraging = NA
  param_grid$true_sd_foraging = NA
  param_grid$model_average_foraging = NA
  param_grid$model_sd_foraging = NA
  param_grid$counts = NA
  param_grid$num_unobserved = NA
}


##### Simulate data ######
# Get task ID from SLURM environment variable
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
params <- param_grid[task_id, ]

# Get landscape from saved file
fq <- readRDS(paste0(sprintf("landscapes/random_field_range10_large/landscape_%03d", params$landscape_id), ".rds"))

# Run simulation
result <- draw_bees_colony_restricted(
  sample_size     = params$sample_size,
  landscape_size  = 1500,
  colonygrid_size = 700,
  trapgrid_size   = 300,
  number_traps    = 25,
  number_colonies = params$colony_density,
  colony_sizes    = rep(100, params$colony_density),
  rho            = params$rho,
  theta           = 0.5,
  resource_landscape = fq,
  nesting_landscape = NULL,
  distance_decay = params$distance_decay
)

# Save output
outfilepath <- sprintf("simulate_data/methods_comparison/data/sim_result_%03d", task_id)
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)
saveRDS(result[[1]], paste(outfilepath, "/yobs.RDS", sep = ""))
saveRDS(result[[2]], paste(outfilepath, "/colonydata.RDS", sep = ""))
saveRDS(result[[3]], paste(outfilepath, "/trapdata.RDS", sep = ""))

#save coloony metrics to param_grid
nonzero = yobs[rowSums(yobs) > 0,]
zero = yobs[rowSums(yobs) ==0,]
param_grid$counts = list(rowSums(nonzero))
param_grid$num_unobserved = nrow(zero)
param_grid$true_average_foraging = mean(result[[2]]$foraging_range)
param_grid$true_sd_foraging - sd(result[[2]]$foraging_range)
saveRDS(param_grid, "simulate_data/methods_comparison/param_grid.rds")
