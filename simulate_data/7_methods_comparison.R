##### Accuracy comparison of different simulation/modeling scenarios #####
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

##### Set Environment #####
#setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging") # server
options(mc.cores = parallel::detectCores())

##### Source functions #####
source("simulate_data/src/GeneralizedSimFunctions.R")


##### Create a grid of parameters to test #####
##### Prepare to run in parallel ! #####
# set parameter combinations
landscape_ids = 1:10
rho <- c(50, 75, 100, 125, 150)
colony_density <- c(1000, 2000, 4000)
sample_sizes <- c(1000)

param_grid <- expand.grid(
  landscape_id = landscape_ids,
  rho = rho,
  colony_density =colony_density,
  sample_size = sample_sizes,
  stringsAsFactors = FALSE
)
saveRDS(param_grid, "exponentiated_quadratic_sim/param_grid.rds")

##### Set up run #####
# Get task ID from SLURM environment variable
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
params <- param_grid[task_id, ]

# Get landscape from saved file
fq <- readRDS(paste0(sprintf("landscapes/random_field_range10/landscape_%03d", params$landscape_id), ".rds"))

# Run simulation
result <- draw_N_bees(
  sample_size     = params$sample_size,
  landscape_size  = 1100,
  colonygrid_size = 700,
  trapgrid_size   = 300,
  number_traps    = 25,
  number_colonies = params$colony_density,
  colony_sizes    = rep(100, params$colony_density),
  rho            = params$rho,
  theta           = 0.5,
  resource_landscape = fq,
  batch_size = 1
)

# Save output
outfilepath <- sprintf("exponentiated_quadratic_sim/data/sim_result_%03d", task_id)
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)
saveRDS(result[[1]], paste(outfilepath, "/yobs.RDS", sep = ""))
saveRDS(result[[2]], paste(outfilepath, "/colonydata.RDS", sep = ""))
saveRDS(result[[3]], paste(outfilepath, "/trapdata.RDS", sep = ""))


# First, read in files for each task id
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
param_grid = readRDS("simulate_data/exponentiated_quadratic_sim/param_grid.rds")
results_path = sprintf("simulate_data/exponentiated_quadratic_sim/data/sim_result_%03d", task_id)
colony_data = readRDS(paste(results_path, "/colonydata.RDS", sep =""))
trap_data = readRDS(paste(results_path, "/trapdata.RDS", sep = ""))
yobs = readRDS(paste(results_path, "/yobs.RDS", sep=""))

# Remove all colonies from which ZERO bees were sampled
colony_data = colony_data[rowSums(yobs) > 0,]
yobs = yobs[rowSums(yobs) > 0,]

# Then, format data for stan
data = list()
data$C = nrow(yobs)
data$K = 25
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$y = yobs
data$lowerbound = 200
data$upperbound = 900
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$rho_center = 100
data$rho_sd_log = 0.5

# Fit Stan model!
stanFit = stan(file = "models/EQmodel.stan",
               data = data, seed = 5838299,
               chains = 4, cores = 4,
               iter = 10000,
               verbose = TRUE)
print("Model complete.")
saveRDS(stanFit, file=paste(results_path,"/stanFit.RDS", sep =""))
print("Model saved.")

