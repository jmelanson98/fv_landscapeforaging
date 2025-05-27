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
    colony_density =colony_density,
    sample_size = sample_sizes,
    stringsAsFactors = FALSE
  )
  saveRDS(param_grid, "simulate_data/methods_comparison/param_grid.rds")
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



##### Fit Models ######

### Format most basic data for Stan
data = list()
data$K = 25
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$lowerbound = 400
data$upperbound = 1100
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$rho_center = 4.5
data$rho_sd_log = 0.5

# Save original yobs
yobs.orig = yobs

### Correctly subset yobs based on approach
if (model_approach == "all"){
  data$y = yobs
  data$C = nrow(yobs)
} else if (model_approach == "singletons"){
  colony_data = colony_data[rowSums(yobs) > 0,]
  yobs = yobs[rowSums(yobs) > 0,]
  data$y = yobs
  data$C = nrow(yobs)
} else if (model_approach == "doubletons"){
  colony_data = colony_data[rowSums(yobs) > 1,]
  yobs = yobs[rowSums(yobs) > 1,]
  data$y = yobs
  data$C = nrow(yobs)
} else if (model_approach == "centroid"){
  colony_data = colony_data[rowSums(yobs) > 1,]
  yobs = yobs[rowSums(yobs) > 1,]
  data$y = yobs
  data$C = nrow(yobs)
} else {
  print("Not a valid approach.")
}



### Fit Stan model if necessary
if (model_approach != "centroid"){
  #select stan model to fit
  stanfile = paste("models/", distance_decay, ".stan", sep = "")
  
  #fit and save model
  stanFit = stan(file = stanfile,
                 data = data, seed = 5838299,
                 chains = 4, cores = 4,
                 iter = 10000,
                 verbose = TRUE)
  print("Model complete.")
  saveRDS(stanFit, file=paste(results_path,"/stanFit.RDS", sep =""))
  print("Model saved.")
}



##### Estimate and plot landscape-level foraging distance ######
if (model_approach == "centroid"){
  get_avg_distance_to_centroid <- function(counts, coords) {
    total <- sum(counts)
    if (total == 0) return(NA)
    
    # calculate weighted centroid
    centroid_x <- sum(counts * coords$trap_x) / total
    centroid_y <- sum(counts * coords$trap_y) / total
    
    # calculate distance of each trap from centroid
    dists <- sqrt((coords$trap_x - centroid_x)^2 + (coords$trap_y - centroid_y)^2)
    
    # calculate weighted average distance to centroid
    mean_dist <- sum(dists * counts) / total
    return(mean_dist)
  }
  
  #apply function to yobs
  avg_dists <- apply(yobs_temp, 1, function(counts) {
    get_avg_distance_to_centroid(counts, trap_data)
  })
  
} else if
