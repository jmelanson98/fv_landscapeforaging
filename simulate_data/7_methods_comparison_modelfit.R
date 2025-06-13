##### Accuracy comparison of different simulation/modeling scenarios #####
##### THIS SCRIPT IS FOR THE MODEL FITTING AND PLOTTING PHASE ######
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
library(ggpubr)
library(filelock)

##### Set Environment #####
#setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging") # server
options(mc.cores = parallel::detectCores())


##### Source functions #####
source("simulate_data/src/GeneralizedSimFunctions.R")


##### Load in param grid and simulated data #####
param_grid = readRDS("simulate_data/methods_comparison/param_grid.rds")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
inputfilepath <- sprintf("simulate_data/methods_comparison/data/sim_result_%03d", task_id)
yobs = readRDS(paste(inputfilepath, "/yobs.RDS", sep = ""))
colony_data = readRDS(paste(inputfilepath, "/colonydata.RDS", sep = ""))
trap_data = readRDS(paste(inputfilepath, "/trapdata.RDS", sep = ""))


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
data$rho_sd = 0.5


### Current params
current_params = param_grid[task_id,]


### Correctly subset yobs based on approach
if (current_params$model_approach == "all"){
  data$y = yobs
  data$C = nrow(yobs)
} else if (current_params$model_approach == "singletons"){
  colony_data = colony_data[rowSums(yobs) > 0,]
  yobs = yobs[rowSums(yobs) > 0,]
  data$y = yobs
  data$C = nrow(yobs)
} else if (current_params$model_approach == "doubletons"){
  colony_data = colony_data[rowSums(yobs) > 1,]
  yobs = yobs[rowSums(yobs) > 1,]
  data$y = yobs
  data$C = nrow(yobs)
} else if (current_params$model_approach == "centroid"){
  colony_data = colony_data[rowSums(yobs) > 1,]
  yobs = yobs[rowSums(yobs) > 1,]
  data$y = yobs
  data$C = nrow(yobs)
} else {
  print("Not a valid approach.")
}



### Fit Stan model if necessary
if (current_params$model_approach != "centroid"){
  #select stan model to fit
  stanfile = paste("models/", current_params$distance_decay, ".stan", sep = "")
  
  #fit and save model
  stanFit = stan(file = stanfile,
                 data = data, seed = 5838299,
                 chains = 4, cores = 4,
                 iter = 10000,
                 verbose = TRUE)
  print("Model complete.")
  saveRDS(stanFit, file=paste(inputfilepath,"/stanFit.RDS", sep =""))
  print("Model saved.")
}



##### Estimate and save landscape-level foraging distance ######
# use lock to make sure multiple tasks don't write to output at the same time
lockfile <- "output.lock"
rds_file <- "simulate_data/methods_comparison/output.rds"


if (current_params$model_approach == "centroid"){
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
  avg_dists <- apply(yobs, 1, function(counts) {
    get_avg_distance_to_centroid(counts, trap_data)
  })
  
  # save in output grid
  # try to acquire the lock (waits up to 60 seconds)
  lock <- lock(lockfile, timeout = 60000)
  
  if (!is.null(lock)) {
    df <- readRDS(rds_file)
    
    df$model_average_foraging[df$task_id == task_id] = mean(avg_dists)
    df$model_sd_foraging[df$task_id == task_id] = sd(avg_dists)
    
    # Save the updated dataframe
    saveRDS(df, rds_file)
    
    # Release lock
    unlock(lock)
  } else {
    stop("Could not acquire lock on output file.")
  }
  
} else {
  colony_data$model_estimate = summary(stanFit, pars = c("colony_dist"))$summary[,1]
  
  # save output
  # try to acquire the lock (waits up to 60 seconds)
  lock <- lock(lockfile, timeout = 60000)
  
  if (!is.null(lock)) {
    df <- readRDS(rds_file)
    
    df$model_average_foraging[df$task_id == task_id] = mean(colony_data$model_estimate)
    df$model_sd_foraging[df$task_id == task_id] = sd(colony_data$model_estimate)
    df$model_mu[df$task_id == task_id] = summary(stanFit, pars = c("mu"))$summary[,1]
    
    # Save the updated dataframe
    tmp <- tempfile()
    saveRDS(df, tmp)
    file.rename(tmp, rds_file)
    
    # Release lock
    unlock(lock)
  } else {
    stop("Could not acquire lock on output file.")
  }
}
