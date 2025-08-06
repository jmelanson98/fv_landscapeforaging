##### How well can we estimate landscape effects? #####
##### THIS SCRIPT IS FOR THE MODEL FITTING AND PLOTTING PHASE ######
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
#library(ggpubr)
library(filelock)

##### Set Environment #####
#setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging") # server
options(mc.cores = parallel::detectCores())


##### Source functions #####
source("simulate_data/src/GeneralizedSimFunctions.R")


##### Load in param grid and simulated data #####
param_grid = readRDS("simulate_data/foraging_distance/methods_comparison/landscape_effects/output_sim.rds")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
params = param_grid[param_grid$task_id == task_id,]
inputfilepath <- sprintf("simulate_data/foraging_distance/methods_comparison/landscape_effects/data/sim_result_%03d", task_id)
yobs = readRDS(paste(inputfilepath, "/yobs.RDS", sep = ""))
colony_data = readRDS(paste(inputfilepath, "/colonydata.RDS", sep = ""))
trap_data = readRDS(paste(inputfilepath, "/trapdata.RDS", sep = ""))
print("Loading landscape metric raster.")
repo_name = paste0("simulate_data/landscapes/random_field_range", params$autocorrelation)
landscape_name = as.character(sprintf("/landscape_%03d", params$landscape_id))
lmq = readRDS(paste0(repo_name, landscape_name, ".rds"))


### Current params
current_params = params

##### Fit Models ######

# Subset only doubletons
yobs_doubleton = yobs[rowSums(yobs) > 1,] 
doubleton_colonies = colony_data[rowSums(yobs) > 1,]

# Prep data list for Stan
data = list()
data$y = yobs_doubleton
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


#select stan model to fit
stanfile = "models/multinomial_landscape.stan"
  
#fit and save model
stanFit = stan(file = stanfile,
               data = data, seed = 5838299,
               chains = 4, cores = 4,
               iter = 6000, warmup = 1000,
               verbose = TRUE)
print("Model complete.")
saveRDS(stanFit, file=paste(inputfilepath,"/stanFit.RDS", sep =""))
print("Model saved.")


#### save stan fit results
colony_data_partial = colony_data[rowSums(yobs) > 1,]
colony_data_partial$model_estimate = summary(stanFit, pars = c("colony_dist"))$summary[,1]
saveRDS(colony_data_partial, paste(inputfilepath, "/colonydata_partial.RDS", sep = ""))

# save the result of *this* simulation only
current_params$model_average_foraging = mean(colony_data$model_estimate)
current_params$model_sd_foraging = sd(colony_data$model_estimate)
current_params$model_mu = summary(stanFit, pars = c("mu"))$summary[,1]
current_params$model_rho = summary(stanFit, pars = c("rho"))$summary[,1]
current_params$model_alpha = summary(stanFit, pars = c("alpha"))$summary[,1]
current_params$model_alpha_min = summary(stanFit, pars = c("alpha"))$summary[,4]
current_params$model_alpha_max = summary(stanFit, pars = c("alpha"))$summary[,8]




##### Calculate sibling-separation distance with mean-center approach ######
# Calculate mean centers for each colony and grab associated data
for (i in 1:nrow(yobs_doubleton)){
  row = yobs_doubleton[i,]
  
  #mean center
  doubleton_colonies$mean_x[i] = sum(row*trap_data$trap_x)/sum(row)
  doubleton_colonies$mean_y[i] = sum(row*trap_data$trap_y)/sum(row)
  
  #mean distance of workers from mean center
  x_dist = sum(abs(trap_data$trap_x - doubleton_colonies$mean_x[i])*row)/sum(row)
  y_dist = sum(abs(trap_data$trap_y - doubleton_colonies$mean_y[i])*row)/sum(row)
  doubleton_colonies$mean_dist[i] = sqrt(x_dist^2 + y_dist^2)
  
  # landscape value at mean center
  doubleton_colonies$center_landscape[i] = lmq[ceiling(doubleton_colonies$mean_x[i]), ceiling(doubleton_colonies$mean_y[i])]
  
}

# set center_landscape values to numeric
doubleton_colonies$center_landscape = as.numeric(doubleton_colonies$center_landscape)

# Linear model of mean_dist ~ landscape
regression = lm(mean_dist ~ center_landscape, data = doubleton_colonies)
summary(regression)
current_params$mean_center_dist = mean(doubleton_colonies$mean_dist)
current_params$mean_center_alpha = summary(regression)$coefficients["center_landscape","Estimate"]
current_params$mean_center_alpha_sd = summary(regression)$coefficients["center_landscape","Std. Error"]

# Save values for the current iteration!
saveRDS(current_params, paste(inputfilepath, "/iteration_output.rds", sep = ""))
result = tibble(current_params)


# use lock to make sure multiple tasks don't write to output at the same time
lockfile <- "output.lock"
rds_file <- "simulate_data/foraging_distance/methods_comparison/landscape_effects/output_fit.rds"

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

