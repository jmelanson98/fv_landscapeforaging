##### Evaluate Pope Stan Fits -- Convergence + Accuracy #####
# Script Initiated: April 27, 2025
# By: Jenna Melanson
# Goals:
### Assess convergence of stanfit objects from 3_testPopeModel.R
### Assess fit accuracy (recreate figs 3 + 4 from Pope & Jha)
### Explore match between simulated data and true data


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
rstan_options(auto_write = TRUE) 
options(mc.cores = parallel::detectCores())

##### Source functions #####
source("simulate_data/0_PopeSimFunctions.R")

##### Load in data #####
param_grid = readRDS("simulate_data/batch_sim1/param_grid.rds")




##### Check for errors in stan fit
# list all .err files
err_dir <- "simulate_data/logs"
err_files <- list.files(err_dir, pattern = "\\.err$", full.names = TRUE)

# Function to check each file
check_err_file <- function(file) {
  lines <- readLines(file, warn = FALSE)

  list(
    file = basename(file),
    unfinished = any(grepl("DUE TO TIME LIMIT", lines, ignore.case = TRUE)),
    has_divergent = any(grepl("divergent", lines, ignore.case = TRUE)),
    has_low_ess = any(grepl("Effective Samples Size", lines, ignore.case = TRUE)),
    has_error = any(grepl("error", lines, ignore.case = TRUE))
  )
}

# Apply to all files
results <- lapply(err_files, check_err_file)

# Convert to a data frame
error_summary <- do.call(rbind, lapply(results, as.data.frame))

# Join with all_sim_df
param_grid$id = as.integer(rownames(param_grid))
error_summary$id = as.integer(sapply(strsplit(error_summary$file, "[_.]"), function(x) x[3]))
error_summary = left_join(error_summary, param_grid, by = "id")

# plot results
samplesize = ggplot(error_summary, aes(x = sample_size, colour = has_divergent)) +
  geom_histogram() +
  theme_minimal() +
  labs(title = "Divergent Transitions x Sample Size")

beta = ggplot(error_summary, aes(x = beta, colour = has_divergent)) +
  geom_histogram() +
  theme_minimal() +
  labs(title = "Divergent Transitions x Simulated Foraging Distance")

# remove any sims that had errors
param_grid = error_summary[error_summary$has_error == FALSE,]

# currently the saved version of these figures are from initial model fits with 
# divergent transitions -- before tightening prior on beta and decreasing step size
# ggsave("figures/simfigs/divergenttransitions_samplesize_informativebeta.jpg",
#        samplesize, width = 1500, height = 1000, units = "px")
# 
# ggsave("figures/simfigs/divergenttransitions_foragingdist_informativebeta.jpg",
#        beta, width = 1500, height = 1000, units = "px")


##### Plot simulated vs model-predicted colony average foraging distance
# e.g., figure 3 from Pope & Jha

#initate a data frame to store values for each simulation
columns = c("samplesize", "beta", "landscape_id", "colony_size_bin", 
            "true_colony_avg", "true_colony_sd", "model_colony_avg", "model_colony_sd")
allsim_colonies = data.frame(matrix(nrow = 2*nrow(param_grid), ncol = length(columns))) 
colnames(allsim_colonies) = columns
allsim_colonies$id = rep(param_grid$id,2)
allsim_colonies$colony_size_bin = c(rep(c("1-3"), nrow(param_grid)), rep(c("4+"), nrow(param_grid)))

for (sim in param_grid$id){
  # load in results for sim
  results_path = sprintf("simulate_data/batch_sim1/data/sim_result_%03d", sim)
  colony_data = colony_data = readRDS(paste(results_path, "/colonydata.RDS", sep =""))
  yobs = readRDS(paste(results_path, "/yobs.RDS", sep=""))
  stanFit = readRDS(paste(results_path, "/stanFitGQ.RDS", sep=""))
  
  # remove zero colonies
  print(class(yobs))
  print(dim(yobs))
  colony_data = colony_data[rowSums(yobs)>0,]
  yobs = yobs[rowSums(yobs) >0,]

  #extract model estimates
  colony_data$model_estimate = summary(stanFit, pars = c("colony_dist"))$summary[,1]
  
  #record observed colony sizes
  colony_data$observed_size = rowSums(yobs)
  
  # create categorical variable for colony size
 if(max(colony_data$observed_size > 3)){
 colony_data$size_bin = cut(colony_data$observed_size, breaks = c(0, 3, max(colony_data$observed_size)), labels = c("1-3", "4+"))
  } else {
  colony_data$size_bin = rep(c("1-3"), nrow(colony_data))
  }
  for (bin in c("1-3", "4+")){
    # add info to summary df
    allsim_colonies$samplesize[allsim_colonies$id == sim] = param_grid$sample_size[param_grid$id == sim]
    allsim_colonies$beta[allsim_colonies$id == sim] = param_grid$beta[param_grid$id == sim]
    allsim_colonies$landscape_id[allsim_colonies$id == sim] = param_grid$landscape_id[param_grid$id == sim]
    allsim_colonies$true_colony_avg[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = mean(colony_data$foraging_range[colony_data$size_bin == bin])
    allsim_colonies$true_colony_sd[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = sd(colony_data$foraging_range[colony_data$size_bin == bin])
    allsim_colonies$model_colony_avg[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = mean(colony_data$model_estimate[colony_data$size_bin == bin])
    allsim_colonies$model_colony_sd[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = sd(colony_data$model_estimate[colony_data$size_bin == bin])
    allsim_colonies$landscape_mean[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = summary(stanFit, pars = c("land_dist"))$summary[,1]
    allsim_colonies$landscape_sd[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = summary(stanFit, pars = c("land_dist"))$summary[,3]
    
    # save betas and thetas
    allsim_colonies$beta_mean[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = summary(stanFit, pars = c("beta"))$summary[,1]
    allsim_colonies$beta_msd[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = summary(stanFit, pars = c("beta"))$summary[,3]
    allsim_colonies$theta_mean[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = summary(stanFit, pars = c("theta"))$summary[,1]
    allsim_colonies$theta_sd[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = summary(stanFit, pars = c("theta"))$summary[,3]
    
  }}

colonyplot = ggplot(data = allsim_colonies, 
       aes(x = true_colony_avg, y = model_colony_avg, color = colony_size_bin)) +
  geom_point() +
  geom_errorbar(aes(ymin=model_colony_avg-model_colony_sd, ymax=model_colony_avg+model_colony_sd)) +
  geom_errorbarh(aes(xmin = true_colony_avg - true_colony_sd, xmax = true_colony_avg + true_colony_sd)) +
  geom_abline(intercept = 0, slope =1) +
  labs(x = "True colony foraging distance", y = "Model estimated colony foraging distance", color = "Number of workers per colony") +
  theme_minimal()

ggsave("figures/simfigs/Popefig3GQ.jpg", colonyplot, width = 4000, height = 1000, units= "px")
write.csv(allsim_colonies, "simulate_data/batch_sim1/forplotting.csv")
