##### Test accuracy of landscape effects #####
# Script Initiated: June 15, 2025
# By: Jenna Melanson
# Goals:
### Check whether models can accurately return estimates of landscape effects on foraging distance;
### assuming that a "true" landscape metric value is known for each 


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
library(tibble)
library(ggpubr)

##### Set Environment #####
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
rstan_options(auto_write = TRUE) 
options(mc.cores = parallel::detectCores())

##### Source functions #####
source("simulate_data/src/GeneralizedSimFunctions.R")

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
  alpha = 3.5,
  resource_landscape = fq,
  nesting_landscape = NULL,
  distance_decay = "exponential"
)

# Save results
saveRDS(result, "simulate_data/methods_comparison/landscape_effects/simdata.rds")
result = readRDS("simulate_data/methods_comparison/landscape_effects/simdata.rds")

# Write outputs to variables
yobs = result[[1]]
colony_data = result[[2]]
trap_data = result[[3]]

# First try for only detected colonies
yobs_detected = yobs[rowSums(yobs > 0),]
colony_data_detected = colony_data[rowSums(yobs > 0),]

# Prep data list for Stan
data = list()
data$y = yobs_detected
data$C = nrow(data$y)
data$K = ncol(data$y)
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$lowerbound = 400
data$upperbound = 1100
data$landscape = colony_data_detected$landscape_metric
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$rho_center = 3.5
data$rho_sd = 0.5

#select stan model to fit
stanfile = paste("models/landscape_exp.stan")

#fit and save model
stanFitLandscape = stan(file = stanfile,
                  data = data, seed = 5838299,
                  chains = 4, cores = 4,
                  iter = 10000,
                  verbose = TRUE)
saveRDS(stanFitAll, file="simulate_data/methods_comparison/observed_vs_unobserved/stanFitAll.RDS")

