##### What is the effect of a better prior on colony location? #####
# Script Initiated: July 16, 2025
# By: Jenna Melanson
# Goals:
### Determine whether soft prior on colony location helps (rather than hard limit)


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

##### Simulate some data #####
# use 2000 background colonies and 1000 observations?
# this distribution won't be totally accurate but it's hard to fit a model with 8000 background colonies

# Get landscape from saved file
fq <- readRDS("simulate_data/landscapes/landscapes/random_field_range10/landscape_001.rds")
fq = terra::rast(fq)

# Run simulation
result <- draw_bees_colony_restricted(
  sample_size     = 1000,
  landscape_size  = 1500,
  colonygrid_size = 700,
  trapgrid_size   = 300,
  number_traps    = 25,
  number_colonies = 2000,
  colony_sizes    = rep(100, 2000),
  rho            = 50,
  theta           = 0.5,
  resource_landscape = fq,
  nesting_landscape = NULL,
  distance_decay = "exponential"
)

# Save results
saveRDS(result, "simulate_data/methods_comparison/observed_vs_unobserved/simdata.rds")
result = readRDS("simulate_data/methods_comparison/observed_vs_unobserved/simdata.rds")

# Write outputs to variables
yobs = result[[1]]
colony_data = result[[2]]
trap_data = result[[3]]

# Make some subsets within and without zero colonies
yobs_detected = yobs[rowSums(yobs) > 0,] 

# Prep data list for Stan
data = list()
data$y = yobs_detected
data$C = nrow(data$y)
data$K = 25
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$lowerbound = 0
data$upperbound = 1500
data$nestrange = 175
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$rho_center = 4.5
data$rho_sd = 0.5


#select stan model to fit
stanfile = paste("models/colonydetectionprior.stan")

#fit and save model
stanFitCPrior = stan(file = stanfile,
                             data = data, seed = 5838299,
                             chains = 4, cores = 4,
                             iter = 2000,
                             verbose = TRUE)
saveRDS(stanFit, file="simulate_data/methods_comparison/observed_vs_unobserved/stanFitAll.RDS")
