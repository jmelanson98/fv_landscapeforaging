##### Incorporating multiple landscapes and timepoints into simulation and modeling workflow #####
# Script Initiated: September 15, 2025
# By: Jenna Melanson

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
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
options(mc.cores = parallel::detectCores())

##### Source functions #####
source("simulate_data/src/GeneralizedSimFunctions.R")

##### Simulate multi-landscape, multi-timepoint data #####
result = draw_multi_landscape_timepoint(sample_size = 1000,
                                        num_landscape = 3,
                                        num_timepoints = 2,
                                        landscape_size = 1500,
                                        trapgrid_size = 300,
                                        number_traps = 25,
                                        number_colonies = 10000,
                                        colony_sizes = rep(100,10000),
                                        rho = 50,
                                        theta = 0.5,
                                        distance_decay = "exponential")
yobs_sum = result[[1]][[1]] + result[[1]][[2]]
floral = result[[4]]
trap_data = result[[3]]

# check distribution of sibship sizes
tbl <- table(rowSums(yobs_sum)[rowSums(yobs_sum) > 0])
prop_tbl <- prop.table(tbl)  # gives proportions
barplot(prop_tbl, ylab = "Proportion", xlab = "Count category")

##### Fit with basic multinomial model #####
yobs = yobs_sum[rowSums(yobs_sum) >0,]

data = list()
data$y = yobs
data$C = nrow(data$y)
data$K = ncol(data$y)
data$lowerbound = 0
data$upper_y = 1500
data$upper_x = 4500
data$floral = rowMeans(floral)
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))


#select stan model to fit
stanfile = "models/multinomial_landscape.stan"

#fit and save model
stanFit = stan(file = stanfile,
               data = data, seed = 5838299,
               chains = 4, cores = 4,
               iter = 6000, warmup = 1000,
               verbose = TRUE)