##### Can a centroid/mean center approach detect landscape effects? #####
# Script Initiated: July 11, 2025
# By: Jenna Melanson
# Goals:
### Check whether centroid-regression can detect effects of landscape on foraging distance


##### Load packages #####
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
library(posterior)
library(terra)
library(lme4)


##### Set Environment #####
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
#setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local

##### Source functions #####
source("simulate_data/src/GeneralizedSimFunctions.R")

##### Load in floral quality landscape #####
fq = readRDS("simulate_data/landscapes/landscapes/random_field_range10/landscape_003.rds")

##### Simulate landscape characteristics "landscape" #####
landscape_char = simulateLandscapeRaster(landscape_size = 1500, resource_range = 400)
saveRDS(landscape_char, "simulate_data/landscapes/landscapes/random_field_range10/landscape_001.rds")
#landscape_char = readRDS("simulate_data/landscapes/landscapes/random_field_range10/landscape_001.rds")

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
  theta           = 0.5,
  alpha = 0.5,
  resource_landscape = fq,
  configuration = landscape_char,
  nesting_landscape = NULL,
  distance_decay = "exponential"
)

# Save results
saveRDS(result, "simulate_data/methods_comparison/landscape_effects/centroids/simdata.rds")
#result = readRDS("simulate_data/methods_comparison/landscape_effects/centroids/simdata.rds")

# Write outputs to variables
yobs = result[[1]]
colony_data = result[[2]]
trap_data = result[[3]]

# Filter yobs + colony data to only doubletons
ydoubleton = yobs[rowSums(yobs) > 1,]
doubleton_colonies = colony_data[rowSums(yobs) > 1,]

# Calculate mean centers for each colony and grab associated data
for (i in 1:nrow(ydoubleton)){
  row = ydoubleton[i,]
  
  #mean center
  doubleton_colonies$mean_x[i] = sum(row*trap_data$trap_x)/sum(row)
  doubleton_colonies$mean_y[i] = sum(row*trap_data$trap_y)/sum(row)
  
  #mean distance of workers from mean center
  x_dist = sum(abs(trap_data$trap_x - doubleton_colonies$mean_x[i])*row)/sum(row)
  y_dist = sum(abs(trap_data$trap_y - doubleton_colonies$mean_y[i])*row)/sum(row)
  doubleton_colonies$mean_dist[i] = sqrt(x_dist^2 + y_dist^2)
  
  # landscape value at mean center
  doubleton_colonies$center_landscape[i] = landscape_char[ceiling(doubleton_colonies$mean_x[i]), ceiling(doubleton_colonies$mean_y[i])]
  
}

# set center_landscape values to numeric? right now they're a list?
doubleton_colonies$center_landscape = as.numeric(doubleton_colonies$center_landscape)

# Linear model of mean_dist ~ landscape
regression = lm(mean_dist ~ center_landscape, data = doubleton_colonies)
summary(regression)


# Lil plot
mean_center_plot = ggplot(doubleton_colonies, aes(x = center_landscape, y = mean_dist)) +
  geom_point() +
  xlab("Landscape characteristic at mean center") +
  ylab("Colony mean foraging distance") +
  theme_minimal() + 
  geom_smooth(method='lm')

ggsave("simulate_data/methods_comparison/landscape_effects/centroids/mean_center_plot", mean_center_plot,
       width = 500, height = 500, units = "px")
