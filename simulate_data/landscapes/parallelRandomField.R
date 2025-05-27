##### Simulate resource landscapes and save #####
# Script Initiated: April 10, 2025
# By: Jenna Melanson
# Goal: self contained script for simulating resoure landscapes on AllianceCan server

#set working directory on server -- change if working locally!!
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging/simulate_data")

##### Load packages #####
library(matrixStats)
library(sp)
library(gstat)
library(raster)
library(rasterVis)
library(parallel)
library(future)
library(furrr)
library(terra)

##### Source Helper Functions #####
source("src/GeneralizedSimFunctions.R")

##### Create directory to save landscapes #####
dir.create("landscapes/random_field_range10_large", recursive = TRUE, showWarnings = FALSE)

##### Set up run #####
# Get task ID from SLURM environment variable
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Simulate landscape
fq <- simulateLandscapeRaster(landscape_size = 1500, resource_range = 10)

# Save landscape
saveRDS(fq, file = paste(sprintf("landscapes/random_field_range10/landscape_%03d", task_id), ".rds", sep = ""))

