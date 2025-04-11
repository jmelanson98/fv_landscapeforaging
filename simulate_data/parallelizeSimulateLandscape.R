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

##### Source Helper Functions #####
source("PopeSimFunctions.R")

##### Prepare to run in parallel ! #####

#simulate landscapes and save them
landscape_ids  = 1:10
for (i in landscape_ids) {
  fq <- simulateLandscape(landscape_size = 1100, resource_range = 10)
  saveRDS(fq, file = paste0("landscapes/random_field_range10/landscape_", i, ".rds"))
}