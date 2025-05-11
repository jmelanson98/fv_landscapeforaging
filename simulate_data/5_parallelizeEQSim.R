##### Run parallel simulations on server using exponentiated quadratic #####
# Script Initiated: May 9, 2025
# By: Jenna Melanson
# Goal: self contained script for running multiple iterations of Pope simulation on AllianceCan server
# in this case, using an exponentiated quadratic to explain patterns of bee visitation

#set working directory on server -- change if working locally!!
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging/simulate_data")

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
library(gstat)

##### Source Helper Functions #####
source("0_EQSimFunctions.R")

##### Prepare to run in parallel ! #####
# set parameter combinations
landscape_ids = 1:10
beta <- c(1/50, 1/75, 1/100, 1/125, 1/150)
colony_density <- c(1000, 2000, 4000)
sample_sizes <- c(1000)

param_grid <- expand.grid(
  landscape_id = landscape_ids,
  rhos = rhos,
  colony_density =colony_density,
  sample_size = sample_sizes,
  stringsAsFactors = FALSE
)
saveRDS(param_grid, "exponentiated_quadratic_sim/param_grid.rds")

##### Set up run #####
# Get task ID from SLURM environment variable
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
params <- param_grid[task_id, ]

# Get landscape from saved file
fq <- readRDS(paste0(sprintf("landscapes/random_field_range10/landscape_%03d", params$landscape_id), ".rds"))

# Run simulation
result <- draw_N_bees(
  sample_size     = params$sample_size,
  landscape_size  = 1100,
  colonygrid_size = 700,
  trapgrid_size   = 300,
  number_traps    = 25,
  number_colonies = params$colony_density,
  colony_sizes    = rep(100, params$colony_density),
  rho            = params$beta,
  theta           = 0.5,
  resource_landscape = fq,
  batch_size = 1
)

# Save output
outfilepath <- sprintf("exponentiated_quadratic_sim/data/sim_result_%03d", task_id)
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)
saveRDS(result[[1]], paste(outfilepath, "/yobs.RDS", sep = ""))
saveRDS(result[[2]], paste(outfilepath, "/colonydata.RDS", sep = ""))
saveRDS(result[[3]], paste(outfilepath, "/trapdata.RDS", sep = ""))


