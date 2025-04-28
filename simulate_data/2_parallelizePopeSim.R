##### Run parallel simulations of Pope model on server #####
# Script Initiated: April 8, 2025
# By: Jenna Melanson
# Goal: self contained script for running multiple iterations of Pope simulation on AllianceCan server

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
source("PopeSimFunctions.R")

##### Prepare to run in parallel ! #####
# set parameter combinations
landscape_ids = 1:10
betas <- c(-1/10, -1/25, -1/50, -1/75, -1/100)
sample_sizes <- c(250, 500, 1000)

param_grid <- expand.grid(
  landscape_id = landscape_ids,
  beta = betas,
  sample_size = sample_sizes,
  stringsAsFactors = FALSE
)
saveRDS(param_grid, "batch_sim1/param_grid.rds")

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
  number_colonies = 1000,
  colony_sizes    = rep(100, 1000),
  beta            = params$beta,
  theta           = 0.5,
  resource_landscape = fq,
  batch_size = 1
)

# Save output
outfilepath <- sprintf("batch_sim1/data/sim_result_%03d", task_id)
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)
saveRDS(result[[1]], paste(outfilepath, "/yobs.RDS", sep = ""))
saveRDS(result[[2]], paste(outfilepath, "/colonydata.RDS", sep = ""))
saveRDS(result[[3]], paste(outfilepath, "/trapdata.RDS", sep = ""))


# for later
# files <- list.files("results", pattern = "sim_result_.*rds", full.names = TRUE)
# results <- lapply(files, readRDS)
# param_grid <- readRDS("param_grid.rds")
# 
# full_results <- bind_cols(param_grid, tibble(result = results))

