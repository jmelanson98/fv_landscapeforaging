##### Run models with and without unobserved colonies #####
# Script Initiated: June 17, 2025
# By: J. Melanson
# Goals:
### Test some methods for parallelizing stan code on the server with cmdstan

##### Load packages #####
library(cmdstanr)
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
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
set_cmdstan_path("/home/melanson/projects/def-ckremen/melanson/cmdstan")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))


##### Load in some simulation data!! #####

## NOTE: this is simulation 74!! this data is for landscape = 4, rho = 50, exponential decay
yobs = readRDS("simulate_data/methods_comparison/data/sim_result_074/yobs.RDS")
colony_data = readRDS("simulate_data/methods_comparison/data/sim_result_074/colonydata.RDS")
trap_data = readRDS("simulate_data/methods_comparison/data/sim_result_074/trapdata.RDS")


# Prep data list for Stan
data = list()
data$y = yobs
data$C = nrow(data$y)
data$K = ncol(data$y)
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$lowerbound = 400
data$upperbound = 1100
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$rho_center = 4.5
data$rho_sd = 0.5

# Table of conditions
model = c("reduce_sum_exponential", "reduce_sum_exponentialGQ")
threads = c(4, 8)
grid = expand.grid(model = model, 
                   threads = threads)
grid$task = 1:nrow(grid)

current = grid[grid$task == task_id,]

print(current$model)
print(current$threads)

#select stan model to fit
exe_file <- paste("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging/models/", current$model, sep = "")
mod <- cmdstan_model(exe_file, compile = FALSE)


# Run sampling
threads_per_chain = current$threads

# add grainsize to data list
grainsize <- max(floor(C / (threads_per_chain * 5)), 1)
data_list$grainsize = grainsize


#fit and save model
print("Starting sampling.")
fit <- mod$sample(
  data = data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = threads_per_chain,
  refresh = 100,
  iter_warmup = 1000,
  iter_sampling = 5000
)

saveRDS(fit, paste("methods_comparison/observed_vs_unobserved/reducesum/", current$thread, current$model, "stanFit.rds"))
