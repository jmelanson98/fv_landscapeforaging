###### Get cumulative colony posteriors with neutral simulation
### Simulate 100 null datasets
### Get cumulative colony posteriors
### January 16, 2026
### J. Melanson

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
library(sf)
library(terra)
library(stringr)

##### Set Environment #####
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/UBC/bombus_project"

setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "/home/melanson/projects/def-ckremen/melanson"
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
source("src/analysis_functions.R")


################################################################################
###############               SIMULATE DATA                 ###################
################################################################################
 
# Load in data / models
if (task_id < 101){
  species = "mixtus"
  fit = readRDS("analysis/kernel_locations/foraging_modelfits/simpleforaging_mixtus_steepgradient.rds")
  sibs1 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
  sibs2 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
  draw_index = task_id
} else {
  species = "impatiens"
  fit = readRDS("analysis/kernel_locations/foraging_modelfits/simpleforaging_impatiens_steepgradient.rds")
  sibs1 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
  sibs2 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
  draw_index = task_id -100
}

effort1 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
effort2 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))
samplepoints = read.csv(paste0(bombus_path, "/raw_data/allsamplepoints.csv"), header = FALSE)

# Get colony counts (real data)
data = prep_stan_simpleforaging_bothyears(sibs1,
                                          sibs2,
                                          effort1,
                                          effort2,
                                          samplepoints)


##### Simulate datasets #####

# Get trap data
traps_m = data[[3]]

# Get rho draws from fitted model
rhodraws = rstan::extract(fit, pars = "rho")$rho
iterrho = rhodraws[draw_index*4000/100]

# Simulate one dataset
sim = draw_simple_true_sites(sample_size = 6000,
                          number_colonies = 18000,
                          colony_sizes = rep(20,18000),
                          trap_data = traps_m,
                          rho = iterrho,
                          distance_decay = "exponentiated_quadratic")

# Draw colonies from simulated dataset to match real data
counts_summary = data[[2]] %>% 
  group_by(stansibkey) %>%
  summarize(n = sum(count)) %>%
  group_by(n) %>%
  summarize(num_per_bin = n())
simcounts = sim %>% 
  group_by(colonyid) %>%
  summarize(n = sum(counts)) %>%
  filter(n %in% counts_summary$n) %>%
  group_by(n) %>%
  tidyr::nest() %>%
  left_join(counts_summary, by = "n")
simIDs_tokeep = map2(simcounts$data, simcounts$num_per_bin,
                          ~ slice_sample(.x, n = .y, replace = FALSE)) %>%
                      bind_rows()
sim_tokeep = sim[sim$colonyid %in% simIDs_tokeep$colonyid,]

if(length(unique(sim_tokeep$colonyid)) < length(unique(data[[2]]$stansibkey))){
  print("Not enough colonies saved... :(")
} else {
  print("Colonies simulated and saved! :)")
}

# Write to file
dir = paste0("analysis/kernel_locations/neutral_data/", species, "_dataset", task_id)
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
write.csv(sim_tokeep, paste0(dir, "/simulation", task_id, ".csv"), row.names = FALSE)



################################################################################
###############                  FIT MODELS                  ###################
################################################################################