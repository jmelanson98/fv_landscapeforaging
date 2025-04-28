##### Evaluate Pope Stan Fits -- Convergence + Accuracy #####
# Script Initiated: April 27, 2025
# By: Jenna Melanson
# Goals:
### Assess convergence of stanfit objects from 3_testPopeModel.R
### Assess fit accuracy (recreate figs 3 + 4 from Pope & Jha)
### Explore match between simulated data and true data


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

##### Set Environment #####
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
#setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging") # server
rstan_options(auto_write = TRUE) 
options(mc.cores = parallel::detectCores())

##### Source functions #####
source("simulate_data/0_PopeSimFunctions.R")

##### Load in data #####
all_sim_df = readRDS("simulate_data/batch_sim1/all_sim_df.RDS")




##### Check for errors in stan fit
# list all .err files
err_dir <- "simulate_data/logs"
err_files <- list.files(err_dir, pattern = "\\.err$", full.names = TRUE)

# Function to check each file
check_err_file <- function(file) {
  lines <- readLines(file, warn = FALSE)

  list(
    file = basename(file),
    unfinished = any(grepl("DUE TO TIME LIMIT", lines, ignore.case = TRUE)),
    has_divergent = any(grepl("divergent", lines, ignore.case = TRUE)),
    has_low_ess = any(grepl("Effective Samples Size", lines, ignore.case = TRUE)),
    has_error = any(grepl("error", lines, ignore.case = TRUE))
  )
}

# Apply to all files
results <- lapply(err_files, check_err_file)

# Convert to a data frame
error_summary <- do.call(rbind, lapply(results, as.data.frame))

# Join with all_sim_df
all_sim_df$id = as.integer(rownames(all_sim_df))
error_summary$id = as.integer(sapply(strsplit(error_summary$file, "[_.]"), function(x) x[3]))
error_summary = left_join(error_summary, all_sim_df, by = "id")
error_summary$true_avg <- sapply(error_summary$true_colony_foraging, function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    NA_real_
  } else {
    mean(unlist(x), na.rm = TRUE)
  }
})


# plot results
samplesize = ggplot(error_summary, aes(x = sample_size, colour = has_divergent)) +
  geom_histogram() +
  theme_minimal() +
  labs(title = "Divergent Transitions x Sample Size")

beta = ggplot(error_summary, aes(x = true_avg*5, colour = has_divergent)) +
  geom_histogram() +
  theme_minimal() +
  labs(title = "Divergent Transitions x Simulated Foraging Distance")

# currently the saved version of these figures are from initial model fits with 
# divergent transitions -- before tightening prior on beta and decreasing step size
# ggsave("figures/simfigs/divergenttransitions_samplesize_informativebeta.jpg",
#        samplesize, width = 1500, height = 1000, units = "px")
# 
# ggsave("figures/simfigs/divergenttransitions_foragingdist_informativebeta.jpg",
#        beta, width = 1500, height = 1000, units = "px")




##### Plot simulated vs model-predicted landscape average foraging distance 
# e.g., figure 4 from Pope & Jha
all_sim_df$mean_true_foraging = sapply(all_sim_df$true_colony_foraging, mean)
fig4 = ggplot(all_sim_df, aes(x = as.numeric(mean_true_foraging), y = as.numeric(average_foraging), color = sample_size)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0)
fig4


##### Plot simulated vs model-predicted colony average foraging distance
# e.g., figure 3 from Pope & Jha

# make a tible with unlisted estimates
long_df <- pmap_dfr(
  .l = as_tibble(all_sim_df) %>% 
  dplyr::select(landscape_id, beta, sample_size, 
                             average_foraging, mean_true_foraging, 
                             true_colony_foraging, colony_foraging_estimates, 
                             colony_sizes),
  function(landscape_id, beta, sample_size,
                              average_foraging, mean_true_foraging, 
                              true_colony_foraging, colony_foraging_estimates, 
                              colony_sizes) {
  tibble(
    landscape_id = landscape_id,
    beta = beta,
    sample_size = sample_size,
    average_foraging = average_foraging,
    mean_true_foraging = mean_true_foraging,
    true_colony_foraging = true_colony_foraging,
    colony_foraging_estimates = colony_foraging_estimates,
    colony_sizes = colony_sizes
  )
})


# plot
partial_df = filter(long_df, colony_sizes <= 10)
fig3 = ggplot(partial_df, aes(x = true_colony_foraging, y = colony_foraging_estimates, color = colony_sizes)) +
  geom_point() +
  geom_abline(slope = 1, intercept =0)
fig3
