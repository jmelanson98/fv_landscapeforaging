##### Fit Pope Stan code to simulated data #####
# Script Initiated: April 10, 2025
# By: Jenna Melanson
# Goals:
### Fit Pope Stan code using a multiple iterations of simulated data
### 10 landscapes x 5 beta values x 3 sampling intensities
### Generate figure analogous to fig 1A + plots of colony posterior distributions


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
#setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging") # server
rstan_options(auto_write = TRUE) 
options(mc.cores = parallel::detectCores())

##### Source functions #####
source("simulate_data/0_PopeSimFunctions.R")


##### Replicate Figures 3 & 4 #####
## this must be run on server because my computer doesn't have the memory, RIP

# First, read in files for each task id
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
param_grid = readRDS("simulate_data/batch_sim1/param_grid.rds")
results_path = sprintf("simulate_data/batch_sim1/data/sim_result_%03d", task_id)
colony_data = readRDS(paste(results_path, "/colonydata.RDS", sep =""))
trap_data = readRDS(paste(results_path, "/trapdata.RDS", sep = ""))
yobs = readRDS(paste(results_path, "/yobs.RDS", sep=""))

# Remove all colonies from which ZERO bees were sampled
colony_data = colony_data[rowSums(yobs) > 0,]
yobs = yobs[rowSums(yobs) > 0,]

# Then, format data for stan
data = list()
data$C = nrow(yobs)
data$K = 25
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$y = yobs
data$lowerbound = 200
data$upperbound = 900
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$priorBe = 0.1

# Function to estimate necessary RAM for model
# estimate_stanfit_size_gb <- function(num_params, num_iters, num_chains, overhead = 2.5) {
#   (num_params * num_iters * num_chains * 8 * overhead) / (1024^3)
# }

# Fit Stan model!
stanFitGQ = stan(file = "models/popeModified.stan",
               data = data, seed = 5838299,
               warmup = 1000, iter = 10000,
               control = list(adapt_delta = 0.999,
                              stepsize = 0.001,
                              max_treedepth = 20),
               chains = 4, cores = 4,
               verbose = TRUE)
print("Model complete.")
saveRDS(stanFitGQ, file=paste(results_path,"/stanFitGQ.RDS", sep =""))
print("Model saved.")

#Plot the posteriors of some colonies
plot_list = list()
numplots = 12
legends = list()

for (i in 1:numplots){
  c_id = colony_data$colonyid[i]
  delta_draws = as.data.frame(rstan::extract(stanFitGQ, pars = "delta")$delta[, c_id,])
  colnames(delta_draws) = c("x","y")
  trap_data$trap_count = yobs[c_id, ]
  
  p = ggplot(delta_draws, aes(x = x, y = y)) +
    geom_density_2d_filled(alpha = 0.8) +
    
    #plot colony location
    geom_point(data = colony_data[c_id,], aes(x = colony_x, y = colony_y), colour = "lightgreen", show.legend = TRUE) +
    
    #plot trap locations / sizes / quality
    geom_point(data = trap_data, aes(x = trap_x, y = trap_y, size = trap_count, colour = fq)) +
    scale_colour_gradient(low = "white", high = "red") +
    scale_size_continuous(limits = c(0,10), range = c(1, 5)) +
    
    #miscellaneous
    ggtitle(paste("Colony", c_id)) +
    labs(title = paste("Colony", c_id), 
         size = "Number of Captures", 
         gradient = "Trap Floral Quality",
         level = "Colony Location & Posterior Distribution") +
    coord_equal() +
    theme_minimal()
  
  # save legend
  g <- ggplotGrob(p)
  legend_index <- which(g$layout$name == "guide-box-right")
  legend <- g$grobs[[legend_index]]
  
  # remove legend from plot
  p <- p + theme(legend.position = "none")
  
  #save plot
  plot_list[[i]] = p
  legends[[1]] = legend
}

fig = grid.arrange(grobs = plot_list, ncol = 3)
fig = grid.arrange(fig, legends[[1]], ncol = 2, widths = c(4,1))
ggsave(paste(results_path, "/colony_posteriorsGQ.jpg", sep = ""), fig, height = 3000, width = 4000, units = "px")
