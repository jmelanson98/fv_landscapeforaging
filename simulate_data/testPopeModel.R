##### Fit Pope Stan code to simulated data #####
# Script Initiated: April 10, 2025
# By: Jenna Melanson
# Goals:
    ### Fit Pope Stan code using a multiple iterations of simulated data
    ### 10 landscapes x 5 beta values x 3 sampling intensities
    ### Generate figure analogous to fig 1A + plots of colony posterior distributions
    ### Generate figures analogous to fig 3 + 4 from Pope & Jha

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
#setwd("~/Documents/UBC/bombus_project/fvlandscape_foraging") # local
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging") # server
rstan_options(auto_write = TRUE) 
options(mc.cores = parallel::detectCores())

##### Source functions #####
source("simulate_data/PopeSimFunctions.R")

# ####Replicate Figure 1A from Pope & Jha #####
# 
# # Make 4 simulations: 
# ## beta = -1/25 and theta = 0.25 (high distance constraint, low floral attraction)
# ## beta = -1/75 and theta = 0.25 (low distance constraint, low floral attraction)
# ## beta = -1/25 and theta = 0.5 (high distance constraint, high floral attraction)
# ## beta = -1/75 and theta = 0.5 (low distance constraint, high floral attraction)
# 
# set.seed(123) # set seed so that same colony locations are drawn each time
# sim1 = draw_N_bees(sample_size = 50,
#                     landscape_size = 1100,
#                     colonygrid_size = 700,
#                     trapgrid_size = 300,
#                     number_traps = 25,
#                     number_colonies = 10,
#                     colony_sizes = rep(100,10),
#                     beta = -1/25,
#                     theta = 0.25,
#                     resource_landscape = landscape,
#                     batch_size = 1)
# 
# set.seed(123)
# sim2 = draw_N_bees(sample_size = 50,
#                    landscape_size = 1100,
#                    colonygrid_size = 700,
#                    trapgrid_size = 300,
#                    number_traps = 25,
#                    number_colonies = 10,
#                    colony_sizes = rep(100,10),
#                    beta = -1/50,
#                    theta = 0.25,
#                    resource_landscape = landscape,
#                    batch_size = 1)
# set.seed(123)
# sim3 = draw_N_bees(sample_size = 50,
#                    landscape_size = 1100,
#                    colonygrid_size = 700,
#                    trapgrid_size = 300,
#                    number_traps = 25,
#                    number_colonies = 10,
#                    colony_sizes = rep(100,10),
#                    beta = -1/25,
#                    theta = 1,
#                    resource_landscape = landscape,
#                    batch_size = 1)
# set.seed(123)
# sim4 = draw_N_bees(sample_size = 50,
#                    landscape_size = 1100,
#                    colonygrid_size = 700,
#                    trapgrid_size = 300,
#                    number_traps = 25,
#                    number_colonies = 10,
#                    colony_sizes = rep(100,10),
#                    beta = -1/50,
#                    theta = 1,
#                    resource_landscape = landscape,
#                    batch_size = 1)
# 
# ##### Plot visitation rates + traps + true colony location #####
# 
# # first, chose a colony to track
# colony_chosen = 3
# #3 is pretty good and so is 7
# 
# # make visitation matrix into a df
# visitationdf1 = as.data.frame(sim1[[4]][[colony_chosen]]) %>%
#   mutate(x = 1:nrow(sim1[[4]][[colony_chosen]])) %>%  # Create y coordinate
#   pivot_longer(-x, names_to = "y", values_to = "value") %>%
#   mutate(y = as.integer(gsub("V", "", y))) %>% # Convert column names to integer x
#   filter(x > 300 & x < 800 & y > 300 & y < 800)
# 
# # plot
# sim1plot = ggplot() +
#   #plot visitation data
#   geom_tile(data = visitationdf1, aes(x = x, y = y, fill = value)) +
#   coord_equal() +
#   scale_fill_continuous(name = "Vistation Rate") +
#   
#   #plot colony location
#   geom_point(data = sim1[[2]][colony_chosen,], aes(x = colony_x, y = colony_y), colour = "lightgreen", show.legend = TRUE) +
#   
#   #plot trap locations / sizes / quality
#   geom_point(data = sim1[[3]], aes(x = trap_x, y = trap_y, size = sim1[[1]][colony_chosen,], colour = fq)) +
#   scale_colour_gradient(low = "lightpink", high = "red") +
#   scale_size_continuous(limits = c(0,20), range = c(1, 8)) +
#   
#   #other
#   labs(title = "Distance Constraint: High, Forage Attractiveness: Low", y = "Northing", x = "", colour = "Trap Forage Quality", size = "Number of captures", fill = "Visitation Rate") +
#   xlim(c(390,650)) +
#   ylim(c(525,725)) +
#   theme_minimal() +
#   theme(legend.position = "none")
# 
# # repeat for the next 3 simulations
# # convert to dataframes
# visitationdf2 = as.data.frame(sim2[[4]][[colony_chosen]]) %>%
#   mutate(x = 1:nrow(sim2[[4]][[colony_chosen]])) %>%  # Create y coordinate
#   pivot_longer(-x, names_to = "y", values_to = "value") %>%
#   mutate(y = as.integer(gsub("V", "", y))) %>% # Convert column names to integer x
#   filter(x > 300 & x < 800 & y > 300 & y < 800)
# 
# visitationdf3 = as.data.frame(sim3[[4]][[colony_chosen]]) %>%
#   mutate(x = 1:nrow(sim3[[4]][[colony_chosen]])) %>%  # Create y coordinate
#   pivot_longer(-x, names_to = "y", values_to = "value") %>%
#   mutate(y = as.integer(gsub("V", "", y))) %>% # Convert column names to integer x
#   filter(x > 300 & x < 800 & y > 300 & y < 800)
# 
# visitationdf4 = as.data.frame(sim4[[4]][[colony_chosen]]) %>%
#   mutate(x = 1:nrow(sim4[[4]][[colony_chosen]])) %>%  # Create y coordinate
#   pivot_longer(-x, names_to = "y", values_to = "value") %>%
#   mutate(y = as.integer(gsub("V", "", y))) %>% # Convert column names to integer x
#   filter(x > 300 & x < 800 & y > 300 & y < 800)
# 
# # plot
# sim2plot = ggplot() +
#   #plot visitation data
#   geom_tile(data = visitationdf2, aes(x = x, y = y, fill = value)) +
#   coord_equal() +
#   scale_fill_continuous(name = "Vistation Rate") +
#   
#   #plot colony location
#   geom_point(data = sim2[[2]][colony_chosen,], aes(x = colony_x, y = colony_y), colour = "lightgreen", show.legend = TRUE) +
#   
#   #plot trap locations / sizes / quality
#   geom_point(data = sim2[[3]], aes(x = trap_x, y = trap_y, size = sim2[[1]][colony_chosen,], colour = fq)) +
#   scale_colour_gradient(low = "lightpink", high = "red") +
#   scale_size_continuous(limits = c(0,20), range = c(1, 8)) +
#   
#   #other
#   labs(title = "Distance Constraint: Low, Forage Attractiveness: Low", y = "", x = "", colour = "Trap Forage Quality", size = "Number of captures", fill = "Visitation Rate") +
#   xlim(c(390,650)) +
#   ylim(c(525,725)) +
#   theme_minimal() +
#   theme(legend.position = "none")
# 
# sim3plot = ggplot() +
#   #plot visitation data
#   geom_tile(data = visitationdf3, aes(x = x, y = y, fill = value)) +
#   coord_equal() +
#   scale_fill_continuous(name = "Vistation Rate") +
#   
#   #plot colony location
#   geom_point(data = sim3[[2]][colony_chosen,], aes(x = colony_x, y = colony_y), colour = "lightgreen", show.legend = TRUE) +
#   
#   #plot trap locations / sizes / quality
#   geom_point(data = sim3[[3]], aes(x = trap_x, y = trap_y, size = sim3[[1]][colony_chosen,], colour = fq)) +
#   scale_colour_gradient(low = "lightpink", high = "red") +
#   scale_size_continuous(limits = c(0,20), range = c(1, 8)) +
#   
#   #other
#   labs(title = "Distance Constraint: High, Forage Attractiveness: High", y = "Northing", x = "Easting", colour = "Trap Forage Quality", size = "Number of captures", fill = "Visitation Rate") +
#   xlim(c(390,650)) +
#   ylim(c(525,725)) +
#   theme_minimal() +
#   theme(legend.position = "none")
# 
# sim4plot = ggplot() +
#   #plot visitation data
#   geom_tile(data = visitationdf4, aes(x = x, y = y, fill = value)) +
#   coord_equal() +
#   scale_fill_continuous(name = "Vistation Rate") +
#   
#   #plot colony location
#   geom_point(data = sim4[[2]][colony_chosen,], aes(x = colony_x, y = colony_y), colour = "lightgreen", show.legend = TRUE) +
#   
#   #plot trap locations / sizes / quality
#   geom_point(data = sim4[[3]], aes(x = trap_x, y = trap_y, size = sim4[[1]][colony_chosen,], colour = fq)) +
#   scale_colour_gradient(low = "lightpink", high = "red") +
#   scale_size_continuous(limits = c(0,20), range = c(1, 8)) +
#   
#   #other
#   labs(title = "Distance Constraint: Low, Forage Attractiveness: High", y = "", x = "Easting", colour = "Trap Forage Quality", size = "Number of captures", fill = "Visitation Rate") +
#   xlim(c(390,650)) +
#   ylim(c(525,725)) +
#   theme_minimal() +
#   theme(legend.position = "right")
# 
# 
# #save legend for plot grid
# g <- ggplotGrob(sim4plot)
# legend_index <- which(g$layout$name == "guide-box-right")
# legend <- g$grobs[[legend_index]]
# 
# # remove legend from plot
# sim4plot <- sim4plot + theme(legend.position = "none")
# 
# # make a grid plot!
# popefig1A = grid.arrange(sim1plot, sim2plot, 
#                          sim3plot, sim4plot,
#                          ncol = 2)
# popefig1A = grid.arrange(popefig1A, legend,
#                          ncol = 2, widths = c(4,1))
# 
# 
# # save plot
# ggsave("figures/simfigs/popeFig1A.jpg", popefig1A,
#        width = 3500, height = 2500, units = "px")



##### Replicate Figures 3 & 4 #####
## this must be run on server because my computer doesn't have the memory, RIP

# First, read in files for each task id
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
param_grid = readRDS("simulate_data/batch_sim1/param_grid.rds")
results_path = sprintf("simulate_data/batch_sim1/data/sim_result_%03d", task_id)
colony_data = readRDS(paste(results_path, "/colonydata.RDS", sep =""))
trap_data = readRDS(paste(results_path, "/trapdata.RDS", sep = ""))
yobs = readRDS(paste(results_path, "/yobs.RDS", sep=""))

# Then, format data for stan
data = list()
data$C = 1000
data$K = 25
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$y = yobs
data$lowerbound = 200
data$upperbound = 900
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 5

# Fit Stan model!
stanFit = stan(file = "models/pope_consgenetics.stan",
                    data = data, seed = 5838299,
                    warmup = 1000, iter = 10000,
                    chains = 4, cores = 4,
                    verbose = TRUE)
sprint("Model complete.")
saveRDS(stanFit, file=paste(results_path,"/stanFit.RDS", sep =""))
print("Model saved.")

#Plot the posteriors of ten colonies
plot_list = list()
poscols = colony_data$colonyid[rowSums(yobs) > 1]
numplots = 12
legends = list()
  
for (i in 1:numplots){
  c_id = poscols[i]
  delta_draws = as.data.frame(rstan::extract(stanFit, pars = "delta")$delta[, c_id,])
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
ggsave(paste(results_path, "/colony_posteriors.jpg", sep = ""), fig, height = 3000, width = 4000, units = "px")


##### Recreate Figure 3 from Pope & Jha #####

# for each iteration, draw:
    ## delta
    ## beta
    ## theta

# this is a 3D array: dim 1 = iteration, dim 2 = colony, dim 3 = x and y coordinates
deltas = rstan::extract(stanFit, pars = "delta")$delta

# these are 1D arrays: one per iteration
betas = rstan::extract(stanFit, pars = "beta")$beta
thetas = rstan::extract(stanFit, pars = "theta")$theta

# number of iterations
t = length(betas)

# observed colony sizes
colony_sizes = rowSums(yobs)

# initialize vector to record estimated foraging distance for each colony
estimated_foraging = c()

for (i in colony_data$colonyid){
  # initialize sum of foraging distance for all iterations
  d_i = 0
  
  for (iter in 1:t){
    #initialize foraging dist + normalizing factor for a single iteration and single colony
    d_it = 0
    udV = 0
    
    # calculate d_i for a single iteration and a single colony
    for (k in trap_data$trapid){
      dist_x = deltas[iter, colony, 1] - trap_data$trap_x[k]
      dist_y = deltas[iter, colony, 2] - trap_data$trap_y[k]
      dist = sqrt(dist_x^2 + dist_y^2)
      
      lambdaik = betas[iter]*dist + thetas[iter]*trap_data$fq[k]
      
      d_it = d_it + dist*lambdaik
      udV = udV + lambdaik
    }
    
    #normalize by big lambda
    d_it = d_it/udV
    
    # sum over t
    d_i = d_i + d_it
  }
  
  # compute average d_i
  d_i = d_i/t
  
  #add to vector
  estimated_foraging = c(estimated_foraging, d_i)
}



# Check if all_sim_df has been created yet; if not, initialize it

if (!file.exists("simulate_data/batch_sim1/all_sim_df.RDS")){
  all_sim_df = param_grid
  all_sim_df$colony_foraging_estimates = rep(NA,150)
  all_sim_df$colony_sizes = rep(NA,150)
  all_sim_df$average_foraging = rep(NA,150)
}

all_sim_df$colony_foraging_estimates[50] = estimated_foraging