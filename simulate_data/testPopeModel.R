##### Test Pope Stan code on simulated data #####
# Script Initiated: April 10, 2025
# By: Jenna Melanson
# Goals:
    ### Test Pope Stan code using a single iteration of simulated data
    ### Test on 10 landscapes x 5 beta values x 3 sampling intensities
    ### Generate figures analagous to figs 1 + 3 from Pope conservation genetics paper

##### Set Environment #####
#setwd("~/Documents/UBC/bombus_project/fvlandscape_foraging") # local
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging") # server
rstan_options(auto_write = TRUE) 
options(mc.cores = parallel::detectCores())

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
stanFit = stan(file = "models/popeModified.stan",
                    data = data, seed = 5838299,
                    warmup = 1000, iter = 2025,
                    chains = 4, cores = 4,
                    verbose = TRUE)
save(stanFit, file=paste(results_path,"stanFit.Rdata", sep =""))


# Run posterior predictive check
draws = rstan::extract(stanFit, pars = "y_rep")$y_rep
dens_plot = pp_check(data$y, yrep = y_rep, fun = "dens_overlay")
ggsave(paste(results_path, "ppc.jpg", sep = ""), dens_plot, height = 1000, width = 1000, units = "px")

#Plot the posteriors of ten colonies
plot_list = list()
for (i in 1:10){
  delta_draws = rstan::extract(stanFit, pars = "delta")$delta[, i,]
  colnames(delta_draws = x,y)
  
  p = ggplot(delta_draws, aes(x = x, y = y)) +
    geom_density_2d_filled(alpha = 0.8) +
    
    #plot colony location
    geom_point(data = colony_data[i,], aes(x = colony_x, y = colony_y), colour = "lightgreen", show.legend = TRUE) +
    
    #plot trap locations / sizes / quality
    geom_point(data = trap_data, aes(x = trap_x, y = trap_y, size = yobs[i,], colour = fq)) +
    scale_colour_gradient(low = "lightpink", high = "red") +
    scale_size_continuous(limits = c(0,20), range = c(1, 8)) +
    
    coord_equal() +
    theme_minimal()
  
  #save plot
  plot_list[[i]] = p
}

fig = grid.arrange(grobs = plot_list, ncol = 2)
ggsave(paste(results_path, "colony_posteriors.jpg", sep = ""), fit, height = 5000, width = 2000, units = "px")


