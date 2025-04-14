##### Test Pope Stan code on simulated data #####
# Script Initiated: April 10, 2025
# By: Jenna Melanson
# Goals:
    ### Test Pope Stan code using a single iteration of simulated data
    ### Test on 10 landscapes x 5 beta values x 3 sampling intensities
    ### Generate figures analagous to figs 1 + 3 from Pope conservation genetics paper

##### Set Environment #####
setwd("~/Documents/UBC/bombus_project/fvlandscape_foraging")

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


##### Generate One Iteration #####
landscape = simulateLandscape(1100, 10)
iter1 = draw_N_bees(sample_size = 1000,
                    landscape_size = 1100,
                    colonygrid_size = 700,
                    trapgrid_size = 300,
                    number_traps = 25,
                    number_colonies = 500,
                    colony_sizes = rep(100,500),
                    beta = -1/50,
                    theta = 0.5,
                    resource_landscape = landscape,
                    batch_size = 1)


##### Format Data #####
data = list()
data$C = 500
data$K = 25
data$trap = as.matrix(cbind(iter1[[3]]$trap_x, iter1[[3]]$trap_y))
data$y = iter1[[1]]
data$lowerbound = 200
data$upperbound = 900
data$floral = iter1[[3]]$fq
data$priorVa = 1
data$priorCo = 5

##### Fit Stan Model ! #####
one_iter_fit = stan(file = "models/popeModified.stan",
                  data = data, seed = 5838299,
                  warmup = 1000, iter = 2025,
                  chains = 4, cores = 4,
                  verbose = TRUE)

summary(one_iter_fit)
traceplot(one_iter_fit)

##### Run Posterior Predictive Check #####
draws = rstan::extract(one_iter_fit, pars = "y_rep")$y_rep
pp_check(data$y, yrep = y_rep, fun = "dens_overlay")

####Replicate Figure 1A from Pope & Jha #####

# Make 4 simulations: 
## beta = -1/25 and theta = 0.25 (high distance constraint, low floral attraction)
## beta = -1/75 and theta = 0.25 (low distance constraint, low floral attraction)
## beta = -1/25 and theta = 0.5 (high distance constraint, high floral attraction)
## beta = -1/75 and theta = 0.5 (low distance constraint, high floral attraction)

set.seed(123) # set seed so that same colony locations are drawn each time
sim1 = draw_N_bees(sample_size = 50,
                    landscape_size = 1100,
                    colonygrid_size = 700,
                    trapgrid_size = 300,
                    number_traps = 25,
                    number_colonies = 10,
                    colony_sizes = rep(100,10),
                    beta = -1/25,
                    theta = 0.25,
                    resource_landscape = landscape,
                    batch_size = 1)

set.seed(123)
sim2 = draw_N_bees(sample_size = 50,
                   landscape_size = 1100,
                   colonygrid_size = 700,
                   trapgrid_size = 300,
                   number_traps = 25,
                   number_colonies = 10,
                   colony_sizes = rep(100,10),
                   beta = -1/50,
                   theta = 0.25,
                   resource_landscape = landscape,
                   batch_size = 1)
set.seed(123)
sim3 = draw_N_bees(sample_size = 50,
                   landscape_size = 1100,
                   colonygrid_size = 700,
                   trapgrid_size = 300,
                   number_traps = 25,
                   number_colonies = 10,
                   colony_sizes = rep(100,10),
                   beta = -1/25,
                   theta = 1,
                   resource_landscape = landscape,
                   batch_size = 1)
set.seed(123)
sim4 = draw_N_bees(sample_size = 50,
                   landscape_size = 1100,
                   colonygrid_size = 700,
                   trapgrid_size = 300,
                   number_traps = 25,
                   number_colonies = 10,
                   colony_sizes = rep(100,10),
                   beta = -1/50,
                   theta = 1,
                   resource_landscape = landscape,
                   batch_size = 1)

##### Plot visitation rates + traps + true colony location #####

# first, chose a colony to track
colony_chosen = 3
#3 is pretty good and so is 7

# make visitation matrix into a df
visitationdf1 = as.data.frame(sim1[[4]][[colony_chosen]]) %>%
  mutate(x = 1:nrow(sim1[[4]][[colony_chosen]])) %>%  # Create y coordinate
  pivot_longer(-x, names_to = "y", values_to = "value") %>%
  mutate(y = as.integer(gsub("V", "", y))) %>% # Convert column names to integer x
  filter(x > 300 & x < 800 & y > 300 & y < 800)

# plot
sim1plot = ggplot() +
  #plot visitation data
  geom_tile(data = visitationdf1, aes(x = x, y = y, fill = value)) +
  coord_equal() +
  scale_fill_continuous(name = "Vistation Rate") +
  
  #plot colony location
  geom_point(data = sim1[[2]][colony_chosen,], aes(x = colony_x, y = colony_y), colour = "lightgreen", show.legend = TRUE) +
  
  #plot trap locations / sizes / quality
  geom_point(data = sim1[[3]], aes(x = trap_x, y = trap_y, size = sim1[[1]][colony_chosen,], colour = fq)) +
  scale_colour_gradient(low = "lightpink", high = "red") +
  scale_size_continuous(limits = c(0,20), range = c(1, 8)) +
  
  #other
  labs(title = "Distance Constraint: High, Forage Attractiveness: Low", y = "Northing", x = "", colour = "Trap Forage Quality", size = "Number of captures", fill = "Visitation Rate") +
  xlim(c(390,650)) +
  ylim(c(525,725)) +
  theme_minimal() +
  theme(legend.position = "none")

# repeat for the next 3 simulations
# convert to dataframes
visitationdf2 = as.data.frame(sim2[[4]][[colony_chosen]]) %>%
  mutate(x = 1:nrow(sim2[[4]][[colony_chosen]])) %>%  # Create y coordinate
  pivot_longer(-x, names_to = "y", values_to = "value") %>%
  mutate(y = as.integer(gsub("V", "", y))) %>% # Convert column names to integer x
  filter(x > 300 & x < 800 & y > 300 & y < 800)

visitationdf3 = as.data.frame(sim3[[4]][[colony_chosen]]) %>%
  mutate(x = 1:nrow(sim3[[4]][[colony_chosen]])) %>%  # Create y coordinate
  pivot_longer(-x, names_to = "y", values_to = "value") %>%
  mutate(y = as.integer(gsub("V", "", y))) %>% # Convert column names to integer x
  filter(x > 300 & x < 800 & y > 300 & y < 800)

visitationdf4 = as.data.frame(sim4[[4]][[colony_chosen]]) %>%
  mutate(x = 1:nrow(sim4[[4]][[colony_chosen]])) %>%  # Create y coordinate
  pivot_longer(-x, names_to = "y", values_to = "value") %>%
  mutate(y = as.integer(gsub("V", "", y))) %>% # Convert column names to integer x
  filter(x > 300 & x < 800 & y > 300 & y < 800)

# plot
sim2plot = ggplot() +
  #plot visitation data
  geom_tile(data = visitationdf2, aes(x = x, y = y, fill = value)) +
  coord_equal() +
  scale_fill_continuous(name = "Vistation Rate") +
  
  #plot colony location
  geom_point(data = sim2[[2]][colony_chosen,], aes(x = colony_x, y = colony_y), colour = "lightgreen", show.legend = TRUE) +
  
  #plot trap locations / sizes / quality
  geom_point(data = sim2[[3]], aes(x = trap_x, y = trap_y, size = sim2[[1]][colony_chosen,], colour = fq)) +
  scale_colour_gradient(low = "lightpink", high = "red") +
  scale_size_continuous(limits = c(0,20), range = c(1, 8)) +
  
  #other
  labs(title = "Distance Constraint: Low, Forage Attractiveness: Low", y = "", x = "", colour = "Trap Forage Quality", size = "Number of captures", fill = "Visitation Rate") +
  xlim(c(390,650)) +
  ylim(c(525,725)) +
  theme_minimal() +
  theme(legend.position = "none")

sim3plot = ggplot() +
  #plot visitation data
  geom_tile(data = visitationdf3, aes(x = x, y = y, fill = value)) +
  coord_equal() +
  scale_fill_continuous(name = "Vistation Rate") +
  
  #plot colony location
  geom_point(data = sim3[[2]][colony_chosen,], aes(x = colony_x, y = colony_y), colour = "lightgreen", show.legend = TRUE) +
  
  #plot trap locations / sizes / quality
  geom_point(data = sim3[[3]], aes(x = trap_x, y = trap_y, size = sim3[[1]][colony_chosen,], colour = fq)) +
  scale_colour_gradient(low = "lightpink", high = "red") +
  scale_size_continuous(limits = c(0,20), range = c(1, 8)) +
  
  #other
  labs(title = "Distance Constraint: High, Forage Attractiveness: High", y = "Northing", x = "Easting", colour = "Trap Forage Quality", size = "Number of captures", fill = "Visitation Rate") +
  xlim(c(390,650)) +
  ylim(c(525,725)) +
  theme_minimal() +
  theme(legend.position = "none")

sim4plot = ggplot() +
  #plot visitation data
  geom_tile(data = visitationdf4, aes(x = x, y = y, fill = value)) +
  coord_equal() +
  scale_fill_continuous(name = "Vistation Rate") +
  
  #plot colony location
  geom_point(data = sim4[[2]][colony_chosen,], aes(x = colony_x, y = colony_y), colour = "lightgreen", show.legend = TRUE) +
  
  #plot trap locations / sizes / quality
  geom_point(data = sim4[[3]], aes(x = trap_x, y = trap_y, size = sim4[[1]][colony_chosen,], colour = fq)) +
  scale_colour_gradient(low = "lightpink", high = "red") +
  scale_size_continuous(limits = c(0,20), range = c(1, 8)) +
  
  #other
  labs(title = "Distance Constraint: Low, Forage Attractiveness: High", y = "", x = "Easting", colour = "Trap Forage Quality", size = "Number of captures", fill = "Visitation Rate") +
  xlim(c(390,650)) +
  ylim(c(525,725)) +
  theme_minimal() +
  theme(legend.position = "right")


#save legend for plot grid
g <- ggplotGrob(sim4plot)
legend_index <- which(g$layout$name == "guide-box-right")
legend <- g$grobs[[legend_index]]

# remove legend from plot
sim4plot <- sim4plot + theme(legend.position = "none")

# make a grid plot!
popefig1A = grid.arrange(sim1plot, sim2plot, 
                         sim3plot, sim4plot,
                         ncol = 2)
popefig1A = grid.arrange(popefig1A, legend,
                         ncol = 2, widths = c(4,1))


# save plot
ggsave("figures/simfigs/popeFig1A.jpg", popefig1A,
       width = 3500, height = 2500, units = "px")
