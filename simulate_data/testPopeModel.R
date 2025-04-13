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
sim1 = draw_N_bees(sample_size = 1000,
                    landscape_size = 1100,
                    colonygrid_size = 700,
                    trapgrid_size = 300,
                    number_traps = 25,
                    number_colonies = 500,
                    colony_sizes = rep(100,500),
                    beta = -1/25,
                    theta = 0.25,
                    resource_landscape = landscape,
                    batch_size = 1)
set.seed(123)
sim2 = draw_N_bees(sample_size = 1000,
                   landscape_size = 1100,
                   colonygrid_size = 700,
                   trapgrid_size = 300,
                   number_traps = 25,
                   number_colonies = 500,
                   colony_sizes = rep(100,500),
                   beta = -1/75,
                   theta = 0.25,
                   resource_landscape = landscape,
                   batch_size = 1)
set.seed(123)
sim3 = draw_N_bees(sample_size = 1000,
                   landscape_size = 1100,
                   colonygrid_size = 700,
                   trapgrid_size = 300,
                   number_traps = 25,
                   number_colonies = 500,
                   colony_sizes = rep(100,500),
                   beta = -1/25,
                   theta = 0.5,
                   resource_landscape = landscape,
                   batch_size = 1)
set.seed(123)
sim4 = draw_N_bees(sample_size = 1000,
                   landscape_size = 1100,
                   colonygrid_size = 700,
                   trapgrid_size = 300,
                   number_traps = 25,
                   number_colonies = 500,
                   colony_sizes = rep(100,500),
                   beta = -1/75,
                   theta = 0.5,
                   resource_landscape = landscape,
                   batch_size = 1)


##### Plot visitation rates + traps + true colony location #####

# first, chose a colony to track
colony_chosen = 1

vraster1 = sim1[[4]]
ggplot(data = iter1[[3]], aes(x = trap_x, y = trap_y, size = iter1[[2]][colony_chosen,])) +
  geom_point() +
  geom_point(data = iter[[2]][colony_chosen], aes(x = colony_x, y = colony_y, colour = "red") +
  geom_raster(vraster1)