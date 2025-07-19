##### Run models with and without unobserved colonies #####
# Script Initiated: June 14, 2025
# By: Jenna Melanson
# Goals:
### Determine whether removing unobserved colonies from simulated data
### has any bearing on model fit, specifically global intercept mu


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
#library(ggpubr)
library(terra)

##### Set Environment #####
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
#setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
rstan_options(auto_write = TRUE) 
options(mc.cores = parallel::detectCores())

##### Source functions #####
source("simulate_data/src/GeneralizedSimFunctions.R")

##### Simulate some data #####
# use 2000 background colonies and 1000 observations?
# this distribution won't be totally accurate but it's hard to fit a model with 8000 background colonies

# Get landscape from saved file
fq <- readRDS("simulate_data/landscapes/random_field_range10/landscape_003.rds")

# Run simulation
#result <- draw_bees_colony_restricted(
#  sample_size     = 1000,
#  landscape_size  = 1500,
#  colonygrid_size = 700,
#  trapgrid_size   = 300,
#  number_traps    = 25,
#  number_colonies = 2000,
#  colony_sizes    = rep(100, 2000),
#  rho            = 50,
#  theta           = 0.5,
#  resource_landscape = fq,
#  nesting_landscape = NULL,
#  distance_decay = "exponential"
#)

# Save results
#saveRDS(result, "simulate_data/foraging_distance/methods_comparison/observed_vs_unobserved/simdata.rds")
result = readRDS("simulate_data/foraging_distance/methods_comparison/observed_vs_unobserved/simdata.rds")

# Write outputs to variables
yobs = result[[1]]
colony_data = result[[2]]
trap_data = result[[3]]

# Make some subsets within and without zero colonies
yobs_detected = yobs[rowSums(yobs) > 0,] 

# Prep data list for Stan
data = list()
data$K = 25
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$lowerbound = 0
data$upperbound = 1500
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$rho_center = 4.5
data$rho_sd = 0.5
data$y = yobs_detected
data$C = nrow(data$y)


#select stan model to fit
stanfile = paste("models/exponential.stan")

#fit and save model
#stanFitObs_widerlimit = stan(file = stanfile,
#               data = data, seed = 5838299,
#               chains = 4, cores = 4,
#               iter = 4000, warmup = 1000,
#               verbose = TRUE)
#saveRDS(stanFitObs_widerlimit, file="simulate_data/foraging_distance/methods_comparison/observed_vs_unobserved/stanFitObsWiderLimit.RDS")
stanFitObs_widerlimit = readRDS("simulate_data/foraging_distance/methods_comparison/observed_vs_unobserved/stanFitObs_widerlimit.RDS")
# when I ran these, the generated quantities block of exponential.stan had an error :(
# so ignore the colony dist estimates, they're all wrong

# okay so initial thoughts...maintaining the zeros makes a HUGE difference
# for all colonies included, rho = 53 
# for only observed colonies, rho = 199
# (true rho = 50)
# this seems too extreme to me true....going to rerun the stanFitObserved
# second run -- same results

# plot colony posteriors
#Plot the posteriors of ten colonies
plot_list = list()
numplots = 6
legends = list()
colony_observed = colony_data[rowSums(yobs) > 0,]

for (i in 1:numplots){
  c_id = colony_observed$colonyid[i]
  delta_draws = as.data.frame(rstan::extract(stanFitObs_widerlimit, pars = "delta")$delta[, c_id,])
  colnames(delta_draws) = c("x","y")
  trap_data$trap_count = yobs_detected[c_id, ]
  
  p = ggplot(delta_draws, aes(x = x, y = y)) +
    geom_density_2d_filled(alpha = 0.8) +
    
    #plot colony location
    geom_point(data = colony_observed[c_id,], aes(x = colony_x, y = colony_y), colour = "lightgreen", show.legend = TRUE) +
    
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
  #g <- ggplotGrob(p)
  #legend_index <- which(g$layout$name == "guide-box-right")
  #legend <- g$grobs[[legend_index]]
  
  # remove legend from plot
  #p <- p + theme(legend.position = "none")
  
  #save plot
  plot_list[[i]] = p
  #legends[[1]] = legend
}

fig = grid.arrange(grobs = plot_list, ncol = 3)
#fig = grid.arrange(fig, legends[[1]], ncol = 2, widths = c(4,1))
ggsave("figures/simfigs/colonyposteriors_widelimits.png", fig, height = 750, width = 1500)


# # check worker distributions
# # prep real data
# allspecs = read.csv("data/siblingships/allsibships_cleaned.csv")
# mixsum = allspecs %>% filter(final_id == "B. mixtus") %>%
#   filter(!is.na(ClusterIndex)) %>%
#   group_by(ClusterIndex) %>%
#   summarize(n = n())
# impsum = allspecs %>% filter(final_id == "B. impatiens") %>%
#   filter(!is.na(ClusterIndex)) %>%
#   group_by(ClusterIndex) %>%
#   summarize(n = n())
# 
# # simulated data
# simsum = data.frame(n = rowSums(yobs_detected))
# 
# # combine
# mix = data.frame(n = mixsum$n)
# imp = data.frame(n = impsum$n)
# df_combined = bind_rows(
#   simsum %>% mutate(group = "sim"),
#   mix %>% mutate(group = "mix"),
#   imp %>% mutate(group = "imp")
# )
# 
# # bin manually
# df_props <- df_combined %>%
#   group_by(group, n) %>%
#   summarise(count = n(), .groups = "drop") %>%
#   group_by(group) %>%
#   mutate(proportion = count / sum(count))
# 
# # plot
# ggplot(df_props, aes(x = n, y = proportion, color = group)) +
#   geom_segment(aes(x = n - 0.5, xend = n + 0.5), linewidth = 0.6) +
#   geom_point(size = 3) +
#   labs(
#     x = "Number of Sibs",
#     y = "Proportion"
#   ) +
#   theme_minimal() +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
