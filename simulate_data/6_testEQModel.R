##### Fit exponentiated quadratic code to simulated data #####
# Script Initiated: May 9, 2025
# By: Jenna Melanson
# Goals:
### Fit exponentiated quadratic using multiple iterations of simulated data
### 10 landscapes x 5 rho values x 3 underlying colony densities
### Plot colony posteriors
### Plot colony foraging distance estimates

stan_rhat
stan_ess
traceplot(exfit, pars = "sigma", inc_warmup = TRUE)


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
source("simulate_data/0_EQSimFunctions.R")


##### Fit Models #####
## this must be run on server because my computer doesn't have the memory, RIP

# First, read in files for each task id
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
param_grid = readRDS("simulate_data/exponentiated_quadratic_sim/param_grid.rds")
results_path = sprintf("simulate_data/exponentiated_quadratic_sim/data/sim_result_%03d", task_id)
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
data$rho_center = 100
data$rho_sd_log = 0.5

# Fit Stan model!
stanFit = stan(file = "models/EQmodel.stan",
               data = data, seed = 5838299,
               chains = 4, cores = 4,
               verbose = TRUE)
print("Model complete.")
saveRDS(stanFit, file=paste(results_path,"/stanFit.RDS", sep =""))
print("Model saved.")




##### Check for errors in stan fit #####
# list all .err files
err_file <- sprintf("simulate_data/logs/stanEQ_%s_%s.out",
        Sys.getenv("SLURM_ARRAY_JOB_ID"),
        Sys.getenv("SLURM_ARRAY_TASK_ID"))
lines <- readLines(err_file, warn = FALSE)
error_check = list(
    unfinished = any(grepl("DUE TO TIME LIMIT", lines, ignore.case = TRUE)),
    has_divergent = any(grepl("divergent", lines, ignore.case = TRUE)),
    has_low_ess = any(grepl("Effective Samples Size", lines, ignore.case = TRUE)),
    has_error = any(grepl("error", lines, ignore.case = TRUE))
  )
if (any(unlist(error_check))) {
  print("Errors detected in log:")
  print(Filter(identity, error_check))
  stop("Aborting task due to detected log errors.")
}


##### Plot the posteriors of some colonies #####
plot_list = list()
numplots = 12
legends = list()

for (i in 1:numplots){
  c_id = colony_data$colonyid[i]
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




##### Save simulated vs model-predicted colony average foraging distance #####

# check if the all sim data frame exists, and if not, create one
if (file.exists("simulate_data/exponentiated_quadratic_sim/all_simulation_summary.csv")){
  allsim_colonies = read.csv("simulate_data/exponentiated_quadratic_sim/all_simulation_summary.csv", row.names = 1)
} else{
  columns = c("colony_density", "rho", "landscape_id", "colony_size_bin", 
              "true_colony_avg", "true_colony_sd", "model_colony_avg", "model_colony_sd")
  allsim_colonies = data.frame(matrix(nrow =0, ncol = length(columns))) 
  colnames(allsim_colonies) = columns
}


# extract model estimates
colony_data$model_estimate = summary(stanFit, pars = c("colony_dist"))$summary[,1]

# record observed colony sizes and create categorical variable
colony_data$observed_size = rowSums(yobs)

# create categorical variable for colony size
if(max(colony_data$observed_size > 3)){
  colony_data$size_bin = cut(colony_data$observed_size, breaks = c(0, 3, max(colony_data$observed_size)), labels = c("1-3", "4+"))
  dim = 2
  bin_vec = c("1-3", "4+")
} else {
  colony_data$size_bin = rep(c("1-3"), nrow(colony_data))
  dim = 1
  bin_vec = c("1-3")
}

# make a temporary df and add to allsim_colonies
columns = c("colony_density", "rho", "landscape_id", "colony_size_bin", 
            "true_colony_avg", "true_colony_sd", "model_colony_avg", "model_colony_sd")
temp_df = data.frame(matrix(nrow = dim, ncol = length(columns))) 
colnames(temp_df) = columns

temp_df$colony_density = param_grid$sample_size[task_id]
temp_df$rho = param_grid$rho[task_id]
temp_df$landscape_id = param_grid$landscape_id[task_id]
temp_df$colony_size_bin = bin_vec


for (bin in bin_vec){
    temp_df$true_colony_avg[temp_df$colony_size_bin == bin] = mean(colony_data$foraging_range[colony_data$size_bin == bin])
    temp_df$true_colony_sd[temp_df$colony_size_bin == bin] = sd(colony_data$foraging_range[colony_data$size_bin == bin])
    temp_df$model_colony_avg[temp_df$colony_size_bin == bin] = mean(colony_data$model_estimate[colony_data$size_bin == bin])
    temp_df$model_colony_sd[temp_df$colony_size_bin == bin] = sd(colony_data$model_estimate[colony_data$size_bin == bin])
    
    # save rhos and thetas
    temp_df$rho_mean[temp_df$colony_size_bin == bin] = summary(stanFit, pars = c("rho"))$summary[,1]
    temp_df$rho_sd[temp_df$colony_size_bin == bin] = summary(stanFit, pars = c("rho"))$summary[,3]
    temp_df$theta_mean[temp_df$colony_size_bin == bin] = summary(stanFit, pars = c("theta"))$summary[,1]
    temp_df$theta_sd[temp_df$colony_size_bin == bin] = summary(stanFit, pars = c("theta"))$summary[,3]
}

allsim_colonies = rbind(allsim_colonies, temp_df)
write.csv(allsim_colonies, "simulate_data/exponentiated_quadratic_sim/all_simulation_summary.csv")





##### Plot figure 3 from Pope & Jha #####
allsim_colonies = read.csv("simulate_data/exponentiated_quadratic_sim/all_simulation_summary.csv")

allsim_summary = allsim_colonies %>% group_by(colony_density, rho, colony_size_bin) %>%
  summarize(true_avg = mean(true_colony_avg, na.rm = TRUE),
            true_sd = sqrt(sum(true_colony_sd^2, na.rm = TRUE)),
            model_avg = mean(model_colony_avg, na.rm = TRUE),
            model_sd = sqrt(sum(model_colony_sd^2, na.rm = TRUE))
  ) %>%
  ungroup()

fig3a = ggplot(as.data.frame(allsim_summary[allsim_summary$colony_density == 1000,]), 
                  aes(x = true_avg, y = model_avg,
                      color = colony_size_bin)) +
  geom_point() +
  geom_errorbar(aes(ymin = model_avg - 2*model_sd, ymax = model_avg + 2*model_sd)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "True Foraging Distance",
       y = "Modelled Foraging Distance",
       color = "Number of Siblings",
       title = "underlying colony density = 1000 colonies") +
  theme_minimal()

fig3b = ggplot(as.data.frame(allsim_summary[allsim_summary$colony_density == 2000,]), 
                 aes(x = true_avg, y = model_avg,
                     color = colony_size_bin)) +
  geom_point() +
  geom_errorbar(aes(ymin = model_avg - 2*model_sd, ymax = model_avg + 2*model_sd)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "True Foraging Distance",
       y = "Modelled Foraging Distance",
       color = "Number of Siblings",
       title = "underlying colony density = 2000 colonies") +
  theme_minimal()

fig3c = ggplot(as.data.frame(allsim_summary[allsim_summary$colony_density == 4000,]), 
                 aes(x = true_avg, y = model_avg,
                     color = colony_size_bin)) +
  geom_point() +
  geom_errorbar(aes(ymin = model_avg - 2*model_sd, ymax = model_avg + 2*model_sd)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "True Foraging Distance",
       y = "Modelled Foraging Distance",
       color = "Number of Siblings",
       title = "underlying colony density = 4000 colonies") +
  theme_minimal()

grid = grid.arrange(fig3a, fig3b, fig3c, ncol = 1)
ggsave("figures/simfigs/exponentiated_quadratic_estimates.jpg", grid, width = 2000, height = 3000, units = 'px')
