##### Accuracy comparison of different simulation/modeling scenarios #####
##### THIS SCRIPT IS FOR THE MODEL FITTING AND PLOTTING PHASE ######
# Script Initiated: May 26, 2025
# By: Jenna Melanson
# Goals:
### use (1) higher sample size, (2) accurate background colony density, (3) larger resource landscape to compare conditions:
##### Exponentiated quadratic (all colonies)
##### Exponential (all colonies)
##### exponentiated quadratic (all nonzero colonies)
##### exponential (all nonzero colonies)
##### exponentiated quadratic (doubleton+ colonies)
##### exponential (doubleton+ colonies)
##### centroid approach on exponential data (doubleton+)
##### centroid approach on exp. quadratic data (doubleton+)


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
library(terra)
library(ggpubr)
library(grid.arrange)


##### Set Environment #####
#setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging") # server
options(mc.cores = parallel::detectCores())


##### Source functions #####
source("simulate_data/src/GeneralizedSimFunctions.R")


##### Load in param grid and simulated data #####
param_grid = readRDS("simulate_data/methods_comparison/param_grid.rds")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
inputfilepath <- sprintf("simulate_data/methods_comparison/data/sim_result_%03d", task_id)
yobs = readRDS(paste(inputfilepath, "/yobs.RDS", sep = ""))
colony_data = readRDS(paste(inputfilepath, "/colonydata.RDS", sep = ""))
trap_data = readRDS(paste(inputfilepath, "/trapdata.RDS", sep = ""))


##### Fit Models ######

### Format most basic data for Stan
data = list()
data$K = 25
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$lowerbound = 400
data$upperbound = 1100
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$rho_center = 4.5
data$rho_sd = 0.5


### Current params
current_params = param_grid[task_id,]


### Correctly subset yobs based on approach
if (current_params$model_approach == "all"){
  data$y = yobs
  data$C = nrow(yobs)
} else if (current_params$model_approach == "singletons"){
  colony_data = colony_data[rowSums(yobs) > 0,]
  yobs = yobs[rowSums(yobs) > 0,]
  data$y = yobs
  data$C = nrow(yobs)
} else if (current_params$model_approach == "doubletons"){
  colony_data = colony_data[rowSums(yobs) > 1,]
  yobs = yobs[rowSums(yobs) > 1,]
  data$y = yobs
  data$C = nrow(yobs)
} else if (current_params$model_approach == "centroid"){
  colony_data = colony_data[rowSums(yobs) > 1,]
  yobs = yobs[rowSums(yobs) > 1,]
  data$y = yobs
  data$C = nrow(yobs)
} else {
  print("Not a valid approach.")
}



### Fit Stan model if necessary
if (current_params$model_approach != "centroid"){
  #select stan model to fit
  stanfile = paste("models/", current_params$distance_decay, ".stan", sep = "")
  
  #fit and save model
  stanFit = stan(file = stanfile,
                 data = data, seed = 5838299,
                 chains = 4, cores = 4,
                 iter = 10000,
                 verbose = TRUE)
  print("Model complete.")
  saveRDS(stanFit, file=paste(inputfilepath,"/stanFit.RDS", sep =""))
  print("Model saved.")
}



##### Estimate and plot landscape-level foraging distance ######

if (current_params$model_approach == "centroid"){
  get_avg_distance_to_centroid <- function(counts, coords) {
    total <- sum(counts)
    if (total == 0) return(NA)
    
    # calculate weighted centroid
    centroid_x <- sum(counts * coords$trap_x) / total
    centroid_y <- sum(counts * coords$trap_y) / total
    
    # calculate distance of each trap from centroid
    dists <- sqrt((coords$trap_x - centroid_x)^2 + (coords$trap_y - centroid_y)^2)
    
    # calculate weighted average distance to centroid
    mean_dist <- sum(dists * counts) / total
    return(mean_dist)
  }
  
  #apply function to yobs
  avg_dists <- apply(yobs, 1, function(counts) {
    get_avg_distance_to_centroid(counts, trap_data)
  })
  
  # save in param grid
  param_grid$model_average_foraging = mean(avg_dists)
  param_grid$model_sd_foraging = sd(avg_dists)
  
} else {
  colony_data$model_estimate = summary(stanFit, pars = c("colony_dist"))$summary[,1]
  param_grid$model_average_foraging = mean(colony_data$model_estimate)
  param_grid$model_sd_foraging = sd(colony_data$model_estimate)
  param_grid$model_mu = summary(stanFit, pars = c("mu"))$summary[,1]
}

saveRDS(param_grid, "simulate_data/methods_comparison/param_grid.rds")


###### Plot some figures #######
summary = param_grid %>% group_by(rho, distance_decay, model_approach) %>%
  summarize(true_avg = mean(true_average_foraging, na.rm = TRUE),
            true_sd = sqrt(sum(true_sd_foraging^2, na.rm = TRUE)),
            model_avg = mean(model_average_foraging, na.rm = TRUE),
            model_sd = sqrt(sum(model_sd_foraging^2, na.rm = TRUE))
  ) %>%
  ungroup()

summary$technique = paste(summary$distance_decay, summary$model_approach, sep = " / ")


#### Accuracy plot
fig = ggplot(summary, aes(x = true_avg, y = model_avg,
                          color = technique)) +
  geom_point() +
  geom_errorbar(aes(ymin = model_avg - 2*model_sd, ymax = model_avg + 2*model_sd)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "True Foraging Distance",
       y = "Modelled Foraging Distance",
       color = "Approach") +
  theme_minimal()
ggsave(ggsave("figures/simfigs/methods_comparison_accuracy.jpg", fig, width = 2000, height = 1000, units = 'px'))

#### Worker distribution plots
# prep real data
allspecs = read.csv("data/siblingships/allsibships_cleaned.csv")
mixsum = allspecs %>% filter(final_id == "B. mixtus") %>%
  filter(!is.na(ClusterIndex)) %>%
  group_by(ClusterIndex) %>%
  summarize(n = n())
impsum = allspecs %>% filter(final_id == "B. impatiens") %>%
  filter(!is.na(ClusterIndex)) %>%
  group_by(ClusterIndex) %>%
  summarize(n = n())

mix_df = data.frame(dataset_id = rep("mix", length(mixsum$n)),
                    value = mixsum$n,
                    type = rep("real", length(mixsum$n)),
                    rho = rep(NA, length(mixsum$n)),
                    distance_decay = rep(NA, length(mixsum$n)))

imp_df = data.frame(dataset_id = rep("imp", length(impsum$n)),
                    value = impsum$n,
                    type = rep("real", length(impsum$n)),
                    rho = rep(NA, length(impsum$n)),
                    distance_decay = rep(NA, length(impsum$n)))

real_df = rbind(mix_df, imp_df)

# Prep simulated data
columns = c("dataset_id", "value", "type", "rho", "distance_decay") 
distribution_df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(distribution_df) = columns


for (sim in rownames(param_grid)){
  sim = as.numeric(sim)
  counts = unlist(param_grid$counts[sim])
  temp_df = data.frame(rep(sim, length(counts)),
                       counts,
                       rep("sim", length(counts)),
                       rep(param_grid$rho[sim], length(counts)),
                       rep(param_grid$distance_decay[sim], length(counts))
  )
  colnames(temp_df) = columns
  distribution_df = rbind(distribution_df, temp_df)
}

distribution_df = distribution_df %>% filter(value > 0)
all_df <- rbind(distribution_df, real_df)



# bin the data for plotting
breaks <- seq(0.5, max(all_df$value)+0.5, by = 1)

all_df <- all_df %>%
  mutate(bin = cut(value, breaks = breaks, include.lowest = TRUE, right = FALSE),
         bin_mid = as.numeric(sub("\\[([0-9]+(?:\\.[0-9]+)?),.*", "\\1", bin)) + 0.5)


# loop over each of 10 conditions (rho and distance_decay)
plotlist = list()
legend = NULL
count = 1

for (rho_select in unique(param_grid$rho)){
  for (decay in unique(param_grid$distance_decay)){
    # calculate proportion in each bin
    binned_props <- all_df %>%
      filter(type == "real" | (rho == rho_select & distance_decay == decay)) %>%
      group_by(dataset_id, type, bin_mid) %>%
      summarize(n = n(), .groups = "drop") %>%
      group_by(dataset_id) %>%
      mutate(prop = n / sum(n)) %>%
      ungroup()
    
    # create confidence interval for simulated data (proportion in each bin based on multiple sims)
    sim_summary <- binned_props %>%
      filter(type == "sim") %>%
      group_by(bin_mid) %>%
      summarize(
        ymin = quantile(prop, 0.10),
        ymax = quantile(prop, 0.90),
        ymed = median(prop),
        .groups = "drop"
      )
    
    plot = ggplot() +
      # simulated confidence intervals
      geom_rect(
        data = sim_summary,
        aes(xmin = bin_mid - 0.5, xmax = bin_mid + 0.5, ymin = ymin, ymax = ymax),
        fill = "darkred",
        alpha = 0.3
      ) +
      
      # simulated median
      geom_segment(
        data = sim_summary,
        aes(x = bin_mid - 0.5, xend = bin_mid + 0.5, y = ymed, yend = ymed),
        color = "darkred",
        linewidth = 1
      ) +
      
      # real data lollipop-histogram
      geom_segment(
        data = binned_props %>% filter(type == "real"),
        aes(x = bin_mid - 0.5, xend = bin_mid + 0.5,
            y = prop, yend = prop, color = dataset_id),
        linewidth = 1.2
      ) +
      
      labs(x = "Number of Siblings", y = "Proportion",
           title = paste("Rho = ", rho_select, ", Decay Method = ", decay, sep = "")) +
      scale_color_brewer(palette = "Paired") +
      theme_minimal()
    
    legend = as_ggplot(get_legend(plot))
    plot = plot + theme(legend.position = "none")
    plotlist[[count]] = plot
    count = count + 1
  }
}

grid = grid.arrange(grobs = plotlist, ncol = 2)
sibdistributions = grid.arrange(grid, legend, ncol = 2, widths = c(10,1))
ggsave("figures/simfigs/methods_comparison_sibship_distribution.jpg", sibdistributions,
       height = 5000, width = 2500, units = "px")
