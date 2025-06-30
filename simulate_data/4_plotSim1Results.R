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
library(tibble)
library(ggpubr)

##### Set Environment #####
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
#setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging") # server
rstan_options(auto_write = TRUE) 
options(mc.cores = parallel::detectCores())
sim_path = "simulate_data/methods_comparison/"

##### Load in data #####
param_grid = readRDS(paste(sim_path, "param_grid.rds", sep = ""))

##### Check for errors in stan fit
# list all .err files
err_dir <- "simulate_data/logs"
err_files <- list.files(err_dir, pattern = "^cmdstan.*\\.err$", full.names = TRUE)
out_files <- list.files(err_dir, pattern = "^cmdstan.*\\.out$", full.names = TRUE)


# Function to check each file
check_err_file <- function(file) {
  lines <- readLines(file, warn = FALSE)

  list(
    file = basename(file),
    unfinished = any(grepl("DUE TO TIME LIMIT", lines, ignore.case = TRUE)),
    has_divergent = any(grepl("divergent", lines, ignore.case = TRUE)),
    has_low_ess = any(grepl("Effective Samples Size", lines, ignore.case = TRUE)),
    has_error = any(grepl("error", lines, ignore.case = TRUE)),
    has_undefined = any(grepl("The following variables have undefined values", lines, ignore.case = TRUE)),
    executed = any(grepl("Execution halted", lines, ignore.case = TRUE)),
    unknown_error = any(grepl("finished unexpectedly!", lines, ignore.case = TRUE))
  )
}

check_out_file <- function(file) {
  lines <- readLines(file, warn = FALSE)
  
  list(
    file = basename(file),
    has_divergent = any(grepl("divergent", lines, ignore.case = TRUE)),
    has_low_ess = any(grepl("Effective Sample Size", lines, ignore.case = TRUE)),
    has_error_out = any(grepl("error", lines, ignore.case = TRUE)),
    rejected_init = any(grepl("Rejecting initial value", lines, ignore.case = TRUE)),
    saved = any(grepl("Model saved", lines, ignore.case = TRUE)),
    finished = any(grepl("Model complete", lines, ignore.case = TRUE))
  )
}


# Apply to all files
err_results <- lapply(err_files, check_err_file)
out_results <- lapply(out_files, check_out_file)

# Convert to a data frame
error_summary <- do.call(rbind, lapply(err_results, as.data.frame))
output_summary <- do.call(rbind, lapply(out_results, as.data.frame))

# Join with param grid
param_grid$id = as.integer(rownames(param_grid))
error_summary$id = as.integer(sapply(strsplit(error_summary$file, "[_.]"), function(x) x[3]))
output_summary$id = as.integer(sapply(strsplit(output_summary$file, "[_.]"), function(x) x[3]))


error_summary = left_join(error_summary, param_grid, by = "id")
summary = left_join(error_summary, output_summary, by = "id")

# plot results
ggplot(summary, aes(x = distance_decay, fill = unknown_error)) +
  geom_bar() +
  theme_minimal()


# remove any sims that had errors
param_grid = error_summary[error_summary$has_error == FALSE,]







##### Plot simulated vs model-predicted colony average foraging distance
# e.g., figure 3 from Pope & Jha

#initate a data frame to store values for each simulation
columns = c("samplesize", "beta", "landscape_id", "colony_size_bin", 
            "true_colony_avg", "true_colony_sd", "model_colony_avg", "model_colony_sd")
allsim_colonies = data.frame(matrix(nrow = 2*nrow(param_grid), ncol = length(columns))) 
colnames(allsim_colonies) = columns
allsim_colonies$id = rep(param_grid$id,2)
allsim_colonies$colony_size_bin = c(rep(c("1-3"), nrow(param_grid)), rep(c("4+"), nrow(param_grid)))

for (sim in param_grid$id){
  # load in results for sim
  results_path = sprintf("simulate_data/batch_sim1/data/sim_result_%03d", sim)
  colony_data = colony_data = readRDS(paste(results_path, "/colonydata.RDS", sep =""))
  yobs = readRDS(paste(results_path, "/yobs.RDS", sep=""))
  stanFit = readRDS(paste(results_path, "/stanFitGQ.RDS", sep=""))
  
  # remove zero colonies
  print(class(yobs))
  print(dim(yobs))
  colony_data = colony_data[rowSums(yobs)>0,]
  yobs = yobs[rowSums(yobs) >0,]

  #extract model estimates
  colony_data$model_estimate = summary(stanFit, pars = c("colony_dist"))$summary[,1]
  
  #record observed colony sizes
  colony_data$observed_size = rowSums(yobs)
  
  # create categorical variable for colony size
 if(max(colony_data$observed_size > 3)){
 colony_data$size_bin = cut(colony_data$observed_size, breaks = c(0, 3, max(colony_data$observed_size)), labels = c("1-3", "4+"))
  } else {
  colony_data$size_bin = rep(c("1-3"), nrow(colony_data))
  }
  for (bin in c("1-3", "4+")){
    # add info to summary df
    allsim_colonies$samplesize[allsim_colonies$id == sim] = param_grid$sample_size[param_grid$id == sim]
    allsim_colonies$beta[allsim_colonies$id == sim] = param_grid$beta[param_grid$id == sim]
    allsim_colonies$landscape_id[allsim_colonies$id == sim] = param_grid$landscape_id[param_grid$id == sim]
    allsim_colonies$true_colony_avg[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = mean(colony_data$foraging_range[colony_data$size_bin == bin])
    allsim_colonies$true_colony_sd[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = sd(colony_data$foraging_range[colony_data$size_bin == bin])
    allsim_colonies$model_colony_avg[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = mean(colony_data$model_estimate[colony_data$size_bin == bin])
    allsim_colonies$model_colony_sd[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = sd(colony_data$model_estimate[colony_data$size_bin == bin])
    allsim_colonies$landscape_mean[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = summary(stanFit, pars = c("land_dist"))$summary[,1]
    allsim_colonies$landscape_sd[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = summary(stanFit, pars = c("land_dist"))$summary[,3]
    
    # save betas and thetas
    allsim_colonies$beta_mean[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = summary(stanFit, pars = c("beta"))$summary[,1]
    allsim_colonies$beta_msd[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = summary(stanFit, pars = c("beta"))$summary[,3]
    allsim_colonies$theta_mean[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = summary(stanFit, pars = c("theta"))$summary[,1]
    allsim_colonies$theta_sd[allsim_colonies$id == sim & allsim_colonies$colony_size_bin == bin] = summary(stanFit, pars = c("theta"))$summary[,3]
    
  }}

colonyplot = ggplot(data = allsim_colonies, 
       aes(x = true_colony_avg, y = model_colony_avg, color = colony_size_bin)) +
  geom_point() +
  geom_errorbar(aes(ymin=model_colony_avg-model_colony_sd, ymax=model_colony_avg+model_colony_sd)) +
  geom_errorbarh(aes(xmin = true_colony_avg - true_colony_sd, xmax = true_colony_avg + true_colony_sd)) +
  geom_abline(intercept = 0, slope =1) +
  labs(x = "True colony foraging distance", y = "Model estimated colony foraging distance", color = "Number of workers per colony") +
  theme_minimal()

ggsave("figures/simfigs/Popefig3GQ.jpg", colonyplot, width = 4000, height = 1000, units= "px")
write.csv(allsim_colonies, "simulate_data/batch_sim1/forplotting.csv")


##### Plot figure 3 from Pope & Jha
allsim_colonies = read.csv("simulate_data/batch_sim1/forplotting.csv")

allsim_summary = allsim_colonies %>% group_by(samplesize, beta, colony_size_bin) %>%
  summarize(true_avg = mean(true_colony_avg, na.rm = TRUE),
            true_sd = sqrt(sum(true_colony_sd^2, na.rm = TRUE)),
            model_avg = mean(model_colony_avg, na.rm = TRUE),
            model_sd = sqrt(sum(model_colony_sd^2, na.rm = TRUE))
            ) %>%
  ungroup()

fig31000 = ggplot(as.data.frame(allsim_summary[allsim_summary$samplesize == 1000,]), 
              aes(x = true_avg, y = model_avg,
                        color = colony_size_bin)) +
  geom_point() +
  geom_errorbar(aes(ymin = model_avg - 2*model_sd, ymax = model_avg + 2*model_sd)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "True Foraging Distance",
       y = "Modelled Foraging Distance",
       color = "Number of Siblings",
       title = "sample size = 1000 individuals") +
  theme_minimal()

fig3500 = ggplot(as.data.frame(allsim_summary[allsim_summary$samplesize == 500,]), 
                  aes(x = true_avg, y = model_avg,
                      color = colony_size_bin)) +
  geom_point() +
  geom_errorbar(aes(ymin = model_avg - 2*model_sd, ymax = model_avg + 2*model_sd)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "True Foraging Distance",
       y = "Modelled Foraging Distance",
       color = "Number of Siblings",
       title = "sample size = 500 individuals") +
  theme_minimal()

fig3250 = ggplot(as.data.frame(allsim_summary[allsim_summary$samplesize == 250,]), 
                 aes(x = true_avg, y = model_avg,
                     color = colony_size_bin)) +
  geom_point() +
  geom_errorbar(aes(ymin = model_avg - 2*model_sd, ymax = model_avg + 2*model_sd)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "True Foraging Distance",
       y = "Modelled Foraging Distance",
       color = "Number of Siblings",
       title = "sample size = 250 individuals") +
  theme_minimal()

grid = grid.arrange(fig31000, fig3500, fig3250, ncol = 1)
ggsave("figures/simfigs/Popefig3GQ.jpg", grid, width = 2000, height = 3000, units = 'px')


# Plot distributions of number of siblings for real and simulated data

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

mix = ggplot(mixsum, aes(x = n)) +
  geom_histogram() +
  labs(x = "Number of siblings",
       y = "Frequency",
       title = "B. mixtus") +
  theme_minimal()

imp = ggplot(impsum, aes(x = n)) +
  geom_histogram() +
  labs(x = "Number of siblings",
       y = "Frequency",
       title = "B. impatiens") +
  theme_minimal()

numberofsibs = grid.arrange(mix, imp, ncol = 2)
ggsave("figures/manuscript_figures/numberofsibs.jpg", numberofsibs,
       width = 2000, height = 500, units = "px")

mix_df = data.frame(dataset_id = rep("mix", length(mixsum$n)),
                    value = mixsum$n,
                    type = rep("real", length(mixsum$n)),
                    rhos = rep(NA, length(mixsum$n)),
                    colony_density = rep(NA, length(mixsum$n)))

imp_df = data.frame(dataset_id = rep("imp", length(impsum$n)),
                    value = impsum$n,
                    type = rep("real", length(impsum$n)),
                    rhos = rep(NA, length(impsum$n)),
                    colony_density = rep(NA, length(impsum$n)))

real_df = rbind(mix_df, imp_df)

# prep simulated data
columns = c("dataset_id", "value", "type", "rhos", "colony_density") 
sim_df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(sim_df) = columns

#reload full paramgrid
sim_path = "simulate_data/exponentiated_quadratic_sim/"
param_grid = readRDS(paste(sim_path, "param_grid.rds", sep = ""))

for (sim in rownames(param_grid)){
  sim = as.numeric(sim)
  results_path = sprintf("simulate_data/exponentiated_quadratic_sim/data/sim_result_%03d", sim)
  yobs = readRDS(paste(results_path, "/yobs.RDS", sep=""))
  
  counts = rowSums(yobs)
  temp_df = data.frame(rep(sim, length(counts)),
                       counts,
                       rep("sim", length(counts)),
                       rep(param_grid$rho[sim], length(counts)),
                       rep(param_grid$colony_density[sim], length(counts))
                       )
  colnames(temp_df) = columns
  sim_df = rbind(sim_df, temp_df)
}

sim_df = sim_df %>% filter(value > 0)
all_df <- rbind(sim_df, real_df)



# bin the data for plotting
breaks <- seq(0.5, max(all_df$value)+0.5, by = 1)

all_df <- all_df %>%
  mutate(bin = cut(value, breaks = breaks, include.lowest = TRUE, right = FALSE),
         bin_mid = as.numeric(sub("\\[([0-9]+(?:\\.[0-9]+)?),.*", "\\1", bin)) + 0.5)


# loop over each of 15 conditions (beta and sample size)
plotlist = list()
legend = NULL
count = 1

for (rho in unique(param_grid$rho)){
  for (colonydensity in unique(param_grid$colony_density)){
    # calculate proportion in each bin
    binned_props <- all_df %>%
      filter(type == "real" | (rhos == rho & colony_density == colonydensity)) %>%
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
    print(rho)
    print(paste("Median = ", sim_summary$ymed, sep = ""))
    
    plot = ggplot() +
      # Simulated confidence intervals
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
           title = paste("Rho = ", rho, ", Colony Density = ", colonydensity, sep = "")) +
      scale_color_brewer(palette = "Paired") +
      theme_minimal()
    
    legend = as_ggplot(get_legend(plot))
    plot = plot + theme(legend.position = "none")
    plotlist[[count]] = plot
    count = count + 1
  }
}

grid = grid.arrange(grobs = plotlist, ncol = 3)
sibdistributions = grid.arrange(grid, legend, ncol = 2, widths = c(10,1))
ggsave("figures/simfigs/sibship_size_distributions.jpg", sibdistributions,
       height = 4000, width = 4000, units = "px")
