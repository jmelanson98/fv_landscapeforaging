##### Plot model comparison reults #####
##### plot number of workers "collected" per observed colony ######
##### plot foraging distance estimates #####
# Script Initiated: June 11, 2025
# By: Jenna Melanson


##### Load packages #####
library(ggplot2)
library(reshape2)
library(parallel)
library(future)
library(furrr)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggpubr)

#### Set working directory ####
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")

###############################################
### Check worker distribution plots
###############################################

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
sim_output = readRDS("simulate_data/methods_comparison/output.rds")

columns = c("dataset_id", "value", "type", "rho", "distance_decay") 
distribution_df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(distribution_df) = columns


for (sim in rownames(sim_output)){
  sim = as.numeric(sim)
  counts = unlist(sim_output$counts[sim])
  temp_df = data.frame(rep(sim, length(counts)),
                       counts,
                       rep("sim", length(counts)),
                       rep(sim_output$rho[sim], length(counts)),
                       rep(sim_output$distance_decay[sim], length(counts))
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

for (rho_select in unique(sim_output$rho)){
  for (decay in unique(sim_output$distance_decay)){
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
       height = 5000, width = 3000, units = "px")




###############################################
### Check error logs for stan model fitting
###############################################
# list all .err files
err_dir <- "simulate_data/logs"
err_files <- list.files(err_dir, pattern = "^stan.*\\.err$", full.names = TRUE)
out_files <- list.files(err_dir, pattern = "^stan.*\\.out$", full.names = TRUE)


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
    OOM = any(grepl("oom_kill", lines, ignore.case = TRUE))
  )
}

check_out_file <- function(file) {
  lines <- readLines(file, warn = FALSE)
  
  list(
    file = basename(file),
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
error_summary$task_id = as.integer(sapply(strsplit(error_summary$file, "[_.]"), function(x) x[3]))
output_summary$task_id = as.integer(sapply(strsplit(output_summary$file, "[_.]"), function(x) x[3]))


error_summary = left_join(error_summary, sim_output, by = "task_id")
summary = left_join(error_summary, output_summary, by = "task_id")

# plot results
ggplot(summary[summary$model_approach == "doubletons",], aes(x = true_average_foraging, y = model_average_foraging, color = has_low_ess)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  xlim(c(0,300)) +
  ylim(c(0,300)) +
  theme_minimal()

# bar plot of errors
ggplot(summary, aes(x = model_approach, fill = finished)) +
  geom_bar() +
  theme_minimal()

# plot distribution of OOM fails in singletons as a function of # colonies
summary$numcolonies = sapply(summary$counts, length)
ggplot(summary[summary$model_approach == "singletons",], aes(x = numcolonies, fill = OOM)) +
  geom_histogram() +
  theme_minimal()


###############################################
### Plot some figures of model accuracy
###############################################

accuracy = ggplot(sim_output, aes(x = true_average_foraging, y = model_average_foraging, color = model_approach)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  #xlim(c(0,300)) +
  #ylim(c(0,300)) +
  theme_minimal()

rhovsmu = ggplot(sim_output, aes(x = model_rho, y = model_mu, color = model_approach)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  #xlim(c(0,300)) +
  #ylim(c(0,300)) +
  theme_minimal()

rhoreturn = ggplot(sim_output, aes(x = rho, y = model_rho, color = model_approach)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  #xlim(c(0,300)) +
  #ylim(c(0,300)) +
  theme_minimal()

summary = sim_output %>% group_by(rho, distance_decay, model_approach) %>%
  summarize(true_avg = mean(true_average_foraging, na.rm = TRUE),
            true_sd = sqrt(sum(true_sd_foraging^2, na.rm = TRUE)),
            model_avg = mean(model_average_foraging, na.rm = TRUE),
            model_sd = sqrt(sum(model_sd_foraging^2, na.rm = TRUE))
  ) %>%
  ungroup()

summary$technique = paste(summary$distance_decay, summary$model_approach, sep = " / ")


#### Accuracy plot

# add a jitter for visibility (x direction)
summary$x_jittered <- summary$true_avg + runif(nrow(summary), min = -2, max = 2)

fig = ggplot(summary, aes(x = x_jittered, y = model_avg,
                          color = technique)) +
  geom_point() +
  geom_errorbar(aes(ymin = model_avg - 2*model_sd, ymax = model_avg + 2*model_sd)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "True Foraging Distance",
       y = "Modelled Foraging Distance",
       color = "Approach") +
  theme_minimal()
ggsave("figures/simfigs/methods_comparison_accuracy.jpg", fig, width = 2000, height = 1000, units = 'px')

