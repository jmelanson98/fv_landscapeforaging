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
       height = 5000, width = 2500, units = "px")




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
ggsave("figures/simfigs/methods_comparison_accuracy.jpg", fig, width = 2000, height = 1000, units = 'px')

