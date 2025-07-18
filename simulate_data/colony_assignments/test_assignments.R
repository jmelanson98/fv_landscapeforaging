################################################################################
## Test COLONY performance on simulated multilocus genotypes
## Started by J Melanson
## July 17, 2025
################################################################################


# Load packages and set environment
rm(list = ls())
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")

# first, load in packages
source('simulate_data/src/GeneralizedSimFunctions.R')
source('colony_assignments/src/colony_assignment_functions.R')
library(dplyr)
library(tidyr)
library(stringr)
library(rcolony)
library(genepop)
library(data.table)
library(cowplot)
library(raster)
library(igraph)


# First simulate some observations of bees at multiple landscapes
result = draw_simple_multi_landscape(sample_size = 2000,
                                   num_landscape = 6,
                                   landscape_size = 1500,
                                   trapgrid_size = 300,
                                   number_traps = 25,
                                   number_colonies = 10000,
                                   colony_sizes = rep(100,10000),
                                   rho = 75,
                                   distance_decay = "exponential")
yobs = result[[1]]
yobs_detected = yobs[rowSums(yobs)>0,]
colony_data = result[[2]]
trap_data = result[[3]]


# Plot a line from each colony to each trap where it was observed

# Get row and column indices for non-zero entries
nonzero_idx <- which(yobs > 0, arr.ind = TRUE)

# Build the segment data frame
segment_df <- data.frame(
  colonyid = nonzero_idx[, "row"],
  trapid = nonzero_idx[, "col"],
  n_captures = yobs[nonzero_idx]
)

# add coordinates using the row/column indices
segment_df = segment_df %>%
  left_join(colony_data, by = "colonyid") %>%
  left_join(trap_data, by = "trapid")

ggplot() + 
  geom_point(data = colony_data, aes(x = colony_x, y = colony_y), color = "darkred", size = 0.5) +
  geom_point(data = trap_data, aes(x = trap_x, y = trap_y), color = "darkblue") +
  geom_segment(data = segment_df, aes(x = colony_x, y = colony_y, xend = trap_x, yend = trap_y),
               color = "darkgreen") +
  theme_minimal()


# plot distributions of number of obs per colony
# check worker distributions
# prep real data
mixsibs22 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mixsibs23 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
mixsibs = rbind(mixsibs22, mixsibs23)
impsibs22 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
impsibs23 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
impsibs = rbind(impsibs22, impsibs23)

mixsum = mixsibs %>%
  group_by(sibship_id) %>%
  summarize(n = n())
impsum = impsibs %>%
   group_by(sibship_id) %>%
   summarize(n = n())
 
 # simulated data
 simsum = data.frame(n = rowSums(yobs_detected))
 
 # combine
 mix = data.frame(n = mixsum$n)
 imp = data.frame(n = impsum$n)
 df_combined = bind_rows(
   simsum %>% mutate(group = "sim"),
   mix %>% mutate(group = "mix"),
   imp %>% mutate(group = "imp")
 )
 
 # bin manually
 df_props <- df_combined %>%
   group_by(group, n) %>%
   summarise(count = n(), .groups = "drop") %>%
   group_by(group) %>%
   mutate(proportion = count / sum(count))
 
 # plot
 ggplot(df_props, aes(x = n, y = proportion, color = group)) +
   geom_segment(aes(x = n - 0.5, xend = n + 0.5), linewidth = 0.6) +
   geom_point(size = 3) +
   labs(
     x = "Number of Sibs",
     y = "Proportion"
   ) +
   theme_minimal() +
   scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

