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
library(matrixStats)


# First simulate some observations of bees at multiple landscapes
result = draw_simple_multi_landscape(sample_size = 2000,
                                   num_landscape = 6,
                                   landscape_size = 1500,
                                   trapgrid_size = 300,
                                   number_traps = 25,
                                   number_colonies = 20000,
                                   colony_sizes = rep(100,15000),
                                   rho = 100,
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


 
# load in allele frequency data
impatiens_alellefreq = read.csv("colony_assignments/Colony2_Linux/impatiens2023_final1.AlleleFreq")
mixtus_alellefreq = read.csv("colony_assignments/Colony2_Linux/mixtus2023_final1.AlleleFreq")
colony_data_detected = colony_data[rowSums(yobs) > 1,]

impatienssample_df = data.frame(individual = NA, truecolony = NA,
                                BT10_1 = NA, BT10_2 = NA,
                                B96_1 = NA, B96_2 = NA,
                                BTMS0059_1 = NA, BTMS0059_2 = NA,
                                BTMS0081_1 = NA, BTMS0081_2 = NA,
                                BTMS0062_1 = NA, BTMS0062_2 = NA,
                                B126_1 = NA, B126_2 = NA,
                                BTERN01_1 = NA, BTERN01_2 = NA,
                                B124_1 = NA, B124_2 = NA,
                                BTMS0057_1 = NA, BTMS0057_2 = NA,
                                BT30_1 = NA, BT30_2 = NA,
                                B10_1 = NA, B10_2 = NA,
                                BTMS0083_1 = NA, BTMS0083_2 = NA
                           )
count = 0

for (colony in 1:nrow(colony_data_detected)){
  
  # make a dataframe to hold genotypes for the focal colony
  numsibs = rowSums(yobs_detected)[colony]
  singlesibship_df = data.frame(individual = (count + 1):(count+numsibs), truecolony = rep(colony_data_detected$colonyid[colony], numsibs),
                                  BT10_1 = NA, BT10_2 = NA,
                                  B96_1 = NA, B96_2 = NA,
                                  BTMS0059_1 = NA, BTMS0059_2 = NA,
                                  BTMS0081_1 = NA, BTMS0081_2 = NA,
                                  BTMS0062_1 = NA, BTMS0062_2 = NA,
                                  B126_1 = NA, B126_2 = NA,
                                  BTERN01_1 = NA, BTERN01_2 = NA,
                                  B124_1 = NA, B124_2 = NA,
                                  BTMS0057_1 = NA, BTMS0057_2 = NA,
                                  BT30_1 = NA, BT30_2 = NA,
                                  B10_1 = NA, B10_2 = NA,
                                  BTMS0083_1 = NA, BTMS0083_2 = NA
  )
  count = max(singlesibship_df$individual)
  
  # for each allele, assign some values to the parents and then to each offspring
  allelecounter = 2
  for (allele in unique(impatiens_alellefreq$MarkerID)){
    # assign parental genotypes
    queen = sample(impatiens_alellefreq$AlleleID[impatiens_alellefreq$MarkerID == allele], 
                   size = 2, replace = TRUE, 
                   prob = impatiens_alellefreq$UpdatedFreq[impatiens_alellefreq$MarkerID == allele])
    male = sample(impatiens_alellefreq$AlleleID[impatiens_alellefreq$MarkerID == allele], 
                   size = 1, replace = TRUE, 
                   prob = impatiens_alellefreq$UpdatedFreq[impatiens_alellefreq$MarkerID == allele])
    
    # assign daughter genotypes
    singlesibship_df[,allelecounter + 1] = rep(male, numsibs)
    singlesibship_df[,allelecounter + 2] = sample(queen, size = numsibs, replace = TRUE, prob = c(0.5, 0.5))
    
    allelecounter = allelecounter + 2
  }
  impatienssample_df=rbind(impatienssample_df, singlesibship_df)
}
 