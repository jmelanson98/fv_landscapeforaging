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
source('simulate_data/src/GenotypeSimFunctions.R')
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
result = draw_simple_multi_landscape(sample_size = 1000,
                                   num_landscape = 6,
                                   landscape_size = 1500,
                                   trapgrid_size = 300,
                                   number_traps = 25,
                                   number_colonies = 10000,
                                   colony_sizes = rep(100,10000),
                                   rho = 100,
                                   distance_decay = "exponential")
yobs = result[[1]]
yobs_detected = yobs[rowSums(yobs)>0,]
colony_data = result[[2]]
colony_data_detected = colony_data[rowSums(yobs) > 0,]
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

# create a set of simulations to perform
paternity_probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)
impGenotypesList = list()
mixGenotypesList = list()

# simulate genotypes for each species and mating condition
for (i in 1:length(paternity_probs)){
  impGenotypesList[[i]] = simulateGenotypes(alleleFreqs = impatiens_alellefreq,
                                            colonyDataDetected = colony_data_detected,
                                            observationMatrix = yobs_detected,
                                            probMultiplePaternity = paternity_probs[i])
  
  mixGenotypesList[[i]] = simulateGenotypes(alleleFreqs = mixtus_alellefreq,
                                            colonyDataDetected = colony_data_detected,
                                            observationMatrix = yobs_detected,
                                            probMultiplePaternity = paternity_probs[i])
}

all_data = list(yobs, yobs_detected, colony_data, colony_data_detected, trap_data, impGenotypesList, mixGenotypesList)
saveRDS(all_data, "simulate_data/colony_assignments/test_effective_paternity/firstsim.RDS")


# write files to .txt for COLONY
write.table(impGenotypesList[[1]][,!colnames(impGenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony/impatiens_pp0.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(impGenotypesList[[2]][,!colnames(impGenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony/impatiens_pp0.2.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(impGenotypesList[[3]][,!colnames(impGenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony/impatiens_pp0.4.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(impGenotypesList[[4]][,!colnames(impGenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony/impatiens_pp0.6.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(impGenotypesList[[5]][,!colnames(impGenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony/impatiens_pp0.8.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(impGenotypesList[[6]][,!colnames(impGenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony/impatiens_pp1.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(mixGenotypesList[[1]][,!colnames(mixGenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony/mixtus_pp0.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(mixGenotypesList[[2]][,!colnames(mixGenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony/mixtus_pp0.2.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(mixGenotypesList[[3]][,!colnames(mixGenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony/mixtus_pp0.4.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(mixGenotypesList[[4]][,!colnames(mixGenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony/mixtus_pp0.6.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(mixGenotypesList[[5]][,!colnames(mixGenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony/mixtus_pp0.8.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(mixGenotypesList[[6]][,!colnames(mixGenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony/mixtus_pp1.txt", sep= ",", col.names = FALSE, row.names = FALSE)

# construct error rates files -- all zeros for now
imp_columns = unique(impatiens_alellefreq$MarkerID)
impatiens_error_rates = as.data.frame(matrix(0, nrow = 4, ncol = length(imp_columns)))
names(impatiens_error_rates) = imp_columns
write.table(impatiens_error_rates, "simulate_data/colony_assignments/test_effective_paternity/for_colony/impatiens_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

mix_columns = unique(mixtus_alellefreq$MarkerID)
mixtus_error_rates = as.data.frame(matrix(0, nrow = 4, ncol = length(mix_columns)))
names(mixtus_error_rates) = mix_columns
write.table(mixtus_error_rates, "simulate_data/colony_assignments/test_effective_paternity/for_colony/mixtus_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

# construct .DAT files for COLONY
# manually change these files to run colony with or without female monogamy -- it's faster than running through all the steps again
# no multiple paternity
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/impatiens_pp0.DAT", delim=",")
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/mixtus_pp0.DAT", delim=",")

 # multiple paternity 0.2
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/impatiens_pp0.2.DAT", delim=",")
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/mixtus_pp0.2.DAT", delim=",")

# multiple paternity 0.4
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/impatiens_pp0.4.DAT", delim=",")
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/mixtus_pp0.4.DAT", delim=",")

# multiple paternity 0.6
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/impatiens_pp0.6.DAT", delim=",")
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/mixtus_pp0.6.DAT", delim=",")

# multiple paternity 0.8
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/impatiens_pp0.8.DAT", delim=",")
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/mixtus_pp0.8.DAT", delim=",")

# multiple paternity 1
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/impatiens_pp1.DAT", delim=",")
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/mixtus_pp1.DAT", delim=",")
