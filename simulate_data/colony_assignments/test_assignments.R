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
library(purrr)
library(stringr)
library(rcolony)
library(genepop)
library(data.table)
library(cowplot)
library(raster)
library(igraph)
library(matrixStats)
library(ggplot2)

################################################################################
## Test different rates of multiple paternity in mixtus and impatiens
################################################################################

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


# Check worker distributions
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


 
# Load in allele frequency data
impatiens_alellefreq = read.csv("colony_assignments/Colony2_Linux/impatiens2023_final1.AlleleFreq")
mixtus_alellefreq = read.csv("colony_assignments/Colony2_Linux/mixtus2023_final1.AlleleFreq")

# Create a set of simulations to perform
paternity_probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)
impGenotypesList = list()
mixGenotypesList = list()

# Simulate genotypes for each species and mating condition
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
# all_data = readRDS("simulate_data/colony_assignments/test_effective_paternity/firstsim.RDS")
# yobs = all_data[[1]]
# yobs_detected = all_data[[2]]
# colony_data = all_data[[3]]
# colony_data_detected = all_data[[4]]
# trap_data = all_data[[5]]
# impGenotypesList = all_data[[6]]
# mixGenotypesList = all_data[[7]]

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

# Write files to .csv to save true colony IDs
write.csv(impGenotypesList[[1]], 
            "simulate_data/colony_assignments/test_effective_paternity/true_data/impatiens_pp0.csv")
write.csv(impGenotypesList[[2]], 
            "simulate_data/colony_assignments/test_effective_paternity/true_data/impatiens_pp0.2.csv")
write.csv(impGenotypesList[[3]], 
            "simulate_data/colony_assignments/test_effective_paternity/true_data/impatiens_pp0.4.csv")
write.csv(impGenotypesList[[4]], 
            "simulate_data/colony_assignments/test_effective_paternity/true_data/impatiens_pp0.6.csv")
write.csv(impGenotypesList[[5]], 
            "simulate_data/colony_assignments/test_effective_paternity/true_data/impatiens_pp0.8.csv")
write.csv(impGenotypesList[[6]], 
            "simulate_data/colony_assignments/test_effective_paternity/true_data/impatiens_pp1.csv")
write.csv(mixGenotypesList[[1]], 
            "simulate_data/colony_assignments/test_effective_paternity/true_data/mixtus_pp0.csv")
write.csv(mixGenotypesList[[2]], 
            "simulate_data/colony_assignments/test_effective_paternity/true_data/mixtus_pp0.2.csv")
write.csv(mixGenotypesList[[3]], 
            "simulate_data/colony_assignments/test_effective_paternity/true_data/mixtus_pp0.4.csv")
write.csv(mixGenotypesList[[4]], 
            "simulate_data/colony_assignments/test_effective_paternity/true_data/mixtus_pp0.6.csv")
write.csv(mixGenotypesList[[5]], 
            "simulate_data/colony_assignments/test_effective_paternity/true_data/mixtus_pp0.8.csv")
write.csv(mixGenotypesList[[6]], 
            "simulate_data/colony_assignments/test_effective_paternity/true_data/mixtus_pp1.csv")



# Construct error rates files -- all zeros for now
imp_columns = unique(impatiens_alellefreq$MarkerID)
impatiens_error_rates = as.data.frame(matrix(0, nrow = 4, ncol = length(imp_columns)))
impatiens_error_rates[1,] = imp_columns
write.table(impatiens_error_rates, "simulate_data/colony_assignments/test_effective_paternity/for_colony/impatiens_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

mix_columns = unique(mixtus_alellefreq$MarkerID)
mixtus_error_rates = as.data.frame(matrix(0, nrow = 4, ncol = length(mix_columns)))
mixtus_error_rates[1,] = mix_columns
write.table(mixtus_error_rates, "simulate_data/colony_assignments/test_effective_paternity/for_colony/mixtus_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Construct .DAT files for COLONY
# manually change these files to run colony with or without female monogamy -- it's faster than running through all the steps again
# no multiple paternity
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/impatiens_pp0.DAT", delim=",")
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/mixtus_pp0.DAT", delim=",")

# multiple paternity 0.2
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/impatiens_pp0.2.DAT", delim=",")
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/mixtus_pp0.2.DAT", delim=",")

# multiple paternity 0.4
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/impatiens_pp0.4.DAT", delim=",")
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/mixtus_pp0.4.DAT", delim=",")

# multiple paternity 0.6
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/impatiens_pp0.6.DAT", delim=",")
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/mixtus_pp0.6.DAT", delim=",")

# multiple paternity 0.8
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/impatiens_pp0.8.DAT", delim=",")
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/mixtus_pp0.8.DAT", delim=",")

# multiple paternity 1
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/impatiens_pp1.DAT", delim=",")
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/mixtus_pp1.DAT", delim=",")



# Load in results from COLONY
files = c("impatiens_pp0", "impatiens_pp0.2", "impatiens_pp0.4", "impatiens_pp0.6", "impatiens_pp0.8", "impatiens_pp1",
         "impatiens_pp0_poly", "impatiens_pp0.2_poly", "impatiens_pp0.4_poly", "impatiens_pp0.6_poly", "impatiens_pp0.8_poly", 
         "impatiens_pp1_poly", "mixtus_pp0", "mixtus_pp0.2", "mixtus_pp0.4", "mixtus_pp0.6", "mixtus_pp0.8", "mixtus_pp1",
         "mixtus_pp0_poly", "mixtus_pp0.2_poly", "mixtus_pp0.4_poly", "mixtus_pp0.6_poly", "mixtus_pp0.8_poly", 
         "mixtus_pp1_poly")

errors = data.frame(test_condition = files,
                  numFP = NA,
                  numFN = NA,
                  numTP = NA,
                  total_real = NA,
                  FPR = NA,
                  FNR = NA)
family_plots = list()
for (i in 1:length(files)){
    filename = paste0("simulate_data/colony_assignments/test_effective_paternity/colony_output_withprior/", files[i], ".BestCluster")
    colony_output = as.data.frame(do.call(rbind, strsplit(trimws(readLines(filename)), "\\s+")[-1]))
    colnames(colony_output) = unlist(strsplit(readLines(filename), "\\s+")[1])
    
    genotypesim = paste(unlist(strsplit(files[i], "_"))[1:2], collapse = "_")
    true_data = read.csv(paste0("simulate_data/colony_assignments/test_effective_paternity/true_data/", genotypesim, ".csv"))
    
    # set probability threshold
    prob_thresh = 0.95
    
    # filter colony outputs
    colony_output = colony_output %>% filter(Probability >= prob_thresh)
    
    # make edge lists
    true_edges = true_data %>%
      group_by(truecolony) %>%
      filter(n() > 1) %>%
      summarise(pairs = combn(individual, 2, simplify = FALSE), .groups = "drop") %>%
      mutate(from = map_chr(pairs, 1),
             to = map_chr(pairs, 2)) %>%
      dplyr::select(from, to)
    
    inferred_edges = colony_output %>%
      group_by(ClusterIndex) %>%
      filter(n() > 1) %>%
      summarise(pairs = combn(OffspringID, 2, simplify = FALSE), .groups = "drop") %>%
      mutate(from = map_chr(pairs, 1),
             to = map_chr(pairs, 2)) %>%
      dplyr::select(from, to)
    
    # combine and classify edges by type (FN, FP, TP)
    # get true positives
    tp_edges = inner_join(true_edges, inferred_edges, by = c("from", "to"))
    tp_edges$type = "TP"
    
    # initialize other types as FN and FP
    true_edges$type = "FN"
    inferred_edges$type = "FP"
    
    # remove true positives from FN and FP dataframes
    fn_edges = anti_join(true_edges, tp_edges, by = c("from", "to"))
    fp_edges = anti_join(inferred_edges, tp_edges, by = c("from", "to"))
    
    # combine all edges
    all_edges = rbind(tp_edges, fn_edges, fp_edges)
    
    # make an igraph object
    graph = graph_from_data_frame(all_edges, directed = FALSE)
    E(graph)$color = recode(E(graph)$type, TP = "black", FP = "red", FN = "purple")
    families = plot(graph, 
         edge.color = E(graph)$color,
         vertex.label = NA,
         vertex.size = 2,
         main = files[i]
    )
    family_plots[[i]] = families
    
    # record FPR and FNR
    errors$FPR[errors$test_condition == files[i]] = nrow(fp_edges) / (nrow(tp_edges) + nrow(fn_edges))
    errors$FNR[errors$test_condition == files[i]] = nrow(fn_edges) / (nrow(tp_edges) + nrow(fn_edges))
    errors$total_real[errors$test_condition == files[i]] = nrow(tp_edges) + nrow(fn_edges)
    errors$numFP[errors$test_condition == files[i]] = nrow(fp_edges)
    errors$numFN[errors$test_condition == files[i]] = nrow(fn_edges)
    errors$numTP[errors$test_condition == files[i]] = nrow(tp_edges)
    
    
    }


# make a results table without poly
errors_subset = errors %>%
  filter(!str_detect(test_condition, "poly"))

ggplot(errors) +
  geom_point(aes(x = test_condition, y = numFP, color = "Number FP")) +
  geom_point(aes(x = test_condition, y = numFN, color = "Number FN")) +
  geom_point(aes(x = test_condition, y = total_real, color = "Total true")) +
  xlab("Simulation and COLONY Conditions") +
  ylab("Number of inferred or true relationships") +
  labs(title = "Sibship inclusion: P = 0.95") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(errors_subset, aes(x = test_condition, y = FPR)) +
  geom_point() +
  xlab("Simulation and COLONY Conditions") +
  ylab(expression(FPR == frac(FP, TP + FN))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(errors_subset, aes(x = test_condition, y = FNR)) +
  geom_point() +
  xlab("Simulation and COLONY Conditions") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))


################################################################################
## Test multiple paternity for higher resolution genetic dataset
################################################################################
# Siblingship size prior greatly improves FPR rate when we run COLONY assuming female monogamy, but 
# if we assume female polygamy we get a whole mess of false sibships
# Check: could we improve this tendency by simulating a fake data set with twice as many loci?

# Load in allele frequency data
impatiens_alellefreq = read.csv("colony_assignments/Colony2_Linux/impatiens2023_final1.AlleleFreq")
mixtus_alellefreq = read.csv("colony_assignments/Colony2_Linux/mixtus2023_final1.AlleleFreq")

# combine frequency data!
impatiens_alellefreq$MarkerID = paste0("impatiens_", impatiens_alellefreq$MarkerID)
mixtus_alellefreq$MarkerID = paste0("mixtus_", mixtus_alellefreq$MarkerID)
combined_allelefreq = rbind(impatiens_alellefreq, mixtus_alellefreq)

# Create a set of simulations to perform
paternity_probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)
GenotypesList = list()

# Simulate genotypes for each species and mating condition
for (i in 1:length(paternity_probs)){
  GenotypesList[[i]] = simulateGenotypes(alleleFreqs = combined_allelefreq,
                                            colonyDataDetected = colony_data_detected,
                                            observationMatrix = yobs_detected,
                                            probMultiplePaternity = paternity_probs[i])
}

saveRDS(GenotypesList, "simulate_data/colony_assignments/test_effective_paternity/augmentedsim.RDS")


# write files to .txt for COLONY
write.table(GenotypesList[[1]][,!colnames(GenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/augmented_pp0.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(GenotypesList[[2]][,!colnames(GenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/augmented_pp0.2.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(GenotypesList[[3]][,!colnames(GenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/augmented_pp0.4.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(GenotypesList[[4]][,!colnames(GenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/augmented_pp0.6.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(GenotypesList[[5]][,!colnames(GenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/augmented_pp0.8.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(GenotypesList[[6]][,!colnames(GenotypesList[[1]]) %in% c("truecolony")], 
            "simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/augmented_pp1.txt", sep= ",", col.names = FALSE, row.names = FALSE)

# Write files to .csv to save true colony IDs
write.csv(GenotypesList[[1]], 
          "simulate_data/colony_assignments/test_effective_paternity/true_data/augmented_pp0.csv")
write.csv(GenotypesList[[2]], 
          "simulate_data/colony_assignments/test_effective_paternity/true_data/augmented_pp0.2.csv")
write.csv(GenotypesList[[3]], 
          "simulate_data/colony_assignments/test_effective_paternity/true_data/augmented_pp0.4.csv")
write.csv(GenotypesList[[4]], 
          "simulate_data/colony_assignments/test_effective_paternity/true_data/augmented_pp0.6.csv")
write.csv(GenotypesList[[5]], 
          "simulate_data/colony_assignments/test_effective_paternity/true_data/augmented_pp0.8.csv")
write.csv(GenotypesList[[6]], 
          "simulate_data/colony_assignments/test_effective_paternity/true_data/augmented_pp1.csv")


# Construct error rates files -- all zeros for now
columns = unique(combined_allelefreq$MarkerID)
all_error_rates = as.data.frame(matrix(0, nrow = 4, ncol = length(columns)))
all_error_rates[1,] = columns
write.table(all_error_rates, "simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/augmented_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Construct .DAT files for COLONY
# manually change these files to run colony with or without female monogamy -- it's faster than running through all the steps again
# no multiple paternity
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                                name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/augmented_pp0.DAT", delim=",")

# multiple paternity 0.2
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                                name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/augmented_pp0.2.DAT", delim=",")

# multiple paternity 0.4
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                                name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/augmented_pp0.4.DAT", delim=",")

# multiple paternity 0.6
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                                name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/augmented_pp0.6.DAT", delim=",")

# multiple paternity 0.8
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                                name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/augmented_pp0.8.DAT", delim=",")

# multiple paternity 1
rcolony::build.colony.automatic(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux", 
                                name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux/augmented_pp1.DAT", delim=",")

# Load in results from COLONY
files = c("augmented_pp0", "augmented_pp0.2", "augmented_pp0.4", "augmented_pp0.6", "augmented_pp0.8", "augmented_pp1",
          "augmented_pp0_poly", "augmented_pp0.2_poly", "augmented_pp0.4_poly", "augmented_pp0.6_poly", "augmented_pp0.8_poly", "augmented_pp1_poly")

errors = data.frame(test_condition = files,
                    numFP = NA,
                    numFN = NA,
                    numTP = NA,
                    total_real = NA,
                    FPR = NA,
                    FNR = NA)
family_plots = list()
for (i in 1:length(files)){
  filename = paste0("simulate_data/colony_assignments/test_effective_paternity/colony_output_augmented/", files[i], ".BestCluster")
  colony_output = as.data.frame(do.call(rbind, strsplit(trimws(readLines(filename)), "\\s+")[-1]))
  colnames(colony_output) = unlist(strsplit(readLines(filename), "\\s+")[1])
  
  genotypesim = paste(unlist(strsplit(files[i], "_"))[1:2], collapse = "_")
  true_data = read.csv(paste0("simulate_data/colony_assignments/test_effective_paternity/true_data/", genotypesim, ".csv"))
  
  # set probability threshold
  prob_thresh = 0.95
  
  # filter colony outputs
  colony_output = colony_output %>% filter(Probability >= prob_thresh)
  
  # make edge lists
  true_edges = true_data %>%
    group_by(truecolony) %>%
    filter(n() > 1) %>%
    summarise(pairs = combn(individual, 2, simplify = FALSE), .groups = "drop") %>%
    mutate(from = map_chr(pairs, 1),
           to = map_chr(pairs, 2)) %>%
    dplyr::select(from, to)
  
  inferred_edges = colony_output %>%
    group_by(ClusterIndex) %>%
    filter(n() > 1) %>%
    summarise(pairs = combn(OffspringID, 2, simplify = FALSE), .groups = "drop") %>%
    mutate(from = map_chr(pairs, 1),
           to = map_chr(pairs, 2)) %>%
    dplyr::select(from, to)
  
  # combine and classify edges by type (FN, FP, TP)
  # get true positives
  tp_edges = inner_join(true_edges, inferred_edges, by = c("from", "to"))
  tp_edges$type = "TP"
  
  # initialize other types as FN and FP
  true_edges$type = "FN"
  inferred_edges$type = "FP"
  
  # remove true positives from FN and FP dataframes
  fn_edges = anti_join(true_edges, tp_edges, by = c("from", "to"))
  fp_edges = anti_join(inferred_edges, tp_edges, by = c("from", "to"))
  
  # combine all edges
  all_edges = rbind(tp_edges, fn_edges, fp_edges)
  
  # make an igraph object
  graph = graph_from_data_frame(all_edges, directed = FALSE)
  E(graph)$color = recode(E(graph)$type, TP = "black", FP = "red", FN = "purple")
  families = plot(graph, 
                  edge.color = E(graph)$color,
                  vertex.label = NA,
                  vertex.size = 2,
                  main = files[i]
  )
  family_plots[[i]] = families
  
  # record FPR and FNR
  errors$FPR[errors$test_condition == files[i]] = nrow(fp_edges) / (nrow(tp_edges) + nrow(fn_edges))
  errors$FNR[errors$test_condition == files[i]] = nrow(fn_edges) / (nrow(tp_edges) + nrow(fn_edges))
  errors$total_real[errors$test_condition == files[i]] = nrow(tp_edges) + nrow(fn_edges)
  errors$numFP[errors$test_condition == files[i]] = nrow(fp_edges)
  errors$numFN[errors$test_condition == files[i]] = nrow(fn_edges)
  errors$numTP[errors$test_condition == files[i]] = nrow(tp_edges)
  
  
}


# make a results table without poly
errors_subset = errors %>%
  filter(!str_detect(test_condition, "poly"))

ggplot(errors) +
  geom_point(aes(x = test_condition, y = numFP, color = "Number FP")) +
  geom_point(aes(x = test_condition, y = numFN, color = "Number FN")) +
  geom_point(aes(x = test_condition, y = total_real, color = "Total true")) +
  xlab("Simulation and COLONY Conditions") +
  ylab("Number of inferred or true relationships") +
  labs(title = "Sibship inclusion: P = 0.95") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(errors_subset, aes(x = test_condition, y = FPR)) +
  geom_point() +
  xlab("Simulation and COLONY Conditions") +
  ylab(expression(FPR == frac(FP, TP + FN))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(errors_subset, aes(x = test_condition, y = FNR)) +
  geom_point() +
  xlab("Simulation and COLONY Conditions") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))


################################################################################
## Test accuracy across datasets of different sizes, with realistic error rates
################################################################################

# read in real error rate files
mixtus_error_rates = data.frame(#c("BT10", 0, 0, 0.01),
                                #c("BTMS0104", 0, 0, 0.01),
                                #c("BTMS0057", 0, 0, 0.01),
                                c("BTMS0086", 0, 0, 0.015),
                                c("BTMS0066", 0, 0, 0.01),
                                c("BTMS0062", 0, 0, 0.015),
                                c("BTMS0136", 0, 0, 0.01),
                                #c("BTERN01", 0, 0, 0.0215),
                                c("BTMS0126", 0, 0, 0.0162),
                                #c("BTMS0059", 0, 0, 0.01),
                                c("BL13", 0, 0, 0.0159),
                                c("BTMS0083", 0, 0, 0.01),
                                c("B126", 0, 0, 0.01))
write.table(mixtus_error_rates, "simulate_data/colony_assignments/test_sample_size/for_colony/mixtus_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
mixtus_errors = read.table("simulate_data/colony_assignments/test_sample_size/for_colony/mixtus_error_rates.txt", sep = ",")
colnames(mixtus_errors) = mixtus_errors[1,]
mixtus_errors = mixtus_errors[-1,]

impatiens_error_rates = data.frame(c("BT10", 0, 0, 0.017),
                                   c("B96", 0, 0, 0.027),
                                   c("BTMS0059", 0, 0, 0.017),
                                   c("BTMS0081", 0, 0, 0.016),
                                   c("BL13", 0, 0, 0.011),
                                   c("BTMS0062", 0, 0, 0.022),
                                   #c("B126", 0, 0, 0.01),
                                   c("BTERN01", 0, 0, 0.028),
                                   c("B124", 0, 0, 0.011),
                                   #c("BTMS0057", 0, 0, 0.017),
                                   c("BT30", 0, 0, 0.01),
                                   c("B10", 0, 0, 0.029),
                                   c("BTMS0083", 0, 0, 0.01)
                                   #c("BTMS0073", 0, 0, 0.01),
                                   #c("BT28", 0, 0, 0.01)
)
write.table(impatiens_error_rates, "simulate_data/colony_assignments/test_sample_size/for_colony/impatiens_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
impatiens_errors = read.table("simulate_data/colony_assignments/test_sample_size/for_colony/impatiens_error_rates.txt", sep = ",")
colnames(impatiens_errors) = impatiens_errors[1,]
impatiens_errors = impatiens_errors[-1,]


# get allele names
mixtus_alleles = as.vector(as.matrix(mixtus_error_rates[1, ]))
impatiens_alleles = as.vector(as.matrix(impatiens_error_rates[1,]))

# load in allele frequency files
impatiens_alellefreq = read.csv("colony_assignments/Colony2/initial/impatiens2023.AlleleFreq")
mixtus_alellefreq = read.csv("colony_assignments/Colony2/initial/mixtus2023.AlleleFreq")

# remove loci not being used
impatiens_alellefreq = impatiens_alellefreq %>% filter(MarkerID %in% impatiens_alleles)
mixtus_alellefreq = mixtus_alellefreq %>% filter(MarkerID %in% mixtus_alleles)


# load in real genotype files and get missing rates
mixtus_scores = read.csv("data/merged_by_year/csvs/mixtus_2023_scores_finalloci.csv")[,-1]
impatiens_scores = read.csv("data/merged_by_year/csvs/impatiens_2023_scores_finalloci.csv")[,-1]

mixtus_missing =  setNames(data.frame(matrix(ncol = length(mixtus_alleles), nrow = 0)), mixtus_alleles)
for (i in 1:length(mixtus_alleles)){
  allele = mixtus_alleles[i]
  mixtus_missing[1,i] = sum(colSums(mixtus_scores[,str_detect(colnames(mixtus_scores), allele)] ==0)) / (2*nrow(mixtus_scores))
}

impatiens_missing =  setNames(data.frame(matrix(ncol = length(impatiens_alleles), nrow = 0)), impatiens_alleles)
for (i in 1:length(impatiens_alleles)){
  allele = impatiens_alleles[i]
  impatiens_missing[1,i] = sum(colSums(impatiens_scores[,str_detect(colnames(impatiens_scores), allele)] ==0)) / (2*nrow(impatiens_scores))
}


# Generate three sibship data sets (e.g., 3 x yobs) of 2000 individuals each
numsims = 5
for (i in 1:numsims){
  result = draw_simple_multi_landscape(sample_size = 2000,
                                       num_landscape = 6,
                                       landscape_size = 1500,
                                       trapgrid_size = 300,
                                       number_traps = 25,
                                       number_colonies = 20000,
                                       colony_sizes = rep(100,20000),
                                       rho = 100,
                                       distance_decay = "exponential")
  sibship_data = list(result[[1]], result[[1]][rowSums(result[[1]])>0,], result[[2]], result[[2]][rowSums(result[[1]])>0,], result[[3]])
  saveRDS(sibship_data, paste0("simulate_data/colony_assignments/test_sample_size/sim", i, ".RDS"))
}
#sibdata = readRDS(paste0("simulate_data/colony_assignments/test_sample_size/sim1.RDS"))

# For each data set, simulate genotypes assuming single paternity (mixtus and impatiens)
impGenotypesList = list()
mixGenotypesList = list()

# Simulate genotypes for each species and mating condition
for (i in 1:numsims){
  sibship_data = readRDS(paste0("simulate_data/colony_assignments/test_sample_size/sim", i, ".RDS"))
  # simulate imp genotypes
  tempImp = simulateGenotypes(alleleFreqs = impatiens_alellefreq,
                           colonyDataDetected = sibship_data[[4]],
                           observationMatrix = sibship_data[[2]],
                           trapData = sibship_data[[5]],
                           probMultiplePaternity = 0)
  print("Imp sim done")
  # induce errors and missingness
  impGenotypesList[[i]] = induceErrors(genotypeDF = tempImp,
                                       errorRates = impatiens_errors,
                                       missingRates = impatiens_missing,
                                       alleleFreqs = impatiens_alellefreq)
  print("imp errors done")
  tempMix = simulateGenotypes(alleleFreqs = mixtus_alellefreq,
                              colonyDataDetected = sibship_data[[4]],
                              observationMatrix = sibship_data[[2]],
                              trapData = sibship_data[[5]],
                              probMultiplePaternity = 0)
  print("mix sim done")
  mixGenotypesList[[i]] = induceErrors(genotypeDF = tempMix,
                                       errorRates = mixtus_errors,
                                       missingRates = mixtus_missing,
                                       alleleFreqs = mixtus_alellefreq)
  print("mix errors done")
}

sibship_genotypes = list(mixGenotypesList, impGenotypesList)
saveRDS(sibship_genotypes, "simulate_data/colony_assignments/test_sample_size/sibship_genotypes.RDS")
# sibship_genotypes = readRDS("simulate_data/colony_assignments/test_sample_size/sibship_genotypes.RDS")
# mixGenotypesList = sibship_genotypes[[1]]
# impGenotypesList = sibship_genotypes[[2]]

# Now, subset each dataset to contain 100%, 80%, 60%, 40%, 20%, 10% of individuals
subsets = c(1, 0.8, 0.6, 0.4, 0.2, 0.1)
for (i in 1:length(mixGenotypesList)){
    genotypes_mix = mixGenotypesList[[i]]
    genotypes_imp = impGenotypesList[[i]]
    
  for (j in 1:length(subsets)){
    rows_keep = nrow(genotypes_mix)*subsets[j]
    
    # make subset for mixtus
    genotypes_subset_mix = genotypes_mix[sample(nrow(genotypes_mix), rows_keep), ]
    
    # make sibship exclusion table
    subset_exclusion_mix = createExclusionTable(genotypes_subset_mix)

    # write txt file for colony
    write.table(genotypes_subset_mix[,!colnames(genotypes_subset_mix) %in% c("truecolony", "landscape_id", "trap_id")], 
                paste0("simulate_data/colony_assignments/test_sample_size/for_colony/mixtus_set", i, "_sub", subsets[j], ".txt"),
                sep= ",", col.names = FALSE, row.names = FALSE)
    # save full csv
    write.csv(genotypes_subset_mix, 
              paste0("simulate_data/colony_assignments/test_sample_size/true_data/mixtus_set", i, "_sub", subsets[j], ".csv"))
    # write sibship exclusion table
    write.table(
      subset_exclusion_mix,
      file = paste0("simulate_data/colony_assignments/test_sample_size/for_colony/mixtus_exclusion", i, "_sub", subsets[j], ".txt"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      na = ""
    )
    
    # make subset for impatiens
    genotypes_subset_imp = genotypes_imp[sample(nrow(genotypes_imp), rows_keep), ]
    
    # make exclusion table
    subset_exclusion_imp = createExclusionTable(genotypes_subset_imp)
    
    # write txt file for colony
    write.table(genotypes_subset_imp[,!colnames(genotypes_subset_imp) %in% c("truecolony", "landscape_id", "trap_id")], 
                paste0("simulate_data/colony_assignments/test_sample_size/for_colony/impatiens_set", i, "_sub", subsets[j], ".txt"),
                sep= ",", col.names = FALSE, row.names = FALSE)
    # save full csv
    write.csv(genotypes_subset_imp, 
              paste0("simulate_data/colony_assignments/test_sample_size/true_data/impatiens_set", i, "_sub", subsets[j], ".csv"))
    # write exclusion table
    write.table(
      subset_exclusion_imp,
      file = paste0("simulate_data/colony_assignments/test_sample_size/for_colony/impatiens_exclusion", i, "_sub", subsets[j], ".txt"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      na = ""
    )
  }
}


# Construct .DAT files for COLONY
mixtus_errors_filepath = "simulate_data/colony_assignments/test_sample_size/for_colony/mixtus_error_rates.txt"
impatiens_errors_filepath = "simulate_data/colony_assignments/test_sample_size/for_colony/impatiens_error_rates.txt"

for (i in 1:nsim){
  for (j in 1:length(subsets)){
    for (k in c("exclusion", "no_exclusion")){
    # get sample size, working directory
    size = nrow(mixGenotypesList[[i]]) * subsets[j]
    workingdir = "/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux"
    
    # get genotype filepaths
    mixtus_genotypes_filepath = paste0("simulate_data/colony_assignments/test_sample_size/for_colony/mixtus_set", i, "_sub", subsets[j], ".txt")
    impatiens_genotypes_filepath = paste0("simulate_data/colony_assignments/test_sample_size/for_colony/impatiens_set", i, "_sub", subsets[j], ".txt")
    
    # get exclusion paths
    if (k == "exclusion"){
      mixtus_exclusion_filepath = paste0("simulate_data/colony_assignments/test_sample_size/for_colony/mixtus_exclusion", i, "_sub", subsets[j], ".txt")
      impatiens_exclusion_filepath = paste0("simulate_data/colony_assignments/test_sample_size/for_colony/impatiens_exclusion", i, "_sub", subsets[j], ".txt")
    } else {
      mixtus_exclusion_filepath = NULL
      impatiens_exclusion_filepath = NULL
    }
    
    #build .DAT for mixtus
    rcolony::build.colony.superauto(wd=workingdir, 
                                    name=paste0("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/test_sample_size/colony_output/mixtus_set", i, "_sub", subsets[j], "_", k, ".DAT"), 
                                    datasetname = paste0("mixtus_set", i, "_sub", subsets[j], "_", k),
                                    delim=",",
                                    sample_size = size,
                                    num_loci = 8,
                                    error_rates_path = mixtus_errors_filepath,
                                    genotypes_path = mixtus_genotypes_filepath,
                                    exclusion_path = mixtus_exclusion_filepath
                                    )
    
    # build .DAT for impatiens
    rcolony::build.colony.superauto(wd=workingdir, 
                                    name=paste0("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/test_sample_size/colony_output/impatiens_set", i, "_sub", subsets[j], "_", k, ".DAT"), 
                                    datasetname = paste0("impatiens_set", i, "_sub", subsets[j], "_", k),
                                    delim=",",
                                    sample_size = size,
                                    num_loci = 11,
                                    error_rates_path = impatiens_errors_filepath,
                                    genotypes_path = impatiens_genotypes_filepath,
                                    exclusion_path = impatiens_exclusion_filepath
    )
    }
  }
}



# Load in results from COLONY
errors = data.frame(count = 1:50,
                    test_condition = NA,
                    exclusion = NA,
                    numFP = NA,
                    numFN = NA,
                    numTP = NA,
                    total_real = NA,
                    FPR = NA,
                    FNR = NA)
family_plots = list()
count = 1
for (i in 1:nsim){
  for (j in 1:length(subsets)){
    for (k in c("exclusion")){
      for (species in c("mixtus", "impatiens")){
      name = paste0(species, "_set", i, "_sub", subsets[j], "_", k)
      print(paste0("Loading ", name))
      filename = paste0("simulate_data/colony_assignments/test_sample_size/colony_output/", name, ".BestCluster")
      colony_output = as.data.frame(do.call(rbind, strsplit(trimws(readLines(filename)), "\\s+")[-1]))
      colnames(colony_output) = unlist(strsplit(readLines(filename), "\\s+")[1])
      true_data = read.csv(paste0("simulate_data/colony_assignments/test_sample_size/true_data/", species, "_set", i, "_sub", subsets[j], ".csv"))
      
      print(paste0("Files loaded for ", name))
      
      # set probability threshold
      prob_thresh = 0.95
      
      # filter colony outputs
      colony_output = colony_output %>% filter(Probability >= prob_thresh)
      
      # make edge lists
      true_edges = true_data %>%
        group_by(truecolony) %>%
        filter(n() > 1) %>%
        summarise(pairs = combn(individual, 2, simplify = FALSE), .groups = "drop") %>%
        mutate(from = map_chr(pairs, 1),
               to = map_chr(pairs, 2)) %>%
        dplyr::select(from, to)
      
      inferred_edges = colony_output %>%
        group_by(ClusterIndex) %>%
        filter(n() > 1) %>%
        summarise(pairs = combn(OffspringID, 2, simplify = FALSE), .groups = "drop") %>%
        mutate(from = map_chr(pairs, 1),
               to = map_chr(pairs, 2)) %>%
        dplyr::select(from, to)
      
      # combine and classify edges by type (FN, FP, TP)
      # get true positives
      tp_edges = inner_join(true_edges, inferred_edges, by = c("from", "to"))
      tp_edges$type = "TP"
      
      # initialize other types as FN and FP
      true_edges$type = "FN"
      inferred_edges$type = "FP"
      
      # remove true positives from FN and FP dataframes
      fn_edges = anti_join(true_edges, tp_edges, by = c("from", "to"))
      fp_edges = anti_join(inferred_edges, tp_edges, by = c("from", "to"))
      
      # combine all edges
      all_edges = rbind(tp_edges, fn_edges, fp_edges)
      
      print(paste0("Edges computed for ", name))
      
      # make an igraph object
      graph = graph_from_data_frame(all_edges, directed = FALSE)
      E(graph)$color = recode(E(graph)$type, TP = "black", FP = "red", FN = "purple")
      families = plot(graph, 
                      edge.color = E(graph)$color,
                      vertex.label = NA,
                      vertex.size = 2,
                      main = name
      )
      family_plots[[count]] = families
      
      # record FPR and FNR
      print(paste0("Start saving data for ", name))
      print(paste0("Count = ", count))
      errors$test_condition[count] = name
      print(paste0("Step 1 ", name))
      errors$FPR[count] = nrow(fp_edges) / (nrow(tp_edges) + nrow(fn_edges))
      errors$FNR[count] = nrow(fn_edges) / (nrow(tp_edges) + nrow(fn_edges))
      print(paste0("Midpoint data save for ", name))
      errors$total_real[count] = nrow(tp_edges) + nrow(fn_edges)
      errors$numFP[count] = nrow(fp_edges)
      errors$numFN[count] = nrow(fn_edges)
      errors$numTP[count] = nrow(tp_edges)
      errors$exclusion[count] = k
      
      print(paste0("Data saved for ", name))
      
      count = count +1
  
      }
    }
  }  
}
errors$sub_value = str_extract(errors$test_condition, "(?<=sub)[0-9.]+")
errors$numbees = 2000*as.numeric(errors$sub_value)

ggplot(errors) +
  geom_point(aes(x = numbees, y = numFP, color = "Number FP")) +
  geom_point(aes(x = numbees, y = numFN, color = "Number FN")) +
  geom_point(aes(x = numbees, y = total_real, color = "Total true")) +
  xlab("Simulation and COLONY Conditions") +
  ylab("Number of inferred or true relationships") +
  labs(title = "Sibship inclusion: P = 0.95") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(errors, aes(x = test_condition, y = FPR)) +
  geom_point() +
  xlab("Simulation and COLONY Conditions") +
  ylab(expression(FPR == frac(FP, TP + FN))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(errors, aes(x = test_condition, y = FNR)) +
  geom_point() +
  xlab("Simulation and COLONY Conditions") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(errors, aes(x = numbees, y = FPR)) +
  geom_point() +
  xlab("Simulation and COLONY Conditions") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  ylim(c(0,0.5))

