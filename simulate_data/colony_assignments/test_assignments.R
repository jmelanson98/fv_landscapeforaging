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
library(gridExtra)
library(grid)

################################################################################
## PRACTICE SIMULATING DATA --- IS IT REALISTIC?
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


 
 ################################################################################
 ## NOW START REAL SIMULATIONS FOR TESTING ASSUMPTIONS
 ################################################################################
 
# Generate five sibship data sets (e.g., 5 x yobs) of 2000 individuals each
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
 saveRDS(sibship_data, paste0("simulate_data/colony_assignments/sim_data/sim", i, ".RDS"))
}

 
# Load in allele frequency data
impatiens_allelefreq = read.csv("colony_assignments/Colony2_Linux/impatiens2023.AlleleFreq")
mixtus_allelefreq = read.csv("colony_assignments/Colony2_Linux/mixtus2023.AlleleFreq")

# Remove loci we are not using
mixtus_to_remove = c("BTMS0104", "BTERN01", "BTMS0059")
impatiens_to_remove = c("BL13", "BTMS0073", "BT28")

impatiens_allelefreq = impatiens_allelefreq %>% filter(!MarkerID %in% impatiens_to_remove)
mixtus_allelefreq = mixtus_allelefreq %>% filter(!MarkerID %in% mixtus_to_remove)


# Create real error rate files
mixtus_error_rates = data.frame(c("BT10", 0, 0, 0.01),
                                #c("BTMS0104", 0, 0, 0.01),
                                c("BTMS0057", 0, 0, 0.01),
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
write.table(mixtus_error_rates, "simulate_data/colony_assignments/sim_data/for_colony/mixtus_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
mixtus_errors = read.table("simulate_data/colony_assignments/sim_data/for_colony/mixtus_error_rates.txt", sep = ",")
colnames(mixtus_errors) = mixtus_errors[1,]
mixtus_errors = mixtus_errors[-1,]

impatiens_error_rates = data.frame(c("BT10", 0, 0, 0.017),
                                   c("B96", 0, 0, 0.027),
                                   c("BTMS0059", 0, 0, 0.017),
                                   c("BTMS0081", 0, 0, 0.016),
                                   #c("BL13", 0, 0, 0.011),
                                   c("BTMS0062", 0, 0, 0.022),
                                   c("B126", 0, 0, 0.01),
                                   c("BTERN01", 0, 0, 0.028),
                                   c("B124", 0, 0, 0.011),
                                   c("BTMS0057", 0, 0, 0.017),
                                   c("BT30", 0, 0, 0.01),
                                   c("B10", 0, 0, 0.029),
                                   c("BTMS0083", 0, 0, 0.01)
                                   #c("BTMS0073", 0, 0, 0.01),
                                   #c("BT28", 0, 0, 0.01)
)
write.table(impatiens_error_rates, "simulate_data/colony_assignments/sim_data/for_colony/impatiens_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
impatiens_errors = read.table("simulate_data/colony_assignments/sim_data/for_colony/impatiens_error_rates.txt", sep = ",")
colnames(impatiens_errors) = impatiens_errors[1,]
impatiens_errors = impatiens_errors[-1,]


# Get missingness rates!
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


# For each data set, simulate genotypes assuming single paternity (mixtus and impatiens)
impGenotypesList = list()
mixGenotypesList = list()

for (i in 1:numsims){
  sibship_data = readRDS(paste0("simulate_data/colony_assignments/sim_data/sim", i, ".RDS"))
  # simulate imp genotypes
  tempImp = simulateGenotypes(alleleFreqs = impatiens_allelefreq,
                              colonyDataDetected = sibship_data[[4]],
                              observationMatrix = sibship_data[[2]],
                              trapData = sibship_data[[5]],
                              probMultiplePaternity = 0)
  # induce errors and missingness
  impGenotypesList[[i]] = induceErrors(genotypeDF = tempImp,
                                       errorRates = impatiens_errors,
                                       missingRates = impatiens_missing,
                                       alleleFreqs = impatiens_allelefreq)
  
  # simulate mix genotypes
  tempMix = simulateGenotypes(alleleFreqs = mixtus_allelefreq,
                              colonyDataDetected = sibship_data[[4]],
                              observationMatrix = sibship_data[[2]],
                              trapData = sibship_data[[5]],
                              probMultiplePaternity = 0)
  # induce errors and missingnes
  mixGenotypesList[[i]] = induceErrors(genotypeDF = tempMix,
                                       errorRates = mixtus_errors,
                                       missingRates = mixtus_missing,
                                       alleleFreqs = mixtus_allelefreq)
}

sibship_genotypes = list(mixGenotypesList, impGenotypesList)
saveRDS(sibship_genotypes, "simulate_data/colony_assignments/test_sample_size/sibship_genotypes.RDS")
sibship_genotypes = readRDS("simulate_data/colony_assignments/test_sample_size/sibship_genotypes.RDS")
mixGenotypesList = sibship_genotypes[[1]]
impGenotypesList = sibship_genotypes[[2]]

# Now, subset each dataset to contain 100%, 80%, 60%, 40%, 20% of individuals
subsets = c(1, 0.8, 0.6, 0.4, 0.2)
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
                paste0("simulate_data/colony_assignments/sim_data/for_colony/mixtus_set", i, "_sub", subsets[j], ".txt"),
                sep= ",", col.names = FALSE, row.names = FALSE)
    # save full csv
    write.csv(genotypes_subset_mix, 
              paste0("simulate_data/colony_assignments/sim_data/true_data/mixtus_set", i, "_sub", subsets[j], ".csv"))
    # write sibship exclusion table
    write.table(
      subset_exclusion_mix,
      file = paste0("simulate_data/colony_assignments/sim_data/for_colony/mixtus_exclusion", i, "_sub", subsets[j], ".txt"),
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
                paste0("simulate_data/colony_assignments/sim_data/for_colony/impatiens_set", i, "_sub", subsets[j], ".txt"),
                sep= ",", col.names = FALSE, row.names = FALSE)
    # save full csv
    write.csv(genotypes_subset_imp, 
              paste0("simulate_data/colony_assignments/sim_data/true_data/impatiens_set", i, "_sub", subsets[j], ".csv"))
    # write exclusion table
    write.table(
      subset_exclusion_imp,
      file = paste0("simulate_data/colony_assignments/sim_data/for_colony/impatiens_exclusion", i, "_sub", subsets[j], ".txt"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      na = ""
    )
  }
}


################################################################################
## Test performance of varying Pinc threshold vs repeat runs
################################################################################
# use 0.6X datasets, e.g., 1200 individuals
nsims = 5
nruns = 5
sub = 0.6
workingdir = "/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/test_repetition"
mixtus_errors_filepath = "simulate_data/colony_assignments/sim_data/for_colony/mixtus_error_rates.txt"
impatiens_errors_filepath = "simulate_data/colony_assignments/sim_data/for_colony/impatiens_error_rates.txt"

for (i in 1:nsims){
  for (j in 1:nruns){
    # get genotype filepaths
    mixtus_genotypes_filepath = paste0("simulate_data/colony_assignments/sim_data/for_colony/mixtus_set", i, "_sub0.6.txt")
    impatiens_genotypes_filepath = paste0("simulate_data/colony_assignments/sim_data/for_colony/impatiens_set", i, "_sub0.6.txt")
    
    # get exclusion paths
    mixtus_exclusion_filepath = paste0("simulate_data/colony_assignments/sim_data/for_colony/mixtus_exclusion", i, "_sub0.6.txt")
    impatiens_exclusion_filepath = paste0("simulate_data/colony_assignments/sim_data/for_colony/impatiens_exclusion", i, "_sub0.6.txt")
    
    
    #build .DAT for mixtus
    rcolony::build.colony.superauto(wd=workingdir, 
                                    name=paste0("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/test_repetition/mixtus_set", i, "_sub0.6_run", j, ".DAT"), 
                                    datasetname = paste0("mixtus_set", i, "_sub0.6_run", j),
                                    delim=",",
                                    sample_size = 1200,
                                    num_loci = 10,
                                    sibship_prior = 1,
                                    error_rates_path = mixtus_errors_filepath,
                                    genotypes_path = mixtus_genotypes_filepath,
                                    exclusion_path = mixtus_exclusion_filepath
    )
    
    # build .DAT for impatiens
    rcolony::build.colony.superauto(wd=workingdir, 
                                    name=paste0("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/test_repetition/impatiens_set", i, "_sub0.6_run", j, ".DAT"), 
                                    datasetname = paste0("impatiens_set", i, "_sub0.6_run", j),
                                    delim=",",
                                    sample_size = 1200,
                                    num_loci = 12,
                                    sibship_prior = 1,
                                    error_rates_path = impatiens_errors_filepath,
                                    genotypes_path = impatiens_genotypes_filepath,
                                    exclusion_path = impatiens_exclusion_filepath
    )
  }
}

# Load in results from COLONY
number_runs = c(1, 5)
prob_threshholds = c(1, 0.995, 0.99, 0.975, 0.95)
errors = data.frame(count = 1:200,
                    numbees = 1200,
                    species = NA,
                    numFP = NA,
                    numFN = NA,
                    numTP = NA,
                    total_real = NA,
                    FPR = NA,
                    FNR = NA,
                    prob_thresh = NA,
                    repetition = NA,
                    data_type = NA)
count = 1
for (i in 1:nsims){
  for (species in c("mixtus", "impatiens")){
    for (prob_thresh in prob_threshholds){
      for (runs in number_runs){
        
        # load in true data
        true_data = read.csv(paste0("simulate_data/colony_assignments/sim_data/true_data/", species, "_set", i, "_sub0.6.csv"))
        
        # make true edge lists
        true_edges = true_data %>%
          group_by(truecolony) %>%
          filter(n() > 1) %>%
          summarise(pairs = combn(individual, 2, simplify = FALSE), .groups = "drop") %>%
          mutate(from = map_chr(pairs, 1),
                 to = map_chr(pairs, 2)) %>%
          dplyr::select(from, to)
        
        # get inferred edges from appropriate COLONY output files
        family_inferred_edges = list()
        dyad_inferred_edges = list()
        for (index in 1:runs){
          name = paste0(species, "_set", i, "_sub0.6_run", index)
          
          # first get family files
          familyfile = paste0("simulate_data/colony_assignments/test_repetition/", name, ".BestCluster")
          family = as.data.frame(do.call(rbind, strsplit(trimws(readLines(familyfile)), "\\s+")[-1]))
          colnames(family) = unlist(strsplit(readLines(familyfile), "\\s+")[1])
          
          # then get dyad files
          dyadfile = paste0("simulate_data/colony_assignments/test_repetition/", name, ".FullSibDyad")
          dyads = read.table(dyadfile, sep = ",")
          colnames(dyads) = c("from", "to", "Probability")
          dyads = dyads[-1,]
          
          # filter based on probability threshhold
          family = family %>% filter(as.numeric(Probability) >= prob_thresh)
          dyads = dyads %>% filter(as.numeric(Probability) >= prob_thresh)
          
          # get edge lists for family method
          family_inferred_edges[[index]] = family %>%
            group_by(ClusterIndex) %>%
            filter(n() > 1) %>%
            summarise(pairs = combn(OffspringID, 2, simplify = FALSE), .groups = "drop") %>%
            mutate(from = map_chr(pairs, 1),
                   to = map_chr(pairs, 2)) %>%
            dplyr::select(from, to)
          
          # get edge lists for dyads method
          dyad_inferred_edges[[index]] = dyads[,1:2]
          
        }
        
        ######################################
        # Identify dyads which are present in all of the runs
        # first pairs list!
        consistent_families = family_inferred_edges[[1]]
        consistent_families$pairs =  apply(family_inferred_edges[[1]], 1, paste, collapse = "-")
        consistent_dyads = dyad_inferred_edges[[1]]
        consistent_dyads$pairs = apply(dyad_inferred_edges[[1]], 1, paste, collapse = "-")
        
        family_keep = consistent_families$pairs
        dyad_keep = consistent_dyads$pairs
        
        # compare to other pairs lists!
        for (index in 1:runs){
          current_family_pairs = apply(family_inferred_edges[[index]], 1, paste, collapse = "-")
          current_dyad_pairs = apply(dyad_inferred_edges[[index]], 1, paste, collapse = "-")
          
          # only keep those which are in both!
          family_keep = family_keep[family_keep %in% current_family_pairs]
          dyad_keep = dyad_keep[dyad_keep %in% current_dyad_pairs]
        }
        
        consistent_families = consistent_families %>% 
          filter(consistent_families$pairs %in% family_keep) %>%
          dplyr::select(-pairs)
        consistent_dyads = consistent_dyads %>% 
          filter(consistent_dyads$pairs %in% dyad_keep) %>%
          dplyr::select(-pairs)
        
        #######################################
        # combine and classify edges (FAMILY)
        # get true positives
        tp_edges_family = inner_join(true_edges, consistent_families, by = c("from", "to"))
        tp_edges_family$type = "TP"
        
        # initialize other types as FN and FP
        true_edges$type = "FN"
        consistent_families$type = "FP"
        
        # remove true positives from FN and FP dataframes
        fn_edges_family = anti_join(true_edges, tp_edges_family, by = c("from", "to"))
        fp_edges_family = anti_join(consistent_families, tp_edges_family, by = c("from", "to"))
        
        # combine all edges
        all_edges_family = rbind(tp_edges_family, fn_edges_family, fp_edges_family)
        
        
        #######################################
        # combine and classify edges (DYADS)
        # get true positives
        tp_edges_dyads = inner_join(true_edges, consistent_dyads, by = c("from", "to"))
        tp_edges_dyads$type = "TP"
        
        # initialize other types as FN and FP
        true_edges$type = "FN"
        consistent_dyads$type = "FP"
        
        # remove true positives from FN and FP dataframes
        fn_edges_dyads = anti_join(true_edges, tp_edges_dyads, by = c("from", "to"))
        fp_edges_dyads = anti_join(consistent_dyads, tp_edges_dyads, by = c("from", "to"))
        
        # combine all edges
        all_edges_dyads = rbind(tp_edges_dyads, fn_edges_dyads, fp_edges_dyads)
        
        ##########################################
        # Record results!
        # first for families
        errors$species[count] = species
        errors$numFP[count] = nrow(fp_edges_family)
        errors$numFN[count] = nrow(fn_edges_family)
        errors$numTP[count] = nrow(tp_edges_family)
        errors$total_real[count] = nrow(tp_edges_family) + nrow(fn_edges_family)
        errors$FPR[count] = nrow(fp_edges_family) / (nrow(tp_edges_family) + nrow(fp_edges_family))
        errors$FNR[count] = nrow(fn_edges_family) / (nrow(tp_edges_family) + nrow(fn_edges_family))
        errors$prob_thresh[count] = prob_thresh
        errors$repetition[count] = runs
        errors$data_type[count] = "family"
        
        # then for dyads
        errors$species[count + 1] = species
        errors$numFP[count + 1] = nrow(fp_edges_dyads)
        errors$numFN[count + 1] = nrow(fn_edges_dyads)
        errors$numTP[count + 1] = nrow(tp_edges_dyads)
        errors$total_real[count + 1] = nrow(tp_edges_dyads) + nrow(fn_edges_dyads)
        errors$FPR[count + 1] = nrow(fp_edges_dyads) / (nrow(tp_edges_dyads) + nrow(fp_edges_dyads))
        errors$FNR[count + 1] = nrow(fn_edges_dyads) / (nrow(tp_edges_dyads) + nrow(fn_edges_dyads))
        errors$prob_thresh[count + 1] = prob_thresh
        errors$repetition[count + 1] = runs
        errors$data_type[count + 1] = "dyads"
        
        
        count = count + 2
      }
    }
  }  
}

mixtusfamilyFPR = ggplot(errors[errors$species == "mixtus" & errors$data_type == "family",], aes(x = prob_thresh, y = FPR, colour = as.factor(repetition))) +
  geom_point() +
  xlab("") +
  ylab(expression(FPR == frac(FP, TP + FP))) +
  labs(colour = "COLONY \n runs") +
  scale_color_manual(
    #labels = c("Excluded", "Not excluded"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.1))
impatiensfamilyFPR = ggplot(errors[errors$species == "impatiens" & errors$data_type == "family",], aes(x = prob_thresh, y = FPR, colour = as.factor(repetition))) +
  geom_point() +
  xlab("Probability threshold") +
  ylab(expression(FPR == frac(FP, TP + FP))) +
  labs(colour = "COLONY \n runs") +
  scale_color_manual(
    #labels = c("Excluded", "Not excluded"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.1))
mixtusdyadsFPR = ggplot(errors[errors$species == "mixtus" & errors$data_type == "dyads",], aes(x = prob_thresh, y = FPR, colour = as.factor(repetition))) +
  geom_point() +
  xlab("") +
  ylab("") +
  labs(colour = "COLONY \n runs") +
  scale_color_manual(
    #labels = c("Excluded", "Not excluded"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.1))
impatiensdyadsFPR = ggplot(errors[errors$species == "impatiens" & errors$data_type == "dyads",], aes(x = prob_thresh, y = FPR, colour = as.factor(repetition))) +
  geom_point() +
  xlab("Probability threshold") +
  ylab("") +
  labs(colour = "COLONY \n runs") +
  scale_color_manual(
    #labels = c("Excluded", "Not excluded"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.1))
mixtusfamilyFNR = ggplot(errors[errors$species == "mixtus" & errors$data_type == "family",], aes(x = prob_thresh, y = FNR, colour = as.factor(repetition))) +
  geom_point() +
  xlab("") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  labs(colour = "COLONY \n runs") +
  scale_color_manual(
    #labels = c("Excluded", "Not excluded"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.1))
impatiensfamilyFNR = ggplot(errors[errors$species == "impatiens" & errors$data_type == "family",], aes(x = prob_thresh, y = FNR, colour = as.factor(repetition))) +
  geom_point() +
  xlab("Probability threshold") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  labs(colour = "COLONY \n runs") +
  scale_color_manual(
    #labels = c("Excluded", "Not excluded"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.1))
mixtusdyadsFNR = ggplot(errors[errors$species == "mixtus" & errors$data_type == "dyads",], aes(x = prob_thresh, y = FNR, colour = as.factor(repetition))) +
  geom_point() +
  xlab("") +
  ylab("") +
  labs(colour = "COLONY \n runs") +
  scale_color_manual(
    #labels = c("Excluded", "Not excluded"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.1))
impatiensdyadsFNR = ggplot(errors[errors$species == "impatiens" & errors$data_type == "dyads",], aes(x = prob_thresh, y = FNR, colour = as.factor(repetition))) +
  geom_point() +
  xlab("Probability threshold") +
  ylab("") +
  labs(colour = "COLONY \n runs") +
  scale_color_manual(
    #labels = c("Excluded", "Not excluded"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.1))

g = ggplotGrob(mixtusfamilyFPR)
legend_index = which(g$layout$name == "guide-box-right")
legend = g$grobs[[legend_index]]

# remove legend from plots
mixtusfamilyFPR = mixtusfamilyFPR + theme(legend.position = "none")
mixtusdyadsFPR = mixtusdyadsFPR + theme(legend.position = "none")
impatiensfamilyFPR = impatiensfamilyFPR + theme(legend.position = "none")
impatiensdyadsFPR = impatiensdyadsFPR + theme(legend.position = "none")
mixtusfamilyFNR = mixtusfamilyFNR + theme(legend.position = "none")
mixtusdyadsFNR = mixtusdyadsFNR + theme(legend.position = "none")
impatiensfamilyFNR = impatiensfamilyFNR + theme(legend.position = "none")
impatiensdyadsFNR = impatiensdyadsFNR + theme(legend.position = "none")

# make some text grobs
imp = textGrob(
  expression(italic("Bombus impatiens")),
  gp = gpar(fontsize = 12, col = "black"), rot = 270
)
mix = textGrob(
  expression(italic("Bombus mixtus")),
  gp = gpar(fontsize = 12, col = "black"), rot = 270
)
family = textGrob(
  "Best Family Clusters",
  gp = gpar(fontsize = 12, col = "black")
)
dyad = textGrob(
  "Full Sibling Dyads",
  gp = gpar(fontsize = 12, col = "black")
)


# arrange and plot
FPR_grid = grid.arrange(family, nullGrob(), dyad, nullGrob(),
                        mixtusfamilyFPR, nullGrob(), mixtusdyadsFPR, mix,
                        impatiensfamilyFPR, nullGrob(), impatiensdyadsFPR, imp, ncol =4, widths = c(8,1,8,1), heights = c(1,8,8))
FPR_grid = grid.arrange(FPR_grid, legend, ncol = 2, widths = c(8,2))
FPR_grid = ggdraw() +
  draw_plot(FPR_grid, 0.07, 0, 0.9, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)", "(D)"), size = 14, 
                  x = c(0, 0.4, 0, 0.4), 
                  y = c(1, 1, 0.53, 0.53))
ggsave("docs/appendix_figures/fpr_repetition.jpg", FPR_grid, height = 1000, width = 3000, units = "px")


FNR_grid = grid.arrange(family, nullGrob(), dyad, nullGrob(),
                        mixtusfamilyFNR, nullGrob(), mixtusdyadsFNR, mix,
                        impatiensfamilyFNR, nullGrob(), impatiensdyadsFNR, imp, ncol =4, widths = c(8,1,8,1), heights = c(1,8,8))
FNR_grid = grid.arrange(FNR_grid, legend, ncol = 2, widths = c(8,2))
FNR_grid = ggdraw() +
  draw_plot(FNR_grid, 0.07, 0, 0.9, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)", "(D)"), size = 14, 
                  x = c(0, 0.4, 0, 0.4), 
                  y = c(1, 1, 0.53, 0.53))
ggsave("docs/appendix_figures/fnr_repetition.jpg", FNR_grid, height = 1000, width = 3000, units = "px")



######################DO SOME EXPLORATION OF DIFFERENT NONCIRCULARITY RESOLUTIONS
##############
#########
#####
###
#

links = data.frame(count = 1:10,
                    species = NA,
                    numfams = NA,
                    numlinks = NA,
                    numFN = NA,
                    numTN = NA
)
links$FNprobs <- vector("list", nrow(links))
links$TNprobs <- vector("list", nrow(links))
count = 1
for (species in c("mixtus", "impatiens")){
  for (i in 1:5){
      print(paste0(species, "_set", i, "_sub0.6.csv"))
      
      # set probability threshold
      prob_thresh = 0.995
      
      # load in true data
      true_data = read.csv(paste0("simulate_data/colony_assignments/sim_data/true_data/", species, "_set", i, "_sub1.csv"))
      
      # make true edge lists
      true_edges = true_data %>%
        group_by(truecolony) %>%
        filter(n() > 1) %>%
        summarise(pairs = combn(individual, 2, simplify = FALSE), .groups = "drop") %>%
        mutate(from = map_chr(pairs, 1),
               to = map_chr(pairs, 2)) %>%
        dplyr::select(from, to)
      true_edges$missing_pair = apply(
        true_edges[, c("from", "to")], 
        1, 
        function(x) paste(sort(x), collapse = "-")
      )
      
      # then get dyad files
      dyadfile = paste0("simulate_data/colony_assignments/test_sample_size/colony_output/", species, "_set", i, "_sub1_exclusion.FullSibDyad")
      dyads = read.table(dyadfile, sep = ",")
      colnames(dyads) = c("from", "to", "Probability")
      dyads = dyads[-1,]
      dyads$Probability = as.numeric(dyads$Probability)
      
      # filter based on probability threshhold
      dyadsfiltered = dyads %>% filter(Probability >= prob_thresh)
      
      # add pair names
      dyads$missing_pair = apply(
        dyads[, c("from", "to")], 
        1, 
        function(x) paste(sort(x), collapse = "-")
      )
      
      # get edge lists for dyads method
      dyad_inferred_edges = dyadsfiltered[,1:2]
      
      # make igraph object
      graph = graph_from_data_frame(dyad_inferred_edges, directed = FALSE)
      components <- decompose(graph)
      
      # remove cliques
      noncirculars = components[!sapply(components, is_clique)]
      noncircular_graph <- do.call(disjoint_union, noncirculars)
      noncircular_comps = decompose(noncircular_graph)
      
      # plot
      if(!is.null(noncircular_graph)){
        #E(noncircular_graph)$color = recode(E(noncircular_graph)$type, TP = "black", FP = "red", FN = "purple")
        families = plot(noncircular_graph)
      } else (print("No noncircularity."))
      
      
      missing_links = lapply(components, function(comp) {
        if (!is_clique(comp)){
          # get all possible pairs of nodes (unordered)
          nodes = V(comp)$name
          node_pairs = t(combn(nodes, 2))
          
          # filter out pairs that already have an edge
          missing <- apply(node_pairs, 1, function(pair) {
            !are_adjacent(comp, pair[1], pair[2])
          })
          
          # make collapsed pair names for eaching missing link
          missing_names = apply(node_pairs[missing, , drop = FALSE], 1, function(pair) {
            paste(sort(pair), collapse = "-")
          })
          
          df = data.frame(
            missing_pair = missing_names,
            stringsAsFactors = FALSE
          )
          return(df)
        }
      })
      
      # Combine all rows into one dataframe
      missing_df = do.call(rbind, missing_links)
      
      # Find false negatives
      false_negatives = missing_df$missing_pair[missing_df$missing_pair %in% true_edges$missing_pair]
      true_negatives = missing_df$missing_pair[!missing_df$missing_pair %in% true_edges$missing_pair]
    
      
      # save some data
      links$numfams[count] = length(noncircular_comps)
      links$numlinks[count] = max(nrow(missing_df),0)
      links$numFN[count] = length(false_negatives)
      links$numTN[count] = length(true_negatives)
      links$FNprobs[count] <- list(dyads$Probability[dyads$missing_pair %in% false_negatives])
      links$TNprobs[count] <- list(dyads$Probability[dyads$missing_pair %in% true_negatives])
      
      
      links$species[count] = species
      
      count = count + 1
    
  }
}

links$FNprobs <- lapply(links$FNprobs, function(x) x[x != 0])
links$TNprobs <- lapply(links$TNprobs, function(x) x[x != 0])

link_mixtus = links[links$species == "mixtus",]
link_impatiens = links[links$species == "impatiens",]
# make histogram
png(filename = "docs/appendix_figures/noncircularity.png", width = 800, height = 600)
hist(unlist(link_impatiens$FNprobs), col = rgb(0, 0, 1, 0.5), xlim = c(0.6,1), ylim = c(0,10),
     main = "", xlab = "Probability", breaks = 30)

hist(unlist(link_impatiens$TNprobs), col = rgb(1, 0, 0, 0.5), add = TRUE, breaks = 30)
abline(v = 0.95, col = "black", lty = 2, lwd = 2)  # lty=2 makes it dashed

legend("topleft", legend = c("FN", "TN"),
       fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))
dev.off()


################################################################################
## Incorporate noncircularity resolution into dyad workflow and compare to families
################################################################################
prob_threshholds = c(1, 0.995, 0.99, 0.975, 0.95)
errors = data.frame(count = 1:100,
                    numbees = 1200,
                    species = NA,
                    numFP = NA,
                    numFN = NA,
                    numTP = NA,
                    total_real = NA,
                    FPR = NA,
                    FNR = NA,
                    prob_thresh = NA,
                    data_type = NA)
count = 1
for (i in 1:nsims){
  for (species in c("mixtus", "impatiens")){
    for (prob_thresh in prob_threshholds){
      
      # load in true data
      true_data = read.csv(paste0("simulate_data/colony_assignments/sim_data/true_data/", species, "_set", i, "_sub0.6.csv"))
      
      # make true edge lists
      true_edges = true_data %>%
        group_by(truecolony) %>%
        filter(n() > 1) %>%
        summarise(pairs = combn(individual, 2, simplify = FALSE), .groups = "drop") %>%
        mutate(from = map_chr(pairs, 1),
               to = map_chr(pairs, 2)) %>%
        dplyr::select(from, to)
      
      # get inferred edges from appropriate COLONY output files
      name = paste0(species, "_set", i, "_sub0.6_run1")
        
      ###################################
      # FIRST BY FAMILY METHOD
      familyfile = paste0("simulate_data/colony_assignments/test_repetition/", name, ".BestCluster")
      family = as.data.frame(do.call(rbind, strsplit(trimws(readLines(familyfile)), "\\s+")[-1]))
      colnames(family) = unlist(strsplit(readLines(familyfile), "\\s+")[1])
      
      # filter based on probability threshhold
      family = family %>% filter(as.numeric(Probability) >= prob_thresh)

      # get edge lists for family method
      family_inferred_edges = family %>%
        group_by(ClusterIndex) %>%
        filter(n() > 1) %>%
        summarise(pairs = combn(OffspringID, 2, simplify = FALSE), .groups = "drop") %>%
        mutate(from = map_chr(pairs, 1),
               to = map_chr(pairs, 2)) %>%
        dplyr::select(from, to)
      
      #####################################
      # THEN BY DYAD METHOD
      dyadfile = paste0("simulate_data/colony_assignments/test_repetition/", name, ".FullSibDyad")
      dyads = read.table(dyadfile, sep = ",")
      colnames(dyads) = c("from", "to", "Probability")
      dyads = dyads[-1,]
      dyads$Probability = as.numeric(dyads$Probability)
      
      # filter based on probability threshhold
      dyadsfiltered = dyads %>% filter(Probability >= prob_thresh)
      
      # add pair names
      dyads$missing_pair = apply(
        dyads[, c("from", "to")], 
        1, 
        function(x) paste(sort(x), collapse = "-")
      )
      
      # get edge lists for dyads method
      dyad_inferred_edges = dyadsfiltered[,1:2]
      
      # make igraph object
      graph = graph_from_data_frame(dyad_inferred_edges, directed = FALSE)
      components <- decompose(graph)
      
      # find missing links
      missing_links = lapply(components, function(comp) {
        if (!is_clique(comp)){
          # get all possible pairs of nodes (unordered)
          nodes = V(comp)$name
          node_pairs = t(combn(nodes, 2))
          
          # make a list of node pairs which are no adjacent in the component
          missing <- apply(node_pairs, 1, function(pair) {
            !are_adjacent(comp, pair[1], pair[2])
          })
          
          # make collapsed pair names for each missing link
          missing_names = apply(node_pairs[missing, , drop = FALSE], 1, function(pair) {
            paste(sort(pair), collapse = "-")
          })
          
          df = data.frame(
            missing_pair = missing_names,
            stringsAsFactors = FALSE
          )
          return(df)
        }
      })
      missing_df = do.call(rbind, missing_links)
      
      # Maintain sib pairs that are highly likely or which are missing with P > 0.95
      dyadsfiltered = dyads %>% filter(Probability >= prob_thresh |
                                         (missing_pair %in% missing_df$missing_pair & Probability >= 0.95))
      # Get new edge list
      dyad_inferred_edges = dyadsfiltered[,1:2]
      
      
      # Resolve remaining noncircularity
      # Make new graph
      graph = graph_from_data_frame(dyad_inferred_edges, directed = FALSE)
      components <- decompose(graph)
      
      # For each component, keep only the largest clique
      filtered_components <- lapply(components, function(comp) {
        if (is_clique(comp)) {
          return(comp)
        } else {
          cliques <- largest_cliques(comp)
          # randomly choose one if there's a tie
          chosen_clique <- sample(cliques, 1)[[1]]
          # induce subgraph on that clique
          return(induced_subgraph(comp, chosen_clique))
        }
      })
      
      # Recombine the components into a single graph
      final_graph = do.call(disjoint_union, filtered_components)
      
      # get final dyad edges
      dyad_inferred_edges <- as.data.frame(as_edgelist(final_graph))
      colnames(dyad_inferred_edges) = c("from", "to")

      
      #######################################
      # combine and classify edges (FAMILY)
      # get true positives
      tp_edges_family = inner_join(true_edges, family_inferred_edges, by = c("from", "to"))
      tp_edges_family$type = "TP"
      
      # initialize other types as FN and FP
      true_edges$type = "FN"
      family_inferred_edges$type = "FP"
      
      # remove true positives from FN and FP dataframes
      fn_edges_family = anti_join(true_edges, tp_edges_family, by = c("from", "to"))
      fp_edges_family = anti_join(family_inferred_edges, tp_edges_family, by = c("from", "to"))
      
      # combine all edges
      all_edges_family = rbind(tp_edges_family, fn_edges_family, fp_edges_family)
      
      
      #######################################
      # combine and classify edges (DYADS)
      # get true positives
      tp_edges_dyads = inner_join(true_edges, dyad_inferred_edges, by = c("from", "to"))
      tp_edges_dyads$type = "TP"
      
      # initialize other types as FN and FP
      true_edges$type = "FN"
      dyad_inferred_edges$type = "FP"
      
      # remove true positives from FN and FP dataframes
      fn_edges_dyads = anti_join(true_edges, tp_edges_dyads, by = c("from", "to"))
      fp_edges_dyads = anti_join(dyad_inferred_edges, tp_edges_dyads, by = c("from", "to"))
      
      # combine all edges
      all_edges_dyads = rbind(tp_edges_dyads, fn_edges_dyads, fp_edges_dyads)
      
      ##########################################
      # Record results!
      # first for families
      errors$species[count] = species
      errors$numFP[count] = nrow(fp_edges_family)
      errors$numFN[count] = nrow(fn_edges_family)
      errors$numTP[count] = nrow(tp_edges_family)
      errors$total_real[count] = nrow(tp_edges_family) + nrow(fn_edges_family)
      errors$FPR[count] = nrow(fp_edges_family) / (nrow(tp_edges_family) + nrow(fp_edges_family))
      errors$FNR[count] = nrow(fn_edges_family) / (nrow(tp_edges_family) + nrow(fn_edges_family))
      errors$prob_thresh[count] = prob_thresh
      errors$data_type[count] = "family"
      
      # then for dyads
      errors$species[count + 1] = species
      errors$numFP[count + 1] = nrow(fp_edges_dyads)
      errors$numFN[count + 1] = nrow(fn_edges_dyads)
      errors$numTP[count + 1] = nrow(tp_edges_dyads)
      errors$total_real[count + 1] = nrow(tp_edges_dyads) + nrow(fn_edges_dyads)
      errors$FPR[count + 1] = nrow(fp_edges_dyads) / (nrow(tp_edges_dyads) + nrow(fp_edges_dyads))
      errors$FNR[count + 1] = nrow(fn_edges_dyads) / (nrow(tp_edges_dyads) + nrow(fn_edges_dyads))
      errors$prob_thresh[count + 1] = prob_thresh
      errors$data_type[count + 1] = "dyads"
      
      
      count = count + 2
    }
  }
}  


mixtusfamilyFPR = ggplot(errors[errors$species == "mixtus" & errors$data_type == "family",], aes(x = prob_thresh, y = FPR, alpha = 0.4)) +
  geom_point() +
  xlab("") +
  ylab(expression(FPR == frac(FP, TP + FP))) +
  labs(colour = "COLONY \n runs") +
  theme_minimal() +
  ylim(c(0,0.1)) +
  theme(legend.position = "none")
impatiensfamilyFPR = ggplot(errors[errors$species == "impatiens" & errors$data_type == "family",], aes(x = prob_thresh, y = FPR, alpha = 0.4)) +
  geom_point() +
  xlab("Probability threshold") +
  ylab(expression(FPR == frac(FP, TP + FP))) +
  labs(colour = "COLONY \n runs") +
  theme_minimal() +
  ylim(c(0,0.1)) +
  theme(legend.position = "none")
mixtusdyadsFPR = ggplot(errors[errors$species == "mixtus" & errors$data_type == "dyads",], aes(x = prob_thresh, y = FPR, alpha = 0.4)) +
  geom_point() +
  xlab("") +
  ylab("") +
  labs(colour = "COLONY \n runs") +
  theme_minimal() +
  ylim(c(0,0.1)) +
  theme(legend.position = "none")
impatiensdyadsFPR = ggplot(errors[errors$species == "impatiens" & errors$data_type == "dyads",], aes(x = prob_thresh, y = FPR, alpha = 0.4)) +
  geom_point() +
  xlab("Probability threshold") +
  ylab("") +
  labs(colour = "COLONY \n runs") +
  theme_minimal() +
  ylim(c(0,0.1)) +
  theme(legend.position = "none")
mixtusfamilyFNR = ggplot(errors[errors$species == "mixtus" & errors$data_type == "family",], aes(x = prob_thresh, y = FNR, alpha = 0.4)) +
  geom_point() +
  xlab("") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  labs(colour = "COLONY \n runs") +
  theme_minimal() +
  ylim(c(0,0.1)) +
  theme(legend.position = "none")
impatiensfamilyFNR = ggplot(errors[errors$species == "impatiens" & errors$data_type == "family",], aes(x = prob_thresh, y = FNR, alpha = 0.4)) +
  geom_point() +
  xlab("Probability threshold") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  labs(colour = "COLONY \n runs") +
  theme_minimal() +
  ylim(c(0,0.1)) +
  theme(legend.position = "none")
mixtusdyadsFNR = ggplot(errors[errors$species == "mixtus" & errors$data_type == "dyads",], aes(x = prob_thresh, y = FNR, alpha = 0.4)) +
  geom_point() +
  xlab("") +
  ylab("") +
  labs(colour = "COLONY \n runs") +
  theme_minimal() +
  ylim(c(0,0.1)) +
  theme(legend.position = "none")
impatiensdyadsFNR = ggplot(errors[errors$species == "impatiens" & errors$data_type == "dyads",], aes(x = prob_thresh, y = FNR, alpha = 0.4)) +
  geom_point() +
  xlab("Probability threshold") +
  ylab("") +
  labs(colour = "COLONY \n runs") +
  theme_minimal() +
  ylim(c(0,0.1)) +
  theme(legend.position = "none")

# make some text grobs
imp = textGrob(
  expression(italic("Bombus impatiens")),
  gp = gpar(fontsize = 12, col = "black"), rot = 270
)
mix = textGrob(
  expression(italic("Bombus mixtus")),
  gp = gpar(fontsize = 12, col = "black"), rot = 270
)
family = textGrob(
  "Best Family Clusters",
  gp = gpar(fontsize = 12, col = "black")
)
dyad = textGrob(
  "Full Sibling Dyads",
  gp = gpar(fontsize = 12, col = "black")
)


# arrange and plot
FPR_grid = grid.arrange(family, nullGrob(), dyad, nullGrob(),
                        mixtusfamilyFPR, nullGrob(), mixtusdyadsFPR, mix,
                        impatiensfamilyFPR, nullGrob(), impatiensdyadsFPR, imp, ncol =4, widths = c(8,1,8,1), heights = c(1,8,8))
FPR_grid = ggdraw() +
  draw_plot(FPR_grid, 0.07, 0, 0.9, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)", "(D)"), size = 14, 
                  x = c(0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.53, 0.53))
ggsave("docs/appendix_figures/fpr_dyadsvfamilies.jpg", FPR_grid, height = 1000, width = 3000, units = "px")


FNR_grid = grid.arrange(family, nullGrob(), dyad, nullGrob(),
                        mixtusfamilyFNR, nullGrob(), mixtusdyadsFNR, mix,
                        impatiensfamilyFNR, nullGrob(), impatiensdyadsFNR, imp, ncol =4, widths = c(8,1,8,1), heights = c(1,8,8))
FNR_grid = ggdraw() +
  draw_plot(FNR_grid, 0.07, 0, 0.9, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)", "(D)"), size = 14, 
                  x = c(0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.53, 0.53))
ggsave("docs/appendix_figures/fpr_dyadsvfamilies.jpg", FNR_grid, height = 1000, width = 3000, units = "px")


################################################################################
## Test accuracy across dataset sizes, using exclusion tables (or not)
################################################################################
# Construct .DAT files for COLONY
mixtus_errors_filepath = "simulate_data/colony_assignments/sim_data/for_colony/mixtus_error_rates.txt"
impatiens_errors_filepath = "simulate_data/colony_assignments/sim_data/for_colony/impatiens_error_rates.txt"

for (i in 1:nsims){
  for (j in 1:length(subsets)){
    for (k in c("exclusion", "no_exclusion")){
    # get sample size, working directory
    size = nrow(mixGenotypesList[[i]]) * subsets[j]
    workingdir = "/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/test_sample_size/colony_output"
    
    # get genotype filepaths
    mixtus_genotypes_filepath = paste0("simulate_data/colony_assignments/sim_data/for_colony/mixtus_set", i, "_sub", subsets[j], ".txt")
    impatiens_genotypes_filepath = paste0("simulate_data/colony_assignments/sim_data/for_colony/impatiens_set", i, "_sub", subsets[j], ".txt")
    
    # get exclusion paths
    if (k == "exclusion"){
      mixtus_exclusion_filepath = paste0("simulate_data/colony_assignments/sim_data/for_colony/mixtus_exclusion", i, "_sub", subsets[j], ".txt")
      impatiens_exclusion_filepath = paste0("simulate_data/colony_assignments/sim_data/for_colony/impatiens_exclusion", i, "_sub", subsets[j], ".txt")
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
                                    num_loci = 10,
                                    sibship_prior = 1,
                                    female_monogamy = 1,
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
                                    num_loci = 12,
                                    sibship_prior = 1,
                                    female_monogamy = 1,
                                    error_rates_path = impatiens_errors_filepath,
                                    genotypes_path = impatiens_genotypes_filepath,
                                    exclusion_path = impatiens_exclusion_filepath
    )
    }
  }
}

# Load in results from COLONY
### Check error rates with dyads...
errors_dyad = data.frame(count = 1:100,
                    test_condition = NA,
                    exclusion = NA,
                    numFP = NA,
                    numFN = NA,
                    numTP = NA,
                    total_real = NA,
                    FPR = NA,
                    FNR = NA)
count = 1
for (i in 1:nsims){
  for (j in 1:length(subsets)){
    for (k in c("exclusion", "no_exclusion")){
      for (species in c("mixtus", "impatiens")){
        name = paste0(species, "_set", i, "_sub", subsets[j], "_", k)
        filename = paste0("simulate_data/colony_assignments/test_sample_size/colony_output/", name, ".FullSibDyad")
        dyads = read.table(filename, sep = ",")
        colnames(dyads) = c("from", "to", "Probability")
        dyads = dyads[-1,]
        dyads$Probability = as.numeric(dyads$Probability)
        true_data = read.csv(paste0("simulate_data/colony_assignments/sim_data/true_data/", species, "_set", i, "_sub", subsets[j], ".csv"))
        
        
        # set probability threshold
        prob_thresh = 0.995

        # make edge lists
        true_edges = true_data %>%
          group_by(truecolony) %>%
          filter(n() > 1) %>%
          summarise(pairs = combn(individual, 2, simplify = FALSE), .groups = "drop") %>%
          mutate(from = map_chr(pairs, 1),
                 to = map_chr(pairs, 2)) %>%
          dplyr::select(from, to)
        

        # filter based on probability threshhold
        dyadsfiltered = dyads %>% filter(Probability >= prob_thresh)
        
        # add pair names
        dyads$missing_pair = apply(
          dyads[, c("from", "to")], 
          1, 
          function(x) paste(sort(x), collapse = "-")
        )
        
        # get edge lists for dyads method
        dyad_inferred_edges = dyadsfiltered[,1:2]
        
        # make igraph object
        graph = graph_from_data_frame(dyad_inferred_edges, directed = FALSE)
        components <- decompose(graph)
        
        # find missing links
        missing_links = lapply(components, function(comp) {
          if (!is_clique(comp)){
            # get all possible pairs of nodes (unordered)
            nodes = V(comp)$name
            node_pairs = t(combn(nodes, 2))
            
            # make a list of node pairs which are no adjacent in the component
            missing <- apply(node_pairs, 1, function(pair) {
              !are_adjacent(comp, pair[1], pair[2])
            })
            
            # make collapsed pair names for each missing link
            missing_names = apply(node_pairs[missing, , drop = FALSE], 1, function(pair) {
              paste(sort(pair), collapse = "-")
            })
            
            df = data.frame(
              missing_pair = missing_names,
              stringsAsFactors = FALSE
            )
            return(df)
          }
        })
        missing_df = do.call(rbind, missing_links)
        
        # Maintain sib pairs that are highly likely or which are missing with P > 0.95
        dyadsfiltered = dyads %>% filter(Probability >= prob_thresh |
                                           (missing_pair %in% missing_df$missing_pair & Probability >= 0.95))
        # Get new edge list
        dyad_inferred_edges = dyadsfiltered[,1:2]
        
        
        # Resolve remaining noncircularity
        # Make new graph
        graph = graph_from_data_frame(dyad_inferred_edges, directed = FALSE)
        components <- decompose(graph)
        
        # For each component, keep only the largest clique
        filtered_components <- lapply(components, function(comp) {
          if (is_clique(comp)) {
            return(comp)
          } else {
            cliques <- largest_cliques(comp)
            # randomly choose one if there's a tie
            chosen_clique <- sample(cliques, 1)[[1]]
            # induce subgraph on that clique
            return(induced_subgraph(comp, chosen_clique))
          }
        })
        
        # Recombine the components into a single graph
        final_graph = do.call(disjoint_union, filtered_components)
        
        # get final dyad edges
        inferred_edges <- as.data.frame(as_edgelist(final_graph))
        colnames(inferred_edges) = c("from", "to")
        
        
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
        
        
        # record FPR and FNR
        errors_dyad$test_condition[count] = name
        errors_dyad$FPR[count] = nrow(fp_edges) / (nrow(tp_edges) + nrow(fp_edges))
        errors_dyad$FNR[count] = nrow(fn_edges) / (nrow(tp_edges) + nrow(fn_edges))
        errors_dyad$total_real[count] = nrow(tp_edges) + nrow(fn_edges)
        errors_dyad$numFP[count] = nrow(fp_edges)
        errors_dyad$numFN[count] = nrow(fn_edges)
        errors_dyad$numTP[count] = nrow(tp_edges)
        errors_dyad$exclusion[count] = k
        errors_dyad$species[count] = species
        
        count = count +1
        
      }
    }
  }  
}
errors_dyad$sub_value = str_extract(errors_dyad$test_condition, "(?<=sub)[0-9.]+")
errors_dyad$numbees = 2000*as.numeric(errors_dyad$sub_value)


###### Make some plots!

mixtus_FPR_dyads = ggplot(errors_dyad[errors_dyad$species=="mixtus",], aes(x = numbees, y = FPR, colour = exclusion)) +
  geom_point() +
  xlab("") +
  ylab(expression(FPR == frac(FP, TP + FP))) +
  labs(colour = "Across-site sibships") +
  scale_color_manual(
    labels = c("Excluded", "Not excluded"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.5))

mixtus_FNR_dyads = ggplot(errors_dyad[errors_dyad$species=="mixtus",], aes(x = numbees, y = FNR, colour = exclusion)) +
  geom_point() +
  xlab("Number of Observed Bees") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  labs(colour = "Across-site sibships") +
  scale_color_manual(
    labels = c("Excluded", "Not excluded"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.1))

impatiens_FPR_dyads = ggplot(errors_dyad[errors_dyad$species=="impatiens",], aes(x = numbees, y = FPR, colour = exclusion)) +
  geom_point() +
  xlab("") +
  ylab("") +
  labs(colour = "Across-site sibships") +
  scale_color_manual(
    labels = c("Excluded", "Not excluded"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.5))

impatiens_FNR_dyads = ggplot(errors_dyad[errors_dyad$species=="impatiens",], aes(x = numbees, y = FNR, colour = exclusion)) +
  geom_point() +
  xlab("Number of Observed Bees") +
  ylab("") +
  labs(colour = "Across-site sibships") +
  scale_color_manual(
    labels = c("Excluded", "Not excluded"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal()  +
  ylim(c(0,0.1))


# Make some grid plots for appendix
#get legend
g = ggplotGrob(mixtus_FPR_dyads)
legend_index = which(g$layout$name == "guide-box-right")
legend = g$grobs[[legend_index]]

# remove legend from plots
mixtus_FPR_dyads = mixtus_FPR_dyads + theme(legend.position = "none")
mixtus_FNR_dyads = mixtus_FNR_dyads + theme(legend.position = "none")
impatiens_FPR_dyads = impatiens_FPR_dyads + theme(legend.position = "none")
impatiens_FNR_dyads = impatiens_FNR_dyads + theme(legend.position = "none")

# make some text grobs
imp = textGrob(
  expression(italic("Bombus impatiens")),
  gp = gpar(fontsize = 12, col = "black")
)
mix = textGrob(
  expression(italic("Bombus mixtus")),
  gp = gpar(fontsize = 12, col = "black")
)
fpr = textGrob(
  "False Positives",
  gp = gpar(fontsize = 12, col = "black"), rot = 270
)
fnr = textGrob(
  "False Negatives",
  gp = gpar(fontsize = 12, col = "black"), rot = 270
)


# arrange and plot
grid = plot_grid(
  mix,        nullGrob(),      imp,        nullGrob(),
  mixtus_FPR_dyads, nullGrob(), impatiens_FPR_dyads, fpr,
  mixtus_FNR_dyads, nullGrob(), impatiens_FNR_dyads, fnr,
  ncol = 4,
  rel_widths  = c(8, 1, 8, 1),
  rel_heights = c(1, 8, 8),
  align = "hv", axis = "tblr"
)
grid = grid.arrange(grid, legend, ncol = 2, widths = c(8,2))
grid = ggdraw() +
  draw_plot(grid, 0.05, 0, 0.9, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)", "(D)"), size = 14, 
                  x = c(0, 0.4, 0, 0.4), 
                  y = c(1, 1, 0.53, 0.53))
ggsave("docs/appendix_figures/excl_size.jpg", grid, height = 1000, width = 3000, units = "px")

# Number of false positives doesn't *decrease,* but total number of relationships *increases* to lower FPR
ggplot(errors_dyad) +
  geom_point(aes(x = numbees, y = numFP, color = "Number FP", shape = exclusion)) +
  geom_point(aes(x = numbees, y = numFN, color = "Number FN", shape = exclusion)) +
  geom_point(aes(x = numbees, y = total_real, color = "Total true", shape = exclusion)) +
  xlab("Number of Bees") +
  ylab("Number of inferred or true relationships") +
  labs(title = "Sibship inclusion: P = 0.995") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))


################################################################################
## Test use of sibship prior
################################################################################
# Construct .DAT files for COLONY
mixtus_errors_filepath = "simulate_data/colony_assignments/test_sample_size/for_colony/mixtus_error_rates.txt"
impatiens_errors_filepath = "simulate_data/colony_assignments/test_sample_size/for_colony/impatiens_error_rates.txt"
workingdir = "/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/Colony2_Linux"

for (i in 1:nsims){
  for (j in 1:length(subsets)){
    for (k in c("exclusion")){
      # get sample size, working directory
      size = nrow(mixGenotypesList[[i]]) * subsets[j]
      
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
                                      name=paste0("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/test_sibprior/mixtus_set", i, "_sub", subsets[j], "_", k, "nosibprior.DAT"), 
                                      datasetname = paste0("mixtus_set", i, "_sub", subsets[j], "_", k, "_nosibprior"),
                                      delim=",",
                                      sample_size = size,
                                      num_loci = 10,
                                      error_rates_path = mixtus_errors_filepath,
                                      genotypes_path = mixtus_genotypes_filepath,
                                      exclusion_path = mixtus_exclusion_filepath
      )
      
      # build .DAT for impatiens
      rcolony::build.colony.superauto(wd=workingdir, 
                                      name=paste0("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/simulate_data/colony_assignments/test_sibprior/impatiens_set", i, "_sub", subsets[j], "_", k, "nosibprior.DAT"), 
                                      datasetname = paste0("impatiens_set", i, "_sub", subsets[j], "_", k, "_nosibprior"),
                                      delim=",",
                                      sample_size = size,
                                      num_loci = 12,
                                      error_rates_path = impatiens_errors_filepath,
                                      genotypes_path = impatiens_genotypes_filepath,
                                      exclusion_path = impatiens_exclusion_filepath
      )
    }
  }
}


# Load in results from COLONY
### Check error rates with dyads...
errors_dyad = data.frame(count = 1:100,
                         numbees = NA,
                         species = NA,
                         numFP = NA,
                         numFN = NA,
                         numTP = NA,
                         total_real = NA,
                         FPR = NA,
                         FNR = NA,
                         sibprior = NA)
count = 1
for (i in 1:nsims){
  for (j in 1:length(subsets)){
    for (species in c("mixtus", "impatiens")){
      name = paste0(species, "_set", i, "_sub", subsets[j], "_exclusion")
      true_data = read.csv(paste0("simulate_data/colony_assignments/sim_data/true_data/", species, "_set", i, "_sub", subsets[j], ".csv"))
      
      # load in no prior data
      noprior = read.table(paste0("simulate_data/colony_assignments/test_sibprior/", name, "_nosibprior.FullSibDyad"), sep = ",")
      colnames(noprior) = c("from", "to", "Probability")
      noprior = noprior[-1,]
      noprior$Probability = as.numeric(noprior$Probability)
      
      # load in with prior data
      prior = read.table(paste0("simulate_data/colony_assignments/test_sample_size/colony_output/", name, ".FullSibDyad"), sep = ",")
      colnames(prior) = c("from", "to", "Probability")
      prior = prior[-1,]
      prior$Probability = as.numeric(prior$Probability)
      
      # set probability threshold
      prob_thresh = 0.995
      
      # filter colony outputs
      nopriorfiltered = noprior %>% filter(as.numeric(Probability) >= prob_thresh)
      priorfiltered = prior %>% filter(as.numeric(Probability) >= prob_thresh)
      
      # make edge lists
      true_edges = true_data %>%
        group_by(truecolony) %>%
        filter(n() > 1) %>%
        summarise(pairs = combn(individual, 2, simplify = FALSE), .groups = "drop") %>%
        mutate(from = map_chr(pairs, 1),
               to = map_chr(pairs, 2)) %>%
        dplyr::select(from, to)
      
      # Resolve noncircularities for no prior
      # add pair names
      noprior$missing_pair = apply(
        noprior[, c("from", "to")], 
        1, 
        function(x) paste(sort(x), collapse = "-")
      )
      
      # get edge lists for dyads method
      dyad_inferred_edges = nopriorfiltered[,1:2]
      
      # make igraph object
      graph = graph_from_data_frame(dyad_inferred_edges, directed = FALSE)
      components <- decompose(graph)
      
      # find missing links
      missing_links = lapply(components, function(comp) {
        if (!is_clique(comp)){
          # get all possible pairs of nodes (unordered)
          nodes = V(comp)$name
          node_pairs = t(combn(nodes, 2))
          
          # make a list of node pairs which are no adjacent in the component
          missing <- apply(node_pairs, 1, function(pair) {
            !are_adjacent(comp, pair[1], pair[2])
          })
          
          # make collapsed pair names for each missing link
          missing_names = apply(node_pairs[missing, , drop = FALSE], 1, function(pair) {
            paste(sort(pair), collapse = "-")
          })
          
          df = data.frame(
            missing_pair = missing_names,
            stringsAsFactors = FALSE
          )
          return(df)
        }
      })
      missing_df = do.call(rbind, missing_links)
      
      # Maintain sib pairs that are highly likely or which are missing with P > 0.95
      dyadsfiltered = noprior %>% filter(Probability >= prob_thresh |
                                         (missing_pair %in% missing_df$missing_pair & Probability >= 0.95))
      # Get new edge list
      dyad_inferred_edges = dyadsfiltered[,1:2]
      
      
      # Resolve remaining noncircularity
      # Make new graph
      graph = graph_from_data_frame(dyad_inferred_edges, directed = FALSE)
      components <- decompose(graph)
      
      # For each component, keep only the largest clique
      filtered_components <- lapply(components, function(comp) {
        if (is_clique(comp)) {
          return(comp)
        } else {
          cliques <- largest_cliques(comp)
          # randomly choose one if there's a tie
          chosen_clique <- sample(cliques, 1)[[1]]
          # induce subgraph on that clique
          return(induced_subgraph(comp, chosen_clique))
        }
      })
      
      # Recombine the components into a single graph
      final_graph = do.call(disjoint_union, filtered_components)
      
      # get final dyad edges
      inferred_noprior <- as.data.frame(as_edgelist(final_graph))
      colnames(inferred_noprior) = c("from", "to")
      
      
      # Resolve noncircularities for with prior
      # add pair names
      prior$missing_pair = apply(
        prior[, c("from", "to")], 
        1, 
        function(x) paste(sort(x), collapse = "-")
      )
      
      # get edge lists for dyads method
      dyad_inferred_edges = priorfiltered[,1:2]
      
      # make igraph object
      graph = graph_from_data_frame(dyad_inferred_edges, directed = FALSE)
      components <- decompose(graph)
      
      # find missing links
      missing_links = lapply(components, function(comp) {
        if (!is_clique(comp)){
          # get all possible pairs of nodes (unordered)
          nodes = V(comp)$name
          node_pairs = t(combn(nodes, 2))
          
          # make a list of node pairs which are no adjacent in the component
          missing <- apply(node_pairs, 1, function(pair) {
            !are_adjacent(comp, pair[1], pair[2])
          })
          
          # make collapsed pair names for each missing link
          missing_names = apply(node_pairs[missing, , drop = FALSE], 1, function(pair) {
            paste(sort(pair), collapse = "-")
          })
          
          df = data.frame(
            missing_pair = missing_names,
            stringsAsFactors = FALSE
          )
          return(df)
        }
      })
      missing_df = do.call(rbind, missing_links)
      
      # Maintain sib pairs that are highly likely or which are missing with P > 0.95
      dyadsfiltered = prior %>% filter(Probability >= prob_thresh |
                                         (missing_pair %in% missing_df$missing_pair & Probability >= 0.95))
      # Get new edge list
      dyad_inferred_edges = dyadsfiltered[,1:2]
      
      
      # Resolve remaining noncircularity
      # Make new graph
      graph = graph_from_data_frame(dyad_inferred_edges, directed = FALSE)
      components <- decompose(graph)
      
      # For each component, keep only the largest clique
      filtered_components <- lapply(components, function(comp) {
        if (is_clique(comp)) {
          return(comp)
        } else {
          cliques <- largest_cliques(comp)
          # randomly choose one if there's a tie
          chosen_clique <- sample(cliques, 1)[[1]]
          # induce subgraph on that clique
          return(induced_subgraph(comp, chosen_clique))
        }
      })
      
      # Recombine the components into a single graph
      final_graph = do.call(disjoint_union, filtered_components)
      
      # get final dyad edges
      inferred_prior <- as.data.frame(as_edgelist(final_graph))
      colnames(inferred_prior) = c("from", "to")
      
      # combine and classify edges by type (FN, FP, TP)
      
      ###################################
      # NO PRIOR
      # get true positives
      tp_edges = inner_join(true_edges, inferred_noprior, by = c("from", "to"))
      tp_edges$type = "TP"
      
      # initialize other types as FN and FP
      true_edges$type = "FN"
      inferred_noprior$type = "FP"
      
      # remove true positives from FN and FP dataframes
      fn_edges = anti_join(true_edges, tp_edges, by = c("from", "to"))
      fp_edges = anti_join(inferred_noprior, tp_edges, by = c("from", "to"))
      
      # combine all edges
      all_edges = rbind(tp_edges, fn_edges, fp_edges)
      
      # record FPR and FNR
      errors_dyad$numbees[count] = 2000*subsets[j]
      errors_dyad$species[count] = species
      errors_dyad$FPR[count] = nrow(fp_edges) / (nrow(tp_edges) + nrow(fp_edges))
      errors_dyad$FNR[count] = nrow(fn_edges) / (nrow(tp_edges) + nrow(fn_edges))
      errors_dyad$total_real[count] = nrow(tp_edges) + nrow(fn_edges)
      errors_dyad$numFP[count] = nrow(fp_edges)
      errors_dyad$numFN[count] = nrow(fn_edges)
      errors_dyad$numTP[count] = nrow(tp_edges)
      errors_dyad$sibprior[count] = "no_prior"
      
      ###################################
      # PRIOR
      # get true positives
      tp_edges = inner_join(true_edges, inferred_prior, by = c("from", "to"))
      tp_edges$type = "TP"
      
      # initialize other types as FN and FP
      true_edges$type = "FN"
      inferred_prior$type = "FP"
      
      # remove true positives from FN and FP dataframes
      fn_edges = anti_join(true_edges, tp_edges, by = c("from", "to"))
      fp_edges = anti_join(inferred_prior, tp_edges, by = c("from", "to"))
      
      # combine all edges
      all_edges = rbind(tp_edges, fn_edges, fp_edges)
      
      # record FPR and FNR
      errors_dyad$numbees[count + 1] = 2000*subsets[j]
      errors_dyad$species[count + 1] = species
      errors_dyad$FPR[count + 1] = nrow(fp_edges) / (nrow(tp_edges) + nrow(fp_edges))
      errors_dyad$FNR[count + 1] = nrow(fn_edges) / (nrow(tp_edges) + nrow(fn_edges))
      errors_dyad$total_real[count + 1] = nrow(tp_edges) + nrow(fn_edges)
      errors_dyad$numFP[count + 1] = nrow(fp_edges)
      errors_dyad$numFN[count + 1] = nrow(fn_edges)
      errors_dyad$numTP[count + 1] = nrow(tp_edges)
      errors_dyad$sibprior[count + 1] = "prior"
      
      count = count +2
    }
  }
}  


###### Make some plots!

mixtus_FPR = ggplot(errors_dyad[errors_dyad$species=="mixtus",], aes(x = numbees, y = FPR, colour = sibprior, alpha = 0.4)) +
  geom_point() +
  xlab("") +
  ylab(expression(FPR == frac(FP, TP + FP))) +
  labs(colour = "Sibship Prior") +
  scale_color_manual(
    labels = c("No prior", "Prior"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.15))

mixtus_FNR = ggplot(errors_dyad[errors_dyad$species=="mixtus",], aes(x = numbees, y = FNR, colour = sibprior, alpha = 0.4)) +
  geom_point() +
  xlab("Number of Observed Bees") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  labs(colour = "Sibship Prior") +
  scale_color_manual(
    labels = c("No prior", "Prior"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.1))

impatiens_FPR = ggplot(errors_dyad[errors_dyad$species=="impatiens",], aes(x = numbees, y = FPR, colour = sibprior, alpha = 0.4)) +
  geom_point() +
  xlab("") +
  ylab("") +
  labs(colour = "Sibship Prior") +
  scale_color_manual(
    labels = c("No prior", "Prior"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.15))

impatiens_FNR = ggplot(errors_dyad[errors_dyad$species=="impatiens",], aes(x = numbees, y = FNR, colour = sibprior, alpha = 0.4)) +
  geom_point() +
  xlab("Number of Observed Bees") +
  ylab("") +
  labs(colour = "Sibship Prior") +
  scale_color_manual(
    labels = c("No prior", "Prior"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal()  +
  ylim(c(0,0.1))


# Make some grid plots for appendix
#get legend
mixtus_FPR_noalpha = mixtus_FPR + guides(alpha = "none")

g = ggplotGrob(mixtus_FPR_noalpha)
legend_index = which(g$layout$name == "guide-box-right")
legend = g$grobs[[legend_index]]

# remove legend from plots
mixtus_FPR = mixtus_FPR + theme(legend.position = "none")
mixtus_FNR = mixtus_FNR + theme(legend.position = "none")
impatiens_FPR = impatiens_FPR + theme(legend.position = "none")
impatiens_FNR = impatiens_FNR + theme(legend.position = "none")

# make some text grobs
imp = textGrob(
  expression(italic("Bombus impatiens")),
  gp = gpar(fontsize = 12, col = "black")
)
mix = textGrob(
  expression(italic("Bombus mixtus")),
  gp = gpar(fontsize = 12, col = "black")
)
fpr = textGrob(
  "False Positives",
  gp = gpar(fontsize = 12, col = "black"), rot = 270
)
fnr = textGrob(
  "False Negatives",
  gp = gpar(fontsize = 12, col = "black"), rot = 270
)


# arrange and plot
grid = plot_grid(
  mix,        nullGrob(),      imp,        nullGrob(),
  mixtus_FPR, nullGrob(), impatiens_FPR, fpr,
  mixtus_FNR, nullGrob(), impatiens_FNR, fnr,
  ncol = 4,
  rel_widths  = c(8, 1, 8, 1),
  rel_heights = c(1, 8, 8),
  align = "hv", axis = "tblr"
)
grid = grid.arrange(grid, legend, ncol = 2, widths = c(8,2))
grid = ggdraw() +
  draw_plot(grid, 0.05, 0, 0.9, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)", "(D)"), size = 14, 
                  x = c(0, 0.4, 0, 0.4), 
                  y = c(1, 1, 0.53, 0.53))
ggsave("docs/appendix_figures/sibprior_dyads.jpg", grid, height = 1000, width = 3000, units = "px")

### Check error rates with families...
errors_family = data.frame(count = 1:100,
                         numbees = NA,
                         species = NA,
                         numFP = NA,
                         numFN = NA,
                         numTP = NA,
                         total_real = NA,
                         FPR = NA,
                         FNR = NA,
                         sibprior = NA)
count = 1
for (i in 1:nsims){
  for (j in 1:length(subsets)){
    for (species in c("mixtus", "impatiens")){
      name = paste0(species, "_set", i, "_sub", subsets[j], "_exclusion")
      true_data = read.csv(paste0("simulate_data/colony_assignments/sim_data/true_data/", species, "_set", i, "_sub", subsets[j], ".csv"))
      
      # load in no prior data
      nopriorfile = paste0("simulate_data/colony_assignments/test_sibprior/", name, "_nosibprior.BestCluster")
      noprior = as.data.frame(do.call(rbind, strsplit(trimws(readLines(nopriorfile)), "\\s+")[-1]))
      colnames(noprior) = unlist(strsplit(readLines(nopriorfile), "\\s+")[1])
      
      # load in with prior data
      priorfile = paste0("simulate_data/colony_assignments/test_sample_size/colony_output/", name, ".BestCluster")
      prior = as.data.frame(do.call(rbind, strsplit(trimws(readLines(priorfile)), "\\s+")[-1]))
      colnames(prior) = unlist(strsplit(readLines(priorfile), "\\s+")[1])
      
      
      # set probability threshold
      prob_thresh = 0.995
      
      # filter colony outputs
      noprior = noprior %>% filter(as.numeric(Probability) >= prob_thresh)
      prior = prior %>% filter(as.numeric(Probability) >= prob_thresh)
      
      # make edge lists
      true_edges = true_data %>%
        group_by(truecolony) %>%
        filter(n() > 1) %>%
        summarise(pairs = combn(individual, 2, simplify = FALSE), .groups = "drop") %>%
        mutate(from = map_chr(pairs, 1),
               to = map_chr(pairs, 2)) %>%
        dplyr::select(from, to)
      
      # get edge lists for family method
      inferred_noprior = noprior %>%
        group_by(ClusterIndex) %>%
        filter(n() > 1) %>%
        summarise(pairs = combn(OffspringID, 2, simplify = FALSE), .groups = "drop") %>%
        mutate(from = map_chr(pairs, 1),
               to = map_chr(pairs, 2)) %>%
        dplyr::select(from, to)
      
      inferred_prior = prior %>%
        group_by(ClusterIndex) %>%
        filter(n() > 1) %>%
        summarise(pairs = combn(OffspringID, 2, simplify = FALSE), .groups = "drop") %>%
        mutate(from = map_chr(pairs, 1),
               to = map_chr(pairs, 2)) %>%
        dplyr::select(from, to)
      
      # combine and classify edges by type (FN, FP, TP)
      
      ###################################
      # NO PRIOR
      # get true positives
      tp_edges = inner_join(true_edges, inferred_noprior, by = c("from", "to"))
      tp_edges$type = "TP"
      
      # initialize other types as FN and FP
      true_edges$type = "FN"
      inferred_noprior$type = "FP"
      
      # remove true positives from FN and FP dataframes
      fn_edges = anti_join(true_edges, tp_edges, by = c("from", "to"))
      fp_edges = anti_join(inferred_noprior, tp_edges, by = c("from", "to"))
      
      # combine all edges
      all_edges = rbind(tp_edges, fn_edges, fp_edges)
      
      # record FPR and FNR
      errors_family$numbees[count] = 2000*subsets[j]
      errors_family$species[count] = species
      errors_family$FPR[count] = nrow(fp_edges) / (nrow(tp_edges) + nrow(fp_edges))
      errors_family$FNR[count] = nrow(fn_edges) / (nrow(tp_edges) + nrow(fn_edges))
      errors_family$total_real[count] = nrow(tp_edges) + nrow(fn_edges)
      errors_family$numFP[count] = nrow(fp_edges)
      errors_family$numFN[count] = nrow(fn_edges)
      errors_family$numTP[count] = nrow(tp_edges)
      errors_family$sibprior[count] = "no_prior"
      
      ###################################
      # PRIOR
      # get true positives
      tp_edges = inner_join(true_edges, inferred_prior, by = c("from", "to"))
      tp_edges$type = "TP"
      
      # initialize other types as FN and FP
      true_edges$type = "FN"
      inferred_prior$type = "FP"
      
      # remove true positives from FN and FP dataframes
      fn_edges = anti_join(true_edges, tp_edges, by = c("from", "to"))
      fp_edges = anti_join(inferred_prior, tp_edges, by = c("from", "to"))
      
      # combine all edges
      all_edges = rbind(tp_edges, fn_edges, fp_edges)
      
      # record FPR and FNR
      errors_family$numbees[count + 1] = 2000*subsets[j]
      errors_family$species[count + 1] = species
      errors_family$FPR[count + 1] = nrow(fp_edges) / (nrow(tp_edges) + nrow(fp_edges))
      errors_family$FNR[count + 1] = nrow(fn_edges) / (nrow(tp_edges) + nrow(fn_edges))
      errors_family$total_real[count + 1] = nrow(tp_edges) + nrow(fn_edges)
      errors_family$numFP[count + 1] = nrow(fp_edges)
      errors_family$numFN[count + 1] = nrow(fn_edges)
      errors_family$numTP[count + 1] = nrow(tp_edges)
      errors_family$sibprior[count + 1] = "prior"
      
      count = count +2
    }
  }
}  


###### Make some plots!

mixtus_FPR = ggplot(errors_family[errors_family$species=="mixtus",], aes(x = numbees, y = FPR, colour = sibprior, alpha = 0.4)) +
  geom_point() +
  xlab("") +
  ylab(expression(FPR == frac(FP, TP + FP))) +
  labs(colour = "Sibship Prior") +
  scale_color_manual(
    labels = c("No prior", "Prior"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.6))

mixtus_FNR = ggplot(errors_family[errors_family$species=="mixtus",], aes(x = numbees, y = FNR, colour = sibprior, alpha = 0.4)) +
  geom_point() +
  xlab("Number of Observed Bees") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  labs(colour = "Sibship Prior") +
  scale_color_manual(
    labels = c("No prior", "Prior"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.12))

impatiens_FPR = ggplot(errors_family[errors_family$species=="impatiens",], aes(x = numbees, y = FPR, colour = sibprior, alpha = 0.4)) +
  geom_point() +
  xlab("") +
  ylab("") +
  labs(colour = "Sibship Prior") +
  scale_color_manual(
    labels = c("No prior", "Prior"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.6))

impatiens_FNR = ggplot(errors_family[errors_family$species=="impatiens",], aes(x = numbees, y = FNR, colour = sibprior, alpha = 0.4)) +
  geom_point() +
  xlab("Number of Observed Bees") +
  ylab("") +
  labs(colour = "Sibship Prior") +
  scale_color_manual(
    labels = c("No prior", "Prior"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() + 
  ylim(c(0,0.12))


# Make some grid plots for appendix
#get legend
mixtus_FPR_noalpha = mixtus_FPR + guides(alpha = "none")

g = ggplotGrob(mixtus_FPR_noalpha)
legend_index = which(g$layout$name == "guide-box-right")
legend = g$grobs[[legend_index]]

# remove legend from plots
mixtus_FPR = mixtus_FPR + theme(legend.position = "none")
mixtus_FNR = mixtus_FNR + theme(legend.position = "none")
impatiens_FPR = impatiens_FPR + theme(legend.position = "none")
impatiens_FNR = impatiens_FNR + theme(legend.position = "none")

# make some text grobs
imp = textGrob(
  expression(italic("Bombus impatiens")),
  gp = gpar(fontsize = 12, col = "black")
)
mix = textGrob(
  expression(italic("Bombus mixtus")),
  gp = gpar(fontsize = 12, col = "black")
)
fpr = textGrob(
  "False Positives",
  gp = gpar(fontsize = 12, col = "black"), rot = 270
)
fnr = textGrob(
  "False Negatives",
  gp = gpar(fontsize = 12, col = "black"), rot = 270
)


# arrange and plot
grid = plot_grid(
  mix,        nullGrob(),      imp,        nullGrob(),
  mixtus_FPR, nullGrob(), impatiens_FPR, fpr,
  mixtus_FNR, nullGrob(), impatiens_FNR, fnr,
  ncol = 4,
  rel_widths  = c(8, 1, 8, 1),
  rel_heights = c(1, 8, 8),
  align = "hv", axis = "tblr"
)
grid = grid.arrange(grid, legend, ncol = 2, widths = c(8,2))
grid = ggdraw() +
  draw_plot(grid, 0.05, 0, 0.9, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)", "(D)"), size = 14, 
                  x = c(0, 0.4, 0, 0.4), 
                  y = c(1, 1, 0.53, 0.53))
ggsave("docs/appendix_figures/sibprior_families.jpg", grid, height = 1000, width = 3000, units = "px")

################################################################################
## Test multiple paternity
################################################################################

# Create a set of simulations to perform
paternity_probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)
impGenotypesList = list()
mixGenotypesList = list()

# Simulate siblingships
result = draw_simple_multi_landscape(sample_size = 1200,
                                     num_landscape = 6,
                                     landscape_size = 1500,
                                     trapgrid_size = 300,
                                     number_traps = 25,
                                     number_colonies = 12000,
                                     colony_sizes = rep(100,12000),
                                     rho = 100,
                                     distance_decay = "exponential")
sibship_data = list(result[[1]], result[[1]][rowSums(result[[1]])>0,], result[[2]], result[[2]][rowSums(result[[1]])>0,], result[[3]])
saveRDS(sibship_data, paste0("simulate_data/colony_assignments/test_effective_paternity/sim1.RDS"))

# Simulate genotypes
for (i in 1:length(paternity_probs)){
  # simulate impatiens genotypes
  impGenotypesList[[i]] = simulateGenotypes(alleleFreqs = impatiens_allelefreq,
                                            colonyDataDetected = sibship_data[[4]],
                                            observationMatrix = sibship_data[[2]],
                                            trapData = sibship_data[[5]],
                                            probMultiplePaternity = paternity_probs[i])

  # write to files
  write.table(impGenotypesList[[i]][,!colnames(impGenotypesList[[i]]) %in% c("truecolony", "landscape_id", "trap_id")], 
              file = paste0("simulate_data/colony_assignments/test_effective_paternity/for_colony/impatiens_pp", paternity_probs[i], ".txt"), 
              sep= ",", col.names = FALSE, row.names = FALSE)
  write.csv(impGenotypesList[[i]], 
            file = paste0("simulate_data/colony_assignments/test_effective_paternity/true_data/impatiens_pp", paternity_probs[i], ".csv"))
  
  # simulate mixtus genotypes
  mixGenotypesList[[i]] = simulateGenotypes(alleleFreqs = mixtus_allelefreq,
                                            colonyDataDetected = sibship_data[[4]],
                                            observationMatrix = sibship_data[[2]],
                                            trapData = sibship_data[[5]],
                                            probMultiplePaternity = paternity_probs[i])

  #write to files
  write.table(mixGenotypesList[[i]][,!colnames(mixGenotypesList[[i]]) %in% c("truecolony", "landscape_id", "trap_id")], 
              file = paste0("simulate_data/colony_assignments/test_effective_paternity/for_colony/mixtus_pp", paternity_probs[i], ".txt"), 
              sep= ",", col.names = FALSE, row.names = FALSE)
  write.csv(mixGenotypesList[[i]], 
            file = paste0("simulate_data/colony_assignments/test_effective_paternity/true_data/mixtus_pp", paternity_probs[i], ".csv"))
}

all_data = list(result[[1]], result[[1]][rowSums(result[[1]])>0,], result[[2]], result[[2]][rowSums(result[[1]])>0,], result[[3]], impGenotypesList, mixGenotypesList)
saveRDS(all_data, "simulate_data/colony_assignments/test_effective_paternity/firstsim.RDS")
# all_data = readRDS("simulate_data/colony_assignments/test_effective_paternity/firstsim.RDS")
# yobs = all_data[[1]]
# yobs_detected = all_data[[2]]
# colony_data = all_data[[3]]
# colony_data_detected = all_data[[4]]
# trap_data = all_data[[5]]
# impGenotypesList = all_data[[6]]
# mixGenotypesList = all_data[[7]]

# Construct error rates files -- all zeros for now
imp_columns = unique(impatiens_allelefreq$MarkerID)
impatiens_error_rates = as.data.frame(matrix(0, nrow = 4, ncol = length(imp_columns)))
impatiens_error_rates[1,] = imp_columns
write.table(impatiens_error_rates, "simulate_data/colony_assignments/test_effective_paternity/for_colony/impatiens_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

mix_columns = unique(mixtus_allelefreq$MarkerID)
mixtus_error_rates = as.data.frame(matrix(0, nrow = 4, ncol = length(mix_columns)))
mixtus_error_rates[1,] = mix_columns
write.table(mixtus_error_rates, "simulate_data/colony_assignments/test_effective_paternity/for_colony/mixtus_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Make and write exclusion tables
# Only need one since all the siblingships are the same (just with different genotypes)
exclusion_table = createExclusionTable(impGenotypesList[[1]])
write.table(
  exclusion_table,
  file = paste0("simulate_data/colony_assignments/test_effective_paternity/for_colony/exclusion_table.txt"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  na = ""
)


# Construct .DAT files for COLONY
mixtus_errors_filepath = "simulate_data/colony_assignments/test_effective_paternity/for_colony/mixtus_error_rates.txt"
impatiens_errors_filepath = "simulate_data/colony_assignments/test_effective_paternity/for_colony/impatiens_error_rates.txt"

for (i in paternity_probs){
  # get sample size, working directory
  size = 1200
  workingdir = paste0(getwd(),"/simulate_data/colony_assignments/test_effective_paternity/colony_output")
  
  # get genotype filepaths
  mixtus_genotypes_filepath = paste0("simulate_data/colony_assignments/test_effective_paternity/for_colony/mixtus_pp", i, ".txt")
  impatiens_genotypes_filepath = paste0("simulate_data/colony_assignments/test_effective_paternity/for_colony/impatiens_pp", i, ".txt")
  
  # get exclusion path
  exclusion_filepath = paste0("simulate_data/colony_assignments/test_effective_paternity/for_colony/exclusion_table.txt")
  
  #build .DAT for mixtus
  # monogamous settings
  rcolony::build.colony.superauto(wd=workingdir, 
                                  name=paste0(getwd(), "/simulate_data/colony_assignments/test_effective_paternity/colony_output/mixtus_pp", i, ".DAT"), 
                                  datasetname = paste0("mixtus_pp", i),
                                  delim=",",
                                  sample_size = size,
                                  sibship_prior = 1,
                                  female_monogamy = 1,
                                  num_loci = 10,
                                  error_rates_path = mixtus_errors_filepath,
                                  genotypes_path = mixtus_genotypes_filepath,
                                  exclusion_path = exclusion_filepath
  )
  
  # polygamous settings
  rcolony::build.colony.superauto(wd=workingdir, 
                                  name=paste0(getwd(),"/simulate_data/colony_assignments/test_effective_paternity/colony_output/mixtus_pp", i, "_poly.DAT"), 
                                  datasetname = paste0("mixtus_pp", i, "_poly"),
                                  delim=",",
                                  sample_size = size,
                                  sibship_prior = 1,
                                  female_monogamy = 0,
                                  num_loci = 10,
                                  error_rates_path = mixtus_errors_filepath,
                                  genotypes_path = mixtus_genotypes_filepath,
                                  exclusion_path = exclusion_filepath
  )
  
  # build .DAT for impatiens
  # monogamous settings
  rcolony::build.colony.superauto(wd=workingdir, 
                                  name=paste0(getwd(),"/simulate_data/colony_assignments/test_effective_paternity/colony_output/impatiens_pp", i, ".DAT"), 
                                  datasetname = paste0("impatiens_pp", i),
                                  delim=",",
                                  sample_size = size,
                                  sibship_prior = 1,
                                  female_monogamy = 1,
                                  num_loci = 12,
                                  error_rates_path = impatiens_errors_filepath,
                                  genotypes_path = impatiens_genotypes_filepath,
                                  exclusion_path = exclusion_filepath
  )
  
  # polygamous settings
  rcolony::build.colony.superauto(wd=workingdir, 
                                  name=paste0(getwd(),"/simulate_data/colony_assignments/test_effective_paternity/colony_output/impatiens_pp", i, "_poly.DAT"), 
                                  datasetname = paste0("impatiens_pp", i, "_poly"),
                                  delim=",",
                                  sample_size = size,
                                  sibship_prior = 1,
                                  female_monogamy = 0,
                                  num_loci = 12,
                                  error_rates_path = impatiens_errors_filepath,
                                  genotypes_path = impatiens_genotypes_filepath,
                                  exclusion_path = exclusion_filepath
  )
}

# Load in results from COLONY
errors_family = data.frame(count = 1:24,
                         numbees = NA,
                         species = NA,
                         numFP = NA,
                         numFN = NA,
                         numTP = NA,
                         total_real = NA,
                         FPR = NA,
                         FNR = NA,
                         multiplepaternityrate = NA,
                         modelspec = NA)
count = 1
for (i in paternity_probs){
  for (j in c("_poly", "")){
    for (species in c("mixtus", "impatiens")){
      name = paste0(species, "_pp", i)
      true_data = read.csv(paste0("simulate_data/colony_assignments/test_effective_paternity/true_data/", name, ".csv"))
      
      # load in sibship data
      filename = paste0("simulate_data/colony_assignments/test_effective_paternity/colony_output/", name, j, ".BestCluster")
      colony_output = as.data.frame(do.call(rbind, strsplit(trimws(readLines(filename)), "\\s+")[-1]))
      colnames(colony_output) = unlist(strsplit(readLines(filename), "\\s+")[1])

      # set probability threshold
      prob_thresh = 1
      
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
      
      # record FPR and FNR
      errors_family$numbees[count] = 1200
      errors_family$species[count] = species
      errors_family$FPR[count] = nrow(fp_edges) / (nrow(tp_edges) + nrow(fp_edges))
      errors_family$FNR[count] = nrow(fn_edges) / (nrow(tp_edges) + nrow(fn_edges))
      errors_family$total_real[count] = nrow(tp_edges) + nrow(fn_edges)
      errors_family$numFP[count] = nrow(fp_edges)
      errors_family$numFN[count] = nrow(fn_edges)
      errors_family$numTP[count] = nrow(tp_edges)
      errors_family$multiplepaternityrate[count] = i
      errors_family$modelspec[count] = j
      
      count = count + 1
    }
  }
}  

errors_family$modelspec[errors_family$modelspec != "_poly"] = "_mono"
errors_family$test_condition = paste0(errors_family$multiplepaternityrate, errors_family$modelspec)

# make some plots
mixtus_paternity_FPR = ggplot(errors_family[errors_family$species=="mixtus",], aes(x = multiplepaternityrate, y = FPR, colour = modelspec)) +
  geom_point() +
  xlab("") +
  ylab(expression(FPR == frac(FP, TP + FP))) +
  labs(colour = "COLONY setting") +
  scale_color_manual(
    labels = c("Queen monogamy", "Queen polygamy"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.8))

mixtus_paternity_FNR = ggplot(errors_family[errors_family$species=="mixtus",], aes(x = multiplepaternityrate, y = FNR, colour = modelspec)) +
  geom_point() +
  xlab("Percentage of colonies with two fathers") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  labs(colour = "COLONY setting") +
  scale_color_manual(
    labels = c("Queen monogamy", "Queen polygamy"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal()  +
  ylim(c(0,0.5))

impatiens_paternity_FPR = ggplot(errors_family[errors_family$species=="impatiens",], aes(x = multiplepaternityrate, y = FPR, colour = modelspec)) +
  geom_point() +
  xlab("") +
  ylab("") +
  labs(colour = "COLONY setting") +
  scale_color_manual(
    labels = c("Queen monogamy", "Queen polygamy"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.8))

impatiens_paternity_FNR = ggplot(errors_family[errors_family$species=="impatiens",], aes(x = multiplepaternityrate, y = FNR, colour = modelspec)) +
  geom_point() +
  xlab("Percentage of colonies with two fathers") +
  ylab("") +
  labs(colour = "COLONY setting") +
  scale_color_manual(
    labels = c("Queen monogamy", "Queen polygamy"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal()  +
  ylim(c(0,0.5))


# Make some grid plots for appendix
#get legend
g = ggplotGrob(mixtus_paternity_FPR)
legend_index = which(g$layout$name == "guide-box-right")
legend = g$grobs[[legend_index]]

# remove legend from plots
mixtus_paternity_FPR = mixtus_paternity_FPR + theme(legend.position = "none")
mixtus_paternity_FNR = mixtus_paternity_FNR + theme(legend.position = "none")
impatiens_paternity_FPR = impatiens_paternity_FPR + theme(legend.position = "none")
impatiens_paternity_FNR = impatiens_paternity_FNR + theme(legend.position = "none")

# make some text grobs
imp = textGrob(
  expression(italic("Bombus impatiens")),
  gp = gpar(fontsize = 12, col = "black")
)
mix = textGrob(
  expression(italic("Bombus mixtus")),
  gp = gpar(fontsize = 12, col = "black")
)
fpr = textGrob(
  "False Positives",
  gp = gpar(fontsize = 12, col = "black"), rot = 270
)
fnr = textGrob(
  "False Negatives",
  gp = gpar(fontsize = 12, col = "black"), rot = 270
)


# arrange and plot
grid = plot_grid(
  mix,        nullGrob(),      imp,        nullGrob(),
  mixtus_paternity_FPR, nullGrob(), impatiens_paternity_FPR, fpr,
  mixtus_paternity_FNR, nullGrob(), impatiens_paternity_FNR, fnr,
  ncol = 4,
  rel_widths  = c(8, 1, 8, 1),
  rel_heights = c(1, 8, 8),
  align = "hv", axis = "tblr"
)
grid = grid.arrange(grid, legend, ncol = 2, widths = c(8,2))
grid = ggdraw() +
  draw_plot(grid, 0.05, 0, 0.9, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)", "(D)"), size = 14, 
                  x = c(0, 0.4, 0, 0.4), 
                  y = c(1, 1, 0.53, 0.53))
ggsave("docs/appendix_figures/multiplepaternity.jpg", grid, height = 1000, width = 3000, units = "px")



################################################################################
## Test multiple paternity for higher resolution genetic dataset
################################################################################
# Siblingship size prior greatly improves FPR rate when we run COLONY assuming female monogamy, but 
# if we assume female polygamy we get a whole mess of false sibships
# Check: could we improve this tendency by simulating a fake data set with twice as many loci?

# Load in allele frequency data
impatiens_allelefreq = read.csv("colony_assignments/Colony2_Linux/impatiens2023.AlleleFreq")
mixtus_allelefreq = read.csv("colony_assignments/Colony2_Linux/mixtus2023.AlleleFreq")

# combine frequency data!
impatiens_allelefreq$MarkerID = paste0("impatiens_", impatiens_allelefreq$MarkerID)
mixtus_allelefreq$MarkerID = paste0("mixtus_", mixtus_allelefreq$MarkerID)
combined_allelefreq = rbind(impatiens_allelefreq, mixtus_allelefreq)

# Create a set of simulations to perform
paternity_probs = c(0, 0.2, 0.4, 0.6, 0.8, 1)
GenotypesList = list()

# Simulate genotypes for each species and mating condition
sibship_data = readRDS(paste0("simulate_data/colony_assignments/test_effective_paternity/sim1.RDS"))

# Simulate genotypes
for (i in 1:length(paternity_probs)){
  # simulate impatiens genotypes
  GenotypesList[[i]] = simulateGenotypes(alleleFreqs = combined_allelefreq,
                                            colonyDataDetected = sibship_data[[4]],
                                            observationMatrix = sibship_data[[2]],
                                            trapData = sibship_data[[5]],
                                            probMultiplePaternity = paternity_probs[i])
  
  # write to files
  write.table(GenotypesList[[i]][,!colnames(GenotypesList[[i]]) %in% c("truecolony", "landscape_id", "trap_id")], 
              file = paste0("simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/augmented_pp", paternity_probs[i], ".txt"), 
              sep= ",", col.names = FALSE, row.names = FALSE)
  write.csv(GenotypesList[[i]], 
            file = paste0("simulate_data/colony_assignments/test_effective_paternity/true_data/augmented_pp", paternity_probs[i], ".csv"))
  
}
saveRDS(GenotypesList, "simulate_data/colony_assignments/test_effective_paternity/augmentedsim.RDS")


# Construct error rates files -- all zeros for now
columns = unique(combined_allelefreq$MarkerID)
all_error_rates = as.data.frame(matrix(0, nrow = 4, ncol = length(columns)))
all_error_rates[1,] = columns
write.table(all_error_rates, "simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/augmented_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
errors_filepath = "simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/augmented_error_rates.txt"

# Make and write exclusion tables
# Only need one since all the siblingships are the same (just with different genotypes)
exclusion_table = createExclusionTable(GenotypesList[[1]])
write.table(
  exclusion_table,
  file = paste0("simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/exclusion_table.txt"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  na = ""
)
exclusion_filepath = paste0("simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/exclusion_table.txt")


# Construct .DAT files for COLONY
for (i in paternity_probs){
  # get sample size, working directory
  size = 1200
  workingdir = paste0(getwd(),"/simulate_data/colony_assignments/test_effective_paternity/colony_output_augmented")
  
  # get genotype filepaths
  genotypes_filepath = paste0("simulate_data/colony_assignments/test_effective_paternity/for_colony_augmented/augmented_pp", i, ".txt")
  
  #build .DAT for mixtus
  # monogamous settings
  rcolony::build.colony.superauto(wd=workingdir, 
                                  name=paste0(getwd(), "/simulate_data/colony_assignments/test_effective_paternity/colony_output_augmented/augmented_pp", i, ".DAT"), 
                                  datasetname = paste0("augmented_pp", i),
                                  delim=",",
                                  sample_size = size,
                                  sibship_prior = 1,
                                  female_monogamy = 1,
                                  num_loci = 22,
                                  error_rates_path = errors_filepath,
                                  genotypes_path = genotypes_filepath,
                                  exclusion_path = exclusion_filepath
  )
  
  # polygamous settings
  rcolony::build.colony.superauto(wd=workingdir, 
                                  name=paste0(getwd(),"/simulate_data/colony_assignments/test_effective_paternity/colony_output_augmented/augmented_pp", i, "_poly.DAT"), 
                                  datasetname = paste0("augmented_pp", i, "_poly"),
                                  delim=",",
                                  sample_size = size,
                                  sibship_prior = 1,
                                  female_monogamy = 0,
                                  num_loci = 22,
                                  error_rates_path = errors_filepath,
                                  genotypes_path = genotypes_filepath,
                                  exclusion_path = exclusion_filepath
  )
}

# Load in results from COLONY
errors_family = data.frame(count = 1:12,
                           numbees = NA,
                           numFP = NA,
                           numFN = NA,
                           numTP = NA,
                           total_real = NA,
                           FPR = NA,
                           FNR = NA,
                           multiplepaternityrate = NA,
                           modelspec = NA)
count = 1
for (i in paternity_probs){
  for (j in c("_poly", "")){
    true_data = read.csv(paste0("simulate_data/colony_assignments/test_effective_paternity/true_data/augmented_pp", i, ".csv"))
    
    # load in sibship data
    filename = paste0("simulate_data/colony_assignments/test_effective_paternity/colony_output_augmented/augmented_pp", i, j, ".BestCluster")
    colony_output = as.data.frame(do.call(rbind, strsplit(trimws(readLines(filename)), "\\s+")[-1]))
    colnames(colony_output) = unlist(strsplit(readLines(filename), "\\s+")[1])
    
    # set probability threshold
    prob_thresh = 1
    
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
    
    # record FPR and FNR
    errors_family$numbees[count] = 1200
    errors_family$FPR[count] = nrow(fp_edges) / (nrow(tp_edges) + nrow(fp_edges))
    errors_family$FNR[count] = nrow(fn_edges) / (nrow(tp_edges) + nrow(fn_edges))
    errors_family$total_real[count] = nrow(tp_edges) + nrow(fn_edges)
    errors_family$numFP[count] = nrow(fp_edges)
    errors_family$numFN[count] = nrow(fn_edges)
    errors_family$numTP[count] = nrow(tp_edges)
    errors_family$multiplepaternityrate[count] = i
    errors_family$modelspec[count] = j
    
    count = count + 1
  }
}

# make some plots
paternity_FPR = ggplot(errors_family, aes(x = multiplepaternityrate, y = FPR, colour = modelspec)) +
  geom_point() +
  xlab("") +
  ylab(expression(FPR == frac(FP, TP + FP))) +
  labs(colour = "COLONY setting") +
  scale_color_manual(
    labels = c("Queen monogamy", "Queen polygamy"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal() +
  ylim(c(0,0.8))

paternity_FNR = ggplot(errors_family, aes(x = multiplepaternityrate, y = FNR, colour = modelspec)) +
  geom_point() +
  xlab("Percentage of colonies with two fathers") +
  ylab(expression(FNR == frac(FN, TP + FN))) +
  labs(colour = "COLONY setting") +
  scale_color_manual(
    labels = c("Queen monogamy", "Queen polygamy"),
    values = c("darkslateblue", "salmon")) +
  theme_minimal()  +
  ylim(c(0,0.5))


# Make some grid plots for appendix
#get legend
g = ggplotGrob(paternity_FPR)
legend_index = which(g$layout$name == "guide-box-right")
legend = g$grobs[[legend_index]]

# remove legend from plots
paternity_FPR = paternity_FPR + theme(legend.position = "none")
paternity_FNR = paternity_FNR + theme(legend.position = "none")


# arrange and plot
grid = plot_grid(paternity_FPR, paternity_FNR,
  ncol = 1
)
grid = grid.arrange(grid, legend, ncol = 2, widths = c(8,2))
grid = ggdraw() +
  draw_plot(grid, 0.05, 0, 0.9, 1) +
  draw_plot_label(c("(A)", "(B)"), size = 14, 
                  x = c(0, 0), 
                  y = c(1, 0.53))
ggsave("docs/appendix_figures/augmentedpaternity.jpg", grid, height = 1000, width = 3000, units = "px")
