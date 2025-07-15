########
## Get colony IDs from microsatellite genotyping results--Bombus impatiens
## Started by J Melanson
## February 13, 2025
########

# prep workspace
rm(list = ls())
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")

# first, load in packages
source('colony_assignments/src/init.R')
source('colony_assignments/src/joinFunctions.R')
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

# next; load in specimen data and allele tables
#load field data
specimenData2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2022specimendata.csv", sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023specimendata.csv", sep = ",", header = T, fill = TRUE))
sampleEffort2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023sampledata.csv", sep = ","))
samplePoints = as.data.frame(read.table("/Users/jenna1/Documents/UBC/bombus_project/raw_data/allsamplepoints.csv", sep = ","))
colnames(samplePoints) = c("sample_pt", "gps","landowner","subsite")

# load allele tables
alleles_plex1_all = as.data.frame(read.table("data/from_geneious/impatiens_plex1.csv", sep = ",", header = T))
alleles_plex2_all = as.data.frame(read.table("data/from_geneious/impatiens_plex2.csv", sep = ",", header = T))

# create barcode id for allele tables
alleles_plex1_all$barcode_id = sapply(alleles_plex1_all$Name, 
                                      function(x){paste(unlist(strsplit(x, "_"))[2:4], collapse = "_")})
alleles_plex1_all = alleles_plex1_all[,!colnames(alleles_plex1_all) %in% c("Name")]
alleles_plex2_all$barcode_id = sapply(alleles_plex2_all$Name, 
                                      function(x){paste(unlist(strsplit(x, "_"))[2:4], collapse = "_")})
alleles_plex2_all = alleles_plex2_all[,!colnames(alleles_plex2_all) %in% c("Name")]

############################
# Cleaning
############################
# at some point during processing, impatiensQueens02 and impatiensDNA14 (all on
# the same plate) were SWAPPED with impatiensQueens01 -- so I'll start by fixing that!

importantRows = c("barcode_id", "plate", "well", "final_id")
all_samples = rbind(specimenData2022[,colnames(specimenData2022) %in% importantRows], 
                    specimenData2023[, colnames(specimenData2023) %in% importantRows])
reversed_plex2 = left_join(all_samples, alleles_plex2_all, by = "barcode_id") %>%
  filter(plate == "impatiensQueens01" | plate == "impatiensQueens02" | plate == "impatiens_DNA_14")
# a few samples missing from here -- 185 rather than 192 -- likely because a few
# queens on that plate were actually MIXTUS, and some others failed genotyping

# make a column for the correct plate of each sample
reversed_plex2$correct_plate = NULL
combos = expand.grid(c("A", "B", "C", "D", "E", "F", "G", "H"),
                     c("1", "2", "3"))
plate14_wells = paste0(combos$Var1, combos$Var2)

for(i in 1:nrow(reversed_plex2)){
  if(reversed_plex2$plate[i] == "impatiensQueens02" | reversed_plex2$plate[i] == "impatiens_DNA_14"){
    reversed_plex2$correct_plate[i] = "impatiensQueens01"
  } else if(reversed_plex2$plate[i] == "impatiensQueens01" & reversed_plex2$well[i] %in% plate14_wells){
        reversed_plex2$correct_plate[i] = "impatiens_DNA_14"
    } else if(reversed_plex2$plate[i] == "impatiensQueens01" & !reversed_plex2$well[i] %in% plate14_wells){
          reversed_plex2$correct_plate[i] = "impatiensQueens02"
      }
  }

# make location field in allele scores and specimen table
reversed_plex2$location = paste(reversed_plex2$correct_plate, reversed_plex2$well, sep = "_")
all_samples$location = paste(all_samples$plate, all_samples$well, sep = "_")

#make lookup df for correct barcodes
all_samples_lookup = all_samples[,colnames(all_samples) %in% c("barcode_id", "location")]
colnames(all_samples_lookup) = c("true_barcode", "location") 
  
# join lookup to reversed samples
reversed_plex2 = left_join(reversed_plex2, all_samples_lookup, by = "location")

#clean up extra columns
to_remove = c("barcode_id", "final_id", "plate", "well", "Name", "correct_plate", "location")
reversed_plex2_clean = reversed_plex2[,!colnames(reversed_plex2) %in% to_remove]
colnames(reversed_plex2_clean)[colnames(reversed_plex2_clean) == "true_barcode"] = "barcode_id"

# reattach the scores with correct barcodes to the plex 2 dataframe
alleles_plex2_okay = alleles_plex2_all[!alleles_plex2_all$barcode_id %in% reversed_plex2$barcode_id,]
alleles_plex2_cleaned = rbind(alleles_plex2_okay, reversed_plex2_clean[,colnames(reversed_plex2_clean) %in% colnames(alleles_plex2_okay)])

#check for duplicates in these dataframes
print("Number of rows in alleles 1 dataframe:")
dim(alleles_plex1_all)[1]
print("Number of unique barcodes in alleles 1 dataframe (should be the same):")
length(unique(alleles_plex1_all$barcode_id))
print("Duplicates in final dataframe:")
duplicated = alleles_plex1_all$barcode_id[duplicated(alleles_plex1_all$barcode_id)]
print(duplicated)

print("Number of rows in alleles 2 dataframe:")
dim(alleles_plex2_cleaned)[1]
print("Number of unique barcodes in alleles 2 dataframe (should be the same):")
length(unique(alleles_plex2_cleaned$barcode_id))
print("Duplicates in final dataframe:")
duplicated = alleles_plex2_cleaned$barcode_id[duplicated(alleles_plex2_cleaned$barcode_id)]
print(duplicated)

#remove duplicated rows
plex1_unique <- alleles_plex1_all[!duplicated(alleles_plex1_all$barcode_id), ]
plex2_unique <- alleles_plex2_cleaned[!duplicated(alleles_plex2_cleaned$barcode_id), ]

#merge plexes
microsat_scores = full_join(plex1_unique, plex2_unique, by = "barcode_id")

#make a specimen subset that include year and notes (to filter 2022 vs 2023, queens, correct species)
colsToKeep = c("barcode_id", "year", "notes", "final_id", "sample_id")
specsubset = rbind(specimenData2022[,colnames(specimenData2022) %in% colsToKeep],
                   specimenData2023[,colnames(specimenData2023) %in% colsToKeep])
specimens_withscores = left_join(microsat_scores, specsubset, by = "barcode_id")

# join to sample effort df to get julian dates (only need for 2023)
specimens_withdates = joinJulianDates(specimens_withscores, sampleEffort2023)

# split by year, including early season 2023 queens with 2022 subset
# julian date cutoff based on script "queen_phenology.R"
impatiens2022 = filter(specimens_withdates, final_id == "B. impatiens") %>%
  filter(year == "2022" | (str_detect(notes, "queen") & julian_date < 175)) %>%
  filter(!str_detect(notes, "male"))
impatiens2023 = filter(specimens_withdates, final_id == "B. impatiens") %>%
  filter(year == "2023") %>%
  filter((!str_detect(notes, "queen")) | julian_date > 175) %>%
  filter(!str_detect(notes, "male"))
                          

#check which bees are unscored
specimens = full_join(microsat_scores, specsubset, by = "barcode_id")
impatiens = filter(specimens, final_id == "B. impatiens" & barcode_id != "N/A")
unscoredimpatiens = impatiens[rowSums(is.na(impatiens)) >14,]

##################################
# Prep genotype data for COLONY
##################################

#remove all columns except barcode and scores
# also remove loci with high error rates
colsToRemove = c("year", 
                 "final_id", 
                 "notes", 
                 "sample_id", 
                 "date", 
                 "julian_date", 
                 "BL15...1", 
                 "BL15...2",
                 "BT28...1", 
                 "BT28...2",
                 "BL13...1", 
                 "BL13...2",
                 "BTMS0073...1", 
                 "BTMS0073...2")
impatiens2022 = impatiens2022[, !colnames(impatiens2022) %in% colsToRemove]
impatiens2023 = impatiens2023[, !colnames(impatiens2023) %in% colsToRemove]

#relocate barcode, remove rows with more than 10 NAs, replace NAs with 0's
impatiens2022 = impatiens2022 %>% relocate(barcode_id)
impatiens2022_forcolony = impatiens2022[rowSums(is.na(impatiens2022)) <= 8,]
impatiens2022_forcolony[is.na(impatiens2022_forcolony)] = 0

impatiens2023 = impatiens2023 %>% relocate(barcode_id)
impatiens2023_forcolony = impatiens2023[rowSums(is.na(impatiens2023)) <= 8,]
impatiens2023_forcolony[is.na(impatiens2023_forcolony)] = 0

#write csvs for upload to colony
column_names = colnames(impatiens2022_forcolony)

write.table(impatiens2022_forcolony, "data/merged_by_year/sib_scores_for_colony/impatiens2022_forcolony.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(impatiens2023_forcolony, "data/merged_by_year/sib_scores_for_colony/impatiens2023_forcolony.txt", sep= ",", col.names = FALSE, row.names = FALSE)


write.csv(impatiens2022_forcolony, "data/merged_by_year/csvs/impatiens_2022_scores.csv")
write.csv(impatiens2023_forcolony, "data/merged_by_year/csvs/impatiens_2023_scores.csv")


###########################################
# Prep microsat error rates for COLONY
###########################################

# first value per column: marker name
# second value per column: marker type (codominant for microsats -> 0)
# third value per column: allelic dropout rate ?? (set to 0)
# fourth value per column: error/mutation rate (empirically derived, see check_microsat_error_rates.R)
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

write.table(impatiens_error_rates, "data/merged_by_year/error_rates/impatiens_error_rates.txt", 
            sep= ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

###########################################
# Prep sibship exclusion data for COLONY
###########################################
#basically -- we want to exclude sibships from different regions, as they are 
#quite far apart and therefore very unlikely to be genuine

# each row is excluded sibships for an individual; 
# column 1 is the focal sibling, column 2 is number of excluded sibs, additional columns are excluded sibs
sites = c("W", "SD", "ED", "NR", "HR", "PM")
excluded_sibships_2022 = list()
excluded_sibships_2023 = list()
sample.names.2022 = impatiens2022_forcolony[,1]
sample.names.2023 = impatiens2023_forcolony[,1]

for (i in 1:6){
  site = sites[i]
  
  # grep names for focal site
  site.names.2022 = sample.names.2022[grep(site, sample.names.2022)]
  site.names.2023 = sample.names.2023[grep(site, sample.names.2023)]
  
  # list of excluded siblings
  excluded.2022 = sample.names.2022[sample.names.2022 %!in% site.names.2022]
  numex_2022 = length(excluded.2022)
  excluded.2023 = sample.names.2023[sample.names.2023 %!in% site.names.2023]
  numex_2023 = length(excluded.2023)
  
  # make a matrix of excluded samples
  excluded.matrix.2022 <- matrix(rep(excluded.2022, times = length(site.names.2022)),
                                 nrow = length(site.names.2022), byrow = TRUE)
  colnames(excluded.matrix.2022) <- paste0("excluded_", seq_len(ncol(excluded.matrix.2022)))
  
  excluded.matrix.2023 <- matrix(rep(excluded.2023, times = length(site.names.2023)),
                                 nrow = length(site.names.2023), byrow = TRUE)
  colnames(excluded.matrix.2023) <- paste0("excluded_", seq_len(ncol(excluded.matrix.2023)))
  
  # create full df
  sitedf2022 = data.frame(focal = site.names.2022,
                          num_exc = numex_2022,
                          excluded.matrix.2022,
                          stringsAsFactors = FALSE
  )
  
  sitedf2023 = data.frame(focal = site.names.2023,
                          num_exc = numex_2023,
                          excluded.matrix.2023,
                          stringsAsFactors = FALSE
  )
  
  # save exclusion dfs to list
  excluded_sibships_2022[[i]] = sitedf2022
  excluded_sibships_2023[[i]] = sitedf2023
}

# make exclusion tables!
sibexclusions_2022 = bind_rows(excluded_sibships_2022)
sibexclusions_2023 = bind_rows(excluded_sibships_2023)

sib2022_reduced = sibexclusions_2022[,!colnames(sibexclusions_2022) %in% c("num_exc")]
sib2023_reduced = sibexclusions_2023[,!colnames(sibexclusions_2023) %in% c("num_exc")]

write.table(
  sib2022_reduced,
  file = "data/merged_by_year/sib_exclusions/impatiens_sibexclusions_2022.txt",
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  na = ""
)

write.table(
  sib2023_reduced,
  file = "data/merged_by_year/sib_exclusions/impatiens_sibexclusions_2023.txt",
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  na = ""
)

########################################
# Generate .DAT file and run colony
########################################
# A few important notes!
# For MacOS users: run colony in Rosetta terminal -- colony2 expects Intel versions of shared libraries (x86_64) not Apple Silicon (ARM 64)

# Rcolony (a wrapper package for creating the .dat input file for COLONY2) is a 
#bit out of date and is missing some important arguments. For this reason I have 
# made some small updates to the package, available at https://github.com/jmelanson98/rcolony
# Forked from the excellent original package at https://github.com/jonesor/rcolony

# Changes include:
# - Modification to the writing of siblingship exclusion table (remove excess padding of 
# space at the end of lines)
# - Addition of several arguments required by the latest version of COLONY2:
# line 8 (after dioecious/monoecious): 0/1 for inbreeding (recommended for dioecious: no inbreeding (0))
# line 11 (after mating systems): 0/1 for clone/duplicate inference (0 = no inference, 1 = yes inference)
# line 12 (after clone inference): sibship size scaling (1=yes, 0=no) --> default yes, but if the maximal full sibship size is small (<20) then a run with alternative (no scaling) is necessary
# - Addition of exclusion threshold (0) for known paternity/maternity in cases wehre the number of known parentages is 0
# If you do multiple runs, make sure to reset the seed for each!! you can do this by editing the text file rather than running
# through the interactive bit every time

# build .DAT files for both datasets
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux/impatiens2022_final1.DAT", delim=",")
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux/impatiens2023_final1.DAT", delim=",")

# navigate to COLONY2 sub folder and run the following in terminal
#NOTE: ensure that mixtus2022.DAT and mixtus2023.DAT are in the same directory at the colony2s.out executable

# cd /Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2
# ./colony2s.out IFN:initial/impatiens2022.DAT
# ./colony2s.out IFN:initial/impatiens.DAT



#####################################
## Load in results from colony
#####################################

list_of_results = list()
for (i in 1:5){
  # read in dyad data from COLONY
  fsDyad2022 = trimws(readLines(paste0("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux/impatiens2022_final", i, ".FullSibDyad")))
  
  #format
  dyads2022 = do.call(rbind, strsplit(fsDyad2022, ","))
  colnames(dyads2022) = dyads2022[1,]
  dyads2022 = as.data.frame(dyads2022[-1,])
  dyads2022$Probability = as.numeric(dyads2022$Probability)
  
  # remove sibpairs below some threshold probability
  bestdyads2022 = dyads2022 %>% filter(Probability ==1)
  
  # add results to list
  list_of_results[[i]] = bestdyads2022
}

######################################
## Retain dyads present in all 5 runs
######################################

# first, collapse pairs into an ordered character string
pairs_df1 <- collapse(list_of_results[[1]])
pairs_df2 <- collapse(list_of_results[[2]])
pairs_df3 <- collapse(list_of_results[[3]])
pairs_df4 <- collapse(list_of_results[[4]])
pairs_df5 <- collapse(list_of_results[[5]])

# check matches of df1 against other dfs
in_df2 <- pairs_df1 %in% pairs_df2
in_df3 <- pairs_df1 %in% pairs_df3
in_df4 <- pairs_df1 %in% pairs_df4
in_df5 <- pairs_df1 %in% pairs_df5

bestdyads2022 = list_of_results[[1]]
bestdyads2022$in2 = in_df2
bestdyads2022$in3 =  in_df3
bestdyads2022$in4 =  in_df4
bestdyads2022$in5 =  in_df5
bestdyads2022$inall = in_df2 & in_df3 & in_df4 & in_df5
bestdyads2022$num_true = 1 + rowSums(bestdyads2022[,colnames(bestdyads2022) %in% c("in2", "in3", "in4", "in5")])

# subset to the consistent dyads
consistentdyad_2022 = bestdyads2022[bestdyads2022$inall ==TRUE,]


###########################################################
### Identify missing links to establish complete cliques
###########################################################

# make igraph object
edges = consistentdyad_2022[, c("OffspringID1", "OffspringID2")]
graph_22 = graph_from_data_frame(d = edges, directed = FALSE)
components <- decompose(graph_22)

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


#########################################################
#### Check if missing pairs are present with > 99% prob
########################################################

list_of_results_99 = list()
for (i in 1:5){
  # read in dyad data from COLONY
  fsDyad2022 = trimws(readLines(paste0("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux/impatiens2022_final", i, ".FullSibDyad")))
  
  #format
  dyads2022 = do.call(rbind, strsplit(fsDyad2022, ","))
  colnames(dyads2022) = dyads2022[1,]
  dyads2022 = as.data.frame(dyads2022[-1,])
  dyads2022$Probability = as.numeric(dyads2022$Probability)
  
  # remove sibpairs below some threshold probability
  bestdyads2022 = dyads2022 %>% filter(Probability > 0.99)
  
  # add results to list
  list_of_results_99[[i]] = bestdyads2022
}

# generate pairs
pairs_df1_99 = collapse(list_of_results_99[[1]])
pairs_df2_99 = collapse(list_of_results_99[[2]])
pairs_df3_99 = collapse(list_of_results_99[[3]])
pairs_df4_99 = collapse(list_of_results_99[[4]])
pairs_df5_99 = collapse(list_of_results_99[[5]])

# check matches of df1 against other dfs
in_df2_99 <- pairs_df1_99 %in% pairs_df2_99
in_df3_99 <- pairs_df1_99 %in% pairs_df3_99
in_df4_99 <- pairs_df1_99 %in% pairs_df4_99
in_df5_99 <- pairs_df1_99 %in% pairs_df5_99

gooddyads2022 = list_of_results_99[[1]]
gooddyads2022$in2 = in_df2_99
gooddyads2022$in3 =  in_df3_99
gooddyads2022$in4 =  in_df4_99
gooddyads2022$in5 =  in_df5_99
gooddyads2022$inall = in_df2_99 & in_df3_99 & in_df4_99 & in_df5_99
gooddyads2022$num_true = 1 + rowSums(gooddyads2022[,colnames(gooddyads2022) %in% c("in2_99", "in3_99", "in4_99", "in5_99")])
gooddyads2022$pairname = collapse(gooddyads2022)

# check which missing links are "good" but not "best"
partial99 = gooddyads2022[gooddyads2022$pairname %in% missing_df$missing_pair,]

# add these to consistent dyad
intermediate_set_2022 = rbind(consistentdyad_2022, partial99[,!colnames(partial99) %in% c("pairname", "num_true")])

###########################################################
### Create graph object and plot families as components
###########################################################

# initialize edge color vector
igraph::E(graph_22)$color <- rep("gray", ecount(graph_22))

# color in the new links as we add them
new_edges = unlist(str_split(partial99$pairname, "-"))
graph_22 <- add_edges(graph_22, new_edges)
new_edge_ids <- (ecount(graph_22) - (length(new_edges)/2) + 1):ecount(graph_22)
igraph::E(graph_22)$color[new_edge_ids] <- "black"

#color in all the links that are still missing
still_missing = missing_df[!missing_df$missing_pair %in% partial99$pairname,]
missing_edges = unlist(str_split(still_missing, "-"))
graph_22_withmissing <- add_edges(graph_22, missing_edges)
new_edge_ids <- (ecount(graph_22_withmissing) - (length(new_edges)/2) + 1):ecount(graph_22_withmissing)
igraph::E(graph_22_withmissing)$color[new_edge_ids] <- "red"

# set size and label
V(graph_22_withmissing)$label = V(graph_22_withmissing)$name
V(graph_22_withmissing)$size <- 4
V(graph_22_withmissing)$label.cex <- 1

# Set consistent layout
#set.seed(123)
layout <- layout_with_fr(graph_22_withmissing)

# save
png("figures/manuscript_figures/impatiens2022_familystructure_showmissinglinks.png", width = 2000, height = 2000)
plot(graph_22_withmissing,
     vertex.label.cex = 1,
     vertex.label.color = "black",
     edge.width = 4,
     vertex.size = 2,
     main = "Impatiens Family Structures in 2022 -- Missing Links in Red")
dev.off()


###########################################################
### Remove links to create circular cliques
###########################################################

# make graph
edges = intermediate_set_2022[, c("OffspringID1", "OffspringID2")]
graph_22 = graph_from_data_frame(d = edges, directed = FALSE)
components <- decompose(graph_22)

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
final_graph <- do.call(disjoint_union, filtered_components)

# color by family
V(final_graph)$family_id = components$membership
num_families <- components$no
family_colors <- rainbow(num_families)
V(graph_22)$color = family_colors[V(graph_22)$family_id]

png("figures/manuscript_figures/impatiens2022_familystructure.png", width = 2000, height = 2000)
plot(final_graph,
     vertex.label = name,
     vertex.label.cex = 0.6,
     vertex.label.color = "black",
     edge.width = 2,
     vertex.size = 2,
     main = "Impatiens Family Structures in 2022")
dev.off()











# assign component membership as a vertex attribute
V(graph_22)$family_id = components$membership

# choose a color for each family
num_families <- components$no
family_colors <- rainbow(num_families)
V(graph_22)$color <- family_colors[V(graph_22)$family_id]
