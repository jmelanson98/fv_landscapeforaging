########
## Get colony IDs from microsatellite genotyping results--Bombus mixtus
## Started by J Melanson
## July 10, 2024
########

# prep workspace
rm(list = ls())
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")

# first, load in packages
source('colony_assignments/src/init.R')
source('colony_assignments/src/colony_assignment_functions.R')
library(dplyr)
library(tidyr)
library(stringr)
library(rcolony)
library(genepop)
library(data.table)
library(cowplot)
library(igraph)


# next; load in specimen data and allele tables
#load specimen field data
specimenData2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2022specimendata.csv", sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023specimendata.csv", sep = ",", header = T, fill = TRUE))
sampleEffort2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023sampledata.csv", sep = ","))
samplePoints = as.data.frame(read.table("/Users/jenna1/Documents/UBC/bombus_project/raw_data/allsamplepoints.csv", sep = ","))
colnames(samplePoints) = c("sample_pt", "gps","landowner","subsite")

# load allele tables
alleles_plex1_all = as.data.frame(read.table("data/from_geneious/mixtus_plex1.csv", sep = ",", header = T))
alleles_plex2_all = as.data.frame(read.table("data/from_geneious/mixtus_plex2.csv", sep = ",", header = T))
extraqueens_plex1 = as.data.frame(read.table("data/from_geneious/extramixtusqueens_plex1_binned.csv", sep = ",", header = T))
extraqueens_plex2 = as.data.frame(read.table("data/from_geneious/extramixtusqueens_plex2_binned.csv", sep = ",", header = T))

# modify specimen table for merging
specimenData2022$location = str_remove_all(paste(specimenData2022$plate, specimenData2022$well), "_| ")
specimenData2023$plate[specimenData2023$plate == "mixtus_queen_DNA"] = "mixtusQueen"
specimenData2023$location = str_remove_all(paste(specimenData2023$plate, specimenData2023$well), "_| ")

# make a dataframe that tells us where each barcode is located (in which plate / well location)
barcodesubset = rbind(specimenData2022[,colnames(specimenData2022) %in% c("barcode_id", "location")],
                   specimenData2023[,colnames(specimenData2023) %in% c("barcode_id", "location")])

# split allele table into originals and reruns
# originals are named by location (plate / well), reruns are named by barcode
alleles_plex1_original = alleles_plex1_all[str_detect(alleles_plex1_all$Name, "mixtus"),]
alleles_plex1_rerun = alleles_plex1_all[!str_detect(alleles_plex1_all$Name, "mixtus"),]

alleles_plex2_original = alleles_plex2_all[str_detect(alleles_plex2_all$Name, "mixtus"),]
alleles_plex2_rerun = alleles_plex2_all[!str_detect(alleles_plex2_all$Name, "mixtus"),]

# for originals, modify allele table "name" value to match specimen data frame "plate" and "well"
alleles_plex1_original = alleles_plex1_original %>% 
  separate(Name, sep = "-", c(NA, "plate", "well"), extra = "merge", fill = "left")                                       
alleles_plex1_original$location = str_remove_all(paste(alleles_plex1_original$plate, alleles_plex1_original$well), "-| ")
alleles_plex1_original = alleles_plex1_original[,!colnames(alleles_plex1_original) %in% c("plate", "well")]

alleles_plex2_original = alleles_plex2_original %>% 
  separate(Name, sep = "-", c(NA, "plate", "well"), extra = "merge", fill = "left")                                         
alleles_plex2_original$location = str_remove_all(paste(alleles_plex2_original$plate, alleles_plex2_original$well), "-| ")
alleles_plex2_original = alleles_plex2_original[,!colnames(alleles_plex2_original) %in% c("plate", "well")]

# merge original dataframes to barcode subset and remove location column
alleles_plex1_wbarcode = left_join(alleles_plex1_original, barcodesubset, by = "location")
alleles_plex1_wbarcode = alleles_plex1_wbarcode[,!colnames(alleles_plex1_wbarcode) %in% c("location")]
alleles_plex1_wbarcode = alleles_plex1_wbarcode %>% 
  filter(!is.na(barcode_id))

alleles_plex2_wbarcode = left_join(alleles_plex2_original, barcodesubset, by = "location")
alleles_plex2_wbarcode = alleles_plex2_wbarcode[,!colnames(alleles_plex2_wbarcode) %in% c("location")]
alleles_plex2_wbarcode = alleles_plex2_wbarcode %>% 
  filter(!is.na(barcode_id))

# modify reruns to contain just barcode
alleles_plex1_rerun$barcode_id = gsub("^.{0,6}", "", alleles_plex1_rerun$Name)
alleles_plex1_rerun$barcode_id = gsub("-","_",alleles_plex1_rerun$barcode_id)
alleles_plex1_rerun = alleles_plex1_rerun[,!colnames(alleles_plex1_rerun) %in% c("Name")]

alleles_plex2_rerun$barcode_id = gsub("^.{0,6}", "", alleles_plex2_rerun$Name)
alleles_plex2_rerun$barcode_id = gsub("-","_",alleles_plex2_rerun$barcode_id)
alleles_plex2_rerun = alleles_plex2_rerun[,!colnames(alleles_plex2_rerun) %in% c("Name")]

extraqueens_plex1$barcode_id = gsub("^.{0,6}", "", extraqueens_plex1$Name)
extraqueens_plex1$barcode_id = gsub("-","_",extraqueens_plex1$barcode_id)
extraqueens_plex1 = extraqueens_plex1[,!colnames(extraqueens_plex1) %in% c("Name")]

extraqueens_plex2$barcode_id = gsub("^.{0,6}", "", extraqueens_plex2$Name)
extraqueens_plex2$barcode_id = gsub("-","_",extraqueens_plex2$barcode_id)
extraqueens_plex2 = extraqueens_plex2[,!colnames(extraqueens_plex2) %in% c("Name")]

#order columns so that we can rbind originals to reruns
correctorderp1 = colnames(alleles_plex1_wbarcode)
correctorderp2 = colnames(alleles_plex2_wbarcode)

alleles_plex1_rerun <- alleles_plex1_rerun[, correctorderp1]
alleles_plex2_rerun <- alleles_plex2_rerun[, correctorderp2]
extraqueens_plex1 <- extraqueens_plex1[, correctorderp1]
extraqueens_plex2 <- extraqueens_plex2[, correctorderp2]

alleles_plex1 = rbind(alleles_plex1_wbarcode, alleles_plex1_rerun, extraqueens_plex1)
alleles_plex2 = rbind(alleles_plex2_wbarcode, alleles_plex2_rerun, extraqueens_plex2)

#check for duplicates in these dataframes
print("Number of rows in alleles 1 dataframe:")
dim(alleles_plex1)[1]
print("Number of unique barcodes in alleles 1 dataframe (should be the same):")
length(unique(alleles_plex1$barcode_id))
print("Duplicates in final dataframe:")
duplicated = alleles_plex1$barcode_id[duplicated(alleles_plex1$barcode_id)]
print(duplicated)

print("Number of rows in alleles 2 dataframe:")
dim(alleles_plex2)[1]
print("Number of unique barcodes in alleles 2 dataframe (should be the same):")
length(unique(alleles_plex2$barcode_id))
print("Duplicates in final dataframe:")
duplicated = alleles_plex2$barcode_id[duplicated(alleles_plex2$barcode_id)]
print(duplicated)


#merge plexes
microsat_scores = full_join(alleles_plex1, alleles_plex2, by = "barcode_id")

# filter by year and caste
# make a specimen dataframe that inclues year and notes
rowsToKeep = c("sample_id", "barcode_id", "year", "notes", "final_id")
specsubset = rbind(specimenData2022[,colnames(specimenData2022) %in% rowsToKeep],
                      specimenData2023[,colnames(specimenData2023) %in% rowsToKeep])
specimens_withscores = left_join(microsat_scores, specsubset, by = "barcode_id")

# join to sample effort df to get julian dates (only need for 2023)
specimens_withdates = joinJulianDates(specimens_withscores, sampleEffort2023)

# split by year, including early season 2023 queens with 2022 subset
# julian date cutoff based on script "queen_phenology.R"
mixtus2022 = filter(specimens_withdates, year == "2022" | (str_detect(notes, "queen") & julian_date < 160))
mixtus2023 = filter(specimens_withdates, year == "2023" & (!str_detect(notes, "queen") | julian_date > 160))


#check which bees are unscored
specimens = full_join(microsat_scores, specsubset, by = "barcode_id")
mixtus = filter(specimens, final_id == "B. mixtus" & barcode_id != "N/A")
unscoredmixtus = mixtus[rowSums(is.na(mixtus)) >10,]

##################################
# Prep genotype data for COLONY
##################################

#remove all columns except barcode and scores
# also remove loci with high error rates and those out of equilibrium
colsToRemove = c("year", 
                 "final_id", 
                 "notes", 
                 "sample_id", 
                 "date", 
                 "julian_date", 
                 "BL15...1", 
                 "BL15...2") 
                 #"BTMS0072...1", 
                 #"BTMS0072...2",
                 #"BTERN01...1",
                 #"BTERN01...2",
                 #"BTMS0104...1", 
                 #"BTMS0104...2",
                 #"BTMS0059...1", 
                 #"BTMS0059...2")
mixtus2022 = mixtus2022[, !colnames(mixtus2022) %in% colsToRemove]
mixtus2023 = mixtus2023[, !colnames(mixtus2023) %in% colsToRemove]

#relocate barcode, remove rows with more than 10 NAs, replace NAs with 0's
mixtus2022 = mixtus2022 %>% relocate(barcode_id)
mixtus2022_forcolony = mixtus2022[rowSums(is.na(mixtus2022)) <= 12,]
mixtus2022_forcolony[is.na(mixtus2022_forcolony)] = 0

mixtus2023 = mixtus2023 %>% relocate(barcode_id)
mixtus2023_forcolony = mixtus2023[rowSums(is.na(mixtus2023)) <= 12,]
mixtus2023_forcolony[is.na(mixtus2023_forcolony)] = 0

#write csvs for upload to colony
column_names = colnames(mixtus2022_forcolony)

write.table(mixtus2022_forcolony, "data/merged_by_year/sib_scores_for_colony/mixtus2022_forcolony.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(mixtus2023_forcolony, "data/merged_by_year/sib_scores_for_colony/mixtus2023_forcolony.txt", sep= ",", col.names = FALSE, row.names = FALSE)


write.csv(mixtus2022_forcolony, "data/merged_by_year/csvs/mixtus_2022_scores.csv")
write.csv(mixtus2023_forcolony, "data/merged_by_year/csvs/mixtus_2023_scores.csv")

###########################################
# Prep microsat error rates for COLONY
###########################################

# first value per column: marker name
# second value per column: marker type (codominant for microsats -> 0)
# third value per column: allelic dropout rate ?? (set to 0)
# fourth value per column: error/mutation rate (empirically derived, see check_microsat_error_rates.R)
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
write.table(mixtus_error_rates, "data/merged_by_year/error_rates/mixtus_error_rates.txt", 
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
sample.names.2022 = mixtus2022_forcolony[,1]
sample.names.2023 = mixtus2023_forcolony[,1]
'%!in%' <- function(x,y)!('%in%'(x,y))

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
  file = "data/merged_by_year/sib_exclusions/mixtus_sibexclusions_2022.txt",
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  na = ""
)

write.table(
  sib2023_reduced,
  file = "data/merged_by_year/sib_exclusions/mixtus_sibexclusions_2023.txt",
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


# build .DAT files for both datasets
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux/mixtus2022_final1.DAT", delim=",")
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/_Linux", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux/mixtus2023_final1.DAT", delim=",")

# navigate to COLONY2 sub folder and run the following in terminal
#NOTE: ensure that mixtus2022.DAT and mixtus2023.DAT are in the same directory at the colony2s.out executable

# cd /Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2
# ./colony2s.out IFN:mixtus2022.DAT
# ./colony2s.out IFN:mixtus2023.DAT



#####################################
## Load in results from colony
#####################################

list_of_results = list()
for (i in 1:5){
  # read in dyad data from COLONY
  fsDyad2022 = trimws(readLines(paste0("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux/mixtus2022_final", i, ".FullSibDyad")))
  
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
  fsDyad2022 = trimws(readLines(paste0("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux/mixtus2022_final", i, ".FullSibDyad")))
  
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
intermediate_set_2022 = rbind(consistentdyad_2022, partial99[,!colnames(partial99) %in% c("pairname")])

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
new_edge_ids <- (ecount(graph_22_withmissing) - (length(missing_edges)/2) + 1):ecount(graph_22_withmissing)
igraph::E(graph_22_withmissing)$color[new_edge_ids] <- "red"

# save
set.seed(123)
layout <- layout_with_fr(graph_22_withmissing)
png("figures/manuscript_figures/mixtus2022_familystructure_showmissinglinks.png", width = 2000, height = 2000)
plot(graph_22_withmissing,
     vertex.label = V(graph_22_withmissing)$name,
     vertex.label.cex = 1,
     vertex.label.color = "black",
     edge.width = 4,
     vertex.size = 2,
     main = "mixtus Family Structures in 2022 -- Missing Links in Red")
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
final_graph = do.call(disjoint_union, filtered_components)
comp_info = components(final_graph)

# color by family
V(final_graph)$family_id = comp_info$membership
num_families <- comp_info$no
family_colors <- rainbow(num_families)
V(final_graph)$color = family_colors[V(final_graph)$family_id]

png("figures/manuscript_figures/mixtus2022_familystructure.png", width = 2000, height = 2000)
plot(final_graph,
     vertex.label = V(final_graph)$name,
     vertex.label.cex = 0.6,
     vertex.label.color = "black",
     edge.width = 2,
     vertex.size = 2,
     main = "mixtus Family Structures in 2022")
dev.off()


#####################################
### Make final output dataframe
#####################################
sibship2022 =  data.frame(
  barcode_id = V(final_graph)$name,
  sibship_id = comp_info$membership
)

# identify barcodes not in a sibship
singleton_barcodes <- mixtus2022_forcolony$barcode_id[!mixtus2022_forcolony$barcode_id %in% sibship2022$barcode_id]

# create singleton rows with unique sibship IDs
singletons2022 <- data.frame(
  barcode_id = singleton_barcodes,
  sibship_id = max(sibship2022$sibship_id) + seq_along(singleton_barcodes)
)

#combine the two dataframes and save
mixtus2022colonies = rbind(sibship2022, singletons2022)
write.csv(mixtus2022colonies, "data/siblingships/mixtus_sibships_2022.csv", row.names = FALSE)


###############################################################################
###############################################################################
# Now repeat all of that for 2023 mixtus.....
###############################################################################
###############################################################################


#####################################
## Load in results from colony
#####################################

list_of_results = list()
for (i in 1:5){
  # read in dyad data from COLONY
  fsDyad2023 = trimws(readLines(paste0("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux/mixtus2023_final", i, ".FullSibDyad")))
  
  #format
  dyads2023 = do.call(rbind, strsplit(fsDyad2023, ","))
  colnames(dyads2023) = dyads2023[1,]
  dyads2023 = as.data.frame(dyads2023[-1,])
  dyads2023$Probability = as.numeric(dyads2023$Probability)
  
  # remove sibpairs below some threshold probability
  bestdyads2023 = dyads2023 %>% filter(Probability == 1)
  
  # add results to list
  list_of_results[[i]] = bestdyads2023
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

bestdyads2023 = list_of_results[[1]]
bestdyads2023$in2 = in_df2
bestdyads2023$in3 =  in_df3
bestdyads2023$in4 =  in_df4
bestdyads2023$in5 =  in_df5
bestdyads2023$inall = in_df2 & in_df3 & in_df4 & in_df5
bestdyads2023$num_true = 1 + rowSums(bestdyads2023[,colnames(bestdyads2023) %in% c("in2", "in3", "in4", "in5")])

# subset to the consistent dyads
consistentdyad_2023 = bestdyads2023[bestdyads2023$inall ==TRUE,]


###########################################################
### Identify missing links to establish complete cliques
###########################################################

# make igraph object
edges = consistentdyad_2023[, c("OffspringID1", "OffspringID2")]
graph_23 = graph_from_data_frame(d = edges, directed = FALSE)
components <- decompose(graph_23)

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
  fsDyad2023 = trimws(readLines(paste0("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2_Linux/mixtus2023_final", i, ".FullSibDyad")))
  
  #format
  dyads2023 = do.call(rbind, strsplit(fsDyad2023, ","))
  colnames(dyads2023) = dyads2023[1,]
  dyads2023 = as.data.frame(dyads2023[-1,])
  dyads2023$Probability = as.numeric(dyads2023$Probability)
  
  # remove sibpairs below some threshold probability
  bestdyads2023 = dyads2023 %>% filter(Probability > 0.99)
  
  # add results to list
  list_of_results_99[[i]] = bestdyads2023
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

gooddyads2023 = list_of_results_99[[1]]
gooddyads2023$in2 = in_df2_99
gooddyads2023$in3 =  in_df3_99
gooddyads2023$in4 =  in_df4_99
gooddyads2023$in5 =  in_df5_99
gooddyads2023$inall = in_df2_99 & in_df3_99 & in_df4_99 & in_df5_99
gooddyads2023$num_true = 1 + rowSums(gooddyads2023[,colnames(gooddyads2023) %in% c("in2_99", "in3_99", "in4_99", "in5_99")])
gooddyads2023$pairname = collapse(gooddyads2023)

# check which missing links are "good" but not "best"
partial99 = gooddyads2023[gooddyads2023$pairname %in% missing_df$missing_pair,]

# add these to consistent dyad
intermediate_set_2023 = rbind(consistentdyad_2023, partial99[,!colnames(partial99) %in% c("pairname")])

###########################################################
### Create graph object and plot families as components
###########################################################

# initialize edge color vector
igraph::E(graph_23)$color <- rep("gray", ecount(graph_23))

# color in the new links as we add them
new_edges = unlist(str_split(partial99$pairname, "-"))
graph_23 <- add_edges(graph_23, new_edges)
new_edge_ids <- (ecount(graph_23) - (length(new_edges)/2) + 1):ecount(graph_23)
igraph::E(graph_23)$color[new_edge_ids] <- "black"

#color in all the links that are still missing
still_missing = missing_df[!missing_df$missing_pair %in% partial99$pairname,]
missing_edges = unlist(str_split(still_missing, "-"))
graph_23_withmissing <- add_edges(graph_23, missing_edges)
new_edge_ids <- (ecount(graph_23_withmissing) - (length(missing_edges)/2) + 1):ecount(graph_23_withmissing)
igraph::E(graph_23_withmissing)$color[new_edge_ids] <- "red"

# save
set.seed(123)
layout <- layout_with_fr(graph_23_withmissing)
png("figures/manuscript_figures/mixtus2023_familystructure_showmissinglinks.png", width = 2000, height = 2000)
plot(graph_23_withmissing,
     vertex.label = V(graph_23_withmissing)$name,
     vertex.label.cex = 0.6,
     vertex.label.color = "black",
     edge.width = 4,
     vertex.size = 2,
     main = "mixtus Family Structures in 2023 -- Missing Links in Red")
dev.off()


###########################################################
### Remove links to create circular cliques
###########################################################

# make graph
edges = intermediate_set_2023[, c("OffspringID1", "OffspringID2")]
graph_23 = graph_from_data_frame(d = edges, directed = FALSE)
components <- decompose(graph_23)

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
comp_info = components(final_graph)

# color by family
V(final_graph)$family_id = comp_info$membership
num_families <- comp_info$no
family_colors <- rainbow(num_families)
V(final_graph)$color = family_colors[V(final_graph)$family_id]

png("figures/manuscript_figures/mixtus2023_familystructure.png", width = 2000, height = 2000)
plot(final_graph,
     vertex.label = V(final_graph)$name,
     vertex.label.cex = 0.6,
     vertex.label.color = "black",
     edge.width = 2,
     vertex.size = 2,
     main = "mixtus Family Structures in 2023")
dev.off()

#####################################
### Make final output dataframe
#####################################
sibship2023 =  data.frame(
  barcode_id = V(final_graph)$name,
  sibship_id = max(mixtus2022colonies$sibship_id) + comp_info$membership
)

# identify barcodes not in a sibship
singleton_barcodes <- mixtus2023_forcolony$barcode_id[!mixtus2023_forcolony$barcode_id %in% sibship2023$barcode_id]

# create singleton rows with unique sibship IDs
singletons2023 <- data.frame(
  barcode_id = singleton_barcodes,
  sibship_id = max(sibship2023$sibship_id) + seq_along(singleton_barcodes)
)

#combine the two dataframes and save
mixtus2023colonies = rbind(sibship2023, singletons2023)
write.csv(mixtus2023colonies, "data/siblingships/mixtus_sibships_2023.csv", row.names = FALSE)
