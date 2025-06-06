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
source('colony_assignments/src/joinFunctions.R')
library(dplyr)
library(tidyr)
library(stringr)
library(rcolony)
library(genepop)
library(data.table)
library(cowplot)
library(raster)


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
# also remove loci with high error rates
colsToRemove = c("year", 
                 "final_id", 
                 "notes", 
                 "sample_id", 
                 "date", 
                 "julian_date", 
                 "BL15...1", 
                 "BL15...2", 
                 "BTMS0072...1", 
                 "BTMS0072...2")
mixtus2022 = mixtus2022[, !colnames(mixtus2022) %in% colsToRemove]
mixtus2023 = mixtus2023[, !colnames(mixtus2023) %in% colsToRemove]

#relocate barcode, remove rows with more than 10 NAs, replace NAs with 0's
mixtus2022 = mixtus2022 %>% relocate(barcode_id)
mixtus2022_forcolony = mixtus2022[rowSums(is.na(mixtus2022)) < 10,]
mixtus2022_forcolony[is.na(mixtus2022_forcolony)] = 0

mixtus2023 = mixtus2023 %>% relocate(barcode_id)
mixtus2023_forcolony = mixtus2023[rowSums(is.na(mixtus2023)) < 10,]
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
                         c("BTMS0104", 0, 0, 0.01),
                         c("BTMS0057", 0, 0, 0.01),
                         c("BTMS0086", 0, 0, 0.015),
                         c("BTMS0066", 0, 0, 0.01),
                         c("BTMS0062", 0, 0, 0.015),
                         c("BTMS0136", 0, 0, 0.01),
                         c("BTERN01", 0, 0, 0.0215),
                         c("BTMS0126", 0, 0, 0.0162),
                         c("BTMS0059", 0, 0, 0.01),
                         c("BL13", 0, 0, 0.0159),
                         c("BTMS0083", 0, 0, 0.01),
                         c("B126", 0, 0, 0.01))
write.table(mixtus_error_rates, "data/merged_by_year/error_rates/mixtus_error_rates.txt", 
            sep= "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

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
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/initial", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/initial/mixtus2022.DAT", delim=",")
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/initial", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/initial/mixtus2023.DAT", delim=",")

# navigate to COLONY2 sub folder and run the following in terminal
#NOTE: ensure that mixtus2022.DAT and mixtus2023.DAT are in the same directory at the colony2s.out executable

# cd /Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2
# ./colony2s.out IFN:mixtus2022.DAT
# ./colony2s.out IFN:mixtus2023.DAT



#####################################
## Load in results from colony
#####################################

# read in and format 2022 data
lines2022 <- trimws(readLines("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/initial/mixtus2022.BestCluster"))
split_lines2022 <- strsplit(lines2022, "\\s+")
bestconfig2022 <- do.call(rbind, split_lines2022)
colnames(bestconfig2022) <- bestconfig2022[1, ]
bestconfig2022 <- as.data.frame(bestconfig2022[-1, ])
bestconfig2022$Probability = as.numeric(bestconfig2022$Probability)

# read in and format 2023 data
lines2023 <- trimws(readLines("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/initial/mixtus2023.BestCluster"))
split_lines2023 <- strsplit(lines2023, "\\s+")
bestconfig2023 <- do.call(rbind, split_lines2023)
colnames(bestconfig2023) <- bestconfig2023[1, ]
bestconfig2023 <- as.data.frame(bestconfig2023[-1, ])
bestconfig2023$Probability = as.numeric(bestconfig2023$Probability)

# remove sibling relationships with p < 0.95
counter = max(as.numeric(bestconfig2022$ClusterIndex)) + 1
for (i in 1:nrow(bestconfig2022)){
  # if the probability of inclusion in the cluster is < 0.95, make a new cluster
  if (bestconfig2022$Probability[i] < 0.95){
    bestconfig2022$ClusterIndex[i] = counter
    counter = counter + 1
  }
}

# repeat for 2023
counter = max(as.numeric(bestconfig2023$ClusterIndex)) + 1
for (i in 1:nrow(bestconfig2023)){
  # if the probability of inclusion in the cluster is < 0.95, make a new cluster
  if (bestconfig2023$Probability[i] < 0.95){
    bestconfig2023$ClusterIndex[i] = counter
    counter = counter + 1
  }
}


#write csvs to file
write.csv(bestconfig2022, "data/siblingships/mix_sibships_preliminary_2022.csv")
write.csv(bestconfig2023, "data/siblingships/mix_sibships_preliminary_2023.csv")


#################################################
# More stringent approach for colony assignments
#################################################

# try with FS dyads instead of FS clusters
#run 1
fsDyad2022_1 = trimws(readLines("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/initial/mixtus_2022_monogamous/mixtus_2022_monogamous.FullSibDyad"))
fsDyad2022_1 <- strsplit(fsDyad2022_1, ",")
fsDyad2022_1 <- do.call(rbind, fsDyad2022_1)
colnames(fsDyad2022_1) <- fsDyad2022_1[1, ]
fsDyad2022_1 <- as.data.frame(fsDyad2022_1[-1, ])
fsDyad2022_1$Probability = as.numeric(fsDyad2022_1$Probability)

# run 2
fsDyad2022_2 = trimws(readLines("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/initial/mixtus_2022_monogamous_2/mixtus_2022_monogamous_2.FullSibDyad"))
fsDyad2022_2 <- strsplit(fsDyad2022_2, ",")
fsDyad2022_2 <- do.call(rbind, fsDyad2022_2)
colnames(fsDyad2022_2) <- fsDyad2022_2[1, ]
fsDyad2022_2 <- as.data.frame(fsDyad2022_2[-1, ])
fsDyad2022_2$Probability = as.numeric(fsDyad2022_2$Probability)

# run 3
fsDyad2022_3 = trimws(readLines("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/initial/mixtus_2022_monogamous_3/mixtus_2022_monogamous_3.FullSibDyad"))
fsDyad2022_3 <- strsplit(fsDyad2022_3, ",")
fsDyad2022_3 <- do.call(rbind, fsDyad2022_3)
colnames(fsDyad2022_3) <- fsDyad2022_3[1, ]
fsDyad2022_3 <- as.data.frame(fsDyad2022_3[-1, ])
fsDyad2022_3$Probability = as.numeric(fsDyad2022_3$Probability)

# remove sibpairs below some threshold probability (p = 0.95)
bestdyads2022_1 = fsDyad2022_1 %>% filter(Probability > 0.99)
bestdyads2022_2 = fsDyad2022_2 %>% filter(Probability > 0.99)
bestdyads2022_3 = fsDyad2022_3 %>% filter(Probability > 0.99)

example = c(bestdyads2022_1[2, c("OffspringID1")], bestdyads2022_1[2, c("OffspringID2")])
comparison = example %in% bestdyads2022_1[,c("OffspringID1", "OffspringID2")]


# retain dyads which are present in all 3 runs
# first, collapse pairs into an ordered character string
collapse <- function(df) {
  apply(df, 1, function(row) paste(sort(row[c("OffspringID1", "OffspringID2")]), collapse = "-"))
}

# generate pairs
pairs_df1 <- collapse(bestdyads2022_1)
pairs_df2 <- collapse(bestdyads2022_2)
pairs_df3 <- collapse(bestdyads2022_3)

# check matches
in_df2 <- pairs_df1 %in% pairs_df2
in_df3 <- pairs_df1 %in% pairs_df3

bestdyads2022_1$in2 = in_df2
bestdyads2022_1$in3 =  in_df3
bestdyads2022_1$inall = in_df2 & in_df3

# subset to the consistent subset
consistentdyad_2022 = bestdyads2022_1[bestdyads2022_1$inall ==TRUE,]


# make an undirected graph of families
graph_22 <- graph_from_data_frame(consistentdyad_2022[, c("OffspringID1", "OffspringID2")], directed = FALSE)

# find connected components (families)
components <- components(graph_22)

# assign family ID to each individual
membership_df <- data.frame(
  OffspringID = names(components$membership),
  Family = components$membership
)

# then, check for non-circularity (e.g., A related to B, B related to C, A not related to C)
flags <- lapply(unique(membership_df$Family), function(fam_id) {
  nodes <- membership_df$OffspringID[membership_df$Family == fam_id]
  subg <- induced_subgraph(graph_22, vids = nodes)
  
  n <- gorder(subg)
  expected_edges <- n * (n - 1) / 2
  actual_edges <- gsize(subg)
  
  data.frame(Family = fam_id, FlagIncomplete = actual_edges < expected_edges)
})

flag_df <- bind_rows(flags)

# add flags to output
df_with_family <- consistentdyad_2022 %>%
  left_join(membership_df, by = c("OffspringID1" = "OffspringID")) %>%
  rename(Family = Family)

df_with_flags <- df_with_family %>%
  left_join(flag_df, by = "Family")

############################################
# Remove siblings before running genepop
############################################
nonsibs2022 = left_join(specimenData2022, bestconfig2022, by = c("barcode_id" ="OffspringID")) %>%
  filter(!is.na(ClusterIndex)) %>%
  group_by(ClusterIndex) %>%
  sample_n(1)

testalleles2022 = mixtus2022_forcolony[mixtus2022_forcolony$barcode_id %in% nonsibs2022$barcode_id,]


nonsibs2023 = left_join(specimenData2023, bestconfig2023, by = c("barcode_id" ="OffspringID")) %>%
  filter(!is.na(ClusterIndex)) %>%
  group_by(ClusterIndex) %>%
  sample_n(1)

testalleles2023 = mixtus2023_forcolony[mixtus2023_forcolony$barcode_id %in% nonsibs2023$barcode_id,]

######################################
# Prep alleles table for genepop
#####################################

#combine columns for each locus
ms_m22 = data.table(testalleles2022)
ms_m22[ms_m22 ==0] <- "000"


i1 <- seq(1, length(ms_m22)-1, 2)
i2 <- seq(2, length(ms_m22)-1, 2)
msgenepop_m22 = ms_m22[, Map(paste,
                         .SD[, i1, with = FALSE], .SD[, i2, with = FALSE], 
                         MoreArgs = list(sep="")), 
                   by = "barcode_id"]
colnames(msgenepop_m22) = c("barcode_id", "BT10", "BTMS0104", "BTMS0057", "BTMS0086",
                        "BTMS0066", "BTMS0062", "BTMS0136","BTERN01", "BTMS0126", 
                        "BTMS0059", "BL13", "BTMS0083", "B126")


#select and join alleles

#sort allele table by barcode ID
msgenepop_m22 <- msgenepop_m22[order(msgenepop_m22$barcode_id),]
msgenepop_m22$site = substring(msgenepop_m22$barcode_id, 1, 1)

haplotypes_m22 <- as.data.frame(paste(msgenepop_m22$site, "_", msgenepop_m22$barcode_id, ","," ", 
                                      msgenepop_m22$BT10," ", msgenepop_m22$BTMS0104, " ",
                                      msgenepop_m22$BTMS0057," ", msgenepop_m22$BTMS0086, " ",
                                      msgenepop_m22$BTMS0066," ", msgenepop_m22$BTMS0062, " ",
                                      msgenepop_m22$BTMS0136," ", msgenepop_m22$BTERN01, " ",
                                      msgenepop_m22$BTMS0126," ", msgenepop_m22$BTMS0059," ",
                                      msgenepop_m22$BL13," ", msgenepop_m22$BTMS0083, " ",
                                      msgenepop_m22$B126,
                                  sep = ""))

#Build the Genepop format
sink("data/genepop/ms_genepop_format_mixtus2022.txt")
cat("Title: B. mixtus workers (2022) \n")
cat("BT10 \n")
cat("BTMS0104 \n")
cat("BTMS0057 \n")
cat("BTMS0086 \n")
cat("BTMS0066 \n")
cat("BTMS0062 \n")
cat("BTMS0136 \n")
cat("BTERN01 \n")
cat("BTMS0126 \n")
cat("BTMS0059 \n")
cat("BL13 \n")
cat("BTMS0083 \n")
cat("B126 \n")
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_m22[substring(haplotypes_m22[,1],1,1) == "E",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_m22[substring(haplotypes_m22[,1],1,1) == "W",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_m22[substring(haplotypes_m22[,1],1,1) == "S",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_m22[substring(haplotypes_m22[,1],1,1) == "N",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_m22[substring(haplotypes_m22[,1],1,1) == "P",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_m22[substring(haplotypes_m22[,1],1,1) == "H",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
sink()



######################################
# Run genepop on samples
#####################################

#test for HWE
test_HW("data/genepop/ms_genepop_format_mixtus2022.txt", which="Proba", "data/genepop/mixtusHWE2022.txt.D")

#test for LD
test_LD("data/genepop/ms_genepop_format_mixtus2022.txt","data/genepop/mixtusLD2022.txt.DIS")



######################################
# Repeat for 2023 mixtus !
#####################################

##### Prepare alleles table for genepop

#combine columns for each locus
ms_m23 = data.table(testalleles2023)
ms_m23[ms_m23 ==0] <- "000"

i1 <- seq(1, length(ms_m23)-1, 2)
i2 <- seq(2, length(ms_m23)-1, 2)
msgenepop_m23 = ms_m23[, Map(paste,
                             .SD[, i1, with = FALSE], .SD[, i2, with = FALSE], 
                             MoreArgs = list(sep="")), 
                       by = "barcode_id"]
colnames(msgenepop_m23) = c("barcode_id", "BT10", "BTMS0104", "BTMS0057", "BTMS0086",
                            "BTMS0066", "BTMS0062", "BTMS0136","BTERN01", "BTMS0126", 
                            "BTMS0059", "BL13", "BTMS0083", "B126")


#select and join alleles

#sort allele table by barcode ID
msgenepop_m23 <- msgenepop_m23[order(msgenepop_m23$barcode_id),]
msgenepop_m23$site = substring(msgenepop_m23$barcode_id, 1, 1)

haplotypes_m23 <- as.data.frame(paste(msgenepop_m23$site, "_", msgenepop_m23$barcode_id, ","," ", 
                                      msgenepop_m23$BT10," ", msgenepop_m23$BTMS0104, " ",
                                      msgenepop_m23$BTMS0057," ", msgenepop_m23$BTMS0086, " ",
                                      msgenepop_m23$BTMS0066," ", msgenepop_m23$BTMS0062, " ",
                                      msgenepop_m23$BTMS0136," ", msgenepop_m23$BTERN01, " ",
                                      msgenepop_m23$BTMS0126," ", msgenepop_m23$BTMS0059," ",
                                      msgenepop_m23$BL13," ", msgenepop_m23$BTMS0083, " ",
                                      msgenepop_m23$B126,
                                      sep = ""))

#Build the Genepop format
sink("data/genepop/ms_genepop_format_mixtus2023.txt")
cat("Title: B. mixtus workers (2023) \n")
cat("BT10 \n")
cat("BTMS0104 \n")
cat("BTMS0057 \n")
cat("BTMS0086 \n")
cat("BTMS0066 \n")
cat("BTMS0062 \n")
cat("BTMS0136 \n")
cat("BTERN01 \n")
cat("BTMS0126 \n")
cat("BTMS0059 \n")
cat("BL13 \n")
cat("BTMS0083 \n")
cat("B126 \n")
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_m23[substring(haplotypes_m23[,1],1,1) == "E",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_m23[substring(haplotypes_m23[,1],1,1) == "W",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_m23[substring(haplotypes_m23[,1],1,1) == "S",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_m23[substring(haplotypes_m23[,1],1,1) == "N",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_m23[substring(haplotypes_m23[,1],1,1) == "P",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_m23[substring(haplotypes_m23[,1],1,1) == "H",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
sink()



######################################
# Run genepop on samples
#####################################

#test for HWE
test_HW("data/genepop/ms_genepop_format_mixtus2023.txt", which="Proba", "data/genepop/mixtusHWE2023.txt.D")

#test for LD
test_LD("data/genepop/ms_genepop_format_mixtus2023.txt","data/genepop/mixtusLD2023.txt.DIS")


######################################
# Auto-check LD results
#####################################

LD2022lines = readLines("data/genepop/mixtusLD2022.txt.DIS")
split_linesLD2022 <- strsplit(LD2022lines, "\\s+")
LD2022 = do.call(rbind, split_linesLD2022[13:482])
colnames(LD2022) = LD2022[1,]
colnames(LD2022)[colnames(LD2022) == 'P-Value'] = "pval"
LD2022 = as.data.frame(LD2022[3:nrow(LD2022),])

LD2023lines = readLines("data/genepop/mixtusLD2023.txt.DIS")
split_linesLD2023 <- strsplit(LD2023lines, "\\s+")
LD2023 = do.call(rbind, split_linesLD2023[13:482])
colnames(LD2023) = LD2023[1,]
colnames(LD2023)[colnames(LD2023) == 'P-Value'] = "pval"
LD2023 = as.data.frame(LD2023[3:nrow(LD2023),])

# check out loci that are below Bonferroni-corrected P-value
alpha_mix = 5.3*10^-5

LD2022_failed = LD2022[as.numeric(LD2022$pval) < alpha_mix,]
LD2023_failed = LD2023[as.numeric(LD2023$pval) < alpha_mix, ]


##########################################
# Re-run COLONY with corrected locus sets
##########################################

##### prep genotype data

#remove all columns except barcode and scores
# also discard loci with high errors or out of HWE/LD
colsToRemove = c("year", 
                 "final_id", 
                 "notes", 
                 "sample_id", 
                 "date", 
                 "julian_date", 
                 "BL15...1", 
                 "BL15...2", 
                 "BTMS0072...1", 
                 "BTMS0072...2",
                 "BTERN01...1",
                 "BTERN01...2")
mixtus2022 = mixtus2022[, !colnames(mixtus2022) %in% colsToRemove]
mixtus2023 = mixtus2023[, !colnames(mixtus2023) %in% colsToRemove]

#relocate barcode, remove rows with more than 8 NAs, replace NAs with 0's
# changing this to let some queens that were sequenced on impatiens plate pass the NA filter...
mixtus2022 = mixtus2022 %>% relocate(barcode_id)
mixtus2022_forcolony = mixtus2022[rowSums(is.na(mixtus2022)) < 12,]
mixtus2022_forcolony[is.na(mixtus2022_forcolony)] = 0

mixtus2023 = mixtus2023 %>% relocate(barcode_id)
mixtus2023_forcolony = mixtus2023[rowSums(is.na(mixtus2023)) < 12,]
mixtus2023_forcolony[is.na(mixtus2023_forcolony)] = 0

#write csvs for upload to colony
column_names = colnames(mixtus2022_forcolony)

write.table(mixtus2022_forcolony, "data/merged_by_year/sib_scores_for_colony/mixtus2022_forcolony_finalloci.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(mixtus2023_forcolony, "data/merged_by_year/sib_scores_for_colony/mixtus2023_forcolony_finalloci.txt", sep= ",", col.names = FALSE, row.names = FALSE)


write.csv(mixtus2022_forcolony, "data/merged_by_year/csvs/mixtus_2022_scores_finalloci.csv")
write.csv(mixtus2023_forcolony, "data/merged_by_year/csvs/mixtus_2023_scores_finalloci.csv")


##### prep microsat error rates

# first value per column: marker name
# second value per column: marker type (codominant for microsats -> 0)
# third value per column: allelic dropout rate ?? (set to 0)
# fourth value per column: error/mutation rate (empirically derived, see check_microsat_error_rates.R)
mixtus_error_rates = data.frame(c("BT10", 0, 0, 0.01),
                                c("BTMS0104", 0, 0, 0.01),
                                c("BTMS0057", 0, 0, 0.01),
                                c("BTMS0086", 0, 0, 0.015),
                                c("BTMS0066", 0, 0, 0.01),
                                c("BTMS0062", 0, 0, 0.015),
                                c("BTMS0136", 0, 0, 0.01),
                                c("BTMS0126", 0, 0, 0.0162),
                                c("BTMS0059", 0, 0, 0.01),
                                c("BL13", 0, 0, 0.0159),
                                c("BTMS0083", 0, 0, 0.01),
                                c("B126", 0, 0, 0.01))
write.table(mixtus_error_rates, "data/merged_by_year/error_rates/mixtus_error_rates_finalloci.txt", sep= ",", col.names = FALSE, row.names = FALSE)

###### No sibship exclusion this time...

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
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/final", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/final/mixtus2022_1.DAT", delim=",")
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/final", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/final/mixtus2023_1.DAT", delim=",")

# navigate to COLONY2 sub folder and run the following in terminal
#NOTE: ensure that mixtus2022.DAT and mixtus2023.DAT are in the same directory at the colony2s.out executable

# cd /Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2
# ./colony2s.out IFN:final/mixtus2022_1.DAT
# ./colony2s.out IFN:final/mixtus2023_1.DAT



#####################################
## Load in results from colony UPDATE FROM HERE
#####################################

# read in and format 2022 data
lines2022 <- trimws(readLines("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/final/mixtus2022_1.BestCluster"))
split_lines2022 <- strsplit(lines2022, "\\s+")
bestconfig2022 <- do.call(rbind, split_lines2022)
colnames(bestconfig2022) <- bestconfig2022[1, ]
bestconfig2022 <- as.data.frame(bestconfig2022[-1, ])
bestconfig2022$Probability = as.numeric(bestconfig2022$Probability)

# read in and format 2023 data
lines2023 <- trimws(readLines("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/final/mixtus2023_1.BestCluster"))
split_lines2023 <- strsplit(lines2023, "\\s+")
bestconfig2023 <- do.call(rbind, split_lines2023)
colnames(bestconfig2023) <- bestconfig2023[1, ]
bestconfig2023 <- as.data.frame(bestconfig2023[-1, ])
bestconfig2023$Probability = as.numeric(bestconfig2023$Probability)

# remove sibling relationships with p < 0.95
counter = max(as.numeric(bestconfig2022$ClusterIndex)) + 1
for (i in 1:nrow(bestconfig2022)){
  # if the probability of inclusion in the cluster is < 0.95, make a new cluster
  if (bestconfig2022$Probability[i] < 0.95){
    bestconfig2022$ClusterIndex[i] = counter
    counter = counter + 1
  }
}

# repeat for 2023
counter = max(as.numeric(bestconfig2023$ClusterIndex)) + 1
for (i in 1:nrow(bestconfig2023)){
  # if the probability of inclusion in the cluster is < 0.95, make a new cluster
  if (bestconfig2023$Probability[i] < 0.95){
    bestconfig2023$ClusterIndex[i] = counter
    counter = counter + 1
  }
}


#write csvs to file
write.csv(bestconfig2022, "data/siblingships/mix_sibships_final_2022.csv")
write.csv(bestconfig2023, "data/siblingships/mix_sibships_final_2023.csv")
