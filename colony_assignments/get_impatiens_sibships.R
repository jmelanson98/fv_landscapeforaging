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
library(dplyr)
library(tidyr)
library(stringr)
library(rcolony)
library(genepop)
library(data.table)
library(cowplot)
library(raster)

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
  filter(year == "2022" | (str_detect(notes, "queen") & julian_date < 175))
impatiens2023 = filter(specimens_withdates, final_id == "B. impatiens") %>%
  filter(year == "2023") %>%
  filter((!str_detect(notes, "queen")) | julian_date > 175)
                          

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
                 "BL15...2")
impatiens2022 = impatiens2022[, !colnames(impatiens2022) %in% colsToRemove]
impatiens2023 = impatiens2023[, !colnames(impatiens2023) %in% colsToRemove]

#relocate barcode, remove rows with more than 10 NAs, replace NAs with 0's
impatiens2022 = impatiens2022 %>% relocate(barcode_id)
impatiens2022_forcolony = impatiens2022[rowSums(is.na(impatiens2022)) < 14,]
impatiens2022_forcolony[is.na(impatiens2022_forcolony)] = 0

impatiens2023 = impatiens2023 %>% relocate(barcode_id)
impatiens2023_forcolony = impatiens2023[rowSums(is.na(impatiens2023)) < 10,]
impatiens2023_forcolony[is.na(impatiens2023_forcolony)] = 0

#write csvs for upload to colony
column_names = colnames(impatiens2022_forcolony)

write.table(impatiens2022_forcolony, "data/merged_by_year/impatiens2022_forcolony.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(impatiens2023_forcolony, "data/merged_by_year/impatiens2023_forcolony.txt", sep= ",", col.names = FALSE, row.names = FALSE)


write.csv(impatiens2022_forcolony, "data/merged_by_year/impatiens_2022_scores.csv")
write.csv(impatiens2023_forcolony, "data/merged_by_year/impatiens_2023_scores.csv")


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
                                c("BL13", 0, 0, 0.011),
                                c("BTMS0062", 0, 0, 0.022),
                                c("B126", 0, 0, 0.01),
                                c("BTERN01", 0, 0, 0.028),
                                c("B124", 0, 0, 0.011),
                                c("BTMS0057", 0, 0, 0.017),
                                c("BT30", 0, 0, 0.01),
                                c("B10", 0, 0, 0.029),
                                c("BTMS0083", 0, 0, 0.01),
                                c("BTMS0073", 0, 0, 0.01),
                                c("BT28", 0, 0, 0.01)
                                )
write.table(impatiens_error_rates, "data/merged_by_year/impatiens_error_rates.txt", sep= ",", col.names = FALSE, row.names = FALSE)

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
  file = "data/merged_by_year/impatiens_sibexclusions_2022.txt",
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  na = ""
)

write.table(
  sib2023_reduced,
  file = "data/merged_by_year/impatiens_sibexclusions_2023.txt",
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
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/impatiens2022.DAT", delim=",")
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/impatiens2023.DAT", delim=",")

# navigate to COLONY2 sub folder and run the following in terminal
#NOTE: ensure that mixtus2022.DAT and mixtus2023.DAT are in the same directory at the colony2s.out executable

# cd /Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2
# ./colony2s.out IFN:impatiens2022.DAT
# ./colony2s.out IFN:mixtusimpatiens.DAT



#####################################
## Load in results from colony
#####################################

# read in and format 2022 data
lines2022 <- trimws(readLines("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/impatiens2022.BestCluster"))
split_lines2022 <- strsplit(lines2022, "\\s+")
bestconfig2022 <- do.call(rbind, split_lines2022)
colnames(bestconfig2022) <- bestconfig2022[1, ]
bestconfig2022 <- as.data.frame(bestconfig2022[-1, ])
bestconfig2022$Probability = as.numeric(bestconfig2022$Probability)

# read in and format 2023 data
lines2023 <- trimws(readLines("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/impatiens2023.BestCluster"))
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
write.csv(bestconfig2022, "data/siblingships/imp_sibships_preliminary_2022.csv")
write.csv(bestconfig2023, "data/siblingships/imp_sibships_preliminary_2023.csv")

############################################
# Remove siblings before running genepop
############################################

# for 2022, join specimen dataframes from both years to include queens
ctkeep = c("site", "round", "sample_pt", "sample_id", "year", "barcode_id", "active_flower", "final_id", "notes", "plate", "well")
full = rbind(specimenData2022[,colnames(specimenData2022) %in% ctkeep],
             specimenData2023[,colnames(specimenData2023) %in% ctkeep])
nonsibs2022 = left_join(full, bestconfig2022, by = c("barcode_id" ="OffspringID")) %>%
  filter(!is.na(ClusterIndex)) %>%
  group_by(ClusterIndex) %>%
  sample_n(1)

testalleles2022 = impatiens2022_forcolony[impatiens2022_forcolony$barcode_id %in% nonsibs2022$barcode_id,]

# repeat for 2023
nonsibs2023 = left_join(specimenData2023, bestconfig2023, by = c("barcode_id" ="OffspringID")) %>%
  filter(!is.na(ClusterIndex)) %>%
  group_by(ClusterIndex) %>%
  sample_n(1)

testalleles2023 = impatiens2023_forcolony[impatiens2023_forcolony$barcode_id %in% nonsibs2023$barcode_id,]

######################################
# Prep alleles table for genepop
#####################################

### first for 2022! 
#combine columns for each locus
ms_i22 = data.table(testalleles2022)
ms_i22[ms_i22 ==0] <- "000"


i1 <- seq(1, length(ms_i22)-1, 2)
i2 <- seq(2, length(ms_i22)-1, 2)
msgenepop_i22 = ms_i22[, Map(paste,
                         .SD[, i1, with = FALSE], .SD[, i2, with = FALSE], 
                         MoreArgs = list(sep="")), 
                   by = "barcode_id"]
colnames(msgenepop_i22) = c("barcode_id", "BT10", "B96", "BTMS0059", "BTMS0081",
                        "BL13", "BTMS0062", "B126", "BTERN01","B124", "BTMS0057", 
                        "BT30", "B10", "BTMS0083", "BTMS0073", "BT28")


#select and join alleles

#sort allele table by barcode ID
msgenepop_i22 <- msgenepop_i22[order(msgenepop_i22$barcode_id),]
msgenepop_i22$site = substring(msgenepop_i22$barcode_id, 1, 1)

haplotypes_i22 <- as.data.frame(paste(msgenepop_i22$site, "_", msgenepop_i22$barcode_id, ","," ", 
                                      msgenepop_i22$BT10," ", msgenepop_i22$B96, " ",
                                      msgenepop_i22$BTMS0059," ", msgenepop_i22$BTMS0081, " ",
                                      msgenepop_i22$BL13," ", msgenepop_i22$BTMS0062, " ",
                                      msgenepop_i22$B126," ", msgenepop_i22$BTERN01, " ",
                                      msgenepop_i22$B124," ", msgenepop_i22$BTMS0057," ",
                                      msgenepop_i22$BT30," ", msgenepop_i22$B10, " ",
                                      msgenepop_i22$BTMS0083," ", msgenepop_i22$BTMS0073," ",
                                      msgenepop_i22$BT28,
                                  sep = ""))

#Build the Genepop format
sink("data/genepop/ms_genepop_format_impatiens2022.txt")
cat("Title: B. impatiens workers (2022) \n")
cat("BT10 \n")
cat("B96 \n")
cat("BTMS0059 \n")
cat("BTMS0081 \n")
cat("BL13 \n")
cat("BTMS0062 \n")
cat("B126 \n")
cat("BTERN01 \n")
cat("B124 \n")
cat("BTMS0057 \n")
cat("BT30 \n")
cat("B10 \n")
cat("BTMS0083 \n")
cat("BTMS0073 \n")
cat("BT28 \n")
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_i22[substring(haplotypes_i22[,1],1,1) == "E",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_i22[substring(haplotypes_i22[,1],1,1) == "W",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_i22[substring(haplotypes_i22[,1],1,1) == "S",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_i22[substring(haplotypes_i22[,1],1,1) == "N",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_i22[substring(haplotypes_i22[,1],1,1) == "P",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_i22[substring(haplotypes_i22[,1],1,1) == "H",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
sink()


# then again, for 2023!
#combine columns for each locus
ms_i23 = data.table(testalleles2023)
ms_i23[ms_i23 ==0] <- "000"


i1 <- seq(1, length(ms_i23)-1, 2)
i2 <- seq(2, length(ms_i23)-1, 2)
msgenepop_i23 = ms_i23[, Map(paste,
                             .SD[, i1, with = FALSE], .SD[, i2, with = FALSE], 
                             MoreArgs = list(sep="")), 
                       by = "barcode_id"]
colnames(msgenepop_i23) = c("barcode_id", "BT10", "B96", "BTMS0059", "BTMS0081",
                            "BL13", "BTMS0062", "B126", "BTERN01","B124", "BTMS0057", 
                            "BT30", "B10", "BTMS0083", "BTMS0073", "BT28")


#select and join alleles

#sort allele table by barcode ID
msgenepop_i23 <- msgenepop_i23[order(msgenepop_i23$barcode_id),]
msgenepop_i23$site = substring(msgenepop_i23$barcode_id, 1, 1)

haplotypes_i23 <- as.data.frame(paste(msgenepop_i23$site, "_", msgenepop_i23$barcode_id, ","," ", 
                                      msgenepop_i23$BT10," ", msgenepop_i23$B96, " ",
                                      msgenepop_i23$BTMS0059," ", msgenepop_i23$BTMS0081, " ",
                                      msgenepop_i23$BL13," ", msgenepop_i23$BTMS0062, " ",
                                      msgenepop_i23$B126," ", msgenepop_i23$BTERN01, " ",
                                      msgenepop_i23$B124," ", msgenepop_i23$BTMS0057," ",
                                      msgenepop_i23$BT30," ", msgenepop_i23$B10, " ",
                                      msgenepop_i23$BTMS0083," ", msgenepop_i23$BTMS0073," ",
                                      msgenepop_i23$BT28,
                                      sep = ""))

#Build the Genepop format
sink("data/genepop/ms_genepop_format_impatiens2023.txt")
cat("Title: B. impatiens workers (2023) \n")
cat("BT10 \n")
cat("B96 \n")
cat("BTMS0059 \n")
cat("BTMS0081 \n")
cat("BL13 \n")
cat("BTMS0062 \n")
cat("B126 \n")
cat("BTERN01 \n")
cat("B124 \n")
cat("BTMS0057 \n")
cat("BT30 \n")
cat("B10 \n")
cat("BTMS0083 \n")
cat("BTMS0073 \n")
cat("BT28 \n")
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_i23[substring(haplotypes_i23[,1],1,1) == "E",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_i23[substring(haplotypes_i23[,1],1,1) == "W",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_i23[substring(haplotypes_i23[,1],1,1) == "S",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_i23[substring(haplotypes_i23[,1],1,1) == "N",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_i23[substring(haplotypes_i23[,1],1,1) == "P",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes_i23[substring(haplotypes_i23[,1],1,1) == "H",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
sink()

######################################
# Run genepop on samples
#####################################

#test for HWE
test_HW("data/genepop/ms_genepop_format_impatiens2022.txt", which="Proba", "data/genepop/impatiensHWE2022.txt.D")
test_HW("data/genepop/ms_genepop_format_impatiens2023.txt", which="Proba", "data/genepop/impatiensHWE2023.txt.D")

#test for LD
test_LD("data/genepop/ms_genepop_format_impatiens2022.txt","data/genepop/impatiensLD2022.txt.DIS")
test_LD("data/genepop/ms_genepop_format_impatiens2023.txt","data/genepop/impatiensLD2023.txt.DIS")


######################################
# Auto-check LD results
#####################################

LD2022lines = readLines("data/genepop/impatiensLD2022.txt.DIS")
split_linesLD2022 <- strsplit(LD2022lines, "\\s+")
LD2022 = do.call(rbind, split_linesLD2022[13:482])
colnames(LD2022) = LD2022[1,]
colnames(LD2022)[colnames(LD2022) == 'P-Value'] = "pval"
LD2022 = as.data.frame(LD2022[3:nrow(LD2022),])

LD2023lines = readLines("data/genepop/impatiensLD2023.txt.DIS")
split_linesLD2023 <- strsplit(LD2023lines, "\\s+")
LD2023 = do.call(rbind, split_linesLD2023[13:482])
colnames(LD2023) = LD2023[1,]
colnames(LD2023)[colnames(LD2023) == 'P-Value'] = "pval"
LD2023 = as.data.frame(LD2023[3:nrow(LD2023),])

# check out loci that are below Bonferroni-corrected P-value
alpha_imp = 3.9*10^-5

LD2022_failed = LD2022[as.numeric(LD2022$pval) < alpha_imp,]
LD2023_failed = LD2023[as.numeric(LD2023$pval) < alpha_imp, ]
