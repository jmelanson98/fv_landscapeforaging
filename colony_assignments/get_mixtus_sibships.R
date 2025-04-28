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

#order columns so that we can rbind originals to reruns
correctorderp1 = colnames(alleles_plex1_wbarcode)
correctorderp2 = colnames(alleles_plex2_wbarcode)

alleles_plex1_rerun <- alleles_plex1_rerun[, correctorderp1]
alleles_plex2_rerun <- alleles_plex2_rerun[, correctorderp2]

alleles_plex1 = rbind(alleles_plex1_wbarcode, alleles_plex1_rerun)
alleles_plex2 = rbind(alleles_plex2_wbarcode, alleles_plex2_rerun)

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

#remove duplicated rows
plex1_unique <- alleles_plex1[!duplicated(alleles_plex1$barcode_id), ]
plex2_unique <- alleles_plex2[!duplicated(alleles_plex2$barcode_id), ]

#merge plexes
microsat_scores = full_join(plex1_unique, plex2_unique, by = "barcode_id")

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

write.table(mixtus2022_forcolony, "data/merged_by_year/mixtus2022_forcolony.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(mixtus2023_forcolony, "data/merged_by_year/mixtus2023_forcolony.txt", sep= ",", col.names = FALSE, row.names = FALSE)


write.csv(mixtus2022_forcolony, "data/merged_by_year/mixtus_2022_scores.csv")
write.csv(mixtus2023_forcolony, "data/merged_by_year/mixtus_2023_scores.csv")


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
write.table(mixtus_error_rates, "data/merged_by_year/mixtus_error_rates.txt", sep= ",", col.names = FALSE, row.names = FALSE)

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
  file = "data/merged_by_year/mixtus_sibexclusions_2022.txt",
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  na = ""
)

write.table(
  sib2023_reduced,
  file = "data/merged_by_year/mixtus_sibexclusions_2023.txt",
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
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/mixtus2022.DAT", delim=",")
rcolony::build.colony.input(wd="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2", 
                            name="/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/mixtus2023.DAT", delim=",")

# navigate to COLONY2 sub folder and run the following in terminal
#NOTE: ensure that mixtus2022.DAT and mixtus2023.DAT are in the same directory at the colony2s.out executable

# cd /Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2
# ./colony2s.out IFN=mixtus2022.DAT
# ./colony2s.out IFN=mixtus2023.DAT



#####################################
## Load in results from colony
#####################################

# read in and format 2022 data
lines2022 <- trimws(readLines("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/mixtus2022.BestCluster"))
split_lines2022 <- strsplit(lines2022, "\\s+")
bestconfig2022 <- do.call(rbind, split_lines2022)
colnames(bestconfig2022) <- bestconfig2022[1, ]
bestconfig2022 <- as.data.frame(bestconfig2022[-1, ])

# read in and format 2023 data
lines2023 <- trimws(readLines("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging/colony_assignments/Colony2/mixtus2023.BestCluster"))
split_lines2023 <- strsplit(lines2023, "\\s+")
bestconfig2023 <- do.call(rbind, split_lines2023)
colnames(bestconfig2023) <- bestconfig2023[1, ]
bestconfig2023 <- as.data.frame(bestconfig2023[-1, ])

# remove sibling relationships with p < 0.95
counter = max(as.numeric(bestconfig2022$ClusterIndex)) + 1
for (i in 1:length(nrow(bestconfig2022))){
  # if the probability of inclusion in the cluster is < 0.95, make a new cluster
  if (bestconfig2022$Probability[i] < 0.95){
    bestconfig2022$ClusterIndex[i] = counter
    counter = counter + 1
  }
}

# repeat for 2023
counter = max(as.numeric(bestconfig2023$ClusterIndex)) + 1
for (i in 1:length(nrow(bestconfig2023))){
  # if the probability of inclusion in the cluster is < 0.95, make a new cluster
  if (bestconfig2023$Probability[i] < 0.95){
    bestconfig2023$ClusterIndex[i] = counter
    counter = counter + 1
  }
}


#write csvs to file
write.csv(bestconfig2022, "data/siblingships/mix_sibships_preliminary_2022.csv")
write.csv(bestconfig2023, "data/siblingships/mix_sibships_preliminary_2023.csv")

############################################
# Remove siblings before running genepop
############################################
nonsibs2022 = left_join(specimenData2022, bestconfig2022, by = c("barcode_id" ="OffspringID")) %>%
  filter(!is.na(ClusterIndex)) %>%
  group_by(ClusterIndex) %>%
  sample_n(1)

testalleles2022 = mixtus2022_forcolony[mixtus2022_forcolony$barcode_id %in% nonsibs2022$barcode_id,]

######################################
# Prep alleles table for genepop
#####################################

#combine columns for each locus
ms_new = data.table(testalleles2022)
ms_new[ms_new ==0] <- "000"


i1 <- seq(1, length(ms_new)-1, 2)
i2 <- seq(2, length(ms_new)-1, 2)
msgenepop = ms_new[, Map(paste,
                         .SD[, i1, with = FALSE], .SD[, i2, with = FALSE], 
                         MoreArgs = list(sep="")), 
                   by = "barcode_id"]
colnames(msgenepop) = c("barcode_id", "BT10", "BTMS0104", "BTMS0057", "BTMS0086",
                        "BTMS0066", "BTMS0062", "BTMS0136","BTERN01", "BTMS0126", 
                        "BTMS0059", "BL13", "BTMS0083", "B126")


#select and join alleles

#sort allele table by barcode ID
msgenepop <- msgenepop[order(msgenepop$barcode_id),]
msgenepop$site = substring(msgenepop$barcode_id, 1, 1)

haplotypes <- as.data.frame(paste(msgenepop$site, "_", msgenepop$barcode_id, ","," ", 
                                  msgenepop$BT10," ", msgenepop$BTMS0104, " ",
                                  msgenepop$BTMS0057," ", msgenepop$BTMS0086, " ",
                                  msgenepop$BTMS0066," ", msgenepop$BTMS0062, " ",
                                  msgenepop$BTMS0136," ", msgenepop$BTERN01, " ",
                                  msgenepop$BTMS0126," ", msgenepop$BTMS0059," ",
                                  msgenepop$BL13," ", msgenepop$BTMS0083, " ",
                                  msgenepop$B126,
                                  sep = ""))

#Build the Genepop format
sink("data/genepop/ms_genepop_format_2022.txt")
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
invisible(apply(as.data.frame(haplotypes[substring(haplotypes[,1],1,1) == "E",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes[substring(haplotypes[,1],1,1) == "W",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes[substring(haplotypes[,1],1,1) == "S",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes[substring(haplotypes[,1],1,1) == "N",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes[substring(haplotypes[,1],1,1) == "P",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
invisible(apply(as.data.frame(haplotypes[substring(haplotypes[,1],1,1) == "H",]), 1,function(x) cat(x,"\n")))
cat("Pop \n")
sink()



######################################
# Run genepop on samples
#####################################

#test for HWE
test_HW("data/genepop/ms_genepop_format_2022.txt", which="Proba", "data/genepop/outputHWE2022.txt.D")

#test for LD
test_LD("data/genepop/ms_genepop_format_2022.txt","data/genepop/outputLD2022.txt.DIS")


# for 2022 samples
### BTERN01 Fis values are Fis-hy (haha I'm so funny)
### so is BL13...
### so is BTMS0083

### LD: BTMS0057 and BTMS0066
####### BTMS0126 and BTMS0059
####### BTMS0136 and BL13
####### BTMS0126 and B126
####### BTMS0066 and BTMS0136 *** across all pops?
####### BTMS0066 and BTMS0059













#######################
# some of this code is not for getting sibships...
#######################

sibships2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/mixtus_sibships.csv", sep = ",", header = T))
sibships2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/mixtus2023siblings.csv", sep = ",", header = T))

sib_long_filtered2022 <- sibships2022 %>% 
  pivot_longer(
    cols = `sib_1`:`sib_4`, 
    names_to = "sib",
    values_to = "barcode_id") %>%
  filter(barcode_id != "",
         prob_inc > 0.95)

sibships2023$fullsibshipindex = sibships2023$fullsibshipindex + max(sibships2022$fullsib_index)
sib_long_filtered2023 <- sibships2023 %>% 
  pivot_longer(
    cols = `sib1`:`sib5`, 
    names_to = "sib",
    values_to = "barcode_id") %>%
  filter(barcode_id != "",
         prob_inc > 0.95)

#combine 2022 and 2023 sibs
colnames(sib_long_filtered2023) = colnames(sib_long_filtered2022)
sib_long_filtered2023$barcode_id = gsub("_", "", sib_long_filtered2023$barcode_id)
allsibs = rbind(sib_long_filtered2022, sib_long_filtered2023)

#clean up mixtus specimen data frame
col_to_keep = c("barcode_id", "final_id", "year", "site", "round", "sample_pt", "sample_id", "active_flower", "notes", "pollen")
all_specimens = rbind(specimenData2022[,colnames(specimenData2022) %in% col_to_keep],
                                   specimenData2023[,colnames(specimenData2023) %in% col_to_keep])
all_mixtus = filter(all_specimens, final_id == "B. mixtus" & barcode_id != "N/A" & barcode_id != "")

#join on compressed barcode, then delete compressed barcode
all_mixtus$compressedbarcode = gsub("_", "", all_mixtus$barcode_id)
mixtus_wsibs = full_join(all_mixtus, allsibs, by= c("compressedbarcode" = "barcode_id"))
mixtus_wsibs$grouptime = paste(mixtus_wsibs$fullsib_index, mixtus_wsibs$round, sep="_")

search_list_sametime_sameplace = filter(mixtus_wsibs, !is.na(fullsib_index)) %>%
  group_by(grouptime) %>%
  filter(n() >= 2) %>%
  filter(n_distinct(sample_pt) == 1) %>%
  slice_sample(n =2)

search_list_sametime_differentplace = filter(mixtus_wsibs, !is.na(fullsib_index)) %>%
  group_by(grouptime) %>%
  filter(n() >= 2) %>%
  filter(n_distinct(sample_pt) > 1) %>%
  slice_sample(n =2)

search_list_differenttime_sampeplace = mixtus_wsibs[!mixtus_wsibs$barcode_id 
                                     %in% search_list_sametime$barcode_id,] %>%
  filter(!is.na(fullsib_index)) %>%
  group_by(fullsib_index) %>%
  filter(n() >= 2) %>%
  filter(n_distinct(round) > 1) %>%
  filter(n_distinct(sample_pt) == 1)

search_list_differenttime_differentplace = mixtus_wsibs[!mixtus_wsibs$barcode_id 
                                                    %in% search_list_sametime$barcode_id,] %>%
  filter(!is.na(fullsib_index)) %>%
  group_by(fullsib_index) %>%
  filter(n() >= 2) %>%
  filter(n_distinct(round) > 1) %>%
  filter(n_distinct(sample_pt) > 1)
  
  

search_list2023 = search_list[search_list$year == 2023,]
search_list %>% 
  group_by(fullsib_index, site) %>%
  summarize(num_locations = n_distinct(sample_pt)) -> summarybypt


length(summarybypt$num_locations[summarybypt$num_locations > 1])
  

unique(search_list2023$sample_pt)
search_list2023 %>%
mutate(period = case_when(round < 20 ~ "before",
                          round > 19 ~ "after")) %>%
  relocate(period, .after = round) %>%
  group_by(sample_pt) %>%
  summarize(n=n_distinct(period)) -> locationsummary

#sampling locations over time
all_mixtus = filter(all_specimens, final_id == "B. mixtus" & 
                      barcode_id != "N/A" & 
                      barcode_id != "" & 
                      notes != "queen" & 
                      notes != "male" &
                      year == 2023)

all_mixtus_period = all_mixtus %>% mutate(period = case_when(round < 16 ~ 1,
                          round > 15 & round < 20 ~ 2,
                          round > 19 & round < 24 ~ 3,
                          round > 23 ~ 4)) %>%
  relocate(period, .after = round)

summary_data <- all_mixtus_period %>%
  group_by(sample_pt, period) %>%
  summarise(num_observations = n(), .groups = 'drop')
complete_data <- summary_data %>%
  complete(sample_pt, period = 1:4, fill = list(num_observations = 0))
complete_data = complete_data[complete_data$period != 1,]
filtered_data <- complete_data %>%
  group_by(sample_pt) %>%
  filter(all(num_observations > 0)) %>%
  ungroup()


ggplot(filtered_data, aes(x = period, y = num_observations, color = sample_pt, group = sample_pt)) +
  geom_line() +                      # Use lines to connect points
  geom_point() +                     # Add points for each observation
  labs(title = "B. mixtus Observations Across Rounds",
       x = "Round",
       y = "Number of Observations") +
  theme_minimal()



write.csv(search_list, "mixtussibs_tofind.csv")

search_list %>%
  group_by(fullsibshipindex) %>%
  summarize(pollen_sum = sum(pollen != "n")) -> out

search_list %>%
  group_by(fullsibshipindex) %>%
  summarize(num_sibs = n()) -> sibgroups


######################
### make sib maps
#####################

#add lat and long to samplePoints
#first, wrangle gps coordinates into shape
samplePoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) tail(x, 1))
samplePoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) head(x, 3)[2])
samplePoints = samplePoints[,c("sample_pt", "subsite", "lat", "long")]

#join with siblings
sibgroups = mixtus_wsibs %>% filter(!is.na(sample_pt)) %>% filter(!is.na(fullsib_index))
sibswcoords = merge(sibgroups, samplePoints, by = "sample_pt", all.x = TRUE)

#split 2022 and 2023
sibs22 = sibswcoords %>% filter(year == "2022")
sibs23 = sibswcoords %>% filter(year=="2023")

#make plot grids for both
source("microsatellitecode/src/makeSibMaps.R")
grids2022 = makeSibMaps(sibs22)
length(unique(sibs22$fullsib_index))

grids2023 = makeSibMaps(sibs23)
length(unique(sibs23$fullsib_index))
