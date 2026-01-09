################################################################################
## Get colony IDs from microsatellite genotyping results--Bombus impatiens
## Started by J Melanson
## February 13, 2025
################################################################################

# prep workspace
rm(list = ls())
setwd("/Users/jenna1/fv_landscapeforaging")
workingdir = getwd()

# first, load in packages
source('src/init.R')
source('src/colony_assignment_functions.R')
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
impatiens2023wq = filter(specimens_withdates, final_id == "B. impatiens") %>%
  filter(year == "2023") %>%
  filter((str_detect(notes, "queen")) & julian_date < 175) %>%
  filter(!str_detect(notes, "male"))

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
impatiens2023wq = impatiens2023wq[, !colnames(impatiens2023wq) %in% colsToRemove]

#relocate barcode, remove rows with more than 6 NAs (e.g., at least 8 loci scored), replace NAs with 0's
impatiens2022 = impatiens2022 %>% relocate(barcode_id)
impatiens2022_forcolony = impatiens2022[rowSums(is.na(impatiens2022)) <= 6,]
impatiens2022_forcolony[is.na(impatiens2022_forcolony)] = 0

impatiens2023 = impatiens2023 %>% relocate(barcode_id)
impatiens2023_forcolony = impatiens2023[rowSums(is.na(impatiens2023)) <= 6,]
impatiens2023_forcolony[is.na(impatiens2023_forcolony)] = 0

impatiens2023wq = impatiens2023wq %>% relocate(barcode_id)
impatiens2023wq_forcolony = impatiens2023wq[rowSums(is.na(impatiens2023wq)) <= 6,]
impatiens2023wq_forcolony[is.na(impatiens2023wq_forcolony)] = 0

#write csvs for upload to colony
column_names = colnames(impatiens2022_forcolony)

write.table(impatiens2022_forcolony, "data/merged_by_year/sib_scores_for_colony/impatiens2022_forcolony.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(impatiens2023_forcolony, "data/merged_by_year/sib_scores_for_colony/impatiens2023_forcolony.txt", sep= ",", col.names = FALSE, row.names = FALSE)
write.table(impatiens2023wq_forcolony, "data/merged_by_year/sib_scores_for_colony/impatiens2023queens_forcolony.txt", sep= ",", col.names = FALSE, row.names = FALSE)

write.csv(impatiens2022_forcolony, "data/merged_by_year/csvs/impatiens_2022_scores.csv")
write.csv(impatiens2023_forcolony, "data/merged_by_year/csvs/impatiens_2023_scores.csv")
write.csv(impatiens2023wq_forcolony, "data/merged_by_year/csvs/impatiens_2023wq_scores.csv")


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
excluded_sibships_2023wq = list()
sample.names.2022 = impatiens2022_forcolony[,1]
sample.names.2023 = impatiens2023_forcolony[,1]
sample.names.2023wq = impatiens2023wq_forcolony[,1]

for (i in 1:6){
  site = sites[i]
  
  # grep names for focal site
  site.names.2022 = sample.names.2022[grep(site, sample.names.2022)]
  site.names.2023 = sample.names.2023[grep(site, sample.names.2023)]
  site.names.2023wq = sample.names.2023wq[grep(site, sample.names.2023wq)]
  
  # list of excluded siblings
  excluded.2022 = sample.names.2022[sample.names.2022 %!in% site.names.2022]
  numex_2022 = length(excluded.2022)
  excluded.2023 = sample.names.2023[sample.names.2023 %!in% site.names.2023]
  numex_2023 = length(excluded.2023)
  excluded.2023wq = sample.names.2023wq[sample.names.2023wq %!in% site.names.2023wq]
  numex_2023wq = length(excluded.2023wq)
  
  # make a matrix of excluded samples
  excluded.matrix.2022 <- matrix(rep(excluded.2022, times = length(site.names.2022)),
                                 nrow = length(site.names.2022), byrow = TRUE)
  colnames(excluded.matrix.2022) <- paste0("excluded_", seq_len(ncol(excluded.matrix.2022)))
  
  excluded.matrix.2023 <- matrix(rep(excluded.2023, times = length(site.names.2023)),
                                 nrow = length(site.names.2023), byrow = TRUE)
  colnames(excluded.matrix.2023) <- paste0("excluded_", seq_len(ncol(excluded.matrix.2023)))
  
  excluded.matrix.2023wq <- matrix(rep(excluded.2023wq, times = length(site.names.2023)),
                                 nrow = length(site.names.2023), byrow = TRUE)
  colnames(excluded.matrix.2023wq) <- paste0("excluded_", seq_len(ncol(excluded.matrix.2023wq)))
  
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
  
  sitedf2023wq = data.frame(focal = site.names.2023,
                          num_exc = numex_2023wq,
                          excluded.matrix.2023wq,
                          stringsAsFactors = FALSE
  )
  
  # save exclusion dfs to list
  excluded_sibships_2022[[i]] = sitedf2022
  excluded_sibships_2023[[i]] = sitedf2023
  excluded_sibships_2023wq[[i]] = sitedf2023wq
}

# make exclusion tables!
sibexclusions_2022 = bind_rows(excluded_sibships_2022)
sibexclusions_2023 = bind_rows(excluded_sibships_2023)
sibexclusions_2023wq = bind_rows(excluded_sibships_2023wq)

sib2022_reduced = sibexclusions_2022[,!colnames(sibexclusions_2022) %in% c("num_exc")]
sib2023_reduced = sibexclusions_2023[,!colnames(sibexclusions_2023) %in% c("num_exc")]
sib2023wq_reduced = sibexclusions_2023wq[,!colnames(sibexclusions_2023wq) %in% c("num_exc")]

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

write.table(
  sib2023wq_reduced,
  file = "data/merged_by_year/sib_exclusions/impatiens_maternalexclusions.txt",
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
# I ran this code using the Linux version of COLONY, on the AllianceCan server

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
# - Addition of exclusion threshold (0) for known paternity/maternity in cases where the number of known parentages is 0

# I also made some updated functions which are a bit more "automatic" so that I don't have to run through all of the prompts. See below.


# build .DAT files for both datasets
rcolony::build.colony.superauto(wd=workingdir, 
                                name=paste0(workingdir, "/colony_assignments/Colony2_Linux/impatiens_2022.DAT"), 
                                datasetname = "impatiens2022",
                                delim=",",
                                sample_size = 1362,
                                num_loci = 12,
                                sibship_prior = 0,
                                female_monogamy = 1,
                                error_rates_path = paste0(workingdir, "/data/merged_by_year/error_rates/impatiens_error_rates.txt"),
                                genotypes_path = paste0(workingdir, "/data/merged_by_year/sib_scores_for_colony/impatiens2022_forcolony.txt"),
                                exclusion_path = paste0(workingdir, "/data/merged_by_year/sib_exclusions/impatiens_sibexclusions_2022.txt")
)

rcolony::build.colony.superauto(wd=workingdir, 
                                name=paste0(workingdir, "/colony_assignments/Colony2_Linux/impatiens_2023.DAT"), 
                                datasetname = "impatiens2023",
                                delim=",",
                                sample_size = 2102,
                                num_loci = 12,
                                sibship_prior = 0,
                                female_monogamy = 1,
                                error_rates_path = paste0(workingdir, "/data/merged_by_year/error_rates/impatiens_error_rates.txt"),
                                genotypes_path = paste0(workingdir, "/data/merged_by_year/sib_scores_for_colony/impatiens2023_forcolony.txt"),
                                exclusion_path = paste0(workingdir, "/data/merged_by_year/sib_exclusions/impatiens_sibexclusions_2023.txt")
)

rcolony::build.colony.input(wd=workingdir, 
                                name=paste0(workingdir, "/colony_assignments/Colony2_Linux/impatiens_2023wq.DAT"),
                                delim=",")

###################################################
## Load in results from colony and get sibships
##################################################

# set probability threshold, colony output filenames, specimen datasets
prob_thresh = 0.995
filenameimp22 = "colony_assignments/Colony2_Linux/impatiens2022.FullSibDyad"
filenameimp23 = "colony_assignments/Colony2_Linux/impatiens2023.FullSibDyad"

# get genotyped specimen dataset
colsToKeep = c("site", "round", "sample_pt", "sample_id", "year", "barcode_id", "active_flower", "notes", "final_id", "pollen")
specsubset = rbind(specimenData2022[,colnames(specimenData2022) %in% colsToKeep],
                   specimenData2023[,colnames(specimenData2023) %in% colsToKeep])
genotyped2022 = specsubset[specsubset$barcode_id %in% impatiens2022_forcolony$barcode_id,]
genotyped2023 = specsubset[specsubset$barcode_id %in% impatiens2023_forcolony$barcode_id,]


# get sibships
impsibs2022 = filterSibships(filename = filenameimp22,
                             prob_thresh = prob_thresh,
                             specdata = genotyped2022)
impsibs2023 = filterSibships(filename = filenameimp23,
                             prob_thresh = prob_thresh,
                             specdata = genotyped2023)

write.csv(impsibs2022, "data/siblingships/impatiens_sibships_2022.csv")
write.csv(impsibs2023, "data/siblingships/impatiens_sibships_2023.csv")
