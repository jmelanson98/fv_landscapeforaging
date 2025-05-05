########
## Map sibling locations and make rough foraging distance estimates (regression approach)
## Started by J Melanson
## May 5, 2025
########

###########################
### Prepare environment
###########################
rm(list = ls())
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")

# first, load in packages
source('colony_assignments/src/init.R')
source('colony_assignments/src/joinFunctions.R')
source('colony_assignments/src/makeSibMaps.R')
library(dplyr)
library(tidyr)
library(stringr)
library(rcolony)
library(genepop)
library(data.table)
library(cowplot)
library(raster)


###########################
### Load in data
###########################
# raw data
specimenData2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2022specimendata.csv", sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023specimendata.csv", sep = ",", header = T, fill = TRUE))
sampleEffort2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2022sampledata.csv", sep = ","))
sampleEffort2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023sampledata.csv", sep = ","))
samplePoints = as.data.frame(read.table("/Users/jenna1/Documents/UBC/bombus_project/raw_data/allsamplepoints.csv", sep = ","))
colnames(samplePoints) = c("sample_pt", "gps","landowner","subsite")

# sibships
mix2022 = read.csv("data/siblingships/mix_sibships_preliminary_2022.csv")
mix2023 = read.csv("data/siblingships/mix_sibships_preliminary_2023.csv")
imp2022 = read.csv("data/siblingships/imp_sibships_preliminary_2022.csv")
imp2023 = read.csv("data/siblingships/imp_sibships_preliminary_2023.csv")


######################
### Make sib maps
#####################

# merge 2022 + 2023 specimens
allspecs = rbind(specimenData2022[, colnames(specimenData2022) %in% c("site", "round", "sample_pt", "sample_id", "year", "barcode_id", "active_flower", "final_id", "notes", "pollen")],
                 specimenData2023[, colnames(specimenData2023) %in% c("site", "round", "sample_pt", "sample_id", "year", "barcode_id", "active_flower", "final_id", "notes", "pollen")]) %>%
  filter(final_id == "B. mixtus" | final_id == "B. impatiens")

#add lat and long to samplePoints
#first, wrangle gps coordinates into shape
samplePoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) tail(x, 1))
samplePoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) head(x, 3)[2])
samplePoints = samplePoints[,c("sample_pt", "lat", "long")]

#join coordinates to specimen data
allspecs = left_join(allspecs, samplePoints, by = "sample_pt")

# make unique cluster IDs
mix2023$ClusterIndex = mix2023$ClusterIndex + max(mix2022$ClusterIndex)
imp2023$ClusterIndex = imp2023$ClusterIndex + max(mix2022$ClusterIndex)

# join cluster IDs to specimen dataframe
tokeep = c("ClusterIndex", "OffspringID")
allsibs = rbind(mix2022[, colnames(mix2022) %in% tokeep],
                mix2023[, colnames(mix2023) %in% tokeep],
                imp2022[, colnames(imp2022) %in% tokeep],
                imp2023[, colnames(imp2023) %in% tokeep])
allspecs = left_join(allspecs, allsibs, by = c("barcode_id" = "OffspringID"))
genotypedspecs = filter(allspecs, !is.na(ClusterIndex))

#split 2022 and 2023
mix = genotypedspecs %>% filter(final_id == "B. mixtus")
imp = genotypedspecs %>% filter(final_id == "B. impatiens")

#make plot grids for both
mixgrids = makeSibMaps(mix)
length(unique(sibs22$fullsib_index))

grids2023 = makeSibMaps(sibs23)
length(unique(sibs23$fullsib_index))
