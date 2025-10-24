## Make dataframe of impatiens presence absence for Jens Ulrich paper
## October 23, 2025


# prep workspace
rm(list = ls())
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")
workingdir = getwd()

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


# next; load in specimen data and sample data
specimenData2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2022specimendata.csv", sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023specimendata.csv", sep = ",", header = T, fill = TRUE))
sampleEffort2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2022sampledata.csv", sep = ","))
sampleEffort2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023sampledata.csv", sep = ","))
samplePoints = as.data.frame(read.table("/Users/jenna1/Documents/UBC/bombus_project/raw_data/allsamplepoints.csv", sep = ","))
colnames(samplePoints) = c("sample_point", "gps","landowner","subsite")

# filter and join specimen tables
colstokeep = c("sample_id", "barcode_id", "final_id", "notes")
specimens = rbind(
  specimenData2022[,colnames(specimenData2022) %in% colstokeep],
  specimenData2023[,colnames(specimenData2023) %in% colstokeep])

impatiens_worker_summary = specimens %>% 
  filter(final_id == "B. impatiens") %>%
  filter(!str_detect(notes, "male")) %>%
  filter(!str_detect(notes, "queen")) %>%
  group_by(sample_id) %>%
  summarize(impatiens_worker_abundance = n())

impatiens_male_summary = specimens %>% 
  filter(final_id == "B. impatiens") %>%
  filter(str_detect(notes, "male")) %>%
  group_by(sample_id) %>%
  summarize(impatiens_male_abundance = n())

impatiens_queen_summary = specimens %>% 
  filter(final_id == "B. impatiens") %>%
  filter(str_detect(notes, "queen")) %>%
  group_by(sample_id) %>%
  summarize(impatiens_queen_abundance = n())

# join sample effort tables
sampledImp2022 = sampleEffort2022 %>%
  filter(sampledImp == TRUE)
sampledImp2022 = sampledImp2022[,!colnames(sampledImp2022) %in% c("sampledImp")]

sampleEffort = rbind(sampledImp2022, sampleEffort2023)

# wrangle gps coordinates into shape
samplePoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) tail(x, 1))
samplePoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) head(x, 3)[2])
samplePoints = samplePoints[,c("sample_point", "subsite", "lat", "long")]
samplePoints$lat = as.numeric(samplePoints$lat)
samplePoints$long = as.numeric(samplePoints$long)
samplePoints = samplePoints[,!colnames(samplePoints) %in% c("subsite")]

# join sample effort with coordinates
sampleEffort = full_join(sampleEffort, samplePoints, by = "sample_point")


# join impatiens abundance
fullDF = full_join(sampleEffort, impatiens_worker_summary, by = "sample_id")
fullDF = full_join(fullDF, impatiens_male_summary, by = "sample_id")
fullDF = full_join(fullDF, impatiens_queen_summary, by = "sample_id")

fullDF$impatiens_worker_abundance[is.na(fullDF$impatiens_worker_abundance)] = 0
fullDF$impatiens_queen_abundance[is.na(fullDF$impatiens_queen_abundance)] = 0
fullDF$impatiens_male_abundance[is.na(fullDF$impatiens_male_abundance)] = 0

# save output
write.csv(fullDF, "data/impatiens_presence_absence.csv")
