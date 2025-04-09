########
## Get error rates using re-run plates
########


# first, load in packages
library(dplyr)
library(tidyr)
library(stringr)
library(rcolony)
library(genepop)
library(data.table)

# next; load in allele tables
alleles_plex1 = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/alleles_plex1.csv", sep = ",", header = T))
alleles_plex1_rerun = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/alleles_plex1_rerun.csv", sep = ",", header = T))

#fix name column for consistency
alleles_plex1$Name = gsub("-","",alleles_plex1$Name)
alleles_plex1_rerun$Name = gsub("-","",alleles_plex1_rerun$Name)
alleles_plex1_rerun$Name = gsub("^.{0,5}", "", alleles_plex1_rerun$Name)
rownames(alleles_plex1) = alleles_plex1$Name
rownames(alleles_plex1_rerun) = alleles_plex1_rerun$Name

#compare to see which cells are the same
matching_alleles = alleles_plex1 == alleles_plex1_rerun



#### repeat for plex 2!

#load in allele tables
alleles_plex2 = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/alleles_plex2.csv", sep = ",", header = T))
alleles_plex2_rerun = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/alleles_plex2_rerun.csv", sep = ",", header = T))

#fix name column for consistency
alleles_plex2$Name = gsub("-","",alleles_plex2$Name)
alleles_plex2_rerun$Name = gsub("-","",alleles_plex2_rerun$Name)
alleles_plex2_rerun$Name = gsub("^.{0,5}", "", alleles_plex2_rerun$Name)
rownames(alleles_plex2) = alleles_plex2$Name
rownames(alleles_plex2_rerun) = alleles_plex2_rerun$Name

#compare to see which cells are the same
matching_alleles2 = alleles_plex2 == alleles_plex2_rerun
