########
## Get colony IDs from microsatellite genotyping results
########
rm(list = ls())

# first, load in packages
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
specimenData2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022specimendata.csv", sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2023specimendata.csv", sep = ",", header = T, fill = TRUE))
samplePoints = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/allsamplepoints.csv", sep = ","))
colnames(samplePoints) = c("sample_pt", "gps","landowner","subsite", "twosub")

#load error tables
plex1_init = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/checkerrors/36plex1.csv", sep = ",", header = T))
plex2_init = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/checkerrors/36plex2.csv", sep = ",", header = T))
plex1_rerun = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/checkerrors/36rerunplex1.csv", sep = ",", header = T))
plex2_rerun = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/checkerrors/36rerunplex2.csv", sep = ",", header = T))

#check error rates

plex1merged = inner_join(plex1_init, plex1_rerun, by = "Name", suffix = c("_init", "_rerun"))
plex2merged = inner_join(plex2_init, plex2_rerun, by = "Name", suffix = c("_init", "_rerun"))

# Function to calculate match percentage per column
compare_metrics <- function(df, metrics) {
  results <- metrics %>%
    map_df(~ {
      col1 <- sym(paste0(.x, "_init"))
      col2 <- sym(paste0(.x, "_rerun"))
      
      valid_cases <- df %>% filter(!is.na(!!col1) & !is.na(!!col2))
      mismatches <- valid_cases %>% filter(!!col1 != !!col2)
      
      matches <- nrow(valid_cases) - nrow(mismatches)
      total <- nrow(valid_cases)
      mismatch_names <- mismatches$Name
      
      tibble(
        Metric = .x, 
        Matches = matches, 
        Total_Comparisons = total, 
        Match_Percentage = ifelse(total > 0, (matches / total) * 100, NA),
        Mismatch_Names = ifelse(length(mismatch_names) > 0, paste(mismatch_names, collapse = ", "), "None")
      )
    })
  return(results)
}

# Identify metric columns (excluding "Name")
metric_columns1 <- setdiff(names(plex1_init), "Name")
metric_columns2 <- setdiff(names(plex2_init), "Name")

# Compute match summary
summary_table1 <- compare_metrics(plex1merged, metric_columns1)
summary_table2 <- compare_metrics(plex2merged, metric_columns2)

print(summary_table2)





#load full allele tables
alleles_plex1_all = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/mixtus_allelestable_plex1.csv", sep = ",", header = T))
alleles_plex2_all = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/mixtus_allelestable_plex2.csv", sep = ",", header = T))



#modify allele table "name" value to match specimen data frame "plate" and "well" (for originals)
alleles_plex1_original %>% separate(Name, sep = "-", c(NA, "plate", "well"), extra = "merge", fill = "left") -> alleles_plex1_original                                         
alleles_plex1_original$location = str_remove_all(paste(alleles_plex1_original$plate, alleles_plex1_original$well), "-| ")
alleles_plex1_original = alleles_plex1_original[,!colnames(alleles_plex1_original) %in% c("plate", "well")]

alleles_plex2_original %>% separate(Name, sep = "-", c(NA, "plate", "well"), extra = "merge", fill = "left") -> alleles_plex2_original                                         
alleles_plex2_original$location = str_remove_all(paste(alleles_plex2_original$plate, alleles_plex2_original$well), "-| ")
alleles_plex2_original = alleles_plex2_original[,!colnames(alleles_plex2_original) %in% c("plate", "well")]

#merge original dataframes to barcodesubset so that all four contain barcode_id (remove location)
alleles_plex1_wbarcode = left_join(alleles_plex1_original, barcodesubset, by = "location")
alleles_plex1_wbarcode = alleles_plex1_wbarcode[,!colnames(alleles_plex1_wbarcode) %in% c("location")]
alleles_plex1_wbarcode %>% filter(!is.na(barcode_id)) -> alleles_plex1_wbarcode

alleles_plex2_wbarcode = left_join(alleles_plex2_original, barcodesubset, by = "location")
alleles_plex2_wbarcode = alleles_plex2_wbarcode[,!colnames(alleles_plex2_wbarcode) %in% c("location")]
alleles_plex2_wbarcode %>% filter(!is.na(barcode_id)) -> alleles_plex2_wbarcode

#modify allele table name value of re-runs to match barcode_id
alleles_plex1_rerun$barcode_id = gsub("^.{0,6}", "", alleles_plex1_rerun$Name)
alleles_plex1_rerun$barcode_id = gsub("-","_",alleles_plex1_rerun$barcode_id)
alleles_plex1_rerun = alleles_plex1_rerun[,!colnames(alleles_plex1_rerun) %in% c("Name")]

alleles_plex2_rerun$barcode_id = gsub("^.{0,6}", "", alleles_plex2_rerun$Name)
alleles_plex2_rerun$barcode_id = gsub("-","_",alleles_plex2_rerun$barcode_id)
alleles_plex2_rerun = alleles_plex2_rerun[,!colnames(alleles_plex2_rerun) %in% c("Name")]

#reorder columns so that we can rbind originals to reruns
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

#make a specimen subset that include year and notes (to filter 2022 vs 2023, plus queens)
rowsToKeep = c("barcode_id", "year", "notes", "final_id")
specsubset = rbind(specimenData2022[,colnames(specimenData2022) %in% rowsToKeep],
                   specimenData2023[,colnames(specimenData2023) %in% rowsToKeep])
specimens_withscores = left_join(microsat_scores, specsubset, by = "barcode_id")

mixtus2022 = filter(specimens_withscores, year == "2022" | str_detect(notes, "queen"))
mixtus2023 = filter(specimens_withscores, year == "2023" & !str_detect(notes, "queen"))


#check which bees are unscored
specimens = full_join(microsat_scores, specsubset, by = "barcode_id")
mixtus = filter(specimens, final_id == "B. mixtus" & barcode_id != "N/A")
unscoredmixtus = mixtus[rowSums(is.na(mixtus)) >10,]


##### 
# Prep genotype data for COLONY
#####

#remove all columns except barcode and scores
colsToRemove = c("year", "final_id", "notes")
mixtus2022 = mixtus2022[, !colnames(mixtus2022) %in% colsToRemove]
mixtus2023 = mixtus2023[, !colnames(mixtus2023) %in% colsToRemove]

#relocate barcode and remove rows with more than 10 NAs, replace NAs with 0's
mixtus2022 %>% relocate(barcode_id) -> mixtus2022
mixtus2022_forcolony = mixtus2022[rowSums(is.na(mixtus2022)) < 10,]
mixtus2022_forcolony[is.na(mixtus2022_forcolony)] <- 0

mixtus2023 %>% relocate(barcode_id) -> mixtus2023
mixtus2023_forcolony = mixtus2023[rowSums(is.na(mixtus2023)) < 10,]
mixtus2023_forcolony[is.na(mixtus2023_forcolony)] <- 0

#write csvs for upload to colony
write.csv(mixtus2022_forcolony, "mixtus_2022_scores.csv")
write.csv(mixtus2023_forcolony, "mixtus_2023_scores.csv")

##### 
# Prep sibship exclusion data for COLONY
#####
#basically -- we want to exclude sibships from different regions, as they are 
#quite far apart and therefore very unlikely to be genuine




####
# run colony
####
#this code is not working; use windows GUI on lab computer
build.colony.input(wd=getwd(), name="Colony2.DAT", delim=",")



####
## load in results from colony
####
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


####
# remove siblings 
####
mixtus_wsibs %>%
  filter(!is.na(fullsibshipindex)) %>%
  group_by(fullsibshipindex) %>%
  sample_n(1) -> nonsibs

testalleles = mixtus2023_forcolony[mixtus2023_forcolony$barcode_id %in% nonsibs$barcode_id,]

####
# prep alleles table for genepop
####

#combine columns for each locus
ms_new = data.table(testalleles)
ms_new[ms_new ==0] <- "000"


i1 <- seq(1, length(ms_new)-1, 2)
i2 <- seq(2, length(ms_new)-1, 2)
msgenepop = ms_new[, Map(paste,
                         .SD[, i1, with = FALSE], .SD[, i2, with = FALSE], 
                         MoreArgs = list(sep="")), 
                   by = "barcode_id"]
colnames(msgenepop) = c("barcode_id", "BT10", "BTMS0104", "BTMS0057", 
                        "BTMS0086", "BTMS0066", "BTMS0062", "BTMS0136",
                        "BTERN01", "BTMS0126", "BTMS0072", "BTMS0059", "BL15",
                        "BL13", "BTMS0083", "B126")


#select and join alleles

#sort allele table by barcode ID
msgenepop <- msgenepop[order(msgenepop$barcode_id),]
msgenepop$site = substring(msgenepop$barcode_id, 1, 1)

haplotypes <- as.data.frame(paste(msgenepop$site, "_", msgenepop$barcode_id, ","," ", 
                                  msgenepop$BT10," ", msgenepop$BTMS0104, " ",
                                  msgenepop$BTMS0057," ", msgenepop$BTMS0086, " ",
                                  msgenepop$BTMS0066," ", msgenepop$BTMS0062, " ",
                                  msgenepop$BTMS0136," ", msgenepop$BTERN01, " ",
                                  msgenepop$BTMS0126," ", msgenepop$BTMS0072, " ",
                                  msgenepop$BTMS0059," ", msgenepop$BL15, " ",
                                  msgenepop$BL13," ", msgenepop$BTMS0083, " ",
                                  msgenepop$B126,
                                  sep = ""))

#Build the Genepop format
sink("ms_genepop_format.txt")
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
cat("BTMS0072 \n")
cat("BTMS0059 \n")
cat("BL15 \n")
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



###
#run genepop on samples
####

#test for HWE
test_HW("ms_genepop_format.txt", which='Proba', 'outputHWE.txt.D')

#test for LD
test_LD("ms_genepop_format.txt","outputLD.txt.DIS")



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
