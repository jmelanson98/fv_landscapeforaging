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


#########################
# For Mixtus
##########################
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



#######################
# For Impatiens
######################
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
