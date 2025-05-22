##############################################################
# Create phenology plots of queen observations to determine
# "breakpoint" between mother queens and daughter queens
# Started by J Melanson
# July 30, 2024
##############################################################

# set up workspace
rm(list=ls())
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")

#load packages
source('colony_assignments/src/init.R')
source('colony_assignments/src/joinFunctions.R')
load.packages()


## **********************************************************
## Load in necessary data
## **********************************************************
#load bombus survey data (2023)
specimenData2023 = filter(
  as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023specimendata.csv"), 
                sep = ",", header = T, fill = TRUE),
  !is.na(round))
samplePoints = as.data.frame(read.table("/Users/jenna1/Documents/UBC/bombus_project/raw_data/allsamplepoints.csv", sep = ","))
colnames(samplePoints) = c("sample_pt", "gps","landowner","subsite")
sampleEffort2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023sampledata.csv", sep = ","))
vegData2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023vegetationdata.csv", sep = ","))


## **********************************************************
## Join julian dates to specimen data
## **********************************************************
specimens = joinJulianDates(specimenData2023, sampleEffort2023)

## **********************************************************
## Plot queen phenology for mixtus queens
## **********************************************************
mixtusQueens = filter(specimens,
                      final_id == "B. mixtus",
                      grepl("queen", notes))

ggplot(mixtusQueens, aes(x = julian_date)) +
  geom_bar() +
  theme_classic()
# cutoff = 160 aka June 9? or round 19 (none captured the entire round)


## **********************************************************
## Plot queen phenology for impatiens
## **********************************************************
impatiensQueens = filter(specimens,
                  final_id == "B. impatiens",
                  grepl("queen", notes))

ggplot(impatiensQueens, aes(x = julian_date)) +
  geom_bar() +
  theme_classic()

# cutoff could also be 160 (June 9) o we could push it back a bit to 175 (June 24)
# this would only add 3 bees to the dataset, but two of those were nest searching
# so probably valid to include as mother queens

# i may have captured some impatiens queens as workers in the late summer,
# so will be good to check for those as we go through genotyping
# reclassifying them as queens could help build up that secondary peak

