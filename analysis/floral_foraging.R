###### Fit foraging distance model to real data!
### With floral effects and ~multiple timepoints~
### December 5, 2025
### J. Melanson


##### Load packages #####
library(rstan)
library(matrixStats)
library(sp)
library(gstat)
library(ggplot2)
library(reshape2)
library(raster)
library(rasterVis)
library(parallel)
library(future)
library(furrr)
library(dplyr)
library(tidyr)
library(gridExtra)
library(tibble)
library(sf)
library(stringr)

# Prep workspace
# local
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/UBC/bombus_project"

# remote
setwd("~/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "~/projects/def-ckremen/melanson/"

# Load in data
mixtus_sibs2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mixtus_sibs2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
impatiens_sibs2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
impatiens_sibs2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
specs2022 = read.csv(paste0(bombus_path, "/raw_data/2022specimendata.csv"))
specs2023 = read.csv(paste0(bombus_path, "/raw_data/2023specimendata.csv"))
effort2022 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
effort2023 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))
veg2022 = read.csv(paste0(bombus_path, "/raw_data/2022vegetationdata.csv"))
veg2023 = read.csv(paste0(bombus_path, "/raw_data/2023vegetationdata.csv"))
samplepoints = read.csv(paste0(bombus_path, "/raw_data/allsamplepoints.csv"), header = FALSE)


# Get colonies x trap x timepoint matrix
cxl = impatiens_sibs2022 %>%
  distinct(site, sibshipID)
sereduced = effort2022[,colnames(effort2022) %in% c("site", "sample_point", "round", "sample_id")]
CKT = left_join(cxl, sereduced, by = "site", relationship = "many-to-many")
CKT$round = as.numeric(CKT$round)

sitekey = data.frame(site = c("ED", "HR", "NR", "PM", "SD", "W"),
                     siteid = 1:6)
trapkey = data.frame(sample_point = sort(unique(CKT$sample_point)),
                     trap_id = 1:length(unique(CKT$sample_point)))
  


# Get floral abundance estimates
beeflowers = unique(c(specs2022$active_flower, specs2023$active_flower))
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]

columnlist = c("sample_id", beeflowers)
floral_df_wide = veg2022[,colnames(veg2022) %in% columnlist]
floral_df_long = floral_df_wide %>%
  pivot_longer(!c(sample_id), names_to = "flower", values_to = "floral_abundance") %>%
  mutate(across(-c(sample_id, flower),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(-c(sample_id, flower),
                ~ 10^(.-1))) %>%
  group_by(sample_id) %>%
  summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE)) %>%
  mutate(across(-c(sample_id),
                ~ log10(.)))

# Get trap coordinates
colnames(samplepoints)[1:2] = c("sample_pt", "coord")
samplepoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplepoints$coord), split = " "), function(x) tail(x, 1))
samplepoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplepoints$coord), split = " "), function(x) head(x, 3)[2])
samplepoints$site <- stringr::str_extract(samplepoints$sample_pt, "^[A-Za-z]{1,2}")
samplepoints = samplepoints[,colnames(samplepoints) %in% c("site", "sample_pt", "lat", "long")]

# convert to meter based crs
traps_sf = st_as_sf(samplepoints, coords = c("long", "lat"), crs = 4326)
traps_sf_m = st_transform(traps_sf, 900913)
traps_m = as.data.frame(cbind(traps_sf_m$sample_pt, traps_sf_m$site, st_coordinates(traps_sf_m)))
colnames(traps_m) = c("sample_point", "site", "trap_x", "trap_y")
traps_m$site = as.factor(traps_m$site)


# Get Julian dates
#add julian date to sample effort data frame
effort2022$date = paste(effort2022$day, effort2022$month, effort2022$year)
effort2022$date = gsub(" ", "", effort2022$date, fixed = TRUE)
effort2022$date <- as.POSIXlt(effort2022$date, format = "%d%b%y")
effort2022$julian_date = effort2022$date$yday

# Get nonzero counts
counts = impatiens_sibs2022 %>%
  filter(notes != "male") %>%
  group_by(round, sample_pt, sibshipID) %>%
  summarize(count=n()) %>%
  ungroup()
colnames(counts) = c("round", "sample_point", "sibshipID", "counts")

# Join!
CKT = CKT %>%
  left_join(sitekey) %>%
  left_join(floral_df_long) %>%
  left_join(traps_m) %>%
  left_join(trapkey) %>%
  left_join(counts) %>%
  left_join(effort2022[,colnames(effort2022) %in% c("sample_id", "julian_date")]) %>%
  arrange(round, sibshipID, sample_point)
CKT$counts[is.na(CKT$counts)] = 0 # these are true zeroes
CKT$floral_abundance[CKT$floral_abundance == -Inf] = 0 #questionable -- not true zeroes
CKT$floral_abundance[is.na(CKT$floral_abundance)] = 0 #questionable -- not true zeroes
CKT$trap_x = as.numeric(CKT$trap_x)/1000
CKT$trap_y = as.numeric(CKT$trap_y)/1000


# Get CT length & starts
#CKT$sibround = paste("T", CKT$round, "C", CKT$sibshipID, sep = "-")
CTlengths = CKT %>%
  group_by(sibshipID, round) %>%
  summarize(length = n()) %>%
  arrange(round, sibshipID)
CTstarts = cumsum(c(1, CTlengths$length))[1:length(CTlengths$length)]


# Fit stan model
data = list(
  C = length(unique(CKT$sibshipID)),
  K = length(unique(CKT$sample_point)),
  O = nrow(CKT),
  #CT = length(unique(CKT$sibshipID)) * length(unique(CKT$round)),
  #CTstarts = CTstarts,
  #CTlengths = CTlengths$length,
  trap_pos = cbind(CKT$trap_x, CKT$trap_y),
  colony_id = CKT$sibshipID,
  trap_id = CKT$trap_id,
  landscape_id = CKT$siteid,
  fq = CKT$floral_abundance,
  doy = CKT$julian_date,
  yobs = CKT$counts,
  lower_x = (min(CKT$trap_x) - 5),
  upper_x = (max(CKT$trap_x) + 5),
  lower_y = (min(CKT$trap_y) - 5),
  upper_y = (max(CKT$trap_y) + 5)
)

stanfile = "models/floral_foraging_zinb.stan"

stanFit = stan(file = stanfile,
               data = data, seed = 5838299,
               chains = 4, cores = 4,
               iter = 2000,
               refresh = 1,
               verbose = TRUE)
saveRDS(stanFit, "analysis/tempFit.rds")
