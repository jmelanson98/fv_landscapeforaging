###### Get cumulative colony posteriors with neutral simulation
### Simulate 100 null datasets
### Get cumulative colony posteriors
### January 16, 2026
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
library(terra)
library(stringr)

##### Set Environment #####
setwd("/Users/jenna1/fv_landscapeforaging")
bombus_path = "/Users/jenna1/Documents/UBC/bombus_project"

setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "/home/melanson/projects/def-ckremen/melanson"
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

source("src/analysis_functions.R")
 
# Load in data / models
mix.fit = readRDS("analysis/kernel_locations/foraging_modelfits/simpleforaging_mixtus_steepgradient.rds")
#imp.fit = "analysis/kernel_locations/foraging_modelfits/simpleforaging_impatiens_steepgradient.rds"
mixtus_sibs2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mixtus_sibs2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
impatiens_sibs2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
impatiens_sibs2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
effort1 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
effort2 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))
samplepoints = read.csv(paste0(bombus_path, "/raw_data/allsamplepoints.csv"), header = FALSE)


##### Simulate datasets #####
# Make trap locations into correct format for calculating distance
colnames(samplepoints)[1:2] = c("sample_pt", "coord")
samplepoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplepoints$coord), split = " "), function(x) tail(x, 1))
samplepoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplepoints$coord), split = " "), function(x) head(x, 3)[2])
samplepoints$site <- stringr::str_extract(samplepoints$sample_pt, "^[A-Za-z]{1,2}")
samplepoints = samplepoints[,colnames(samplepoints) %in% c("site", "sample_pt", "lat", "long")]

# convert to meter based crs
traps_sf = st_as_sf(samplepoints, coords = c("long", "lat"), crs = 4326)
traps_sf_m = st_transform(traps_sf, 900913)
traps_m = as.data.frame(cbind(traps_sf_m$sample_pt, traps_sf_m$site, st_coordinates(traps_sf_m)))
colnames(traps_m) = c("sample_pt", "site", "trap_x", "trap_y")
traps_m$site = as.factor(traps_m$site)

# calculate sample effort per trap and year
effort_sum1 = effort1 %>%
  filter(round %in% sibships1$round | round %in% sibships2$round) %>%
  group_by(sample_point, year) %>%
  summarize(total_effort = 5*n())
effort_sum2 = effort2 %>%
  filter(round %in% sibships1$round | round %in% sibships2$round) %>%
  group_by(sample_point, year) %>%
  summarize(total_effort = 5*n())
effort_sum = rbind(effort_sum1, effort_sum2)
traps_m = traps_m %>%
  left_join(effort_sum, by = c("sample_pt" = "sample_point")) %>%
  left_join(site_keys, by = "site") %>%
  arrange(site_id)
traps_m = traps_m[!is.na(traps_m$total_effort),]
traps_m$trap_x = as.numeric(traps_m$trap_x)/1000
traps_m$trap_y = as.numeric(traps_m$trap_y)/1000
traps_m$trap_id = 1:nrow(traps_m)


# Get rho draws from fitted model
rhodraws_mix = rstan::extract(mix.fit, pars = "rho")$rho
#rhodraws_imp = rstan::extract(mix.fit, pars = "rho")$rho
mix_rho = rhodraws_mix[task_id*4000/100]
#imp_rho = rhodraws_imp[task_id*4000/100]

# Simulate one dataset
mixsim = draw_simple_true_sites(sample_size = 2000,
                          number_colonies = 11502,
                          colony_sizes = rep(20,11502),
                          trap_data = traps_m,
                          rho = mix_rho,
                          distance_decay = "exponentiated_quadratic")
# impsim = draw_simple_true_sites(sample_size = 4000,
#                                 number_colonies = 8004,
#                                 colony_sizes = rep(20,8004),
#                                 trap_data = traps_m,
#                                 rho = imp_rho,
#                                 distance_decay = "exponentiated_quadratic")

# Get colony counts (real data)
samplepoints = read.csv(paste0(bombus_path, "/raw_data/allsamplepoints.csv"), header = FALSE)
mixtus_data = prep_stan_simpleforaging_bothyears(mixtus_sibs2022,
                                                 mixtus_sibs2023,
                                                 effort1,
                                                 effort2,
                                                 samplepoints)
impatiens_data = prep_stan_simpleforaging_bothyears(impatiens_sibs2022,
                                                    impatiens_sibs2023,
                                                    effort1,
                                                    effort2,
                                                    samplepoints)
# Draw colonies from simulated dataset to match real data
mixcounts_summary = mixtus_data[[2]] %>% 
  group_by(stansibkey) %>%
  summarize(n = sum(count)) %>%
  group_by(n) %>%
  summarize(num_per_bin = n())
mixsimcounts = mixsim %>% 
  group_by(colonyid) %>%
  summarize(n = sum(counts)) %>%
  group_by(n) %>%
  summarize(numperbin = n())
  
  
  filter(n %in% mixcounts_summary$n) %>%
  group_by(n) %>%
  tidyr::nest() %>%
  left_join(mixcounts_summary, by = "n")
mixsimIDs_tokeep = map2(mixsimcounts$data, mixsimcounts$num_per_bin,
                          ~ slice_sample(.x, n = .y, replace = FALSE)) %>%
                      bind_rows()
mixsim_tokeep = mixsim[mixsim$colonyid %in% mixsimIDs_tokeep$colonyid,]


# impcounts_summary = imptus_data[[2]] %>% 
#   group_by(stansibkey) %>%
#   summarize(n = sum(count)) %>%
#   group_by(n) %>%
#   summarize(num_per_bin = n())
# impsimcounts = impsim %>% 
#   group_by(colonyid) %>%
#   summarize(n = sum(counts)) %>%
#   filter(n %in% impcounts_summary$n) %>%
#   group_by(n) %>%
#   tidyr::nest() %>%
#   left_join(impcounts_summary, by = "n")
# impsimIDs_tokeep = map2(impsimcounts$data, impsimcounts$num_per_bin,
#                         ~ slice_sample(.x, n = .y, replace = FALSE)) %>%
#   bind_rows()
# impsim_tokeep = impsim[impsim$colonyid %in% impsimIDs_tokeep$colonyid,]


# Write to file
mixdir = paste0("analysis/kernel_locations/neutral_data/mixdataset", task_id)
dir.create(mixdir, recursive = TRUE, showWarnings = FALSE)
write.csv(mixsim_tokeep, paste0(mixdir, "/simulation", task_id, ".csv"), row.names = FALSE)

# impdir = paste0("analysis/kernel_locations/neutral_data/impdataset", task_id)
# dir.create(impdir, recursive = TRUE, showWarnings = FALSE)
# write.csv(impsim_tokeep, paste0(impdir, "/simulation", task_id, ".csv"), row.names = FALSE)
