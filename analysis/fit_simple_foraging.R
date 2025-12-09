###### Fit foraging distance model to real data!
### November 14, 2025
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

##### Set Environment #####
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/UBC/bombus_project"

setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "/home/melanson/projects/def-ckremen/melanson"
source("src/analysis_functions.R")


# Load in data
mixtus_sibs2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mixtus_sibs2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
impatiens_sibs2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
impatiens_sibs2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
specs2022 = read.csv(paste0(bombus_path, "/raw_data/2022specimendata.csv"))
specs2023 = read.csv(paste0(bombus_path, "/raw_data/2023specimendata.csv"))
effort2022 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
effort2023 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))
samplepoints = read.csv(paste0(bombus_path, "/raw_data/allsamplepoints.csv"), header = FALSE)



# Prepare data for Stan (with function)
mixtus_data = prep_stan_simpleforaging_bothyears(mixtus_sibs2022,
                                                 mixtus_sibs2023,
                                                 effort2022,
                                                 effort2023,
                                                 samplepoints)
impatiens_data = prep_stan_simpleforaging_bothyears(impatiens_sibs2022,
                                                 impatiens_sibs2023,
                                                 effort2022,
                                                 effort2023,
                                                 samplepoints)


# Get task ID from slurm manager
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
data_list = list(mixtus_data[[1]], impatiens_data[[1]], mixtus_data[[1]], impatiens_data[[1]])
data = data_list[[task_id]]

if (task_id %in% c(1,2)){
  #select stan model to fit
  stanfile = "models/simple_multinomial.stan"
  
  #fit and save model
  stanFit = stan(file = stanfile,
                 data = data, seed = 5838299,
                 chains = 4, cores = 4,
                 control = list(max_treedepth = 15),
                 iter = 4000,
                 verbose = TRUE)
  saveRDS(stanFit, paste0("analysis/foraging_modelfits/simpleforaging_multinomial", task_id, ".rds"))
} else{
  # remove unnecessary data
  data$starts = NULL
  data$lengths = NULL
  data$yn = NULL
  
  #select stan model to fit
  stanfile = "models/simple_ZINB.stan"
  
  #fit and save model
  stanFit = stan(file = stanfile,
                 data = data, seed = 5838299,
                 chains = 4, cores = 4,
                 control = list(max_treedepth = 15),
                 iter = 4000,
                 verbose = TRUE)
  saveRDS(stanFit, paste0("analysis/foraging_modelfits/simpleforaging_ZINB", task_id, ".rds"))
}


#############################################
# Extract values and makes some plots
#############################################

# mixtus 2022
# stanFitm22 = readRDS("analysis/foraging_modelfits/foragingmodel_1.rds")
# summarym22 = rstan::summary(stanFitm22)$summary
# 
# rhom22 = c(summarym22["rho", "2.5%"],
#            summarym22["rho", "mean"],
#            summarym22["rho", "97.5%"])
# rhom22 = rhom22*1000
# 
# # mixtus 2023
# stanFitm23 = readRDS("analysis/foraging_modelfits/foragingmodel_2.rds")
# summarym23 = rstan::summary(stanFitm23)$summary
# 
# rhom23 = c(summarym23["rho", "2.5%"],
#            summarym23["rho", "mean"],
#            summarym23["rho", "97.5%"])
# rhom23 = rhom23*1000
# 
# 
# # impatiens 2022
# stanFiti22 = readRDS("analysis/foraging_modelfits/foragingmodel_3.rds")
# summaryi22 = rstan::summary(stanFiti22)$summary
# 
# rhoi22 = c(summaryi22["rho", "2.5%"],
#            summaryi22["rho", "mean"],
#            summaryi22["rho", "97.5%"])
# rhoi22 = rhoi22*1000
# 
# # impatiens 2023
# stanFiti23 = readRDS("analysis/foraging_modelfits/foragingmodel_4.rds")
# summaryi23 = rstan::summary(stanFiti23)$summary
# 
# rhoi23 = c(summaryi23["rho", "2.5%"],
#            summaryi23["rho", "mean"],
#            summaryi23["rho", "97.5%"])
# rhoi23 = rhoi23*1000
# 
# 
# # Make plot of rhos for each species/year!
# postm22 = as.data.frame(stanFitm22)
# postm23 = as.data.frame(stanFitm23)
# posti22 = as.data.frame(stanFiti22)
# posti23 = as.data.frame(stanFiti23)
# 
# combinedpost =  bind_rows(
#   data.frame(model = "B. mixtus 2022", rho = postm22$rho),
#   data.frame(model = "B. mixtus 2023", rho = postm23$rho),
#   data.frame(model = "B. impatiens 2022", rho = posti22$rho),
#   data.frame(model = "B. impatiens 2023", rho = posti23$rho)
# )
# 
# ggplot(combinedpost, aes(x = rho, fill = model)) +
#   geom_histogram(bins = 30) +
#   facet_wrap(~ model) +
#   theme_minimal()
# 
# ggplot(combinedpost, aes(x = rho, fill = model)) +
#   scale_fill_manual(values = c(medium_gold, lm_gold, faded_strong, faded_green)) +
#   stat_halfeye(point_interval = median_qi, .width = c(0.95)) +
#   facet_wrap(~model, ncol = 1) +
#   theme_minimal()
# 
# year_difference = bind_rows(
#   data.frame(species = "mixtus", rho_diff = postm22$rho - postm23$rho),
#   data.frame(species = "impatiens", rho_diff = posti22$rho - posti23$rho)
# )
# 
# ggplot(year_difference, aes(x = rho_diff, fill = species)) +
#   geom_histogram(alpha = 0.4) +
#   theme_minimal()
# 
# mix_diff = year_difference[year_difference$species == "mixtus",]
# bayesp = sum(mix_diff$rho_diff > 0)/nrow(mix_diff)
# 
# species_difference = bind_rows(
#   data.frame(year = "2022", rho_diff = posti22$rho - postm22$rho),
#   data.frame(year = "2023", rho_diff = posti23$rho - postm23$rho)
# )
# 
# ggplot(species_difference, aes(x = rho_diff, fill = year)) +
#   geom_histogram(alpha = 0.4) +
#   theme_minimal()
