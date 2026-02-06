###### Test impact of floral quality inclusion on model performance
# With simulated datasets
### January 24, 2026
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
library(purrr)


##### Set Environment #####
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/UBC/bombus_project"

setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "/home/melanson/projects/def-ckremen/melanson"
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
source("src/analysis_functions.R")


################################################################################
###############               SIMULATE DATA                 ###################
################################################################################

# # # Load in data
# sibs1 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
# sibs2 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
# effort1 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
# effort2 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))
# samplepoints = read.csv(paste0(bombus_path, "/raw_data/allsamplepoints.csv"), header = FALSE)
# 
# # Get colony counts (real data)
# data = prep_stan_simpleforaging_bothyears(sibs1,
#                                           sibs2,
#                                           effort1,
#                                           effort2,
#                                           samplepoints)
# 
# # Get trap data
# traps_m = data[[3]]
# 
# ##### Simulate datasets #####
# 
# for (i in 1:3){
#   for (rho in c(0.25,0.375, 0.5, 0.625, 0.75)){
#     sim = draw_simple_true_sites(sample_size = 2000,
#                                  number_colonies = 11700,
#                                  colony_sizes = rep(20,11700),
#                                  trap_data = traps_m,
#                                  rho = rho,
#                                  theta = 0.5,
#                                  distance_decay = "exponentiated_quadratic")
#     sim = sim %>% group_by(colonyid) %>% filter(sum(counts)>0)
#     write.csv(sim, paste0("analysis/kernel_locations/modeltest/simulation", i, "_rho", rho, ".csv"))
# 
# 
#   }
# }

################################################################################
###############                  FIT MODELS                  ###################
################################################################################
# make params table
rho = c(0.25,0.375, 0.5, 0.625, 0.75)
sim = c(1,2,3)
Rmax = c(1.23, 2.68, 4.14)
stanfile = c("models/simple_multinomial_uniformdisprior.stan",
          "models/simple_multinomial_uniformdisprior_floral.stan",
          "models/poisson_uniformdisprior_floral.stan")
inclusion = c("all", "doubletons")
params = expand.grid(rho = rho, sim = sim, Rmax = Rmax, stanfile = stanfile, inclusion = inclusion)

stanfilekey = data.frame(stanfile = c("models/simple_multinomial_uniformdisprior.stan",
                           "models/simple_multinomial_uniformdisprior_floral.stan",
                           "models/poisson_uniformdisprior_floral.stan"),
                         key = c("multinomial_notheta", "multinomial_theta", "poisson_theta"))
params = left_join(params, stanfilekey)

# load in dataset
dir = paste0("analysis/kernel_locations/modeltest/simulation",params$sim[task_id], "_rho", params$rho[task_id], ".csv")
dataset = read.csv(dir)

# filter out singletons for one condition
if (params$inclusion[task_id] == "doubletons"){
  dataset = dataset %>% group_by(colonyid) %>% filter(sum(counts)>1)
}

stan_data = prep_stan_simpleforaging_neutraldata(dataset)
stan_data$Rmax = params$Rmax[task_id]
stan_data$steepness = 10

if (params$key[task_id] != "multinomial_notheta"){
  stan_data$fq = dataset$fq}

if(params$key[task_id] == "poisson_theta"){
  stan_data$starts = NULL
  stan_data$lengths = NULL
}

#fit and save model
fit = stan(file = as.character(params$stanfile[task_id]),
                  data = stan_data, seed = 5838299,
                  chains = 4, cores = 4,
                  iter = 2000,
                  verbose = TRUE)
saveRDS(fit, 
        paste0("analysis/kernel_locations/modeltest/", 
               params$key[task_id], 
               "_rho", params$rho[task_id], 
               "_Rmax", params$Rmax[task_id], 
               "_sim", params$sim[task_id],
               "_", params$inclusion[task_id], ".rds"))


# get model summary statistics
params$rho2.5 = NA
params$rho50 = NA
params$rho97.5 = NA
params$theta2.5 = NA
params$theta50 = NA
params$theta97.5 = NA

for (i in 1:nrow(params)){
  fit = readRDS(paste0("analysis/kernel_locations/modeltest/", 
                       params$key[i], 
                       "_rho", params$rho[i], 
                       "_Rmax", params$Rmax[i], 
                       "_sim", params$sim[i],
                       "_", params$inclusion[i], ".rds"))
  summar = summary(fit)$summary
  params$rho2.5[i] = summar["rho", "2.5%"]
  params$rho50[i] = summar["rho", "50%"]
  params$rho97.5[i] = summar["rho", "97.5%"]
  if(params$key[i] != "multinomial_notheta"){
    params$theta2.5[i] = summar["theta", "2.5%"]
    params$theta50[i] = summar["theta", "50%"]
    params$theta97.5[i] = summar["theta", "97.5%"]
  }
  print(paste0("Done with ", i))
}
write.csv(params, "analysis/kernel_locations/modeltest/params.csv")


# plot comparisons
result = read.csv("analysis/kernel_locations/modeltest/params.csv")
result$Rmax = as.factor(result$Rmax)


pd = position_dodge(width = 0.03)

ggplot(result[result$key != "poisson_theta",], aes(x = rho, y = rho50, color = Rmax, group = sim)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(position = pd, size = 1) +
  geom_errorbar(aes(ymin = rho2.5, ymax = rho97.5),
                position = pd,
                width = 0) +
  facet_grid(key~inclusion) +
  theme_bw()
