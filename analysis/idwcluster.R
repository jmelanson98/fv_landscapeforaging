### Calculate idw on cluster


# Load packages
library(terra)
library(sf)
library(raster)
library(data.table)
library(dbscan)
library(tidyr)
library(dplyr)
library(bayesplot)
library(ggplot2)

# Set up workspace
setwd("~/projects/def-ckremen/melanson/fv_landscapeforaging/")
bombus_path = "~/projects/def-ckremen/melanson/"
source("src/analysis_functions.R")


# Load data
fv_points = st_read(paste0(bombus_path, "landscape/fvbombus/fvbombus_points.shp"))
st_crs(fv_points) = 900913
semi_file = paste0(bombus_path, "/landscape/rasters/seminat.tif")
semi = rast(semi_file)

#############################################
# Get RHO values for each species x year
#############################################
# mixtus 2022
stanFitm22 = readRDS("analysis/foraging_modelfits/foragingmodel_1.rds")
summarym22 = rstan::summary(stanFitm22)$summary

rhom22 = c(summarym22["rho", "2.5%"],
           summarym22["rho", "mean"],
           summarym22["rho", "97.5%"])
rhom22 = rhom22*1000

# mixtus 2023
stanFitm23 = readRDS("analysis/foraging_modelfits/foragingmodel_2.rds")
summarym23 = rstan::summary(stanFitm23)$summary

rhom23 = c(summarym23["rho", "2.5%"],
           summarym23["rho", "mean"],
           summarym23["rho", "97.5%"])
rhom23 = rhom23*1000


# impatiens 2022
stanFiti22 = readRDS("analysis/foraging_modelfits/foragingmodel_3.rds")
summaryi22 = rstan::summary(stanFiti22)$summary

rhoi22 = c(summaryi22["rho", "2.5%"],
           summaryi22["rho", "mean"],
           summaryi22["rho", "97.5%"])
rhoi22 = rhoi22*1000

# impatiens 2023
stanFiti23 = readRDS("analysis/foraging_modelfits/foragingmodel_4.rds")
summaryi23 = rstan::summary(stanFiti23)$summary

rhoi23 = c(summaryi23["rho", "2.5%"],
           summaryi23["rho", "mean"],
           summaryi23["rho", "97.5%"])
rhoi23 = rhoi23*1000


#############################################
# Calculate IDW seminatural area per point
#############################################
buffer = 1500

# make columns
fv_points$idwSN_mix22_low = NA
fv_points$idwSN_mix22_mean = NA
fv_points$idwSN_mix22_high = NA
fv_points$idwSN_mix23_low = NA
fv_points$idwSN_mix23_mean = NA
fv_points$idwSN_mix23_high = NA
fv_points$idwSN_imp22_low = NA
fv_points$idwSN_imp22_mean = NA
fv_points$idwSN_imp22_high = NA
fv_points$idwSN_imp23_low = NA
fv_points$idwSN_imp23_mean = NA
fv_points$idwSN_imp23_high = NA

condition_vector = c("idwSN_mix22_low", "idwSN_mix22_mean", "idwSN_mix22_max",
                   "idwSN_mix23_low", "idwSN_mix23_mean", "idwSN_mix22_max",
                   "idwSN_imp22_low", "idwSN_imp22_mean", "idwSN_imp22_max",
                   "idwSN_imp23_low", "idwSN_imp23_mean", "idwSN_imp23_max")
allrho = c(rhom22, rhom23, rhoi22, rhoi23)

# get task id
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

for (i in 1:nrow(fv_points)) {
  fv_points[[condition_vector[task_id]]][i] = compute_idw_area(fv_points[i,], semi, buffer, allrho[task_id]
  if (i %% 10 == 0) cat("Processed point", i, "of", nrow(fv_points), "\n")
}

saved = read.csv("analysis/calculated_metrics.csv")
updated = left_join(saved, as.data.frame(fv_points)[ c("site_id", condition_vector[task_id])], by = c("sample_pt" = "site_id"))
write.csv(updated, "analysis/calculated_metrics.csv")