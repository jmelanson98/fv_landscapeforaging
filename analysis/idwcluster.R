### Calculate idw on cluster


# Load packages
library(terra)
library(sf)
library(raster)
library(data.table)
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

#############################################
# Get RHO values for each species
#############################################

# using multinomial models for now...

# mixtus
mixtusFit = readRDS("analysis/foraging_modelfits/simpleforaging_multinomial1.rds")
summarym = rstan::summary(mixtusFit)$summary

rhom = c(summarym["rho", "2.5%"],
         summarym["rho", "mean"],
         summarym["rho", "97.5%"])
rhom = rhom*1000


# impatiens
impatiensFit = readRDS("analysis/foraging_modelfits/simpleforaging_multinomial2.rds")
summaryi = rstan::summary(impatiensFit)$summary

rhoi = c(summaryi["rho", "2.5%"],
         summaryi["rho", "mean"],
         summaryi["rho", "97.5%"])
rhoi = rhoi*1000



#############################################
# Calculate IDW seminatural area per point
#############################################
buffer = 1500

# for mixtus
fv_points$idwSN_mix = NA
fv_points$idwHAY_mix = NA
fv_points$idwIND_mix = NA
fv_points$idwURB_mix = NA
fv_points$idwANN_mix = NA
fv_points$idwBLU_mix = NA
fv_points$idwPER_mix = NA
fv_points$idwFAL_mix = NA
fv_points$idwSN_imp = NA
fv_points$idwHAY_imp = NA
fv_points$idwIND_imp = NA
fv_points$idwURB_imp = NA
fv_points$idwANN_imp = NA
fv_points$idwBLU_imp = NA
fv_points$idwPER_imp = NA
fv_points$idwFAL_imp = NA


condition_vector = c("idwSN_mix", "idwHAY_mix", "idwIND_mix", "idwURB_mix", "idwANN_mix", "idwBLU_mix", "idwPER_mix", "idwFAL_mix",
                     "idwSN_imp", "idwHAY_imp", "idwIND_imp", "idwURB_imp", "idwANN_imp", "idwBLU_imp", "idwPER_imp", "idwFAL_imp")
filepaths = c("/landscape/rasters/seminat.tif",
              "/landscape/rasters/hay.tif",
              "/landscape/rasters/industrial.tif",
              "/landscape/rasters/urban.tif",
              "/landscape/rasters/annual.tif",
              "/landscape/rasters/blueberry.tif",
              "/landscape/rasters/perennial.tif",
              "/landscape/rasters/fallow.tif",
              "/landscape/rasters/seminat.tif",
              "/landscape/rasters/hay.tif",
              "/landscape/rasters/industrial.tif",
              "/landscape/rasters/urban.tif",
              "/landscape/rasters/annual.tif",
              "/landscape/rasters/blueberry.tif",
              "/landscape/rasters/perennial.tif",
              "/landscape/rasters/fallow.tif"
              )


# get task id
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

landscape = rast(paste0(bombus_path, filepaths[[task_id]]))
crs(landscape) = "epsg:900913"

for (i in 1:nrow(fv_points)) {
  if (task_id %in% 1:8){
    rho = rhom[2]
  } else {
    rho = rhoi[2]
  }
  print(rho)
  
  fv_points[[condition_vector[task_id]]][i] = compute_idw_area(fv_points[i,], landscape, buffer, rho)
  if (i %% 10 == 0) cat("Processed point", i, "of", nrow(fv_points), "\n")
}

saved = read.csv("analysis/calculated_metrics.csv")
updated = left_join(saved, as.data.frame(fv_points)[ c("site_id", condition_vector[task_id])], by = c("sample_pt" = "site_id"))
write.csv(updated, "analysis/calculated_metrics.csv", row.names = FALSE)
