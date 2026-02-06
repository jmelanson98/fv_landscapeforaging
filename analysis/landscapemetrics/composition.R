##### Calculate landscape composition at 1500m
### Feb 5, 2026
### J. Melanson


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
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"
setwd("~/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "~/projects/def-ckremen/melanson/"
source("src/analysis_functions.R")

#####################################################
# Load and prep landcover data
#####################################################
landscape_raster = raster(paste0(bombus_path, "landscape/rasters/FValley_lc_1res.tif"))
fv_points = st_read(paste0(bombus_path, "landscape/fvbombus/fvbombus_points.shp"))

#Change CRS to meters
fv_points = st_transform(fv_points, 32610)
landscape_raster = rast(landscape_raster)
crs(landscape_raster) = "EPSG:3857"
landscape_raster = project(landscape_raster, terra::crs(fv_points))


#####################################################
# Check composition of all landcover classes
#####################################################
composition_df = sample_lsm(landscape_raster,
                            y = fv_points,
                            plot_id = fv_points$site_id,
                            shape = "circle",
                            size = 1500,
                            return_raster = TRUE,
                            what = "lsm_c_ca")
write.csv(composition_df, "analysis/landscapemetrics/composition1500.csv")





