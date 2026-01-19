##### Calculate distance weighted landscape metrics
### Nov 18, 2025
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
setwd("/Users/jenna1/fv_landscapeforaging")
bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"
#bombus_path = "~/projects/def-ckremen/melanson/"
source("src/analysis_functions.R")

#####################################################
# Load and prep landcover data
#####################################################
landscape_raster = raster(paste0(bombus_path, "landscape/rasters/FValley_lc_1res.tif"))
fv_points = st_read(paste0(bombus_path, "landscape/fvbombus/fvbombus_points.shp"))
landcover = read.csv(paste0(bombus_path, "landscape/landcover.csv"))

#Change CRS to meters
crs(landscape_raster) = 900913
st_crs(fv_points) = 900913

# Make sure it's a spatrast object
landscape_raster = rast(landscape_raster)


#####################################################
# Get seminatural raster
#####################################################

# Change land cover raster to 1/0 (seminatural/other)
# Combine hedgerows (blackberry + trees), grassy margins, forest, wetland
vals = c(2, 6, 8, 10, 16)
semi = ifel(landscape_raster %in% vals, 1L, 0L,
             filename = paste0(bombus_path, "landscape/rasters/seminat.tif"),
             overwrite = TRUE,
             wopt = list(datatype="INT1U"))
semi_file = paste0(bombus_path, "/landscape/rasters/seminat.tif")
semi = rast(semi_file)

#####################################################
# Get annual crops raster
#####################################################

# Change land cover raster to 1/0 (seminatural/other)
# combine annual crops + polyculture
vals = c(1, 13)
annual = ifel(landscape_raster %in% vals, 1L, 0L,
            filename = paste0(bombus_path, "landscape/rasters/annual.tif"),
            overwrite = TRUE,
            wopt = list(datatype="INT1U"))
annual_file = paste0(bombus_path, "/landscape/rasters/annual.tif")
annual = rast(annual_file)


#####################################################
# Get hay and pasture raster
#####################################################

# Change land cover raster to 1/0 (seminatural/other)
vals = c(7, 11)
hay = ifel(landscape_raster %in% vals, 1L, 0L,
            filename = paste0(bombus_path, "landscape/rasters/hay.tif"),
            overwrite = TRUE,
            wopt = list(datatype="INT1U"))
hay_file = paste0(bombus_path, "/landscape/rasters/hay.tif")
hay = rast(hay_file)


#####################################################
# Get urban area raster
#####################################################

# Change land cover raster to 1/0 (seminatural/other)
vals = c(14)
urban = ifel(landscape_raster %in% vals, 1L, 0L,
            filename = paste0(bombus_path, "landscape/rasters/urban.tif"),
            overwrite = TRUE,
            wopt = list(datatype="INT1U"))
urban_file = paste0(bombus_path, "/landscape/rasters/urban.tif")
urban = rast(urban_file)

#####################################################
# Get industrial raster
#####################################################

# Change land cover raster to 1/0 (seminatural/other)
vals = c(9)
industrial = ifel(landscape_raster %in% vals, 1L, 0L,
            filename = paste0(bombus_path, "landscape/rasters/industrial.tif"),
            overwrite = TRUE,
            wopt = list(datatype="INT1U"))
industrial_file = paste0(bombus_path, "/landscape/rasters/industrial.tif")
industrial = rast(industrial_file)


#####################################################
# Get blueberry raster
#####################################################

# Change land cover raster to 1/0 (seminatural/other)
vals = c(3)
blueberry = ifel(landscape_raster %in% vals, 1L, 0L,
            filename = paste0(bombus_path, "landscape/rasters/blueberry.tif"),
            overwrite = TRUE,
            wopt = list(datatype="INT1U"))
blueberry_file = paste0(bombus_path, "/landscape/rasters/blueberry.tif")
blueberry = rast(blueberry_file)


#####################################################
# Get other perennial raster
#####################################################

# Change land cover raster to 1/0 (seminatural/other)
# combine cranberry + perennial
vals = c(4, 12)
perennial = ifel(landscape_raster %in% vals, 1L, 0L,
            filename = paste0(bombus_path, "landscape/rasters/perennial.tif"),
            overwrite = TRUE,
            wopt = list(datatype="INT1U"))
perennial_file = paste0(bombus_path, "/landscape/rasters/perennial.tif")
perennial = rast(semi_file)


#####################################################
# Get fallow raster
#####################################################

# Change land cover raster to 1/0 (seminatural/other)
vals = c(5)
fallow = ifel(landscape_raster %in% vals, 1L, 0L,
            filename = paste0(bombus_path, "landscape/rasters/fallow.tif"),
            overwrite = TRUE,
            wopt = list(datatype="INT1U"))
fallow_file = paste0(bombus_path, "/landscape/rasters/fallow.tif")
fallow = rast(fallow_file)



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
# see idwcluster.R --- script for running on server

#############################################
# Calculate IJI per point
#############################################
foraging_80_factor = sqrt(-2*log(1-0.8))

# for mixtus
iji_mix_mean = calculateIJI(landcover.raster = landscape_raster, 
                              site.shapefile = fv_points, 
                              landcover.classification = landcover, 
                              buffer.sizes = rhom[2]*foraging_80_factor)
colnames(iji_mix_mean) = c("sample_pt", "landscape_iji_mix_mean")

# for impatiens
iji_imp_mean = calculateIJI(landcover.raster = landscape_raster, 
                                          site.shapefile = fv_points, 
                                          landcover.classification = landcover, 
                                          buffer.sizes = rhoi[2]*foraging_80_factor)
colnames(iji_imp_mean) = c("sample_pt", "landscape_iji_imp_mean")

# for queen models
iji_250 = calculateIJI(landcover.raster = landscape_raster, 
                             site.shapefile = fv_points, 
                             landcover.classification = landcover, 
                             buffer.sizes = c(250))
iji_1000 = calculateIJI(landcover.raster = landscape_raster, 
                             site.shapefile = fv_points, 
                             landcover.classification = landcover, 
                             buffer.sizes = c(1000))



#############################################
# Calculate SN per point
#############################################
sn_250 = calculateSN(landcover.raster = landscape_raster, 
                       site.shapefile = fv_points, 
                       landcover.classification = landcover, 
                       buffer.sizes = c(250))
sn_1000 = calculateSN(landcover.raster = landscape_raster, 
                        site.shapefile = fv_points, 
                        landcover.classification = landcover, 
                        buffer.sizes = c(1000))


#############################################
# Add to landscape metrics file
#############################################
saved = read.csv("analysis/calculated_metrics.csv")
updated = saved %>%
  left_join(iji_mix_mean) %>%
  left_join(iji_imp_mean) %>%
  left_join(iji_1000) %>%
  left_join(iji_250) %>%
  left_join(sn_250) %>%
  left_join(sn_1000)

write.csv(updated, "analysis/calculated_metrics.csv", row.names = FALSE)
