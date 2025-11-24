##### Calculate distance weighted landscape metrics
### Nov 18, 2025
### J. Melanson


# Load packages
library(terra)
library(sf)
library(raster)
library(data.table)
library(dbscan)

# Set up workspace
setwd("/Users/jenna1/fv_landscapeforaging")
bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"
source("src/analysis_functions.R")

#####################################################
# Get seminatural raster
#####################################################

# Load in raster
landscape_raster = raster(paste0(bombus_path, "landscape/rasters/FValley_lc_1res.tif"))
fv_points = st_read(paste0(bombus_path, "landscape/fvbombus/fvbombus_points.shp"))
landcover = read.csv(paste0(bombus_path, "landscape/landcover.csv"))

#Change CRS to meters
crs(landscape_raster) = 900913
st_crs(fv_points) = 900913

# Change land cover raster to 1/0 (seminatural/other)
landscape_raster = rast(landscape_raster)
vals = c(2, 6, 8, 10, 16)
semi = ifel(landscape_raster %in% vals, 1L, 0L,
             filename = paste0(bombus_path, "landscape/rasters/seminat.tif"),
             overwrite = TRUE,
             wopt = list(datatype="INT1U"))
semi_file = paste0(bombus_path, "/landscape/rasters/seminat.tif")
semi = rast(semi_file)


#############################################
# Get RHO values for each species x year
#############################################
# impatiens 2022
stanFit = readRDS("analysis/foraging_modelfits/foragingmodel_3.rds")
summary = summary(stanFit)$summary

rho = c(summary["rho", "2.5%"],
        summary["rho", "mean"],
        summary["rho", "97.5%"])
rho = rho*1000


#############################################
# Calculate IDW seminatural area per point
#############################################
buffer = 1500

# for impatiens 2022
fv_points$idwSN_imp22_low = NA
fv_points$idwSN_imp22_mean = NA
fv_points$idwSN_imp22_high = NA

for (i in 1:nrow(fv_points)) {
    fv_points$idwSN_imp22_low[i] <- compute_idw_area(fv_points[i,], semi, buffer, rho[1])
    fv_points$idwSN_imp22_mean[i] <- compute_idw_area(fv_points[i,], semi, buffer, rho[1])
    fv_points$idwSN_imp22_mean[i] <- compute_idw_area(fv_points[i,], semi, buffer, rho[1])
    if (i %% 10 == 0) cat("Processed point", i, "of", nrow(fv_points), "\n")
}

write.csv(fv_points, "analysis/calculated_metrics.csv")


#############################################
# Calculate IJI per point
#############################################

# for impatiens 2022
foraging_80_factor = sqrt(-2*log(1-0.8))
iji_imp22 = calculateLandscapeMetrics(landcover.raster = landscape_raster, 
                                      site.shapefile = fv_points, 
                                      landcover.classification = landcover, 
                                      buffer.sizes = rho*foraging_80_factor)


