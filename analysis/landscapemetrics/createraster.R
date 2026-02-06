##### Create raster from shapefile
## From T Kelly
## Updated Feb 5, 2025 by J Melanson

# Load packages and source data
library(sf)
source("~/projects/def-ckremen/melanson/landscape/tyler_code/Load_Packages.R")
source("src/analysis_functions.R")
load.packages()


# Load the landscape data
lc = read_sf("~/projects/def-ckremen/melanson/landscape/delta_3000m/delta_landcover3000m.shp")
lc = lc %>% st_set_crs(900913)
lc = st_transform(lc, 32610)


raster_100 <- create.raster(lc, "lc_type", 100)
rasterVis::levelplot(raster_100)
print("res 100 raster done!")

start.time <- Sys.time()
raster_2 <- create.raster(lc, lc_types = "lc_type", resolution = 2)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken


raster_data <- raster(extent(sf_object), resolution = resolution)

#save the raster with the correct resolution
writeRaster(raster_2, paste0("~/projects/def-ckremen/melanson/landscape/full_rasters/FValley_", resolution, "res.tif"))
