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

###########################
# Get seminatural raster
###########################

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





#############################################
# Calculate IDW seminatural area per point
#############################################

semi_file = paste0(bombus_path, "landscape/rasters/seminat.tif")
semi = rast(semi_file)

# function to compute IDW natural area for a single point

compute_idw_area = function(points, raster, buffer, rho) {
  # create buffer around point
  buf = vect(st_buffer(points, dist = buffer))
  
  # extract values and coordinates for points inside buffer
  vals_coords = terra::extract(raster, buf, cells = TRUE, xy = TRUE)
  vals = vals_coords[, "FValley_lc_1res"]
  coords = vals_coords[, c("x", "y")]
  
  # get focal point coords
  pt_coords = st_coordinates(points)
  
  # calculate vector of distances from focal point
  d = sqrt((coords[,1] - pt_coords[1])^2 + (coords[,2] - pt_coords[2])^2)
  
  # compute weights
  w = exp(-0.5 * (d/rho)^2)
  
  # weighted sum
  pixel_area = res(raster)[1] * res(raster)[2]
  dw_area = sum(vals * w, na.rm = TRUE) * pixel_area
  
  return(dw_area)
  
}


# loop over all points
fv_points$idwSemiNat = NA
buffer = 1500
rho = 480

for (i in 1:nrow(fv_points)) {
  fv_points$idwSemiNat[i] <- compute_idw_area(fv_points[i,], semi, buffer, rho)
  if (i %% 10 == 0) cat("Processed point", i, "of", nrow(fv_points), "\n")
}
