##### Calculate distance weighted landscape metrics
### Nov 18, 2025
### J. Melanson


# Load packages
library(terra)
library(sf)
library(raster)

# Set up workspace
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")
bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"

###########################
# Get seminatural raster
###########################

# Load in raster
landscape_raster = raster(paste(bombus_path, "landscape/rasters/FValley_lc_1res.tif", sep = ""))
fv_points = read_sf(paste(bombus_path, "landscape/fvbombus/fvbombus_points.shp", sep = ""))
landcover = read.csv(paste(bombus_path, "landscape/landcover.csv", sep = ""))

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
semi = raster(paste(bombus_path, "landscape/rasters/seminat.tif", sep = ""))
semi = rast(semi)



#######################################
# Get coordinates of seminatural cells
#######################################
 
# get coordinates of all semi-natural cells
semi_cells = which(values(semi) == 1)
semi_coords = xyFromCell(semi, semi_cells)

# convert points to simple matrix
pts_coords = st_coordinates(fv_points)

# define inverse distance function
inv_weight = function(dist, rho) exp(-0.5*(dist/rho)^2)

# compute weighted sum per point
results = vector("numeric", nrow(pts_coords))

for (i in seq_len(nrow(pts_coords))) {
  pt <- pts_coords[i, ]
  # compute distances to all semi cells
  d <- sqrt((semi_coords[,1] - pt[1])^2 + (semi_coords[,2] - pt[2])^2)
  # keep only cells within 2 km
  keep <- d <= 2000
  if (sum(keep) == 0) { results[i] <- 0; next }
  d_keep <- d[keep]
  results[i] <- sum(inv_weight(d_keep))  # or multiply by 1 for each cell, or fractional area if needed
}

# attach results to sf points
fv_points$weighted_semi <- results
