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


# Load in raster
landscape_raster = raster(paste(bombus_path, "landscape/rasters/FValley_lc_1res.tif", sep = ""))
fv_points = read_sf(paste(bombus_path, "landscape/fvbombus/fvbombus_points.shp", sep = ""))
landcover = read_csv(paste(bombus_path, "landscape/landcover.csv", sep = ""))

#Change CRS to meters
crs(landscape_raster) = 900913
st_crs(fv_points) = 900913

# Change land cover raster to 1/0 (seminatural/other)
semi = (landscape_raster == 3 | landscape_raster == | landscape_raster == )


# Build moving window raster
radius = 1000
rho = 600

res_xy = res(landscape_raster)
if (!all(abs(res_xy - res_xy[1]) < .Machine$double.eps)) {
  warning("raster has non-square cells; code assumes roughly equal x/y resolution")
}
res_cell = res_xy[1]
cells = ceiling(radius / res_cell)
grid = expand.grid(i = -cells:cells, j = -cells:cells)
dists = sqrt((grid$i * res_cell)^2 + (grid$j * res_cell)^2)

# add weights
weights_vec = numeric(length(dists))
# set center to zero and use neighbouring cells only
center_idx = which(dists == 0)
dists[center_idx] = NA   # or tiny value?
weights_vec = ifelse(dists <= radius, exp(-dists^2/(2*rho^2)), 0)

