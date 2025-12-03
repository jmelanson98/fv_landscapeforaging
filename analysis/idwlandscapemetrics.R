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


# Make plot of rhos for each species/year!
postm22 = as.data.frame(stanFitm22)
postm23 = as.data.frame(stanFitm23)
posti22 = as.data.frame(stanFiti22)
posti23 = as.data.frame(stanFiti23)

combinedpost =  bind_rows(
  data.frame(model = "B. mixtus 2022", rho = postm22$rho),
  data.frame(model = "B. mixtus 2023", rho = postm23$rho),
  data.frame(model = "B. impatiens 2022", rho = posti22$rho),
  data.frame(model = "B. impatiens 2023", rho = posti23$rho)
)

ggplot(combinedpost, aes(x = rho, fill = model)) +
  geom_histogram(alpha = 0.4) +
  theme_minimal()

year_difference = bind_rows(
  data.frame(species = "mixtus", rho_diff = postm22$rho - postm23$rho),
  data.frame(species = "impatiens", rho_diff = posti22$rho - posti23$rho)
)

ggplot(year_difference, aes(x = rho_diff, fill = species)) +
  geom_histogram(alpha = 0.4) +
  theme_minimal()

mix_diff = year_difference[year_difference$species == "mixtus",]
bayesp = sum(mix_diff$rho_diff > 0)/nrow(mix_diff)

species_difference = bind_rows(
  data.frame(year = "2022", rho_diff = posti22$rho - postm22$rho),
  data.frame(year = "2023", rho_diff = posti23$rho - postm23$rho)
)

ggplot(species_difference, aes(x = rho_diff, fill = year)) +
  geom_histogram(alpha = 0.4) +
  theme_minimal()

#############################################
# Calculate IDW seminatural area per point
#############################################
buffer = 1500

# for mixtus 2022
fv_points$idwSN_mix22_low = NA
fv_points$idwSN_mix22_mean = NA
fv_points$idwSN_mix22_high = NA

for (i in 1:nrow(fv_points)) {
  fv_points$idwSN_mix22_low[i] = compute_idw_area(fv_points[i,], semi, buffer, rhom22[1])
  fv_points$idwSN_mix22_mean[i] = compute_idw_area(fv_points[i,], semi, buffer, rhom22[2])
  fv_points$idwSN_mix22_high[i] = compute_idw_area(fv_points[i,], semi, buffer, rhom22[3])
  if (i %% 10 == 0) cat("Processed point", i, "of", nrow(fv_points), "\n")
}

# for mixtus 2023
fv_points$idwSN_mix23_low = NA
fv_points$idwSN_mix23_mean = NA
fv_points$idwSN_mix23_high = NA

for (i in 1:nrow(fv_points)) {
  fv_points$idwSN_mix23_low[i] = compute_idw_area(fv_points[i,], semi, buffer, rhom23[1])
  fv_points$idwSN_mix23_mean[i] = compute_idw_area(fv_points[i,], semi, buffer, rhom23[2])
  fv_points$idwSN_mix23_high[i] = compute_idw_area(fv_points[i,], semi, buffer, rhom23[3])
  if (i %% 10 == 0) cat("Processed point", i, "of", nrow(fv_points), "\n")
}

# for impatiens 2022
fv_points$idwSN_imp22_low = NA
fv_points$idwSN_imp22_mean = NA
fv_points$idwSN_imp22_high = NA

for (i in 1:nrow(fv_points)) {
    fv_points$idwSN_imp22_low[i] = compute_idw_area(fv_points[i,], semi, buffer, rhoi22[1])
    fv_points$idwSN_imp22_mean[i] = compute_idw_area(fv_points[i,], semi, buffer, rhoi22[2])
    fv_points$idwSN_imp22_high[i] = compute_idw_area(fv_points[i,], semi, buffer, rhoi22[3])
    if (i %% 10 == 0) cat("Processed point", i, "of", nrow(fv_points), "\n")
}

# for impatiens 2023
fv_points$idwSN_imp23_low = NA
fv_points$idwSN_imp23_mean = NA
fv_points$idwSN_imp23_high = NA

for (i in 1:nrow(fv_points)) {
  fv_points$idwSN_imp23_low[i] = compute_idw_area(fv_points[i,], semi, buffer, rhoi23[1])
  fv_points$idwSN_imp23_mean[i] = compute_idw_area(fv_points[i,], semi, buffer, rhoi23[2])
  fv_points$idwSN_imp23_high[i] = compute_idw_area(fv_points[i,], semi, buffer, rhoi23[3])
  if (i %% 10 == 0) cat("Processed point", i, "of", nrow(fv_points), "\n")
}


#############################################
# Calculate IJI per point
#############################################
foraging_80_factor = sqrt(-2*log(1-0.8))

# for mixtus 2022
iji_mix22_min = calculateIJI(landcover.raster = landscape_raster, 
                             site.shapefile = fv_points, 
                             landcover.classification = landcover, 
                             buffer.sizes = rhom22[1]*foraging_80_factor)
iji_mix22_mean = calculateIJI(landcover.raster = landscape_raster, 
                              site.shapefile = fv_points, 
                              landcover.classification = landcover, 
                              buffer.sizes = rhom22[2]*foraging_80_factor)
iji_mix22_max = calculateIJI(landcover.raster = landscape_raster, 
                             site.shapefile = fv_points, 
                             landcover.classification = landcover, 
                             buffer.sizes = rhom22[3]*foraging_80_factor)

# for mixtus 2023
iji_mix23_min = calculateIJI(landcover.raster = landscape_raster, 
                             site.shapefile = fv_points, 
                             landcover.classification = landcover, 
                             buffer.sizes = rhom23[1]*foraging_80_factor)
iji_mix23_mean = calculateIJI(landcover.raster = landscape_raster, 
                              site.shapefile = fv_points, 
                              landcover.classification = landcover, 
                              buffer.sizes = rhom23[2]*foraging_80_factor)
iji_mix23_max = calculateIJI(landcover.raster = landscape_raster, 
                             site.shapefile = fv_points, 
                             landcover.classification = landcover, 
                             buffer.sizes = rhom23[3]*foraging_80_factor)



# for impatiens 2022
iji_imp22_min = calculateIJI(landcover.raster = landscape_raster, 
                                      site.shapefile = fv_points, 
                                      landcover.classification = landcover, 
                                      buffer.sizes = rhoi22[1]*foraging_80_factor)
iji_imp22_mean = calculateIJI(landcover.raster = landscape_raster, 
                                          site.shapefile = fv_points, 
                                          landcover.classification = landcover, 
                                          buffer.sizes = rhoi22[2]*foraging_80_factor)
iji_imp22_max = calculateIJI(landcover.raster = landscape_raster, 
                                          site.shapefile = fv_points, 
                                          landcover.classification = landcover, 
                                          buffer.sizes = rhoi22[3]*foraging_80_factor)

# for impatiens 2023
iji_imp23_min = calculateIJI(landcover.raster = landscape_raster, 
                             site.shapefile = fv_points, 
                             landcover.classification = landcover, 
                             buffer.sizes = rhoi23[1]*foraging_80_factor)
iji_imp23_mean = calculateIJI(landcover.raster = landscape_raster, 
                              site.shapefile = fv_points, 
                              landcover.classification = landcover, 
                              buffer.sizes = rhoi23[2]*foraging_80_factor)
iji_imp23_max = calculateIJI(landcover.raster = landscape_raster, 
                             site.shapefile = fv_points, 
                             landcover.classification = landcover, 
                             buffer.sizes = rhoi23[3]*foraging_80_factor)


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
colnames(saved)[colnames(saved) == "X"] = "sample_pt"
updated = saved %>%
  left_join(iji_mix22_min) %>%
  left_join(iji_mix22_mean) %>%
  left_join(iji_mix22_max) %>%
  left_join(iji_mix23_min) %>%
  left_join(iji_mix23_mean) %>%
  left_join(iji_mix23_max) %>%
  left_join(iji_imp22_min) %>%
  left_join(iji_imp22_mean) %>%
  left_join(iji_imp22_max) %>%
  left_join(iji_imp23_min) %>%
  left_join(iji_imp23_mean) %>%
  left_join(iji_imp23_max) %>%
  left_join(iji_1000) %>%
  left_join(iji_250) %>%
  left_join(sn_250) %>%
  left_join(sn_1000)

write.csv(updated, "analysis/calculated_metrics.csv")
