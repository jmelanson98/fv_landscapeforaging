###### Functions used in analysis
#### October 29, 2025
#### J. Melanson


###################################################
# Functions for calculating landscape metrics
##################################################

calculateIJI <- function(landcover.raster, 
                         site.shapefile, 
                         landcover.classification, 
                         buffer.sizes){
  #this function will:
  #### - calculate area of different land classes in 500m buffer around each sample_pt
  #### - calculate landscape shdi in 500m buffer around each sample_pt
  #### - return a dataframe containing prop_blueberry, prop_edge, and shdi for each sample_pt
  
  # load packages
  library(landscapemetrics)
  
  #set and change CRS
  crs(landcover.raster) = 900913
  st_crs(site.shapefile) = 900913
  
  # initiate landscape.dat dataframe
  landscape.dat = as.data.frame(c(site.shapefile$site_id))
  colnames(landscape.dat) = c("sample_pt")
  
  
  for (buffer in buffer.sizes){
    # #area of each class in each buffer zone
    # res1_lsm_c_ca <- sample_lsm(landcover.raster,
    #                             y = site.shapefile,
    #                             plot_id = site.shapefile$site_id,
    #                             shape = "circle",
    #                             size = buffer,
    #                             return_raster = TRUE,
    #                             what = "lsm_c_ca")
    # 
    # #set class IDs
    # res1_lsm_c_ca = merge(res1_lsm_c_ca, landcover.classification[,c("lc_type", "class")], by = "class")
    # colnames(res1_lsm_c_ca)[colnames(res1_lsm_c_ca) == "plot_id"] = "sample_pt"
    # 
    # # get buffer areas
    # res1_lsm_c_ca %>% 
    #   group_by(sample_pt) %>%
    #   summarise(buffer_area=sum(value)) -> buffer_areas
    # res1_lsm_c_ca = merge(res1_lsm_c_ca, buffer_areas, by = "sample_pt")
    
    # interspersion and juxtaposition
    res1_lsm_l_iji <- sample_lsm(landcover.raster,
                                  y = site.shapefile,
                                  plot_id = site.shapefile$site_id,
                                  shape = "circle",
                                  size = buffer,
                                  return_raster = TRUE,
                                  what = "lsm_l_iji")
    res1_iji = res1_lsm_l_iji[,c("plot_id", "value")]
    buffer_iji = paste("landscape_iji_", buffer, sep = "")
    colnames(res1_iji) = c("sample_pt", buffer_iji)
    
    
    # #calculate landscape seminatural area (hedgerow + blackberry + low edge + forest + wetland)
    # res1_seminat = filter(res1_lsm_c_ca, lc_type == "blackberry" | lc_type == "lowedge" | lc_type == "hedgerow" | lc_type == "forest" | lc_type == "wetland")
    # res1_seminat$split_prop = res1_seminat$value/res1_seminat$buffer_area
    # res1_seminat %>%
    #   group_by(sample_pt) %>%
    #   summarize(prop_seminat = sum(split_prop)) -> res1_seminat_summary
    # buffer_edge = paste("prop_seminat_", buffer, sep = "")
    # colnames(res1_seminat_summary) = c("sample_pt", buffer_edge)
    
    #add prop_seminat and iji for given buffer size to landscape.dat
    #landscape.dat = left_join(landscape.dat, res1_seminat_summary, by = "sample_pt")
    landscape.dat = left_join(landscape.dat, res1_iji, by = "sample_pt")
    
    #update
    print(paste("Buffer size", buffer, "done.", sep = " "))
    
  }
  
  return(landscape.dat)
}



calculateSN <- function(landcover.raster, 
                         site.shapefile, 
                         landcover.classification, 
                         buffer.sizes){
  #this function will:
  #### - calculate area of different land classes in 500m buffer around each sample_pt
  #### - calculate landscape shdi in 500m buffer around each sample_pt
  #### - return a dataframe containing prop_blueberry, prop_edge, and shdi for each sample_pt
  
  # load packages
  library(landscapemetrics)
  
  #set and change CRS
  crs(landcover.raster) = 900913
  st_crs(site.shapefile) = 900913
  
  # initiate landscape.dat dataframe
  landscape.dat = as.data.frame(c(site.shapefile$site_id))
  colnames(landscape.dat) = c("sample_pt")
  
  
  for (buffer in buffer.sizes){
    #area of each class in each buffer zone
    res1_lsm_c_ca <- sample_lsm(landcover.raster,
                                y = site.shapefile,
                                plot_id = site.shapefile$site_id,
                                shape = "circle",
                                size = buffer,
                                return_raster = TRUE,
                                what = "lsm_c_ca")

    #set class IDs
    res1_lsm_c_ca = merge(res1_lsm_c_ca, landcover.classification[,c("lc_type", "class")], by = "class")
    colnames(res1_lsm_c_ca)[colnames(res1_lsm_c_ca) == "plot_id"] = "sample_pt"

    # get buffer areas
    res1_lsm_c_ca %>%
      group_by(sample_pt) %>%
      summarise(buffer_area=sum(value)) -> buffer_areas
    res1_lsm_c_ca = merge(res1_lsm_c_ca, buffer_areas, by = "sample_pt")
    
    # #calculate landscape seminatural area (hedgerow + blackberry + low edge + forest + wetland)
    res1_seminat = filter(res1_lsm_c_ca, lc_type == "blackberry" | lc_type == "lowedge" | lc_type == "hedgerow" | lc_type == "forest" | lc_type == "wetland")
    res1_seminat$split_prop = res1_seminat$value/res1_seminat$buffer_area
    res1_seminat %>%
      group_by(sample_pt) %>%
      summarize(prop_seminat = sum(split_prop)) -> res1_seminat_summary
    buffer_edge = paste("prop_seminat_", buffer, sep = "")
    colnames(res1_seminat_summary) = c("sample_pt", buffer_edge)
    
    #add prop_seminat for given buffer size to landscape.dat
    landscape.dat = left_join(landscape.dat, res1_seminat_summary, by = "sample_pt")
    
    #update
    print(paste("Buffer size", buffer, "done.", sep = " "))
    
  }
  
  return(landscape.dat)
}



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





#############################################
# Functions for foraging distance analysis
#############################################

# function to prep data for stan
prep_stan_simpleforaging = function(sibships,
                                    specimens,
                                    effort,
                                    samplepoints){
  
  # Create site ids
  sibships$site = as.factor(sibships$site)
  site_keys = data.frame(site = levels(sibships$site),
                         site_id = 1:length(unique(sibships$site)))
  
  # Make trap locations into correct format for calculating distance
  colnames(samplepoints)[1:2] = c("sample_pt", "coord")
  samplepoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplepoints$coord), split = " "), function(x) tail(x, 1))
  samplepoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplepoints$coord), split = " "), function(x) head(x, 3)[2])
  samplepoints$site <- stringr::str_extract(samplepoints$sample_pt, "^[A-Za-z]{1,2}")
  samplepoints = samplepoints[,colnames(samplepoints) %in% c("site", "sample_pt", "lat", "long")]
  
  # convert to meter based crs
  traps_sf = st_as_sf(samplepoints, coords = c("long", "lat"), crs = 4326)
  traps_sf_m = st_transform(traps_sf, 900913)
  traps_m = as.data.frame(cbind(traps_sf_m$sample_pt, traps_sf_m$site, st_coordinates(traps_sf_m)))
  colnames(traps_m) = c("sample_pt", "site", "trap_x", "trap_y")
  traps_m$site = as.factor(traps_m$site)
  
  # calculate sample effort per trap (2023)
  effort_sum = effort %>%
    group_by(sample_point) %>%
    summarize(total_effort = 5*n())
  traps_m = traps_m %>%
    left_join(effort_sum[,colnames(effort_sum) %in% c("sample_point", "total_effort")], by = c("sample_pt" = "sample_point")) %>%
    left_join(site_keys, by = "site") %>%
    arrange(site_id)
  
  # traps_m is a dataframe/matrix of all trap coordinates
  traps_m = traps_m[!is.na(traps_m$total_effort),]
  traps_m$trap_x = as.numeric(traps_m$trap_x)
  traps_m$trap_y = as.numeric(traps_m$trap_y)
  
  # number of traps per landscape
  traps_n = traps_m %>% 
    group_by(site_id) %>%
    summarize(num_traps = n()) %>%
    arrange(site_id)
  traps_start = cumsum(c(1, traps_n$num_traps))[1:length(traps_n$num_traps)]
  
  # Get counts of bees in traps
  counts = sibships %>%
    filter(notes != "male") %>%
    group_by(site, sample_pt, sibshipID) %>%
    summarize(count=n()) %>%
    ungroup()
  
  filled_counts = counts %>%
    select(-sample_pt) %>%
    distinct(site, sibshipID) %>%
    inner_join(traps_m, by = "site", relationship = "many-to-many") %>%
    left_join(counts, by = c("site", "sample_pt", "sibshipID")) %>%
    mutate(count = replace_na(count, 0)) %>%
    arrange(sibshipID)
  filled_counts$trap_x = as.numeric(filled_counts$trap_x)
  filled_counts$trap_y = as.numeric(filled_counts$trap_y)
  
  # Get landscape id for each colony
  colony_land = sibships %>%
    filter(notes != "male") %>%
    distinct(site, sibshipID) %>%
    left_join(site_keys, by = "site") %>%
    arrange(sibshipID)
  
  # Get number of observations per siblingship
  yobs_n = filled_counts %>% 
    group_by(sibshipID) %>%
    summarize(num_obs = n())
  
  # Get the start indices for observations of each siblingship
  yobs_start = cumsum(c(1, yobs_n$num_obs))[1:length(yobs_n$num_obs)]
  
  # Make data list for stan
  # DISTANCES IN KM, NOT METERS
  stan_data <- list(
    C = length(colony_land$site_id),
    L = length(traps_n$site_id),
    total_traps = nrow(traps_m),
    total_obs = nrow(filled_counts),
    trap_pos = cbind(traps_m$trap_x/1000, traps_m$trap_y/1000), # matrix total_traps x 2, ordered by site, sample_pt
    sample_effort = traps_m$total_effort,
    traps_start = traps_start, # int[L]
    traps_n = traps_n$num_traps, # int[L]
    colony_land = colony_land$site_id, # int[C], ordered by sibshipID
    y_flat = filled_counts$count, # filled counts is ordered by site, then by sib ID, then by trap
    y_start = yobs_start,
    y_n = yobs_n$num_obs,
    lower_x = (min(traps_m$trap_x) - 5000)/1000,
    upper_x = (max(traps_m$trap_x) + 5000)/1000,
    lower_y = (min(traps_m$trap_y) - 5000)/1000,
    upper_y = (max(traps_m$trap_y) + 5000)/1000
  )
  
  out = list(stan_data, filled_counts, traps_m)
  
  return(out)
}


# function to prep data for stan --- both years together
prep_stan_simpleforaging_bothyears = function(sibships1,
                                              sibships2,
                                              effort1,
                                              effort2,
                                              samplepoints){
  
  # Adjust sibship IDs
  sibships2$sibshipID = sibships2$sibshipID + max(sibships1$sibshipID)
  
  # Create site ids
  sibships1$site = as.factor(sibships1$site)
  sibships2$site = as.factor(sibships2$site)
  site_keys = data.frame(site = levels(sibships1$site),
                         site_id = 1:length(unique(sibships1$site)))
  
  # Make trap locations into correct format for calculating distance
  colnames(samplepoints)[1:2] = c("sample_pt", "coord")
  samplepoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplepoints$coord), split = " "), function(x) tail(x, 1))
  samplepoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplepoints$coord), split = " "), function(x) head(x, 3)[2])
  samplepoints$site <- stringr::str_extract(samplepoints$sample_pt, "^[A-Za-z]{1,2}")
  samplepoints = samplepoints[,colnames(samplepoints) %in% c("site", "sample_pt", "lat", "long")]
  
  # convert to meter based crs
  traps_sf = st_as_sf(samplepoints, coords = c("long", "lat"), crs = 4326)
  traps_sf_m = st_transform(traps_sf, 900913)
  traps_m = as.data.frame(cbind(traps_sf_m$sample_pt, traps_sf_m$site, st_coordinates(traps_sf_m)))
  colnames(traps_m) = c("sample_pt", "site", "trap_x", "trap_y")
  traps_m$site = as.factor(traps_m$site)
  
  # calculate sample effort per trap and year
  effort_sum1 = effort1 %>%
    filter(round %in% sibships1$round) %>%
    group_by(sample_point, year) %>%
    summarize(total_effort = 5*n())
  effort_sum2 = effort2 %>%
    filter(round %in% sibships2$round) %>%
    group_by(sample_point, year) %>%
    summarize(total_effort = 5*n())
  effort_sum = rbind(effort_sum1, effort_sum2)
  traps_m = traps_m %>%
    left_join(effort_sum, by = c("sample_pt" = "sample_point")) %>%
    left_join(site_keys, by = "site") %>%
    arrange(site_id)
  traps_m = traps_m[!is.na(traps_m$total_effort),]
  traps_m$trap_x = as.numeric(traps_m$trap_x)/1000
  traps_m$trap_y = as.numeric(traps_m$trap_y)/1000
  traps_m$trap_id = 1:nrow(traps_m)
  
  # Make sibships x sites df
  cxl1 = sibships1 %>%
    distinct(site, sibshipID, year)
  cxl2 = sibships2 %>%
    distinct(site, sibshipID, year)
  cxl = rbind(cxl1, cxl2)
  
  # Join sibships x sites with traps
  cxl = cxl %>% full_join(traps_m, relationship = "many-to-many")
  
  # Get counts of bees in traps
  counts1 = sibships1 %>%
    filter(notes != "male") %>%
    group_by(site, sample_pt, year, sibshipID) %>%
    summarize(count=n()) %>%
    ungroup()
  counts2 = sibships2 %>%
    filter(notes != "male") %>%
    group_by(site, sample_pt, year, sibshipID) %>%
    summarize(count=n()) %>%
    ungroup()
  counts = rbind(counts1, counts2)
  
  
  filled_counts = cxl %>%
    full_join(counts) %>%
    mutate(count = replace_na(count, 0)) %>%
    arrange(sibshipID)
  filled_counts$yn = ifelse(filled_counts$count > 0, 1, 0)
  filled_counts$pointyear = paste0(filled_counts$year, filled_counts$sample_pt)
  
  # Get the start indices for observations of each siblingship
  lengths = filled_counts %>%
    group_by(sibshipID) %>%
    summarize(length = n()) %>%
    arrange(sibshipID)
  starts = cumsum(c(1, lengths$length))[1:length(lengths$length)]
  
  # Make data list for stan
  # DISTANCES IN KM, NOT METERS
  stan_data <- list(
    C = length(unique(filled_counts$sibshipID)),
    K = length(unique(filled_counts$pointyear)),
    O = nrow(filled_counts),
    starts = starts,
    lengths = lengths$length,
    trap_pos = cbind(filled_counts$trap_x, filled_counts$trap_y),
    colony_id = filled_counts$sibshipID,
    trap_id = filled_counts$trap_id,
    sample_effort = filled_counts$total_effort,
    y_obs = filled_counts$count,
    yn = filled_counts$yn,
    lower_x = (min(filled_counts$trap_x) - 5),
    upper_x = (max(filled_counts$trap_x) + 5),
    lower_y = (min(filled_counts$trap_y) - 5),
    upper_y = (max(filled_counts$trap_y) + 5)
  )
  
  out = list(stan_data, filled_counts)
  
  return(out)
}