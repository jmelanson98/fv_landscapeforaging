makeSibMaps <- function(sibswcoords){
  
  # Load raster files
  raster_path = "/Users/jenna1/Documents/UBC/bombus_project/landscape/rasters/zipped clipped rasters"
  westham_basemap <- raster(paste(raster_path, "/westham_masked.tif", sep = ""))
  sd_basemap <- raster(paste(raster_path, "/sd_masked.tif", sep = ""))
  ed_basemap <- raster(paste(raster_path, "/ed_masked.tif", sep = ""))
  nr_basemap <- raster(paste(raster_path, "/nr_masked.tif", sep = ""))
  hr_basemap <- raster(paste(raster_path, "/hr_masked.tif", sep = ""))
  pm_basemap <- raster(paste(raster_path, "/pm_masked.tif", sep = ""))
  
  
  #jitter sibs w coords
  sibswcoords = sibswcoords %>% group_by(ClusterIndex) %>% filter(n() > 1) %>% ungroup()
  points_jittered <- sibswcoords %>%
    mutate(
      lat = as.numeric(lat) + runif(nrow(sibswcoords), -0.0005, 0.0005),
      long = as.numeric(long) + runif(nrow(sibswcoords), -0.0005, 0.0005)
    )
  
  #make input dataframe into an sf object
  site_sf <- st_as_sf(points_jittered, coords = c("long", "lat"), crs = 4326)
  
  #just bees from the site in question (sf version)
  westham_sf = filter(site_sf, site == "W")
  sd_sf = filter(site_sf, site == "SD")
  ed_sf = filter(site_sf, site == "ED")
  nr_sf = filter(site_sf, site == "NR")
  hr_sf = filter(site_sf, site == "HR")
  pm_sf = filter(site_sf, site == "PM")
  
  #just bees from the site in question (df version)
  westham_df = filter(points_jittered, site == "W")
  sd_df = filter(points_jittered, site == "SD")
  ed_df = filter(points_jittered, site == "ED")
  nr_df = filter(points_jittered, site == "NR")
  hr_df = filter(points_jittered, site == "HR")
  pm_df = filter(points_jittered, site == "PM")
  
  #downsample the rasters because they're too hefty
  westham_downsampled = projectRaster(aggregate(westham_basemap, fact = 10), crs = "+proj=longlat +datum=WGS84")
  sd_downsampled = projectRaster(aggregate(sd_basemap, fact = 10), crs = "+proj=longlat +datum=WGS84")
  ed_downsampled = projectRaster(aggregate(ed_basemap, fact = 10), crs = "+proj=longlat +datum=WGS84")
  nr_downsampled = projectRaster(aggregate(nr_basemap, fact = 10), crs = "+proj=longlat +datum=WGS84")
  hr_downsampled = projectRaster(aggregate(hr_basemap, fact = 10), crs = "+proj=longlat +datum=WGS84")
  pm_downsampled = projectRaster(aggregate(pm_basemap, fact = 10), crs = "+proj=longlat +datum=WGS84")

 
  
  westham_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(westham_downsampled, xy = TRUE), aes(x = x, y = y, fill = westham_masked)) +
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +
    geom_sf(data = westham_sf, aes(color = ClusterIndex)) +
    scale_color_viridis(option = "C") +
    geom_path(data = westham_df, aes(x = long, y = lat, group = ClusterIndex, color = ClusterIndex)) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude", title = "Westham Island") +
    theme(legend.position = "none",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.4)

  sd_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(sd_downsampled, xy = TRUE), aes(x = x, y = y, fill = sd_masked)) +
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +
    geom_sf(data = sd_sf, aes(color = ClusterIndex)) +
    scale_color_viridis(option = "C") +
    geom_path(data = sd_df, aes(x = long, y = lat, group = ClusterIndex, color = ClusterIndex)) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude", title = "South Delta") +
    theme(legend.position = "none",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.4)
  
  
  ed_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(ed_downsampled, xy = TRUE), aes(x = x, y = y, fill = ed_masked)) +
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +
    geom_sf(data = ed_sf, aes(color = ClusterIndex)) +
    scale_color_viridis(option = "C") +
    geom_path(data = ed_df, aes(x = long, y = lat, group = ClusterIndex, color = ClusterIndex)) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude", title = "East Delta") +
    theme(legend.position = "none",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.4)
  
  
  nr_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(nr_downsampled, xy = TRUE), aes(x = x, y = y, fill = nr_masked)) +
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +
    geom_sf(data = nr_sf, aes(color = ClusterIndex)) +
    scale_color_viridis(option = "C") +
    geom_path(data = nr_df, aes(x = long, y = lat, group = ClusterIndex, color = ClusterIndex)) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude", title = "Nicomekl River") +
    theme(legend.position = "none",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.4)
  
  
  hr_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(hr_downsampled, xy = TRUE), aes(x = x, y = y, fill = hr_masked)) +
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +
    geom_sf(data = hr_sf, aes(color = ClusterIndex)) +
    scale_color_viridis(option = "C") +
    geom_path(data = hr_df, aes(x = long, y = lat, group = ClusterIndex, color = ClusterIndex)) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude", title = "Harvie Road") +
    theme(legend.position = "none",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.4)
  
  
  pm_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(pm_downsampled, xy = TRUE), aes(x = x, y = y, fill = pm_masked)) +
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +
    geom_sf(data = pm_sf, aes(color = ClusterIndex)) +
    scale_color_viridis(option = "C") +
    geom_path(data = pm_df, aes(x = long, y = lat, group = ClusterIndex, color = ClusterIndex)) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude", title = "Pitt Meadows") +
    theme(legend.position = "none",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none", size = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.3)

  
  #make each plot into a grob
  g1 = ggplotGrob(westham_plot)
  g2 = ggplotGrob(sd_plot)
  g3 = ggplotGrob(ed_plot)
  g4 = ggplotGrob(nr_plot)
  g5 = ggplotGrob(hr_plot)
  g6 = ggplotGrob(pm_plot)
  
  #make spacers
  blank_plot <- ggplot() + 
    theme_void() + 
    theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
  b = ggplotGrob(blank_plot)
  
  #make grid plot of maps
  grid = plot_grid(
    g1, g2, g3,
    b, b, b,
    g5, g4, g6,
    ncol = 3,
    nrow = 3,
    rel_heights = c(1, 0.05, 1),
    rel_widths = c(1, 1, 1),
    align = "hv",
    axis = "tb"
  )
  
  return(grid)
  #final_plot <- ggdraw() +
  #  draw_plot(interactiongrid, 0.015, 0, 1, 1) +
  #  draw_plot_label(c("a", "b", "c"), 
  #                  x = c(0, 0, 0), 
   #                 y = c(0.97, 0.63, 0.3))
  #y = c(0.97, 0.73, 0.48))
  #print(final_plot)
  
  
}