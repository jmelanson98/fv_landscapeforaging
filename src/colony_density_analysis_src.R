###### Colony density analysis source functions
#### October 29, 2025
#### J. Melanson


calculateLandscapeMetrics <- function(landcover.raster, site.shapefile, landcover.classification, buffer.sizes = c(500)){
  #this function will:
  #### - calculate area of different land classes in 500m buffer around each sample_pt
  #### - calculate landscape shdi in 500m buffer around each sample_pt
  #### - return a dataframe containing prop_blueberry, prop_edge, and shdi for each sample_pt
  
  #set and change CRS
  crs(landcover.raster) <- 900913
  st_crs(site.shapefile) <- 900913
  
  # initiate landscape.dat dataframe
  landscape.dat = as.data.frame(c(fv_points2022$site_id))
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
    colnames(res1_lsm_c_ca)[colnames(res1_lsm_c_ca) == "plot_id"] <- "sample_pt"
    
    #shannon diversity
    res1_lsm_l_shdi <- sample_lsm(landcover.raster,
                                  y = site.shapefile,
                                  plot_id = site.shapefile$site_id,
                                  shape = "circle",
                                  size = buffer,
                                  return_raster = TRUE,
                                  what = "lsm_l_shdi")
    res1_shdi = res1_lsm_l_shdi[,c("plot_id", "value")]
    buffer_shdi = paste("landscape_shdi_", buffer, sep = "")
    colnames(res1_shdi) = c("sample_pt", buffer_shdi)
    
    #calculate % landscape blueberry
    res1_lsm_c_ca %>% 
      group_by(sample_pt) %>%
      summarise(buffer_area=sum(value)) -> buffer_areas
    res1_lsm_c_ca = merge(res1_lsm_c_ca, buffer_areas, by = "sample_pt")
    res1_blueberry = filter(res1_lsm_c_ca, lc_type == "blueberry")
    res1_blueberry$prop_blueberry = res1_blueberry$value/res1_blueberry$buffer_area
    res1_blueberry_sub = res1_blueberry[,c("sample_pt", "prop_blueberry")]
    buffer_blueberry = paste("prop_blueberry_", buffer, sep = "")
    colnames(res1_blueberry_sub) = c("sample_pt", buffer_blueberry)
    
    
    #calculate landscape edge area (hedgerow + blackberry + low edge)
    res1_edge = filter(res1_lsm_c_ca, lc_type == "blackberry" | lc_type == "lowedge" | lc_type == "hedgerow")
    res1_edge$split_prop = res1_edge$value/res1_edge$buffer_area
    res1_edge %>%
      group_by(sample_pt) %>%
      summarize(prop_edge = sum(split_prop)) -> res1_edge_summary
    buffer_edge = paste("prop_edge_", buffer, sep = "")
    colnames(res1_edge_summary) = c("sample_pt", buffer_edge)
    
    #add prop_blueberry, prop_edge, and shdi for given buffer size to landscape.dat
    landscape.dat = left_join(landscape.dat, res1_blueberry_sub, by = "sample_pt")
    landscape.dat = left_join(landscape.dat, res1_edge_summary, by = "sample_pt")
    landscape.dat = left_join(landscape.dat, res1_shdi, by = "sample_pt")
    
    #update
    print(paste("Buffer size", buffer, "done.", sep = " "))
    
  }
  
  return(landscape.dat)
}