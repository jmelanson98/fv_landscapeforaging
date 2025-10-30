###### Colony density analysis source functions
#### October 29, 2025
#### J. Melanson


calculateLandscapeMetrics <- function(landcover.raster, 
                                      site.shapefile, 
                                      landcover.classification, 
                                      buffer.sizes = c(500)){
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
    
    
    #calculate landscape seminatural area (hedgerow + blackberry + low edge + forest + wetland)
    res1_seminat = filter(res1_lsm_c_ca, lc_type == "blackberry" | lc_type == "lowedge" | lc_type == "hedgerow" | lc_type == "forest" | lc_type == "wetland")
    res1_seminat$split_prop = res1_seminat$value/res1_seminat$buffer_area
    res1_seminat %>%
      group_by(sample_pt) %>%
      summarize(prop_seminat = sum(split_prop)) -> res1_seminat_summary
    buffer_edge = paste("prop_seminat_", buffer, sep = "")
    colnames(res1_seminat_summary) = c("sample_pt", buffer_edge)
    
    #add prop_seminat and iji for given buffer size to landscape.dat
    landscape.dat = left_join(landscape.dat, res1_seminat_summary, by = "sample_pt")
    landscape.dat = left_join(landscape.dat, res1_iji, by = "sample_pt")
    
    #update
    print(paste("Buffer size", buffer, "done.", sep = " "))
    
  }
  
  return(list(landscape.dat, res1_lsm_c_ca))
}