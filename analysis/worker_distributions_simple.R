########
## Map sibling locations and make rough foraging distance estimates (regression approach)
## Started by J Melanson
## May 5, 2025
########

######################
### make sib maps
#####################

#add lat and long to samplePoints
#first, wrangle gps coordinates into shape
samplePoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) tail(x, 1))
samplePoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) head(x, 3)[2])
samplePoints = samplePoints[,c("sample_pt", "subsite", "lat", "long")]

#join with siblings
sibgroups = mixtus_wsibs %>% filter(!is.na(sample_pt)) %>% filter(!is.na(fullsib_index))
sibswcoords = merge(sibgroups, samplePoints, by = "sample_pt", all.x = TRUE)

#split 2022 and 2023
sibs22 = sibswcoords %>% filter(year == "2022")
sibs23 = sibswcoords %>% filter(year=="2023")

#make plot grids for both
source("microsatellitecode/src/makeSibMaps.R")
grids2022 = makeSibMaps(sibs22)
length(unique(sibs22$fullsib_index))

grids2023 = makeSibMaps(sibs23)
length(unique(sibs23$fullsib_index))
