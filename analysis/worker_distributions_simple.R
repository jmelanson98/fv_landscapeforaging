########
## Map sibling locations and make rough foraging distance estimates (regression approach)
## Started by J Melanson
## May 5, 2025
########

###########################
### Prepare environment
###########################
rm(list = ls())
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")

# first, load in packages
source('colony_assignments/src/init.R')
source('colony_assignments/src/joinFunctions.R')
source('colony_assignments/src/makeSibMaps.R')
library(dplyr)
library(tidyr)
library(stringr)
library(rcolony)
library(genepop)
library(data.table)
library(cowplot)
library(raster)
library(sf)
library(ggplot2)
library(viridis)
library(ggspatial)
library(geodist)


###########################
### Load in data
###########################
# raw data
specimenData2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2022specimendata.csv", sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023specimendata.csv", sep = ",", header = T, fill = TRUE))
sampleEffort2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2022sampledata.csv", sep = ","))
sampleEffort2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/bombus_project/raw_data/2023sampledata.csv", sep = ","))
samplePoints = as.data.frame(read.table("/Users/jenna1/Documents/UBC/bombus_project/raw_data/allsamplepoints.csv", sep = ","))
colnames(samplePoints) = c("sample_pt", "gps","landowner","subsite")

# sibships
mix2022 = read.csv("data/siblingships/mix_sibships_final_2022.csv")
mix2023 = read.csv("data/siblingships/mix_sibships_final_2023.csv")
imp2022 = read.csv("data/siblingships/imp_sibships_final_2022.csv")
imp2023 = read.csv("data/siblingships/imp_sibships_preliminary_2023.csv")


######################
### Make sib maps
######################

# merge 2022 + 2023 specimens
allspecs = rbind(specimenData2022[, colnames(specimenData2022) %in% c("site", "round", "sample_pt", "sample_id", "year", "barcode_id", "active_flower", "final_id", "notes", "pollen", "plate")],
                 specimenData2023[, colnames(specimenData2023) %in% c("site", "round", "sample_pt", "sample_id", "year", "barcode_id", "active_flower", "final_id", "notes", "pollen", "plate")]) %>%
  filter(final_id == "B. mixtus" | final_id == "B. impatiens")

#add lat and long to samplePoints
#first, wrangle gps coordinates into shape
samplePoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) tail(x, 1))
samplePoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) head(x, 3)[2])
samplePoints = samplePoints[,c("sample_pt", "lat", "long")]

#join coordinates to specimen data
allspecs = left_join(allspecs, samplePoints, by = "sample_pt")

# make unique cluster IDs
mix2023$ClusterIndex = mix2023$ClusterIndex + max(mix2022$ClusterIndex)
imp2023$ClusterIndex = imp2023$ClusterIndex + max(imp2022$ClusterIndex)

# join cluster IDs to specimen dataframe
tokeep = c("ClusterIndex", "OffspringID")
allsibs = rbind(mix2022[, colnames(mix2022) %in% tokeep],
                mix2023[, colnames(mix2023) %in% tokeep],
                imp2022[, colnames(imp2022) %in% tokeep],
                imp2023[, colnames(imp2023) %in% tokeep])
allspecs = left_join(allspecs, allsibs, by = c("barcode_id" = "OffspringID"))
genotypedspecs = filter(allspecs, !is.na(ClusterIndex))

#split species
mix = genotypedspecs %>% filter(final_id == "B. mixtus")
imp = genotypedspecs %>% filter(final_id == "B. impatiens")

#make plot grids for both
nonsingletonmix = mix %>% 
  group_by(ClusterIndex) %>%
  filter(n() > 1)
length(unique(nonsingletonmix$ClusterIndex))
mixgrids = makeSibMaps(mix)
ggsave("figures/manuscript_figures/mixtusmap_bothyears.jpg", mixgrids,
       width = 4000, height = 2500, units = "px")

nonsingletonimp = imp %>% 
  group_by(ClusterIndex) %>%
  filter(n() > 1)
length(unique(nonsingletonimp$ClusterIndex))
summary = imp %>%
  group_by(ClusterIndex) %>%
  summarize(n = n())
ggplot(summary, aes(x = n))+
  geom_histogram()
impgrids = makeSibMaps(imp)
ggsave("figures/manuscript_figures/impatiensmap_bothyears.jpg", impgrids,
       width = 4000, height = 2500, units = "px")

##########################################################
### Calculate average pairwise distance between siblings
##########################################################

# initiate vector to save pairwise distance values
mix_centroid = c()
imp_centroid = c()
points_centroids = c()

# loop through sib pairs (or triplets, or quadruplets...) and record distances
for(sibship in unique(nonsingletonmix$ClusterIndex)){
  sibs = nonsingletonmix[nonsingletonmix$ClusterIndex == sibship,colnames(nonsingletonmix) %in% c("barcode_id", "lat", "long")]
  sibs = st_as_sf(sibs, coords = c("long", "lat"), crs = 4326)
  centroid = st_centroid(st_union(sibs))
  dists <- as.numeric(st_distance(sibs, centroid))
  mix_centroid = c(mix_centroid, dists)
}

ggplot(as.data.frame(mix_centroid), aes(x = mix_centroid)) +
  geom_histogram() +
  labs(title = "B. mixtus", y = "count", x = "distance from sibship centroid (m)") +
  theme_minimal()


for(sibship in unique(nonsingletonimp$ClusterIndex)){
  sibs = nonsingletonimp[nonsingletonimp$ClusterIndex == sibship,colnames(nonsingletonimp) %in% c("barcode_id", "lat", "long")]
  sibs = st_as_sf(sibs, coords = c("long", "lat"), crs = 4326)
  centroid = st_centroid(st_union(sibs))
  dists <- as.numeric(st_distance(sibs, centroid))
  imp_centroid = c(imp_centroid, dists)
}

ggplot(as.data.frame(imp_centroid), aes(x = imp_centroid)) +
  geom_histogram() +
  labs(title = "B. impatiens", y = "count", x = "distance from sibship centroid (m)") +
  theme_minimal()


# for all trap locations
df_sf <- st_as_sf(samplePoints, coords = c("long", "lat"), crs = 4326)
dist_matrix <- st_distance(df_sf)
dist_matrix <- as.matrix(dist_matrix)

# for trap locations restricted by site
sf_df <- samplePoints %>%
  mutate(group_letter = substr(sample_pt, 1, 1)) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)

# Function to compute pairwise distances within each group
group_dists <- sf_df %>%
  group_split(group_letter) %>%
  map_dfr(function(group) {
    if (nrow(group) < 2) return(NULL)
    
    combn(nrow(group), 2, simplify = FALSE) %>%
      map_dfr(~{
        idx1 <- .x[1]
        idx2 <- .x[2]
        dist <- st_distance(group[idx1, ], group[idx2, ]) |> as.numeric()
        tibble(
          pt1 = group$sample_pt[idx1],
          pt2 = group$sample_pt[idx2],
          group_letter = group$group_letter[idx1],
          distance_m = dist
        )
      })
  })


write.csv(allspecs, "data/siblingships/allsibships_cleaned.csv")



### check number of sib pairs in same vs different sites
sum = nonsingletonmix %>% group_by(ClusterIndex) %>% summarize(n = n(), 
                                                               num_sites = length(unique(site)),
                                                               plate1 = plate[1],
                                                               plate2 = plate[2])
differentsites = sum %>% filter(num_sites > 1)
plates = as.data.frame(c(differentsites$plate1, differentsites$plate2))
colnames(plates) = c("plates")
platesum_falsesibs = plates %>% group_by(plates) %>% summarize(nfalse=n())

samesite = sum %>% filter(num_sites == 1)
plates = as.data.frame(c(samesite$plate1, samesite$plate2))
colnames(plates) = c("plates")
platesum_truesibs = plates %>% group_by(plates) %>% summarize(ntrue=n())

platesum = full_join(platesum_falsesibs, platesum_truesibs)
