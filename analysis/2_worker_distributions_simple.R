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
source('src/init.R')
source('src/makeSibMaps.R')
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
bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"

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
mix2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mix2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
imp2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
imp2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")

##########################################################
### Calculate pairwise distance between traps
##########################################################

# Get trap data
fv_points = st_read(paste0(bombus_path, "landscape/fvbombus/fvbombus_points.shp"))
traps_sf_m = st_transform(fv_points, 32610)
traps_m = data.frame(sample_point = traps_sf_m$site_id,
                     trap_x = st_coordinates(traps_sf_m)[,1],
                      trap_y = st_coordinates(traps_sf_m)[,2])

sampleEffort2022 = left_join(sampleEffort2022, traps_m, relationship = "many-to-one")
sampleEffort2023 = left_join(sampleEffort2023, traps_m, relationship = "many-to-one")

# get pairwise distance matrix 2022
dmat22 = as.matrix(dist(sampleEffort2022[, c("trap_x", "trap_y")]))
idx22 = upper.tri(dmat22, diag = TRUE)

pairwise_distances22 = data.frame(trap_1 = sampleEffort2022$sample_id[row(dmat22)[idx22]],
                                 trap_2 = sampleEffort2022$sample_id[col(dmat22)[idx22]],
                                 site_1 = sampleEffort2022$site[row(dmat22)[idx22]],
                                 site_2 = sampleEffort2022$site[col(dmat22)[idx22]],
                                 distance = dmat22[idx22])

pairwise_distances22 = pairwise_distances22 %>% filter(site_1 == site_2)
pairwise_distances22$group = "Trap pairs"

# get pairwise distance matrix 2023
dmat23 = as.matrix(dist(sampleEffort2023[, c("trap_x", "trap_y")]))
idx23 = upper.tri(dmat23, diag = TRUE)

pairwise_distances23 = data.frame(trap_1 = sampleEffort2023$sample_id[row(dmat23)[idx23]],
                                  trap_2 = sampleEffort2023$sample_id[col(dmat23)[idx23]],
                                  site_1 = sampleEffort2023$site[row(dmat23)[idx23]],
                                  site_2 = sampleEffort2023$site[col(dmat23)[idx23]],
                                  distance = dmat23[idx23])

pairwise_distances23 = pairwise_distances23 %>% filter(site_1 == site_2)
pairwise_distances23$group = "Trap pairs"



##########################################################
### Calculate pairwise distance between siblings
##########################################################

mix2022_filtered = mix2022 %>% 
  left_join(traps_m, by = c("sample_pt" = "sample_point")) %>%
  filter(!str_detect(notes, "queen")) %>%
  filter(!str_detect(notes, "male")) %>%
  group_by(sibshipID) %>%
  filter(n() > 1)
dmat_msibs22 = as.matrix(dist(mix2022_filtered[, c("trap_x", "trap_y")]))
idx = upper.tri(dmat_msibs22, diag = FALSE)

sib_distances_m22 = data.frame(sib1 = mix2022_filtered$barcode_id[row(dmat_msibs22)[idx]],
                           sib2 = mix2022_filtered$barcode_id[col(dmat_msibs22)[idx]],
                          sibship_1 = mix2022_filtered$sibshipID[row(dmat_msibs22)[idx]],
                          sibship_2 = mix2022_filtered$sibshipID[col(dmat_msibs22)[idx]],
                          distance = dmat_msibs22[idx])

sib_distances_m22 = sib_distances_m22 %>% filter(sibship_1 == sibship_2)
sib_distances_m22$group = "B. mixtus sib-pairs"

mix2023_filtered = mix2023 %>% 
  left_join(traps_m, by = c("sample_pt" = "sample_point")) %>%
  filter(!str_detect(notes, "queen")) %>%
  filter(!str_detect(notes, "male")) %>%
  group_by(sibshipID) %>%
  filter(n() > 1)
dmat_msibs23 = as.matrix(dist(mix2023_filtered[, c("trap_x", "trap_y")]))
idx = upper.tri(dmat_msibs23, diag = FALSE)

sib_distances_m23 = data.frame(sib1 = mix2023_filtered$barcode_id[row(dmat_msibs23)[idx]],
                               sib2 = mix2023_filtered$barcode_id[col(dmat_msibs23)[idx]],
                               sibship_1 = mix2023_filtered$sibshipID[row(dmat_msibs23)[idx]],
                               sibship_2 = mix2023_filtered$sibshipID[col(dmat_msibs23)[idx]],
                               distance = dmat_msibs23[idx])

sib_distances_m23 = sib_distances_m23 %>% filter(sibship_1 == sibship_2)
sib_distances_m23$group = "B. mixtus sib-pairs"

imp2022_filtered = imp2022 %>% 
  left_join(traps_m, by = c("sample_pt" = "sample_point")) %>%
  filter(!str_detect(notes, "queen")) %>%
  filter(!str_detect(notes, "male")) %>%
  group_by(sibshipID) %>%
  filter(n() > 1)
dmat_isibs22 = as.matrix(dist(imp2022_filtered[, c("trap_x", "trap_y")]))
idx = upper.tri(dmat_isibs22, diag = FALSE)

sib_distances_i22 = data.frame(sib1 = imp2022_filtered$barcode_id[row(dmat_isibs22)[idx]],
                               sib2 = imp2022_filtered$barcode_id[col(dmat_isibs22)[idx]],
                               sibship_1 = imp2022_filtered$sibshipID[row(dmat_isibs22)[idx]],
                               sibship_2 = imp2022_filtered$sibshipID[col(dmat_isibs22)[idx]],
                               distance = dmat_isibs22[idx])

sib_distances_i22 = sib_distances_i22 %>% filter(sibship_1 == sibship_2)
sib_distances_i22$group = "B. impatiens sib-pairs"

imp2023_filtered = imp2023 %>% 
  left_join(traps_m, by = c("sample_pt" = "sample_point")) %>%
  filter(!str_detect(notes, "queen")) %>%
  filter(!str_detect(notes, "male")) %>%
  group_by(sibshipID) %>%
  filter(n() > 1)
dmat_isibs23 = as.matrix(dist(imp2023_filtered[, c("trap_x", "trap_y")]))
idx = upper.tri(dmat_isibs23, diag = FALSE)

sib_distances_i23 = data.frame(sib1 = imp2023_filtered$barcode_id[row(dmat_isibs23)[idx]],
                               sib2 = imp2023_filtered$barcode_id[col(dmat_isibs23)[idx]],
                               sibship_1 = imp2023_filtered$sibshipID[row(dmat_isibs23)[idx]],
                               sibship_2 = imp2023_filtered$sibshipID[col(dmat_isibs23)[idx]],
                               distance = dmat_isibs23[idx])

sib_distances_i23 = sib_distances_i23 %>% filter(sibship_1 == sibship_2)
sib_distances_i23$group = "B. impatiens sib-pairs"



##########################################################
### Make figure of pairwise distances
##########################################################

df2022 = rbind(pairwise_distances22[,c("distance", "group")], 
               sib_distances_i22[,c("distance", "group")], 
               sib_distances_m22[,c("distance", "group")])
df2022$year = 2022
df2023 = rbind(pairwise_distances23[,c("distance", "group")], 
               sib_distances_i23[,c("distance", "group")], 
               sib_distances_m23[,c("distance", "group")])
df2023$year = 2023
full_df = rbind(df2022, df2023)
full_df$distance = full_df$distance/2

ggplot(full_df, aes(x = distance)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 75,
                 color = "black") +
  facet_grid(group~year) +
  ylab("Density") +
  xlab("Distance (meters)") +
  theme_bw()

# Assign data to bins
breaks = 25 + 100*(0:13)

full_df$distance[full_df$distance < 25] = 25
full_df$dist_bins = cut(full_df$distance, 
                       breaks = breaks, 
                       right = FALSE, 
                       include.lowest = TRUE)
summary = full_df %>%
  group_by(group, dist_bins) %>%
  summarize(count = n())
sibs = summary %>% filter(group != "Trap pairs")
traps = summary %>% filter(group == "Trap pairs")
colnames(traps) = c("group", "dist_bins", "trap_count")
sibs = left_join(sibs, traps[,c("dist_bins", "trap_count")])
sibs$modified_count = sibs$count / sibs$trap_count

sibs = sibs %>%
  group_by(group) %>%
  mutate(prop = modified_count/sum(modified_count))

# Now get expected values based on 1/r^3 decay
int_r <- function(a, b) {
    log(b) - log(a)
  }
  
total = int_r(25, 1225)
perbin =  data.frame(r_start = breaks[-length(breaks)], r_end = breaks[-1])
perbin$integral = mapply(int_r, perbin$r_start, perbin$r_end)
perbin$proportion = perbin$integral / total
perbin$dist_bins = as.factor(paste0("[", perbin$r_start, ",", perbin$r_end, ")"))
perbin$dist_bins <- as.factor(
  paste0(
    "[",
    ifelse(perbin$r_start >= 1000,
           format(perbin$r_start, scientific = TRUE, digits = 3),
           perbin$r_start),
    ",",
    ifelse(perbin$r_end >= 1000,
           format(perbin$r_end, scientific = TRUE, digits = 3),
           perbin$r_end),
    ")"
  )
)

ggplot() +
  geom_point(data = sibs, aes(x = dist_bins, y = prop, color = group)) +
  geom_point(data = perbin, aes(x = dist_bins, y = proportion)) +
  theme_bw()

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
