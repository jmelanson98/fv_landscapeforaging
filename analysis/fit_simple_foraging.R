###### Fit foraging distance model to real data!
### November 14, 2025
### J. Melanson

##### Load packages #####
library(rstan)
library(matrixStats)
library(sp)
library(gstat)
library(ggplot2)
library(reshape2)
library(raster)
library(rasterVis)
library(parallel)
library(future)
library(furrr)
library(dplyr)
library(tidyr)
library(gridExtra)
library(tibble)
library(sf)

##### Set Environment #####
setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging") # local
options(mc.cores = parallel::detectCores())
bombus_path = "/Users/jenna1/Documents/UBC/bombus_project"


# Load in data
mixtus_sibs2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mixtus_sibs2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
impatiens_sibs2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
impatiens_sibs2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
specs2022 = read.csv(paste0(bombus_path, "/raw_data/2022specimendata.csv"))
specs2023 = read.csv(paste0(bombus_path, "/raw_data/2023specimendata.csv"))
effort2022 = read.csv(paste0(bombus_path, "/raw_data/2022sampledata.csv"))
effort2023 = read.csv(paste0(bombus_path, "/raw_data/2023sampledata.csv"))
veg2022 = read.csv(paste0(bombus_path, "/raw_data/2022vegetationdata.csv"))
veg2023 = read.csv(paste0(bombus_path, "/raw_data/2023vegetationdata.csv"))
samplepoints = read.csv(paste0(bombus_path, "/raw_data/allsamplepoints.csv"), header = FALSE)

#impatiens_sibs2023 = read.csv(paste0(bombus_path, "/impatiens_sibships_2023.csv"))

# Make trap locations into correct format for calculating distance
colnames(samplepoints)[1:2] = c("sample_pt", "coord")
samplepoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplepoints$coord), split = " "), function(x) tail(x, 1))
samplepoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplepoints$coord), split = " "), function(x) head(x, 3)[2])
samplepoints$site <- stringr::str_extract(samplepoints$sample_pt, "^[A-Za-z]{1,2}")
samplepoints = samplepoints[,colnames(samplepoints) %in% c("site", "sample_pt", "lat", "long")]

# convert to meter based crs
traps_sf = st_as_sf(samplepoints, coords = c("long", "lat"), crs = 4326)
traps_sf_m = st_transform(traps_sf, 3005)  # BC Albers
traps_m = as.data.frame(cbind(traps_sf_m$sample_pt, traps_sf_m$site, st_coordinates(traps_sf_m)))
colnames(traps_m) = c("sample_pt", "site", "trap_x", "trap_y")
traps_m$site = as.factor(traps_m$site)

# calculate sample effort per trap (2023)
effort_sum23 = effort2023 %>%
  group_by(sample_point) %>%
  summarize(total_effort = 5*n())
traps_m = left_join(traps_m, effort_sum23[,colnames(effort_sum23) %in% c("sample_point", "total_effort")], by = c("sample_pt" = "sample_point"))

# Prepare data for Stan
# traps_m is a dataframe/matrix of all trap coordinates
traps_m_2023 = traps_m[!is.na(traps_m$total_effort),]
traps_m_2023$trap_x = as.numeric(traps_m_2023$trap_x)
traps_m_2023$trap_y = as.numeric(traps_m_2023$trap_y)

# number of traps per landscape
traps_n_2023 = traps_m_2023 %>% 
  group_by(site) %>%
  summarize(num_traps = n())
traps_start_2023 = cumsum(c(1, traps_n_2023$num_traps))[1:length(traps_n_2023$num_traps)]

# Get counts of bees in traps
counts = impatiens_sibs2023 %>%
  filter(notes != "male") %>%
  group_by(site, sample_pt, sibshipID) %>%
  summarize(count=n()) %>%
  ungroup()

filled_counts = counts %>%
  select(-sample_pt) %>%
  distinct(site, sibshipID) %>%
  inner_join(traps_m_2023, by = "site", relationship = "many-to-many") %>%
  left_join(counts, by = c("site", "sample_pt", "sibshipID")) %>%
  mutate(count = replace_na(count, 0))
filled_counts$trap_x = as.numeric(filled_counts$trap_x)
filled_counts$trap_y = as.numeric(filled_counts$trap_y)

# Get landscape id for each colony
site_keys = data.frame(site = unique(impatiens_sibs2023$site),
               site_id = 1:6)
colony_land = impatiens_sibs2023 %>%
  filter(notes != "male") %>%
  distinct(site, sibshipID) %>%
  left_join(site_keys, by = "site") %>%
  arrange(sibshipID)

# Get number of observations per siblingship
yobs_n_2023 = filled_counts %>% 
  group_by(sibshipID) %>%
  summarize(num_obs = n())

# Get the start indices for observations of each siblingship
yobs_start_2023 = cumsum(c(1, yobs_n_2023$num_obs))[1:length(yobs_n_2023$num_obs)]

# Get upper and lower coordinate limits for each landscape
site_centroids = filled_counts %>%
  group_by(site) %>%
  summarise(center_x = mean(trap_x),
            center_y = mean(trap_y))

# Make data list for stan
# DISTANCES IN KM, NOT METERS
stan_data <- list(
  C = length(colony_land$site_id),
  L = length(traps_n_2023$site),
  total_traps = nrow(traps_m_2023),
  total_obs = nrow(filled_counts),
  trap_pos = cbind(traps_m_2023$trap_x/1000, traps_m_2023$trap_y/1000), # matrix total_traps x 2, ordered by site, sample_pt
  traps_start = traps_start_2023, # int[L]
  traps_n = traps_n_2023$num_traps, # int[L]
  colony_land = colony_land$site_id, # int[C], ordered by sibshipID
  y_flat = filled_counts$count, # filled counts is ordered by site, then by sib ID, then by trap
  y_start = yobs_start_2023,
  y_n = yobs_n_2023$num_obs,
  lower_x = (min(traps_m_2023$trap_x) - 2000)/1000,
  upper_x = (max(traps_m_2023$trap_x) + 2000)/1000,
  lower_y = (min(traps_m_2023$trap_y) - 2000)/1000,
  upper_y = (max(traps_m_2023$trap_y) + 2000)/1000,
  site_centroids = cbind(site_centroids$center_x/1000,site_centroids$center_y/1000)
  
)


#select stan model to fit
stanfile = "models/simple_multinomial.stan"

#fit and save model
stanFit_ht = stan(file = stanfile,
               data = stan_data, seed = 5838299,
               chains = 4, cores = 4,
               control = list(max_treedepth = 15),
               init = "0",
               verbose = TRUE)
saveRDS(stanFit, "analysis/foragingmodel_badtreedepth.RDS")


#Plot the posteriors of some colonies
plot_list = list()
numplots = 12
legends = list()

for (i in 1:numplots){
  delta_draws = cbind(x = rstan::extract(stanFit, pars = "delta_x")$delta[, i],
                      y = rstan::extract(stanFit, pars = "delta_y")$delta[, i])
  traps_temp = full_join(traps_m_2023, filled_counts[filled_counts$sibshipID ==i,])
  traps_temp$count[is.na(traps_temp$count)] = 0
  
  p = ggplot(delta_draws, aes(x = x, y = y)) +
    geom_density_2d_filled(alpha = 0.8) +
    
    #plot trap locations / sizes / quality
    geom_point(data = traps_temp, aes(x = trap_x, y = trap_y, size = count, colour = "red")) +
    scale_size_continuous(limits = c(0,10), range = c(1, 5)) +
    
    #miscellaneous
    labs(title = paste("Colony", i), 
         size = "Number of Captures",
         level = "Colony Posterior") +
    guides(colour = "none") +
    coord_equal() +
    theme_bw()
  
  # save legend
  g <- ggplotGrob(p)
  legend_index <- which(g$layout$name == "guide-box-right")
  legend <- g$grobs[[legend_index]]
  
  # remove legend from plot
  p <- p + theme(legend.position = "none")
  
  #save plot
  plot_list[[i]] = p
  legends[[1]] = legend
}

fig = grid.arrange(grobs = plot_list, ncol = 3)
fig = grid.arrange(fig, legends[[1]], ncol = 2, widths = c(4,1))
ggsave(paste(results_path, "/colony_posteriorsGQ.jpg", sep = ""), fig, height = 3000, width = 4000, units = "px")
