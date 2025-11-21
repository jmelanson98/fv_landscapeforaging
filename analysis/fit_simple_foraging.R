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
samplepoints = read.csv(paste0(bombus_path, "/raw_data/allsamplepoints.csv"), header = FALSE)


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
  mutate(count = replace_na(count, 0)) %>%
  arrange(sibshipID)
filled_counts$trap_x = as.numeric(filled_counts$trap_x)
filled_counts$trap_y = as.numeric(filled_counts$trap_y)

# Get landscape id for each colony
site_keys = data.frame(site = unique(impatiens_sibs2023$site),
               site_id = 1:length(unique(impatiens_sibs2023$site)))
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

# Make data list for stan
# DISTANCES IN KM, NOT METERS
stan_data <- list(
  C = length(colony_land$site_id),
  L = length(traps_n_2023$site),
  total_traps = nrow(traps_m_2023),
  total_obs = nrow(filled_counts),
  trap_pos = cbind(traps_m_2023$trap_x/1000, traps_m_2023$trap_y/1000), # matrix total_traps x 2, ordered by site, sample_pt
  sample_effort = traps_m_2023$total_effort,
  traps_start = array(traps_start_2023, dim = 1), # int[L]
  traps_n = array(traps_n_2023$num_traps, dim=1), # int[L]
  colony_land = array(colony_land$site_id, dim = 1), # int[C], ordered by sibshipID
  y_flat = filled_counts$count, # filled counts is ordered by site, then by sib ID, then by trap
  y_start = array(yobs_start_2023, dim = 1),
  y_n = array(yobs_n_2023$num_obs, dim = 1),
  lower_x = (min(traps_m_2023$trap_x) - 5000)/1000,
  upper_x = (max(traps_m_2023$trap_x) + 5000)/1000,
  lower_y = (min(traps_m_2023$trap_y) - 5000)/1000,
  upper_y = (max(traps_m_2023$trap_y) + 5000)/1000
)


#select stan model to fit
stanfile = "models/simple_multinomial.stan"

#fit and save model
stanFitED = stan(file = stanfile,
               data = stan_data, seed = 5838299,
               chains = 4, cores = 4,
               control = list(max_treedepth = 15),
               iter = 10000,
               verbose = TRUE)
saveRDS(stanFit, "analysis/foragingmodel_nosampleeffort.RDS")


#Plot the posteriors of some colonies
plot_list = list()
numplots = 3
legends = list()

for (i in 1:numplots){
  delta_draws = cbind(x = rstan::extract(stanFitED, pars = "delta_x")$delta[, i],
                      y = rstan::extract(stanFitED, pars = "delta_y")$delta[, i])
  traps_temp = full_join(traps_m_2023, filled_counts[filled_counts$sibshipID ==i,])
  traps_temp$count[is.na(traps_temp$count)] = 0
  traps_temp$trap_x = traps_temp$trap_x/1000
  traps_temp$trap_y = traps_temp$trap_y/1000
  
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
    #xlim(c(1220, 1225)) +
    #ylim(c(455,460)) +
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
