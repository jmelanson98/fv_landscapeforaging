####### Model colony densities at transect level
### December 2 2025
### J Melanson

setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")

# Load packages
source('src/GeneralizedSimFunctions.R')
source('src/GenotypeSimFunctions.R')
source('src/colony_assignment_functions.R')
source('src/analysis_functions.R')
library(tidyverse)
library(scales)
library(fitdistrplus)
library(rstan)
library(truncdist)
library(matrixStats)
library(bayesplot)
library(cowplot)
library(sf)
library(grid)

# Set color scheme for plots
faded_pale = "#D2E4D4"
faded_light = "#B6D3B8"
faded_medium = "#609F65"
faded_strong = "#4E8353"
faded_green = "#355938"
faded_dark = "#1B2D1C"


light_gold = "#F7EAC0"
lm_gold = "#F2DC97"
medium_gold = "#ECCF6F"
gold = "#E2B41D"
dark_gold = "#B99318"
darker_gold = "#907313"

###################################################
### Load in data
###################################################
bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"
mix2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mix2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
imp2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
imp2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
specimenData2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022specimendata.csv"), sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023specimendata.csv"), sep = ",", header = T))
vegData2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022vegetationdata.csv"), sep = ",", header = T))
vegData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023vegetationdata.csv"), sep = ",", header = T))
sampleData2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022sampledata.csv"), sep = ",", header = T))
sampleData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023sampledata.csv"), sep = ",", header = T))
landscape_metrics = read.csv("analysis/calculated_metrics.csv")

rename = c("X", "sample_pt", "site_id", "site_name", "type", "site_type", "transect", 
             "latitude", "longitude", "project_co", "subsite", "geometry",
             "iji_mix22_min", "iji_mix22_mean", "iji_mix22_max", "iji_mix23_min", 
             "iji_mix23_mean", "iji_mix23_max", "iji_imp22_min", "iji_imp22_mean", 
             "iji_imp22_max", "iji_imp23_min", "iji_imp23_mean", "iji_imp23_max",
             "iji_1000", "iji_250", "sn_250", "sn_1000"
             )
keep = colnames(landscape_metrics)[29:ncol(landscape_metrics)]
newcols = c(rename, keep)
colnames(landscape_metrics) = newcols

###################################################
### Prep data for model
###################################################

# First get colony number per transect
pertransect = mix2022 %>%
  group_by(sample_pt) %>%
  summarize(num_colonies = length(unique(sibshipID)))

# Get landscape metrics
subset = c("sample_pt", "idwSN_mix22_mean", "iji_mix22_mean")
pertransect = pertransect %>%
  full_join(landscape_metrics[,colnames(landscape_metrics) %in% subset])
pertransect$num_colonies[is.na(pertransect$num_colonies)] = 0

# Get sample effort
effort = sampleData2022 %>%
  group_by(sample_point) %>%
  summarize(effort = 5*n())
pertransect = pertransect %>%
  full_join(effort, by = c("sample_pt" = "sample_point"))

# Get average floral abundance
beeflowers = unique(c(specimenData2022$active_flower, specimenData2023$active_flower))
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]

columnlist = c("sample_point", "sample_id", beeflowers)
floral_df_wide = vegData2022[,colnames(vegData2022) %in% columnlist]
floral_df_long = floral_df_wide %>%
  pivot_longer(!c(sample_point, sample_id), names_to = "flower", values_to = "floral_abundance") %>%
  mutate(across(-c(sample_point, sample_id, flower),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(-c(sample_point, sample_id, flower),
                ~ 10^(.-1))) %>%
  group_by(sample_point, sample_id) %>%
  summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE)) %>%
  group_by(sample_point) %>%
  summarize(floral_abundance = mean(floral_abundance, na.rm = TRUE)) %>%
  mutate(across(-c(sample_point),
                ~ log10(.)))

pertransect = pertransect %>%
  full_join(floral_df_long, by = c("sample_pt" = "sample_point"))

# remove unvisited transects
pertransect = pertransect[!is.na(pertransect$effort),]

# clean up single inf value
pertransect$floral_abundance[pertransect$floral_abundance == -Inf] = 0

# convert inverse weighted m^2 to inverse weighted ha
pertransect$idwSN_mix22_mean = pertransect$idwSN_mix22_mean/10000


###################################################
### Fit stan model
###################################################
stanfile = paste("models/cabund_transects.stan")

data = list()
data$T = nrow(pertransect)
data$idw = pertransect$idwSN_mix22_mean
data$iji = pertransect$iji_mix22_mean
data$effort = pertransect$effort
data$fq = pertransect$floral_abundance
data$yobs = pertransect$num_colonies

fit = stan(file = stanfile,
           data = data, seed = 5838299,
           chains = 4, cores = 4,
           verbose = TRUE)
plot(fit, pars = c("beta_idw", "beta_iji", "beta_eff"), ci_level = 0.8, outer_level = 0.95)



###################################################
### Prep data for model -- mixtus 2023
###################################################

# First get colony number per transect
pertransect = mix2023 %>%
  group_by(sample_pt) %>%
  summarize(num_colonies = length(unique(sibshipID)))

# Get landscape metrics
subset = c("sample_pt", "idwSN_mix23_mean", "iji_mix23_mean")
pertransect = pertransect %>%
  full_join(landscape_metrics[,colnames(landscape_metrics) %in% subset])
pertransect$num_colonies[is.na(pertransect$num_colonies)] = 0

# Get sample effort
effort = sampleData2023 %>%
  group_by(sample_point) %>%
  summarize(effort = 5*n())
pertransect = pertransect %>%
  full_join(effort, by = c("sample_pt" = "sample_point"))

# Get average floral abundance
beeflowers = unique(c(specimenData2022$active_flower, specimenData2023$active_flower))
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]

columnlist = c("sample_point", "sample_id", beeflowers)
floral_df_wide = vegData2023[,colnames(vegData2023) %in% columnlist]
floral_df_long = floral_df_wide %>%
  pivot_longer(!c(sample_point, sample_id), names_to = "flower", values_to = "floral_abundance") %>%
  mutate(across(-c(sample_point, sample_id, flower),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(-c(sample_point, sample_id, flower),
                ~ 10^(.-1))) %>%
  group_by(sample_point, sample_id) %>%
  summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE)) %>%
  group_by(sample_point) %>%
  summarize(floral_abundance = mean(floral_abundance, na.rm = TRUE)) %>%
  mutate(across(-c(sample_point),
                ~ log10(.)))

pertransect = pertransect %>%
  full_join(floral_df_long, by = c("sample_pt" = "sample_point"))

# remove unvisited transects
pertransect = pertransect[!is.na(pertransect$effort),]

# clean up single inf value
pertransect$floral_abundance[pertransect$floral_abundance == -Inf] = 0
pertransect = pertransect[!is.na(pertransect$floral_abundance),]

# convert inverse weighted m^2 to inverse weighted ha
#pertransect$idwSN_mix23_mean = pertransect$idwSN_mix23_mean/10000


###################################################
### Fit stan model
###################################################
stanfile = paste("models/cabund_transects.stan")

data = list()
data$T = nrow(pertransect)
data$idw = pertransect$idwSN_mix23_mean
data$iji = pertransect$iji_mix23_mean
data$effort = pertransect$effort
data$fq = pertransect$floral_abundance
data$yobs = pertransect$num_colonies

fit23 = stan(file = stanfile,
           data = data, seed = 5838299,
           chains = 4, cores = 4,
           verbose = TRUE)
plot(fit23)
plot(fit23, pars = c("beta_idw", "beta_iji"), ci_level = 0.8, outer_level = 0.95)



###################################################
### Prep data for model -- impatiens 2022
###################################################

# First get colony number per transect
pertransect = imp2022 %>%
  group_by(sample_pt) %>%
  summarize(num_colonies = length(unique(sibshipID)))

# Get landscape metrics
subset = c("sample_pt", "idwSN_imp22_mean", "iji_imp22_mean")
pertransect = pertransect %>%
  full_join(landscape_metrics[,colnames(landscape_metrics) %in% subset])
pertransect$num_colonies[is.na(pertransect$num_colonies)] = 0

# Get sample effort
effort = sampleData2022 %>%
  group_by(sample_point) %>%
  summarize(effort = 5*n())
pertransect = pertransect %>%
  full_join(effort, by = c("sample_pt" = "sample_point"))

# Get average floral abundance
beeflowers = unique(c(specimenData2022$active_flower, specimenData2023$active_flower))
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]

columnlist = c("sample_point", "sample_id", beeflowers)
floral_df_wide = vegData2022[,colnames(vegData2022) %in% columnlist]
floral_df_long = floral_df_wide %>%
  pivot_longer(!c(sample_point, sample_id), names_to = "flower", values_to = "floral_abundance") %>%
  mutate(across(-c(sample_point, sample_id, flower),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(-c(sample_point, sample_id, flower),
                ~ 10^(.-1))) %>%
  group_by(sample_point, sample_id) %>%
  summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE)) %>%
  group_by(sample_point) %>%
  summarize(floral_abundance = mean(floral_abundance, na.rm = TRUE)) %>%
  mutate(across(-c(sample_point),
                ~ log10(.)))

pertransect = pertransect %>%
  full_join(floral_df_long, by = c("sample_pt" = "sample_point"))

# remove unvisited transects
pertransect = pertransect[!is.na(pertransect$effort),]

# clean up single inf value
pertransect$floral_abundance[pertransect$floral_abundance == -Inf] = 0
pertransect = pertransect[!is.na(pertransect$floral_abundance),]

# convert inverse weighted m^2 to inverse weighted ha
pertransect$idwSN_imp22_mean = pertransect$idwSN_imp22_mean/10000


###################################################
### Fit stan model
###################################################
stanfile = paste("models/cabund_transects.stan")

data = list()
data$T = nrow(pertransect)
data$idw = pertransect$idwSN_imp22_mean
data$iji = pertransect$iji_imp22_mean
data$effort = pertransect$effort
data$fq = pertransect$floral_abundance
data$yobs = pertransect$num_colonies

fiti22 = stan(file = stanfile,
           data = data, seed = 5838299,
           chains = 4, cores = 4,
           verbose = TRUE)
plot(fiti22, pars = c("beta_idw", "beta_iji"), ci_level = 0.8, outer_level = 0.95)


###################################################
### Prep data for model -- impatiens 2023
###################################################

# First get colony number per transect
pertransect = imp2023 %>%
  group_by(sample_pt) %>%
  summarize(num_colonies = length(unique(sibshipID)))

# Get landscape metrics
subset = c("sample_pt", "idwSN_imp23_mean", "iji_imp23_mean")
pertransect = pertransect %>%
  full_join(landscape_metrics[,colnames(landscape_metrics) %in% subset])
pertransect$num_colonies[is.na(pertransect$num_colonies)] = 0

# Get sample effort
effort = sampleData2023 %>%
  group_by(sample_point) %>%
  summarize(effort = 5*n())
pertransect = pertransect %>%
  full_join(effort, by = c("sample_pt" = "sample_point"))

# Get average floral abundance
beeflowers = unique(c(specimenData2022$active_flower, specimenData2023$active_flower))
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]

columnlist = c("sample_point", "sample_id", beeflowers)
floral_df_wide = vegData2023[,colnames(vegData2023) %in% columnlist]
floral_df_long = floral_df_wide %>%
  pivot_longer(!c(sample_point, sample_id), names_to = "flower", values_to = "floral_abundance") %>%
  mutate(across(-c(sample_point, sample_id, flower),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(-c(sample_point, sample_id, flower),
                ~ 10^(.-1))) %>%
  group_by(sample_point, sample_id) %>%
  summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE)) %>%
  group_by(sample_point) %>%
  summarize(floral_abundance = mean(floral_abundance, na.rm = TRUE)) %>%
  mutate(across(-c(sample_point),
                ~ log10(.)))

pertransect = pertransect %>%
  full_join(floral_df_long, by = c("sample_pt" = "sample_point"))

# remove unvisited transects
pertransect = pertransect[!is.na(pertransect$effort),]

# clean up single inf value
pertransect$floral_abundance[pertransect$floral_abundance == -Inf] = 0
pertransect = pertransect[!is.na(pertransect$floral_abundance),]

# convert inverse weighted m^2 to inverse weighted ha
pertransect$idwSN_imp23_mean = pertransect$idwSN_imp23_mean/10000


###################################################
### Fit stan model
###################################################
stanfile = paste("models/cabund_transects.stan")

data = list()
data$T = nrow(pertransect)
data$idw = pertransect$idwSN_imp23_mean
data$iji = pertransect$iji_imp23_mean
data$effort = pertransect$effort
data$fq = pertransect$floral_abundance
data$yobs = pertransect$num_colonies

fiti23 = stan(file = stanfile,
             data = data, seed = 5838299,
             chains = 4, cores = 4,
             verbose = TRUE)
plot(fiti23)
plot(fiti23, pars = c("beta_idw", "beta_iji"), ci_level = 0.8, outer_level = 0.95)
plot(fiti23, pars = c("beta_fq"), ci_level = 0.8, outer_level = 0.95)




# Plot all together
models = list(
  mixtus2022  = fit,
  mixtus2023  = fit23,
  impatiens2022 = fiti22,
  impatiens2023 = fiti23
)
params = c("beta_iji", "beta_idw")

df_long = lapply(names(models), function(nm) {
  d <- as_draws_df(models[[nm]])
  d <- posterior::subset_draws(d, variable = params)
  as.data.frame(d) %>% 
    mutate(model = nm)
}) %>% 
  bind_rows()

df_long = df_long %>%
  pivot_longer(cols = all_of(params), names_to = "parameter", values_to = "value")

params_plot = ggplot(df_long, aes(x = value, y = parameter, fill = parameter)) +
  scale_fill_manual(values = c(faded_strong, lm_gold)) +
  stat_halfeye(point_interval = median_qi, .width = c(0.95)) +
  facet_wrap(~model, ncol = 1) +
  theme_minimal()
ggsave("figures/manuscript_figures/colony_density_per_transect.jpg",
       params_plot, units = "px", height = 3000, width = 2000)
