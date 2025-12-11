####### Model colony densities at transect level
### December 2 2025
### J Melanson

setwd("/Users/jenna1/fv_landscapeforaging")

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

###################################################
### Prep data for model -- mixtus
###################################################
mix2023$sibshipID = mix2023$sibshipID + max(mix2022$sibshipID)
mixsibs = rbind(mix2022, mix2023)

# First get colony number per transect
pertransectm = mixsibs %>%
  group_by(sample_pt, year) %>%
  summarize(num_colonies = length(unique(sibshipID)),
            num_individuals = n())


# Get sample effort
effort22 = sampleData2022 %>%
  group_by(sample_point, year) %>%
  summarize(effort = 5*n())
effort23 = sampleData2023 %>%
  group_by(sample_point, year) %>%
  summarize(effort = 5*n())
effort = rbind(effort22, effort23)
pertransectm = pertransectm %>%
  full_join(effort, by = c("sample_pt" = "sample_point", "year"))

# Add transect_id
transectkey = data.frame(sample_pt = unique(pertransectm$sample_pt),
                         transect_id = 1:length(unique(pertransectm$sample_pt)))
pertransectm = left_join(pertransectm, transectkey)

# Get landscape metrics
subset = c("sample_pt", "idwSN_mix", "landscape_iji_mix_mean")
pertransectm = pertransectm %>%
  full_join(landscape_metrics[,colnames(landscape_metrics) %in% subset])
pertransectm$num_colonies[is.na(pertransectm$num_colonies)] = 0
pertransectm$num_individuals[is.na(pertransectm$num_individuals)] = 0

# Get average floral abundance
beeflowers = unique(c(specimenData2022$active_flower, specimenData2023$active_flower))
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]

columnlist = c("sample_point", "sample_id", beeflowers)
floral_df_wide = vegData2022[,colnames(vegData2022) %in% columnlist]
floral_df_long22 = floral_df_wide %>%
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
floral_df_long22$year = 2022

floral_df_wide = vegData2023[,colnames(vegData2023) %in% columnlist]
floral_df_long23 = floral_df_wide %>%
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
floral_df_long23$year = 2023


floral_df_long = rbind(floral_df_long22, floral_df_long23)

pertransectm = pertransectm %>%
  full_join(floral_df_long, by = c("sample_pt" = "sample_point", "year"))

# remove unvisited transects
pertransectm = pertransectm[!is.na(pertransectm$effort),]

# clean up single inf value
pertransectm$floral_abundance[pertransectm$floral_abundance == -Inf] = 0
pertransectm$floral_abundance[is.na(pertransectm$floral_abundance)] = 0

# convert inverse weighted m^2 to inverse weighted ha
pertransectm$idwSN_mix = pertransectm$idwSN_mix/10000

# change year to 0/1
pertransectm$year = ifelse(pertransectm$year == 2022, 0, 1)

###################################################
### Fit stan model
###################################################
stanfile = paste("models/cabund_transects.stan")

data = list()
data$T = length(unique(pertransectm$sample_pt))
data$O = nrow(pertransectm)
data$idw = pertransectm$idwSN_mix
data$iji = pertransectm$landscape_iji_mix_mean
data$effort = pertransectm$effort
data$fq = pertransectm$floral_abundance
data$yobs = pertransectm$num_colonies
data$year = pertransectm$year
data$transect_id = pertransectm$transect_id

fit = stan(file = stanfile,
           data = data, seed = 5838299,
           chains = 4, cores = 4,
           verbose = TRUE)
plot(fit, pars = c("beta_idw", "beta_iji"), ci_level = 0.8, outer_level = 0.95)
saveRDS(fit, "analysis/colony_modelfits/mixtusfit.rds")
mixtusfit = readRDS("analysis/colony_modelfits/mixtusfit.rds")


###################################################
### Prep data for model -- impatiens
###################################################

imp2023$sibshipID = imp2023$sibshipID + max(imp2022$sibshipID)
impsibs = rbind(imp2022, imp2023)

# First get colony number per transect
pertransecti = impsibs %>%
  group_by(sample_pt, year) %>%
  summarize(num_colonies = length(unique(sibshipID)),
            num_individuals = n())


# Get sample effort
effort22 = sampleData2022 %>%
  group_by(sample_point, year) %>%
  summarize(effort = 5*n())
effort23 = sampleData2023 %>%
  group_by(sample_point, year) %>%
  summarize(effort = 5*n())
effort = rbind(effort22, effort23)
pertransecti = pertransecti %>%
  full_join(effort, by = c("sample_pt" = "sample_point", "year"))

# Add transect_id
transectkey = data.frame(sample_pt = unique(pertransecti$sample_pt),
                         transect_id = 1:length(unique(pertransecti$sample_pt)))
pertransecti = left_join(pertransecti, transectkey)

# Get landscape metrics
subset = c("sample_pt", "idwSN_imp", "landscape_iji_imp_mean")
pertransecti = pertransecti %>%
  full_join(landscape_metrics[,colnames(landscape_metrics) %in% subset])
pertransecti$num_colonies[is.na(pertransecti$num_colonies)] = 0
pertransecti$num_individuals[is.na(pertransecti$num_individuals)] = 0

# Get average floral abundance
beeflowers = unique(c(specimenData2022$active_flower, specimenData2023$active_flower))
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]

columnlist = c("sample_point", "sample_id", beeflowers)
floral_df_wide = vegData2022[,colnames(vegData2022) %in% columnlist]
floral_df_long22 = floral_df_wide %>%
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
floral_df_long22$year = 2022

floral_df_wide = vegData2023[,colnames(vegData2023) %in% columnlist]
floral_df_long23 = floral_df_wide %>%
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
floral_df_long23$year = 2023


floral_df_long = rbind(floral_df_long22, floral_df_long23)

pertransecti = pertransecti %>%
  full_join(floral_df_long, by = c("sample_pt" = "sample_point", "year"))

# remove unvisited transects
pertransecti = pertransecti[!is.na(pertransecti$effort),]

# clean up single inf value
pertransecti$floral_abundance[pertransecti$floral_abundance == -Inf] = 0
pertransecti$floral_abundance[is.na(pertransecti$floral_abundance)] = 0

# convert inverse weighted m^2 to inverse weighted ha
pertransecti$idwSN_imp = pertransecti$idwSN_imp/10000

# change year to 0/1
pertransecti$year = ifelse(pertransecti$year == 2022, 0, 1)

###################################################
### Fit stan model
###################################################
stanfile = paste("models/cabund_transects.stan")

data = list()
data$T = length(unique(pertransecti$sample_pt))
data$O = nrow(pertransecti)
data$idw = pertransecti$idwSN_imp
data$iji = pertransecti$landscape_iji_imp_mean
data$effort = pertransecti$effort
data$fq = pertransecti$floral_abundance
data$yobs = pertransecti$num_colonies
data$year = pertransecti$year
data$transect_id = pertransecti$transect_id

fit = stan(file = stanfile,
           data = data, seed = 5838299,
           chains = 4, cores = 4,
           verbose = TRUE)
plot(fit, pars = c("beta_idw", "beta_iji", "beta_yr"), ci_level = 0.8, outer_level = 0.95)
saveRDS(fit, "analysis/colony_modelfits/impatiensfit.rds")
impatiensfit = readRDS("analysis/colony_modelfits/impatiensfit.rds")





# Plot all together
models = list(
  mixtus= mixtusfit,
  impatiens = impatiensfit
)
params = c("beta_iji", "beta_idw")

df_long = lapply(names(models), function(nm) {
  d = as_draws_df(models[[nm]])
  d = posterior::subset_draws(d, variable = params)
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








################################################################################
### NOW MAKE INDIVIDUAL ABUNDANCE MODELS INSTEAD OF 
################################################################################

###################################################
### Fit stan models
###################################################
stanfile = paste("models/cabund_transects.stan")

data = list()
data$T = nrow(pertransectm22)
data$idw = pertransectm22$idwSN_mix22_mean
data$iji = pertransectm22$iji_mix22_mean
data$effort = pertransectm22$effort
data$fq = pertransectm22$floral_abundance
data$yobs = pertransectm22$num_individuals

fit_ind_m22 = stan(file = stanfile,
              data = data, seed = 5838299,
              chains = 4, cores = 4,
              verbose = TRUE)
plot(fit_ind_m22)
plot(fit_ind_m22, pars = c("beta_idw", "beta_iji"), ci_level = 0.8, outer_level = 0.95)
plot(fit_ind_m22, pars = c("beta_fq"), ci_level = 0.8, outer_level = 0.95)


data = list()
data$T = nrow(pertransectm23)
data$idw = pertransectm23$idwSN_mix23_mean
data$iji = pertransectm23$iji_mix23_mean
data$effort = pertransectm23$effort
data$fq = pertransectm23$floral_abundance
data$yobs = pertransectm23$num_individuals

fit_ind_m23 = stan(file = stanfile,
                   data = data, seed = 5838299,
                   chains = 4, cores = 4,
                   verbose = TRUE)
plot(fit_ind_m23)
plot(fit_ind_m23, pars = c("beta_idw", "beta_iji"), ci_level = 0.8, outer_level = 0.95)
plot(fit_ind_m23, pars = c("beta_fq"), ci_level = 0.8, outer_level = 0.95)


data = list()
data$T = nrow(pertransecti22)
data$idw = pertransecti22$idwSN_imp22_mean
data$iji = pertransecti22$iji_imp22_mean
data$effort = pertransecti22$effort
data$fq = pertransecti22$floral_abundance
data$yobs = pertransecti22$num_individuals

fit_ind_i22 = stan(file = stanfile,
                   data = data, seed = 5838299,
                   chains = 4, cores = 4,
                   verbose = TRUE)
plot(fit_ind_i22)
plot(fit_ind_i22, pars = c("beta_idw", "beta_iji"), ci_level = 0.8, outer_level = 0.95)
plot(fit_ind_i22, pars = c("beta_fq"), ci_level = 0.8, outer_level = 0.95)


data = list()
data$T = nrow(pertransecti23)
data$idw = pertransecti23$idwSN_imp23_mean
data$iji = pertransecti23$iji_imp23_mean
data$effort = pertransecti23$effort
data$fq = pertransecti23$floral_abundance
data$yobs = pertransecti23$num_individuals

fit_ind_i23 = stan(file = stanfile,
                   data = data, seed = 5838299,
                   chains = 4, cores = 4,
                   verbose = TRUE)
plot(fit_ind_i23)
plot(fit_ind_i23, pars = c("beta_idw", "beta_iji"), ci_level = 0.8, outer_level = 0.95)
plot(fit_ind_i23, pars = c("beta_fq"), ci_level = 0.8, outer_level = 0.95)




# Plot all together
models = list(
  mixtus2022  = fit_ind_m22,
  mixtus2023  = fit_ind_m23,
  impatiens2022 = fit_ind_i22,
  impatiens2023 = fit_ind_i23
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
ggsave("figures/manuscript_figures/individuals_per_transect.jpg",
       params_plot, units = "px", height = 3000, width = 2000)




###################################################
### Use hier.part for multiple landcover classes
###################################################

data(amphipod
     )
