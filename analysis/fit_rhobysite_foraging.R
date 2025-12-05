###### Fit foraging distance model to real data!
### With different rho values per landscape!
### December 2, 2025
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
#setwd("/Users/jenna1/fv_landscapeforaging")
#bombus_path = "/Users/jenna1/Documents/UBC/bombus_project"

setwd("/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging")
bombus_path = "/home/melanson/projects/def-ckremen/melanson"
source("src/analysis_functions.R")


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



# Prepare data for Stan (with function)
mix22_stan = prep_stan_simpleforaging(mixtus_sibs2022,
                                      specs2022,
                                      effort2022,
                                      samplepoints)
mix23_stan = prep_stan_simpleforaging(mixtus_sibs2023,
                                      specs2023,
                                      effort2023,
                                      samplepoints)
imp22_stan = prep_stan_simpleforaging(impatiens_sibs2022,
                                      specs2022,
                                      effort2022,
                                      samplepoints)
imp23_stan = prep_stan_simpleforaging(impatiens_sibs2023,
                                      specs2023,
                                      effort2023,
                                      samplepoints)



#select stan model to fit
stanfile = "models/rhobylandscape.stan"


# Get task ID from slurm manager
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
data_list = list(mix22_stan[[1]], mix23_stan[[1]], imp22_stan[[1]], imp23_stan[[1]])
data = data_list[[task_id]]

#fit and save model
stanFit = stan(file = stanfile,
               data = data, seed = 5838299,
               chains = 4, cores = 4,
               control = list(max_treedepth = 15),
               iter = 4000,
               verbose = TRUE)
saveRDS(stanFit, paste0("analysis/foraging_modelfits/foragingmodel_", task_id, "_rhobysite.rds"))



#############################################
# Extract values and make some plots
#############################################

# mixtus 2022
stanFitm22 = readRDS("analysis/foraging_modelfits/foragingmodel_1_rhobysite.rds")
summarym22 = rstan::summary(stanFitm22)$summary

# mixtus 2023
stanFitm23 = readRDS("analysis/foraging_modelfits/foragingmodel_2_rhobysite.rds")
summarym23 = rstan::summary(stanFitm23)$summary

# impatiens 2022
stanFiti22 = readRDS("analysis/foraging_modelfits/foragingmodel_3_rhobysite.rds")
summaryi22 = rstan::summary(stanFiti22)$summary

# impatiens 2023
stanFiti23 = readRDS("analysis/foraging_modelfits/foragingmodel_4_rhobysite.rds")
summaryi23 = rstan::summary(stanFiti23)$summary

# Plot all together
models = list(
  mixtus2022  = stanFitm22,
  mixtus2023  = stanFitm23,
  impatiens2022 = stanFiti22,
  impatiens2023 = stanFiti23
)
params = c("rho[1]", "rho[2]", "rho[3]", "rho[4]", "rho[5]", "rho[6]")

draws_df = lapply(names(models), function(nm) {
  d = as_draws_df(models[[nm]])
  d = posterior::subset_draws(d, variable = params)
  as.data.frame(d) %>% 
    mutate(model = nm)
}) %>% 
  bind_rows()

df_long = draws_df %>%
  pivot_longer(cols = all_of(params), names_to = "parameter", values_to = "value")

site_keys = data.frame(site = c("W", "SD", "ED", "NR", "HR", "PM"),
                       siteids = c(1, 2, 3, 4, 5, 6))

params_plot2 = ggplot(df_long, aes(x = value, y = parameter, fill = parameter)) +
  stat_halfeye(point_interval = median_qi, .width = c(0.95)) +
  scale_fill_discrete(labels = site_keys$site) +
  facet_wrap(~model, ncol = 1) +
  xlab("foraging length scale") +
  ylab("posterior samples") +
  theme_bw() +
  vline_at(0.5) +
  theme(legend.position = c(0.8, 0.8),
        panel.grid = element_blank(), 
        strip.background = element_blank(),
        axis.text = element_text(size = 10, color = "grey20"), 
        axis.title = element_text(size = 11, color = "grey20"),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        panel.border= element_blank(),  
        axis.ticks = element_line(color = "grey30", linewidth = 0.3),
        axis.line = element_line(color = "grey30", linewidth = 0.3))


ggsave("figures/manuscript_figures/foraging_rhobysite.jpg",
       params_plot2, units = "px", height = 3000, width = 2000)

