###### Fit foraging distance model to real data!
### Get estimates of cumulative colony posteriors
### Get individual colony foraging kernels
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
library(grid)
library(gridExtra)
library(tibble)
library(sf)
library(terra)
library(stringr)

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

# Get task ID from slurm manager
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Create params table
Rmax = c(1.23, 2.68, 4.14)
species = c("mixtus", "impatiens")
model = c("normal", "uniform")
params = expand.grid(Rmax = Rmax, species = species, model = model)
params$task_id = 1:nrow(params)

#################################################################
# Prepare data for Stan
#################################################################

if (params$species[task_id] == "mixtus"){
  data = prep_stan_simpleforaging_bothyears(mixtus_sibs2022,
                                            mixtus_sibs2023,
                                            effort2022,
                                            effort2023,
                                            samplepoints)
} else if (params$species[task_id] == "impatiens"){
  data = prep_stan_simpleforaging_bothyears(impatiens_sibs2022,
                                            impatiens_sibs2023,
                                            effort2022,
                                            effort2023,
                                            samplepoints)
}

stan_data = data[[1]]
CKT = data[[2]]


#################################################################
# Think about different likelihood distributions for distance
#################################################################
# # probability density for log1p_exp
# x = seq(from = 0, to = 3, by = 0.01)
# Rmax = 2
# steepness = 100
# flog1p = function(x) {
#   return(exp(-log(1+exp((x-Rmax)*steepness))))
# }
# likelihood = flog1p(x)
# defintegral = integrate(flog1p, lower = 0, upper = Inf)
# prob_density = likelihood/Rmax
# plot(x, likelihood)
# plot(x,prob_density)
# 
# # likelihood for exponential
# loglik2 = -exp(steepness*(x-Rmax))
# lik2 = exp(loglik2)
# plot(x, lik2)
# 
# # likelihood for half normal
# library(fdrtool)
# prob_density = dhalfnorm(x, theta=sqrt(pi/2)/(Rmax/3), log = FALSE)
# plot(x, prob_density)
# hn = function(x) {
#   return(dhalfnorm(x, theta=sqrt(pi/2)/(Rmax/3), log = FALSE))
# }
# defintegral = integrate(hn, lower = 0, upper = Inf)

#################################################################
# Fit models
#################################################################

if (params$model[task_id] == "normal"){
  stanfile = "models/simple_multinomial_normaldisprior.stan"
  stan_data$Rmax = params$Rmax[task_id]
  fit = stan(file = stanfile,
             data = stan_data, seed = 5838299,
             chains = 4, cores = 4,
             iter = 2000,
             verbose = TRUE)
  modelname = paste0(params$species[task_id], "_normal_Rmax", params$Rmax[task_id])
  saveRDS(fit, paste0("analysis/kernel_locations/foraging_modelfits/", modelname, ".rds"))
} else if (params$model[task_id] == "uniform") {
  stanfile = "models/simple_multinomial_uniformdisprior.stan"
  stan_data$Rmax = params$Rmax[task_id]
  stan_data$steepness = 100
  fit = stan(file = stanfile,
             data = stan_data, seed = 5838299,
             chains = 4, cores = 4,
             iter = 2000,
             verbose = TRUE)
  modelname = paste0(params$species[task_id], "_uniform_Rmax", params$Rmax[task_id])
  saveRDS(fit, paste0("analysis/kernel_locations/foraging_modelfits/", modelname,".rds"))
}


#############################################
# Extract values and makes some plots
#############################################
# mix.fit = readRDS("analysis/foraging_modelfits/simpleforaging_multinomial1.rds")
# imp.fit = readRDS("analysis/foraging_modelfits/simpleforaging_multinomial2.rds")
# 
# # mixtus
# summarym = rstan::summary(mix.fit)$summary
# 
# rhom = summarym["rho", "mean"]
# 
# # impatiens
# summaryi = rstan::summary(imp.fit)$summary
# 
# rhoi = summaryi["rho", "mean"]
# 
# # Make plot of rhos for each species/year!
# postm = as.data.frame(mix.fit)
# posti = as.data.frame(imp.fit)
# 
# combinedpost =  bind_rows(
#   data.frame(model = "B. mixtus", rho = postm$rho),
#   data.frame(model = "B. impatiens", rho = posti$rho)
# )
# 
# posteriors = ggplot(combinedpost, aes(x = rho, fill = model, color = model)) +
#   scale_fill_manual(values = c(lm_gold, faded_strong)) +
#   scale_color_manual(values = c(dark_gold, faded_dark)) + 
#   geom_histogram() +
#   facet_wrap(~model, ncol = 1) +
#   ylab("Posterior draws") +
#   xlab(expression(rho)) +
#   theme(legend.position = "none") +
#   guides(fill = "none", color = "none") +
#   theme_minimal() +
#   theme(
#     panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
#     panel.grid = element_line()  # optional — minimal removes minor grids
#   )
# posteriors
# 
# x = seq(0, 2, by = 0.0001)
# mix_y = exp(-0.5*(x/rhom)^2)
# imp_y = exp(-0.5*(x/rhoi)^2)
# df = as.data.frame(cbind(x, mix_y, imp_y))
# 
# r50m = round(rhom*sqrt(-2*log(0.5)), 2)
# r50i = round(rhoi*sqrt(-2*log(0.5)), 2)
# r99m = round(rhom*sqrt(-2*log(0.01)), 2)
# r99i = round(rhoi*sqrt(-2*log(0.01)), 2)
# 
# foraging_decay = ggplot(df) +
#   geom_point(aes(x = x, y = mix_y), color = faded_dark) +
#   geom_point(aes(x = x, y = imp_y), color = dark_gold) +
#   vline_at(v = r50m, color = faded_dark, linetype = "dashed") +
#   vline_at(v = r50i, color = dark_gold, linetype = "dashed") +
#   vline_at(v = r99m, color = faded_dark, linetype = "dashed") +
#   vline_at(v = r99i, color = dark_gold, linetype = "dashed") +
#   annotate("text", x = 0.49, y = 0.95, label = expression(R[50*","*mix] == 0.54 ~ km), parse = TRUE, color = faded_dark, angle = 90, size = 3) +
#   annotate("text", x = 0.63, y = 0.95, label = expression(R[50*","*imp] == 0.67 ~ km), parse = TRUE, color = dark_gold, angle = 90, size = 3) +
#   annotate("text", x = 1.35, y = 0.95, label = expression(R[99*","*mix] == 1.39 ~ km), parse = TRUE, color = faded_dark, angle = 90, size = 3) +
#   annotate("text", x = 1.66, y = 0.95, label = expression(R[99*","*imp] == 1.71 ~ km), parse = TRUE, color = dark_gold, angle = 90, size = 3) +
#   ylab("Relative visitation") +
#   xlab("Distance from nest (km)") +
#   theme_bw()


#Plot the posteriors of some colonies
#colonies = c(5, 7, 11, 18)

delta_draws1 = as.data.frame(cbind(x = rstan::extract(fit, pars = "delta_x")$delta[, 1024],
                      y = rstan::extract(fit, pars = "delta_y")$delta[, 1024]))
delta_draws1$colony = "Colony 1024"
delta_draws2 = as.data.frame(cbind(x = rstan::extract(fit, pars = "delta_x")$delta[, 1373],
                     y = rstan::extract(fit, pars = "delta_y")$delta[, 1373]))
delta_draws2$colony = "Colony 1373"
delta_draws3 = as.data.frame(cbind(x = rstan::extract(fit, pars = "delta_x")$delta[, 1399],
                     y = rstan::extract(fit, pars = "delta_y")$delta[, 1399]))
delta_draws3$colony = "Colony 1399"

plot1 = ggplot(delta_draws1, aes(x = x, y = y)) +
    geom_density_2d_filled(alpha = 0.8, bins = 9) +
    scale_fill_brewer(palette = "Greens") +

    #plot trap locations / sizes / quality
    geom_point(data = CKT[CKT$stansibkey ==1024,], aes(x = trap_x, y = trap_y, size = count), colour = "black") +
    scale_size_continuous(limits = c(0,10), range = c(1, 5)) +
    xlim(c(min(CKT[CKT$stansibkey ==1024,]$trap_x)-1, max(CKT[CKT$stansibkey ==1024,]$trap_x) + 1)) +
    ylim(c(min(CKT[CKT$stansibkey ==1024,]$trap_y)-1, max(CKT[CKT$stansibkey ==1024,]$trap_y) + 1)) +
           
    #miscellaneous
    labs(x = "",
         y = "Latitude",
         title = "Colony 1024") +
    guides(fill = "none", colour = "none", size = "none") +
    theme(legend.position = "none") +
    coord_equal() +
    theme_bw()

plot2 = ggplot(delta_draws2, aes(x = x, y = y)) +
  geom_density_2d_filled(alpha = 0.8, bins = 9) +
  scale_fill_brewer(palette = "Greens") +
  
  #plot trap locations / sizes / quality
  geom_point(data = CKT[CKT$stansibkey ==1373,], aes(x = trap_x, y = trap_y, size = count), colour = "black") +
  scale_size_continuous(limits = c(0,10), range = c(1, 5)) +
  xlim(c(min(CKT[CKT$stansibkey ==1373,]$trap_x)-1, max(CKT[CKT$stansibkey ==1373,]$trap_x) + 1)) +
  ylim(c(min(CKT[CKT$stansibkey ==1373,]$trap_y)-1, max(CKT[CKT$stansibkey ==1373,]$trap_y) + 1)) +
  
  #miscellaneous
  labs(x = "Longitude",
       y = "",
       title = "Colony 1373") +
  guides(fill = "none", colour = "none", size = "none") +
  theme(legend.position = "none") +
  coord_equal() +
  theme_bw()

plot3 = ggplot(delta_draws3, aes(x = x, y = y)) +
  geom_density_2d_filled(alpha = 0.8, bins = 9) +
  scale_fill_brewer(palette = "Greens") +
  
  #plot trap locations / sizes / quality
  geom_point(data = CKT[CKT$stansibkey ==1399,], aes(x = trap_x, y = trap_y, size = count), colour = "black") +
  scale_size_continuous(limits = c(0,10), range = c(1, 5)) +
  xlim(c(min(CKT[CKT$stansibkey ==1399,]$trap_x)-1, max(CKT[CKT$stansibkey ==1399,]$trap_x) + 1)) +
  ylim(c(min(CKT[CKT$stansibkey ==1399,]$trap_y)-1, max(CKT[CKT$stansibkey ==1399,]$trap_y) + 1)) +
  
  #miscellaneous
  labs(x = "",
       y = "",
       title = "Colony 1399") +
  guides(fill = "none", colour = "none", size = "none") +
  theme(legend.position = "none") +
  coord_equal() +
  theme_bw()

#### Make a grid plot
grid = grid.arrange(plot1, nullGrob(), plot2, nullGrob(), plot3,
                    ncol = 5, widths = c(10, 1, 10,1,10))

ggsave(paste0("figures/colonyposteriors", modelname, ".jpg"), grid, height = 3000, width = 3200, units = "px")



#################################################################
# Compute raster layer of cumulative colony posteriors
#################################################################

# # get posterior draws
# all_delta_mix = as.data.frame(cbind(lon = 1000*unlist(rstan::extract(mix.fit, pars = "delta_x")),
#                     lat = 1000*unlist(rstan::extract(mix.fit, pars = "delta_y"))))
# all_delta_imp = as.data.frame(cbind(lon = 1000*unlist(rstan::extract(imp.fit, pars = "delta_x")),
#                                     lat = 1000*unlist(rstan::extract(imp.fit, pars = "delta_y"))))
# 
# 
# # convert posterior draws to spatvector
# deltavector_mix = terra::vect(all_delta_mix[,], geom=c("lon", "lat"), crs=crs(landscape_raster), keepgeom=FALSE)
# deltavector_imp = terra::vect(all_delta_imp[,], geom=c("lon", "lat"), crs=crs(landscape_raster), keepgeom=FALSE)
# 
# 
# # create and fill raster
# rmix_empty = rast(
#   xmin = xmin(deltavector_mix),
#   xmax = xmax(deltavector_mix),
#   ymin = ymin(deltavector_mix),
#   ymax = ymax(deltavector_mix),
#   resolution = 30,
#   crs = crs(deltavector_mix))
# rimp_empty = rast(
#   xmin = xmin(deltavector_imp),
#   xmax = xmax(deltavector_imp),
#   ymin = ymin(deltavector_imp),
#   ymax = ymax(deltavector_imp),
#   resolution = 30,
#   crs = crs(deltavector_imp))
# values(rmix_empty) = 0
# value(rimp_empty) = 0
# deltavector_mix$count = 1
# deltavector_imp$count = 1
# 
# rmix = rasterize(
#   deltavector_mix,
#   rmix_empty,                 
#   field = "count",
#   fun = "sum",       
#   background = 0)
# plot(rmix)
# 
# 
# rimp = rasterize(
#   deltavector_imp,
#   rimp_empty,                 
#   field = "count",
#   fun = "sum",       
#   background = 0)
# plot(rimp)


#################################################################
# Check out log-predictive density and perform LOOCV
#################################################################

mixnormmin = readRDS("~/fv_landscapeforaging/analysis/kernel_locations/foraging_modelfits/mixtus_normal_Rmax1.23.rds")
mixnormmed = readRDS("~/fv_landscapeforaging/analysis/kernel_locations/foraging_modelfits/mixtus_normal_Rmax2.68.rds")
mixnormmax = readRDS("~/fv_landscapeforaging/analysis/kernel_locations/foraging_modelfits/mixtus_normal_Rmax4.14.rds")

log_lik_draws1 = rstan::extract(mixnormmin, pars = "loglik")$loglik[,]
log_lik_1 = data.frame(iter = 1:4000,
                           log_lik = rowSums(log_lik_draws),
                       model = "min")
log_lik_draws2 = rstan::extract(mixnormmed, pars = "loglik")$loglik[,]
log_lik_2 = data.frame(iter = 1:4000,
                       log_lik = rowSums(log_lik_draws),
                       model = "med")
log_lik_draws3 = rstan::extract(mixnormmax, pars = "loglik")$loglik[,]
log_lik_3 = data.frame(iter = 1:4000,
                       log_lik = rowSums(log_lik_draws),
                       model = "max")

log_lik_df = rbind(log_lik_1, log_lik_2, log_lik_3)


ggplot(log_lik_df, aes(x = iter, y = log_lik, color = model)) +
  geom_line() +
  theme_bw()

min = loo(mixnormmin, pars = "loglik")
med = loo(mixnormmed, pars = "loglik")
max = loo(mixnormmax, pars = "loglik")
