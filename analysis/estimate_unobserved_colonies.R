####### Fit truncated Poisson to estimate number of unobserved colonies per site
### Nov 4 2025
### J Melanson
### based on tutorial: https://freerangestats.info/blog/2018/03/20/truncated-poisson

setwd("/Users/jenna1/Documents/UBC/bombus_project/fv_landscapeforaging")

# Load packages
source('src/GeneralizedSimFunctions.R')
source('src/GenotypeSimFunctions.R')
source('src/colony_assignment_functions.R')
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


############################################
### A bit of practice with simulated data
############################################

# Simulate some truncated poisson data
set.seed(321)
a = rpois(1000, 1.3)

# truncated version of data:
b = a[ a > 0]

# graphic:
data_frame(value = c(a, b),
           variable = rep(c("Original data", "Truncated so only observations of 1 or more show up"), c(length(a), length(b)))) %>%
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 1, colour = "white") +
  facet_wrap(~variable, ncol = 1) +
  ggtitle("Comparing full and truncated datasets from a Poisson distribution") +
  labs(y = "Number of observations")

# fitting a model to original works well:
mean(a)
fitdistr(a, "Poisson")

# but obviously not if naively done to the truncated version:
mean(b)
fitdistr(b, "Poisson")


# Maximum likelihood approach
# from fitdistrplus & truncdist
# truncdist helps generate pdfs/values from truncated random variables

# adapted from https://stackoverflow.com/questions/16947799/fitting-a-lognormal-distribution-to-truncated-data-in-r
# probability density function
dtruncated_poisson <- function(x, lambda) {
  dtrunc(x, "pois", a = 0.5, b = Inf, lambda)
}

# cumulative distribution function
ptruncated_poisson <- function(x, lambda) {
  ptrunc(x, "pois", a = 0.5, b = Inf, lambda)
}

fitdistrplus::fitdist(b, "truncated_poisson", start = list(lambda = 0.5)) 
# get some weird warnings on this but it seems to work



# Now the Bayesian approach!

#select stan model to fit
stanfile = paste("models/trunc_pois.stan")

# prep data list for Stan
data = list()
data$x = b
data$lower_limit = 1
data$n = length(b)
data$lambda_mu = mean(b)
data$lambda_sigma = 1


#fit and save model
truncfit = stan(file = stanfile,
                            data = data, seed = 5838299,
                            chains = 4, cores = 4,
                            verbose = TRUE)


# visualize credible interval
plot(truncfit) + 
  ggtitle("Credibility interval for lambda, estimated by Stan from truncated data",
          "(Correct value is 1.3)") +
  labs(y = "Estimated parameters") +
  theme_minimal()



###################################################
### Now apply to simulated Bombus data
###################################################

fit_truncated_models = function(sample_size,
                                number_colonies,
                                plot_width
                                ){
  
  # set color for plots
  # Set color scheme for plots
  faded_light = "#B6D3B8"
  faded_green = "#355938"
  faded_dark = "#1B2D1C"
  dark_green = "#105215"
  light_gold = "#F7EAC0"
  gold = "#E2B41D"
  dark_gold = "#B99318"
  
  # simulate some data
  density1 = draw_simple_multi_landscape(sample_size = sample_size,
                                       num_landscape = 1,
                                       landscape_size = 1000,
                                       trapgrid_size = 300,
                                       number_traps = 25,
                                       number_colonies = number_colonies,
                                       colony_sizes = rep(100,number_colonies),
                                       rho = 75,
                                       distance_decay = "exponential")
  
  # prepare data for stan
  yobs1 = density1[[1]]
  yobs_detected1 = yobs1[rowSums(yobs1)>0,]
  colony_data1 = density1[[2]]
  colony_data1$observed_size = rowSums(yobs1)
  detected_colonies1 = colony_data1[colony_data1$observed_size > 0,]
  detected_colonies1$landscape_id = 1
  trap_data1 = density1[[3]]
  
  
  # what is the *true* number of "accessible" colonies?
  dist99 = quantile(detected_colonies1$min_dist, 0.99)
  
  # how many total colonies are within that distance?
  accessible_colonies99 = colony_data1[colony_data1$min_dist < dist99,]
  
  # make a buffer around the traps to visualize this distance on plots
  traps_sf = st_as_sf(trap_data1, coords = c("trap_x", "trap_y"), crs = NA)
  buffer99 = st_union(traps_sf) |>
    st_buffer(dist = dist99) |>
    st_cast("MULTILINESTRING")
  
  
  # load two stan files
  stanpoisson = paste("models/trunc_poisson.stan")
  stannegbin = paste("models/trunc_negbin.stan")
  
  # prep data lists for stan
  
  # first for poisson
  data_pois = list()
  data_pois$num_sites = 1
  data_pois$num_colonies = nrow(detected_colonies1)
  data_pois$lower_limit = 1
  data_pois$sib_counts = detected_colonies1$observed_size
  data_pois$site_ids = detected_colonies1$landscape_id
  data_pois$lambda_mu = mean(detected_colonies1$observed_size)
  data_pois$lambda_sigma = 1
  
  # then for negative binomial
  data_nb = list()
  data_nb$num_sites = 1
  data_nb$num_colonies = nrow(detected_colonies1)
  data_nb$lower_limit = 1
  data_nb$sib_counts = detected_colonies1$observed_size
  data_nb$site_ids = detected_colonies1$landscape_id
  data_nb$alpha_mu = mean(detected_colonies1$observed_size)
  data_nb$alpha_sigma = 1
  
  
  # fit stan models!
  poissonfit1 = stan(file = stanpoisson,
                  data = data_pois, seed = 5838299,
                  chains = 4, cores = 4,
                  verbose = TRUE)
  negbinfit1 = stan(file = stannegbin,
                    data = data_nb, seed = 5838299,
                    chains = 4, cores = 4,
                    verbose = TRUE)
  
  # plot results for poisson
  pois_extracted = as.data.frame(poissonfit1)
  
  tcoloniesplot = ggplot(pois_extracted, aes(x = total_colonies)) +
    geom_histogram(bins = 100, fill = faded_light, color = faded_green) +
    geom_vline(xintercept = nrow(accessible_colonies99), color = "red", linetype = "dashed") +
    xlab("Total number of colonies") +
    ylab("Posterior samples") +
    xlim(c(0, plot_width)) +
    theme_bw() +
    theme(legend.position = 'none', 
          panel.grid = element_blank(), 
          strip.background = element_blank(),
          axis.text = element_text(size = 10, color = "grey20"), 
          axis.title = element_text(size = 11, color = "grey20"),
          plot.margin = margin(t = 20, r = 5, b = 20, l = 5),
          panel.border= element_blank(),  
          axis.ticks = element_line(color = "grey30", linewidth = 0.3),
          axis.line = element_line(color = "grey30", linewidth = 0.3))
    
  
  ocoloniesplot = ggplot() +
    geom_point(data = colony_data1, aes(x = colony_x, y = colony_y), colour = light_gold, size = 0.2) +
    geom_point(data = detected_colonies1, aes(x = colony_x, y = colony_y), colour = dark_gold, size = 0.2) +
    geom_point(data = trap_data1, aes(x = trap_x, y = trap_y), colour = faded_dark, size = 1) +
    geom_sf(data = buffer99,
            fill = NA,
            colour = dark_gold,
            linewidth = 0.7) +
    scale_x_continuous(breaks = seq(0, 1000, 250)) +
    ylab("Northing") +
    xlab("Easting") +
    theme_bw() +
    theme(legend.position = 'none', 
          panel.grid = element_blank(), 
          strip.background = element_blank(),
          axis.text = element_text(size = 10, color = "grey20"), 
          axis.title = element_text(size = 11, color = "grey20"),
          plot.margin = margin(t = 20, b = 5, l = 0, r = 0),
          panel.border=element_rect(color = "grey30", linewidth = 0.3),  
          axis.ticks = element_line(color = "grey30", linewidth = 0.3))
  
  filtered_colonies = colony_data1[(colony_data1$observed_size > 0 | colony_data1$min_dist < dist99),]
  histplot= ggplot(filtered_colonies, aes(x = observed_size)) +
    geom_histogram(aes(fill = ifelse(observed_size == 0, "unobserved", "observed"))) +
    scale_fill_manual(values = c("unobserved" = light_gold, "observed" = dark_gold)) +
    ylab("Count of colonies") +
    xlab("Count of workers per colony") +
    theme_bw() +
    theme(
          legend.position = c(0.7, 0.7),
          legend.title = element_blank(),
          panel.grid = element_blank(), 
          strip.background = element_blank(),
          axis.text = element_text(size = 10, color = "grey20"), 
          axis.title = element_text(size = 11, color = "grey20"),
          plot.margin = margin(t = 5, b = 20, l = 10, r = 10),
          panel.border=element_rect(color = "grey30", linewidth = 0.3),  
          axis.ticks = element_line(color = "grey30", linewidth = 0.3))
  
  
  
  # plot results for negative binomial
  nb_extracted = as.data.frame(negbinfit1)
  
  tcoloniesplot_nb = ggplot(nb_extracted, aes(x = total_colonies)) +
    geom_histogram(bins = 100, fill = faded_light, color = faded_green) +
    geom_vline(xintercept = nrow(accessible_colonies99), color = "red", linetype = "dashed") +
    xlab("Total number of colonies") +
    ylab("Posterior samples") +
    xlim(c(0, plot_width)) +
    theme_bw() +
    theme(legend.position = 'none', 
          panel.grid = element_blank(), 
          strip.background = element_blank(),
          axis.text = element_text(size = 10, color = "grey20"), 
          axis.title = element_text(size = 11, color = "grey20"),
          plot.margin = margin(t = 20, r = 5, b = 20, l = 5),
          panel.border= element_blank(),  
          axis.ticks = element_line(color = "grey30", linewidth = 0.3),
          axis.line = element_line(color = "grey30", linewidth = 0.3))
  
  # final = ggdraw(tcoloniesplot) +
  #   draw_plot(ocoloniesplot, x = 0.42, y = 0.56, width = 0.18, height = 0.36) +
  #   draw_plot(histplot, x = 0.42, y = 0.12, width = 0.18, height = 0.36) +
  #   draw_plot(tcoloniesplot_nb, x = 0.6, y = 0, width = 0.4, height = 1)
  
  middle = plot_grid(ocoloniesplot, histplot, 
                     ncol = 1,
                     rel_heights = c(3,2))
  final1 = plot_grid(tcoloniesplot, middle, tcoloniesplot_nb,
                    ncol = 3,
                    rel_widths = c(1, 0.8, 1))
  return(final1)
}

panel1 = fit_truncated_models(sample_size = 200,
                              number_colonies = 1000,
                              plot_width = 3000)
panel2 = fit_truncated_models(sample_size = 200,
                              number_colonies = 500,
                              plot_width = 1500)
panel3 = fit_truncated_models(sample_size = 200,
                              number_colonies = 250,
                              plot_width = 750)
final = plot_grid(panel1, nullGrob(), panel2, nullGrob(), panel3, 
                  ncol = 1, labels = c('A', '', 'B', '', 'C'),
                  rel_heights = c(1,0.1,1,0.1,1))

ggsave("docs/appendix_figures/fit_truncated_simdata.jpg", final,
       units = "px", height = 5000, width = 3500)

###################################################
### With real data!
###################################################

# start with mixtus 2022
# load in sibships and spec data
bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"
mix2022 = read.csv("data/siblingships/mixtus_sibships_2022.csv")
mix2023 = read.csv("data/siblingships/mixtus_sibships_2023.csv")
imp2022 = read.csv("data/siblingships/impatiens_sibships_2022.csv")
imp2023 = read.csv("data/siblingships/impatiens_sibships_2023.csv")
specimenData2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022specimendata.csv"), sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023specimendata.csv"), sep = ",", header = T))

fit_truncated_nb_real = function(siblingships,
                              specimendata,
                              focal_year){
  # change site code to numeric for stan
  site_keys = data.frame(site = c("W", "SD", "ED", "NR", "HR", "PM"),
                         siteids = c(1, 2, 3, 4, 5, 6))
  

  # filter to proper samples
  # for mixtus 2022
  siblingships = filter(siblingships, 
                   notes != "male",
                   year == focal_year)
  siblingships = left_join(siblingships, specimendata)
  siblingships = left_join(siblingships, site_keys, by = "site")
  
  # reformat for model
  colonies_by_landscape = siblingships %>%
    group_by(sibshipID, siteids) %>%
    summarize(count = n())
  total_colonies = colonies_by_landscape %>%
    group_by(siteids) %>%
    summarize(total_colonies = n())
  
  
  #select stan model to fit
  stanfile = paste("models/trunc_negbin.stan")
  
  # prep data for stan
  data = list()
  data$num_sites = 6
  data$num_colonies = nrow(colonies_by_landscape)
  data$lower_limit = 1
  data$sib_counts = colonies_by_landscape$count
  data$site_ids = colonies_by_landscape$siteids
  data$alpha_mu = mean(colonies_by_landscape$count)
  data$alpha_sigma = 1
  data$col_per_site = total_colonies$total_colonies
  
  
  # fit stan model!
  fit = stan(file = stanfile,
                  data = data, seed = 5838299,
                  chains = 4, cores = 4,
                  verbose = TRUE)
  return(fit)
}


fit_truncated_pois_real = function(siblingships,
                              specimendata,
                              focal_year){
  # change site code to numeric for stan
  site_keys = data.frame(site = c("W", "SD", "ED", "NR", "HR", "PM"),
                         siteids = c(1, 2, 3, 4, 5, 6))
  
  
  # filter to proper samples
  # for mixtus 2022
  siblingships = filter(siblingships, 
                        notes != "male",
                        year == focal_year)
  siblingships = left_join(siblingships, specimendata)
  siblingships = left_join(siblingships, site_keys, by = "site")
  
  # reformat for model
  colonies_by_landscape = siblingships %>%
    group_by(sibshipID, siteids) %>%
    summarize(count = n())
  total_colonies = colonies_by_landscape %>%
    group_by(siteids) %>%
    summarize(total_colonies = n())
  
  
  #select stan model to fit
  stanfile = paste("models/trunc_poisson.stan")
  
  # prep data for stan
  data = list()
  data$num_sites = 6
  data$num_colonies = nrow(colonies_by_landscape)
  data$lower_limit = 1
  data$sib_counts = colonies_by_landscape$count
  data$site_ids = colonies_by_landscape$siteids
  data$lambda_mu = mean(colonies_by_landscape$count)
  data$lambda_sigma = 1
  data$col_per_site = total_colonies$total_colonies
  
  
  # fit stan model!
  fit = stan(file = stanfile,
             data = data, seed = 5838299,
             chains = 4, cores = 4,
             verbose = TRUE)
  return(fit)
}



# fit negative binomial to all datasets
mix22fit = fit_truncated_nb_real(mix2022, specimenData2022, focal_year = 2022)
saveRDS(mix22fit, "analysis/mixtrunc22fit.RDS")
mix23fit = fit_truncated_nb_real(mix2023, specimenData2023, focal_year = 2023)
saveRDS(mix23fit, "analysis/mixtrunc23fit.RDS")
imp22fit = fit_truncated_nb_real(imp2022, specimenData2022, focal_year = 2022)
saveRDS(imp22fit, "analysis/imptrunc22fit_nb.RDS")
imp23fit = fit_truncated_nb_real(imp2023, specimenData2023, focal_year = 2023)
saveRDS(imp23fit, "analysis/imptrunc23fit_nb.RDS")

# plot
mix22_extracted = as.data.frame(mix22fit)
mix23_extracted = as.data.frame(mix23fit)
imp22_extracted = as.data.frame(imp22fit)
imp23_extracted = as.data.frame(imp23fit)

site_keys = data.frame(site = c("W", "SD", "ED", "NR", "HR", "PM"),
                       siteids = c(1, 2, 3, 4, 5, 6))



# mixtus 2022
mix22_tc = mix22_extracted %>%
  pivot_longer(cols = c(`total_colonies[1]`,
                        `total_colonies[2]`,
                        `total_colonies[3]`,
                        `total_colonies[4]`,
                        `total_colonies[5]`,
                        `total_colonies[6]`),
               names_to = "site",
               values_to = "total_colonies")

mix22_colonies = ggplot(mix22_tc, aes(x = total_colonies, fill = site)) +
  geom_density(alpha = 0.5) +
  scale_fill_discrete(labels = site_keys$site) +
  xlab("Total number of colonies") +
  ylab("Posterior samples") +
  xlim(c(0, 5000)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8),
        panel.grid = element_blank(), 
        strip.background = element_blank(),
        axis.text = element_text(size = 10, color = "grey20"), 
        axis.title = element_text(size = 11, color = "grey20"),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        panel.border= element_blank(),  
        axis.ticks = element_line(color = "grey30", linewidth = 0.3),
        axis.line = element_line(color = "grey30", linewidth = 0.3))
ggsave("figures/manuscript_figures/mixtus2022truncnb.jpg",
       mix22_colonies, units = "px", width = 2000, height = 1000)


# mixtus 2023
mix23_tc = mix23_extracted %>%
  pivot_longer(cols = c(`total_colonies[1]`,
                        `total_colonies[2]`,
                        `total_colonies[3]`,
                        `total_colonies[4]`,
                        `total_colonies[5]`,
                        `total_colonies[6]`),
               names_to = "site",
               values_to = "total_colonies")

mix23_colonies = ggplot(mix23_tc, aes(x = total_colonies, fill = site)) +
  geom_density(alpha = 0.5) +
  scale_fill_discrete(labels = site_keys$site) +
  xlab("Total number of colonies") +
  ylab("Posterior samples") +
  xlim(c(0, 5000)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8),
        panel.grid = element_blank(), 
        strip.background = element_blank(),
        axis.text = element_text(size = 10, color = "grey20"), 
        axis.title = element_text(size = 11, color = "grey20"),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        panel.border= element_blank(),  
        axis.ticks = element_line(color = "grey30", linewidth = 0.3),
        axis.line = element_line(color = "grey30", linewidth = 0.3))
ggsave("figures/manuscript_figures/mixtus2023truncnb.jpg",
       mix23_colonies, units = "px", width = 2000, height = 1000)


# impatiens 2022
imp22_tc = imp22_extracted %>%
  pivot_longer(cols = c(`total_colonies[1]`,
                        `total_colonies[2]`,
                        `total_colonies[3]`,
                        `total_colonies[4]`,
                        `total_colonies[5]`,
                        `total_colonies[6]`),
               names_to = "site",
               values_to = "total_colonies")

imp22_colonies = ggplot(imp22_tc, aes(x = total_colonies, fill = site)) +
  geom_density(alpha = 0.5) +
  scale_fill_discrete(labels = site_keys$site) +
  xlab("Total number of colonies") +
  ylab("Posterior samples") +
  xlim(c(0, 5000)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8),
        panel.grid = element_blank(), 
        strip.background = element_blank(),
        axis.text = element_text(size = 10, color = "grey20"), 
        axis.title = element_text(size = 11, color = "grey20"),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        panel.border= element_blank(),  
        axis.ticks = element_line(color = "grey30", linewidth = 0.3),
        axis.line = element_line(color = "grey30", linewidth = 0.3))
ggsave("figures/manuscript_figures/impatiens2022truncnb.jpg",
       imp22_colonies, units = "px", width = 2000, height = 1000)


# impatiens 2023
imp23_tc = imp23_extracted %>%
  pivot_longer(cols = c(`total_colonies[1]`,
                        `total_colonies[2]`,
                        `total_colonies[3]`,
                        `total_colonies[4]`,
                        `total_colonies[5]`,
                        `total_colonies[6]`),
               names_to = "site",
               values_to = "total_colonies")

imp23_colonies = ggplot(imp23_tc, aes(x = total_colonies, fill = site)) +
  geom_density(alpha = 0.5) +
  scale_fill_discrete(labels = site_keys$site) +
  xlab("Total number of colonies") +
  ylab("Posterior samples") +
  xlim(c(0, 5000)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8),
        panel.grid = element_blank(), 
        strip.background = element_blank(),
        axis.text = element_text(size = 10, color = "grey20"), 
        axis.title = element_text(size = 11, color = "grey20"),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        panel.border= element_blank(),  
        axis.ticks = element_line(color = "grey30", linewidth = 0.3),
        axis.line = element_line(color = "grey30", linewidth = 0.3))
ggsave("figures/manuscript_figures/impatiens2023truncnb.jpg",
       imp23_colonies, units = "px", width = 2000, height = 1000)


# fit truncated poisson
mix22fit_pois = fit_truncated_pois_real(mix2022, specimenData2022, focal_year = 2022)
saveRDS(mix22fit_pois, "analysis/mixtrunc22fit_pois.RDS")
mix23fit_pois = fit_truncated_pois_real(mix2023, specimenData2023, focal_year = 2023)
saveRDS(mix23fit_pois, "analysis/mixtrunc23fit_pois.RDS")
imp22fit_pois = fit_truncated_pois_real(imp2022, specimenData2022, focal_year = 2022)
saveRDS(imp22fit_pois, "analysis/imptrunc22fit_pois.RDS")
imp23fit_pois = fit_truncated_pois_real(imp2023, specimenData2023, focal_year = 2023)
saveRDS(imp23fit_pois, "analysis/imptrunc23fit_pois.RDS")


# extract data
mix22_extracted_pois = as.data.frame(mix22fit_pois)
mix23_extracted_pois = as.data.frame(mix23fit_pois)
imp22_extracted_pois = as.data.frame(imp22fit_pois)
imp23_extracted_pois = as.data.frame(imp23fit_pois)


# mixtus 2022
mix22_tc_pois = mix22_extracted_pois %>%
  pivot_longer(cols = c(`total_colonies[1]`,
                        `total_colonies[2]`,
                        `total_colonies[3]`,
                        `total_colonies[4]`,
                        `total_colonies[5]`,
                        `total_colonies[6]`),
               names_to = "site",
               values_to = "total_colonies")

mix22_colonies_pois = ggplot(mix22_tc_pois, aes(x = total_colonies, fill = site)) +
  geom_density(alpha = 0.5) +
  scale_fill_discrete(labels = site_keys$site) +
  xlab("Total number of colonies") +
  ylab("Posterior samples") +
  xlim(c(0, 5000)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8),
        panel.grid = element_blank(), 
        strip.background = element_blank(),
        axis.text = element_text(size = 10, color = "grey20"), 
        axis.title = element_text(size = 11, color = "grey20"),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        panel.border= element_blank(),  
        axis.ticks = element_line(color = "grey30", linewidth = 0.3),
        axis.line = element_line(color = "grey30", linewidth = 0.3))
ggsave("figures/manuscript_figures/mixtus2022truncpois.jpg",
       mix22_colonies_pois, units = "px", width = 2000, height = 1000)


# mixtus 2023
mix23_tc_pois = mix23_extracted_pois %>%
  pivot_longer(cols = c(`total_colonies[1]`,
                        `total_colonies[2]`,
                        `total_colonies[3]`,
                        `total_colonies[4]`,
                        `total_colonies[5]`,
                        `total_colonies[6]`),
               names_to = "site",
               values_to = "total_colonies")

mix23_colonies_pois = ggplot(mix23_tc_pois, aes(x = total_colonies, fill = site)) +
  geom_density(alpha = 0.5) +
  scale_fill_discrete(labels = site_keys$site) +
  xlab("Total number of colonies") +
  ylab("Posterior samples") +
  xlim(c(0, 5000)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8),
        panel.grid = element_blank(), 
        strip.background = element_blank(),
        axis.text = element_text(size = 10, color = "grey20"), 
        axis.title = element_text(size = 11, color = "grey20"),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        panel.border= element_blank(),  
        axis.ticks = element_line(color = "grey30", linewidth = 0.3),
        axis.line = element_line(color = "grey30", linewidth = 0.3))
ggsave("figures/manuscript_figures/mixtus2023truncpois.jpg",
       mix23_colonies_pois, units = "px", width = 2000, height = 1000)


# impatiens 2022
imp22_tc_pois = imp22_extracted_pois %>%
  pivot_longer(cols = c(`total_colonies[1]`,
                        `total_colonies[2]`,
                        `total_colonies[3]`,
                        `total_colonies[4]`,
                        `total_colonies[5]`,
                        `total_colonies[6]`),
               names_to = "site",
               values_to = "total_colonies")

imp22_colonies_pois = ggplot(imp22_tc_pois, aes(x = total_colonies, fill = site)) +
  geom_density(alpha = 0.5) +
  scale_fill_discrete(labels = site_keys$site) +
  xlab("Total number of colonies") +
  ylab("Posterior samples") +
  xlim(c(0, 5000)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8),
        panel.grid = element_blank(), 
        strip.background = element_blank(),
        axis.text = element_text(size = 10, color = "grey20"), 
        axis.title = element_text(size = 11, color = "grey20"),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        panel.border= element_blank(),  
        axis.ticks = element_line(color = "grey30", linewidth = 0.3),
        axis.line = element_line(color = "grey30", linewidth = 0.3))
ggsave("figures/manuscript_figures/impatiens2022truncpois.jpg",
       imp22_colonies_pois, units = "px", width = 2000, height = 1000)


# impatiens 2023
imp23_tc_pois = imp23_extracted_pois %>%
  pivot_longer(cols = c(`total_colonies[1]`,
                        `total_colonies[2]`,
                        `total_colonies[3]`,
                        `total_colonies[4]`,
                        `total_colonies[5]`,
                        `total_colonies[6]`),
               names_to = "site",
               values_to = "total_colonies")

imp23_colonies_pois = ggplot(imp23_tc_pois, aes(x = total_colonies, fill = site)) +
  geom_density(alpha = 0.5) +
  scale_fill_discrete(labels = site_keys$site) +
  xlab("Total number of colonies") +
  ylab("Posterior samples") +
  xlim(c(0, 5000)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8),
        panel.grid = element_blank(), 
        strip.background = element_blank(),
        axis.text = element_text(size = 10, color = "grey20"), 
        axis.title = element_text(size = 11, color = "grey20"),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        panel.border= element_blank(),  
        axis.ticks = element_line(color = "grey30", linewidth = 0.3),
        axis.line = element_line(color = "grey30", linewidth = 0.3))
ggsave("figures/manuscript_figures/impatiens2023truncpois.jpg",
       imp23_colonies_pois, units = "px", width = 2000, height = 1000)

