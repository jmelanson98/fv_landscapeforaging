####### Model queen abundance at transect level
### December 11 2025
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
library(glmmTMB)
library(performance)
library(pbapply)
library(brms)
library(posterior)


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
specimenData2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022specimendata.csv"), sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023specimendata.csv"), sep = ",", header = T))
vegData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023vegetationdata.csv"), sep = ",", header = T))
sampleData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023sampledata.csv"), sep = ",", header = T))
landscape_metrics = read.csv("analysis/calculated_metrics.csv")


###################################################
### Prep data for model
###################################################

# Determine julian_date cutoff for "spring queens"
#add julian date to sample effort data frame
sampleData2023$date = paste(sampleData2023$day, sampleData2023$month, sampleData2023$year)
sampleData2023$date = gsub(" ", "", sampleData2023$date, fixed = TRUE)
sampleData2023$date <- as.POSIXlt(sampleData2023$date, format = "%d%b%y")
sampleData2023$julian_date = sampleData2023$date$yday

specimenData2023 = left_join(specimenData2023, sampleData2023[c("sample_id", "julian_date")])
queens2023 = filter(specimenData2023, str_detect(notes, "queen")) %>%
  filter(final_id == "B. mixtus" | final_id == "B. impatiens") %>%
  group_by(final_id, julian_date) %>%
  summarize(n=n())
ggplot(queens2023, aes(x = julian_date, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap(~final_id, ncol = 2)

springqueens = filter(specimenData2023, str_detect(notes, "queen")) %>%
  filter(final_id == "B. mixtus" | final_id == "B. impatiens") %>%
  #mutate(nest_searching = active_flower == "nest searching") %>%
  filter(round <= 20) %>% # up to and including all of round 20
  group_by(sample_id, final_id) %>%
  summarize(num_queens = n())
#springqueens$nest_searching = ifelse(springqueens$nest_searching, 1, 0)
springsampling = filter(sampleData2023, round <= 20)
samplespp = expand.grid(sample_id = springsampling$sample_id, 
                        final_id = unique(springqueens$final_id))
                        #nest_searching = c(0,1))

sprqueensampling = samplespp %>%
  left_join(springqueens) %>%
  left_join(springsampling[c("sample_id", "sample_point", "julian_date")])
sprqueensampling$num_queens[is.na(sprqueensampling$num_queens)] = 0

# Get landscape variables
landscape_vars = c("sample_pt", "prop_seminat_250", "prop_seminat_1000", "landscape_iji_1000", "landscape_iji_250")

sprqueensampling = sprqueensampling %>%
  left_join(landscape_metrics[landscape_vars], by = c("sample_point" = "sample_pt"))
sprqueensampling$prop_seminat_250[is.na(sprqueensampling$prop_seminat_250)] = 0
sprqueensampling$landscape_iji_250[is.na(sprqueensampling$landscape_iji_250)] = min(sprqueensampling$landscape_iji_250, na.rm = TRUE)
sprqueensampling$landscape_iji_1000 = sprqueensampling$landscape_iji_1000/10
sprqueensampling$landscape_iji_250 = sprqueensampling$landscape_iji_250/10

# Get floral abundance
beeflowers = unique(c(specimenData2022$active_flower, specimenData2023$active_flower))
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]

columnlist = c("sample_id", beeflowers)
vegData2023$TRRE[vegData2023$sample_id == "18_PM51"] = 2
vegData2023$VICR[vegData2023$sample_id == "18_SD28"] = 2
floral_df_wide = vegData2023[vegData2023$sample_id %in% springsampling$sample_id,colnames(vegData2023) %in% columnlist]
floral_df_long = floral_df_wide %>%
  pivot_longer(!c(sample_id), names_to = "flower", values_to = "floral_abundance") %>%
  mutate(across(-c(sample_id, flower),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(-c(sample_id, flower),
                ~ 10^(.-1))) %>%
  group_by(sample_id) %>%
  summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE)) %>%
  mutate(across(-c(sample_id),
                ~ log10(.)))

# clean up bad values
floral_df_long$floral_abundance[floral_df_long$floral_abundance == -Inf] = 0
floral_df_long$floral_abundance[is.na(floral_df_long$floral_abundance)] = 0

# Join
sprqueensampling = sprqueensampling %>%
  left_join(floral_df_long)

# Split
mixsprqueen = sprqueensampling[sprqueensampling$final_id == "B. mixtus",]
mixsprqueen$julian_date = scale(mixsprqueen$julian_date)
impsprqueen = sprqueensampling[sprqueensampling$final_id == "B. impatiens",]
impsprqueen$julian_date = scale(impsprqueen$julian_date)

# Fit model
formula = formula(num_queens ~
                        prop_seminat_1000 + landscape_iji_1000 + floral_abundance +
                        prop_seminat_250 + julian_date + I(julian_date^2) + (1|sample_point)
)

bf <- bf(formula, family="zero_inflated_poisson")
fit.mix.queen = brm(bf, mixsprqueen,
              cores=4, chains = 4,
              iter = 2000,
              control = list(max_treedepth = 15),
              verbose = TRUE
)


fit.imp.queen = brm(bf, impsprqueen,
                    cores=4, chains = 4,
                    iter = 2000,
                    control = list(max_treedepth = 15),
                    verbose = TRUE
)

p = plot(pp_check(fit))
marginalR2m = performance::r2(fit)
marginalR2i = performance::r2(fit.imp.queen)
fullmodelR2m = bayes_R2(fit)
fullmodelR2i = bayes_R2(fit.imp.queen)


# Nest searching vs foraging?
queens = specimenData2023 %>%
  filter(str_detect(notes, "queen")) %>%
  filter(final_id == "B. mixtus" | final_id == "B. impatiens") %>%
  mutate(nest_searching = active_flower == "nest searching")
table(queens$final_id, queens$nest_searching)
fisher.test(queens$final_id, queens$nest_searching)


# Make model of behaviours
queens$nest_searching = ifelse(queens$nest_searching, 1, 0)
queens = queens %>%
  left_join(landscape_metrics) %>%
  left_join(floral_df_long) %>%
  left_join(sampleData2023[c("sample_id", "julian_date")])

formula = formula(nest_searching ~
                    prop_seminat_1000 + landscape_iji_1000 + floral_abundance +
                    prop_seminat_250 + julian_date + I(julian_date^2) + (1|sample_pt)
)

bf <- bf(formula, family="bernoulli")
fit.mix.behaviour = brm(bf, queens[queens$final_id == "B. mixtus",],
                    cores=4, chains = 4,
                    iter = 4000,
                    control = list(max_treedepth = 15),
                    verbose = TRUE
)

fit.imp.behaviour = brm(bf, queens[queens$final_id == "B. impatiens",],
                        cores=4, chains = 4,
                        iter = 4000,
                        control = list(max_treedepth = 15),
                        verbose = TRUE
)

imp = as_draws_df(fit.imp.behaviour)




###############################################
# LINEAGE TURNOVER
##############################################
summary = impatiens_sibs2022 %>%
  group_by(sibshipID) %>%
  summarize(num_workers = sum(!str_detect(notes, "queen")),
            num_queens = sum(str_detect(notes, "queen"))) %>%
  filter(num_workers > 0)

springqueens_summary = impatiens_sibs2022 %>% 
  filter(str_detect(notes, "queen")) %>%
  group_by(sibshipID) %>%
  summarize(numqueens = n(),
            numsites = length(unique(sample_pt)))


data = list()
data$C = nrow(summary)
data$y = summary$num_queens

stanfile = "models/lineage.stan"

stanFit = stan(file = stanfile,
               data = data, seed = 5838299,
               chains = 4, cores = 4,
               iter = 2000,
               verbose = TRUE)
