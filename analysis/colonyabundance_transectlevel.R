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
subset = c("sample_pt", "landscape_iji_mix_mean", "idwSN_mix", "idwHAY_mix", "idwIND_mix", "idwURB_mix", "idwANN_mix", "idwBLU_mix", "idwPER_mix", "idwFAL_mix")
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
landscape_mix = c("idwSN_mix", "idwHAY_mix", "idwIND_mix", "idwURB_mix", "idwANN_mix", "idwBLU_mix", "idwPER_mix", "idwFAL_mix")
for (i in landscape_mix){
  pertransectm[i] = pertransectm[i]/10000
}

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
subset = c("sample_pt", "landscape_iji_imp_mean", "idwSN_imp", "idwHAY_imp", "idwIND_imp", "idwURB_imp", "idwANN_imp", "idwBLU_imp", "idwPER_imp", "idwFAL_imp")
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
# convert inverse weighted m^2 to inverse weighted ha
landscape_imp = c("idwSN_imp", "idwHAY_imp", "idwIND_imp", "idwURB_imp", "idwANN_imp", "idwBLU_imp", "idwPER_imp", "idwFAL_imp")
for (i in landscape_imp){
  pertransecti[i] = pertransecti[i]/10000
}

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



###################################################
### Use hier.part for multiple landcover classes
###################################################



# mixtus

# fuck hier.part, we ball
all_subsets_mix = unlist(lapply(0:length(landscape_mix),
                             function(k) combn(landscape_mix, k, simplify=FALSE)),
                      recursive=FALSE)


fit_models_mix =  function(varset) {
  # varset: character vector of predictor names
  if (length(varset) == 0) {
    rhs <- "1"
  } else {
    rhs <- paste(varset, collapse = " + ")
  }
  form = as.formula(paste0("num_colonies ~ ", rhs, " + offset(log(effort)) + (1|sample_pt)"))
  print(form)
  fit = glmmTMB(form, family = nbinom2, data = pertransectm)
  r2 = performance::r2(fit)
  # performance::r2 returns tibble with r2_marginal and r2_conditional
  return(r2$R2_marginal)   # use marginal R2 for partitioning fixed effects
}

mixmetrics = pbapply::pblapply(all_subsets_mix, fit_models_mix)




# impatiens

# fuck hier.part, we ball
all_subsets_imp = unlist(lapply(0:length(landscape_imp),
                                function(k) combn(landscape_imp, k, simplify=FALSE)),
                         recursive=FALSE)
fit_models_imp =  function(varset) {
  # varset: character vector of predictor names
  if (length(varset) == 0) {
    rhs <- "1"
  } else {
    rhs <- paste(varset, collapse = " + ")
  }
  form = as.formula(paste0("num_colonies ~ ", rhs, " + offset(log(effort)) + (1|sample_pt)"))
  print(form)
  fit = glmmTMB(form, family = nbinom2, data = pertransecti)
  r2 = performance::r2(fit)
  # performance::r2 returns tibble with r2_marginal and r2_conditional
  return(r2$R2_marginal)   # use marginal R2 for partitioning fixed effects
}
impmetrics = pbapply::pblapply(all_subsets_imp, fit_models_imp)





# store mapping
subset_names = sapply(all_subsets_imp, function(s) if(length(s)==0) "(Intercept)" else paste(s, collapse = "+"))
res_df = data.frame(subset = subset_names, R2 = unlist(impmetrics), stringsAsFactors = FALSE)
R2_lookup = setNames(res_df$R2, res_df$subset)

# Function to compute delta_j(S) = I(S) - I(S\\{j})
delta_j_S <- function(S_vec, j) {
  S_name = if (length(S_vec)==0) "(Intercept)" else paste(S_vec, collapse="+")
  S_minus_j = setdiff(S_vec, j)
  minus_name = if (length(S_minus_j)==0) "(Intercept)" else paste(S_minus_j, collapse="+")
  I_S = R2_lookup[[S_name]]
  I_Sm = R2_lookup[[minus_name]]
  return(I_S - I_Sm)
}

# For each predictor j, average deltas across all S that contain j
indep_contrib = numeric(length(landscape_imp))
names(indep_contrib) = landscape_imp
for (j in landscape_imp) {
  subsets_with_j = Filter(function(s) j %in% s, all_subsets_imp)
  deltas = vapply(subsets_with_j, function(S) delta_j_S(S, j), numeric(1))
  indep_contrib[j] <- mean(deltas, na.rm = TRUE)   # hier.part uses average over 2^(p-1) subsets
}

# Normalize or present absolute contributions
indep_contrib
sum(indep_contrib)  # should equal total explained R2 (approximately, up to numeric issues)


xcan = pertransectm[landscape_mix]
comb.data = data.frame(xcan[, current.comb])
colnames(comb.data) <- colnames(xcan)[current.comb]
data <- data.frame(y, comb.data)
depv <- names(data)[1]
n.comb <- dim(comb.data)[2]
xs <- vector("character", n.comb)
for (i in 1:(n.comb)) xs[i] <- paste(names(comb.data)[i], "+",
                                         sep = "")
xs[n.comb+1] = "log(effort)"
xss <- paste(xs, collapse = " ", sep = "")
formu <- stats::formula(paste(depv, "~", xss, sep = ""))



###################################################
### Fit models with all predictors in brms
###################################################
formula.mix = formula(num_colonies ~
                        idwSN_mix + idwHAY_mix + idwIND_mix + idwURB_mix +
                        idwANN_mix + idwBLU_mix + idwPER_mix + idwFAL_mix +
                        landscape_iji_mix_mean + year + offset(log(effort)) +
                        (1|sample_pt)
)

bf.mix <- bf(formula.mix, family="negbinomial")
fit.mix = brm(bf.mix, pertransectm,
                      cores=4, chains = 4,
              verbose = TRUE
)
p = plot(pp_check(fit.mix))
p + xlim(0, 50)
fullmodelR2 = bayes_R2(fit.mix)
marginalR2 = performance::r2(fit.mix, component = "fixed")


# impatiens
formula.imp = formula(num_colonies ~
                        idwSN_imp + idwHAY_imp + idwIND_imp + idwURB_imp +
                        idwANN_imp + idwBLU_imp + idwPER_imp + idwFAL_imp +
                        landscape_iji_imp_mean + year + offset(log(effort)) +
                        (1|sample_pt)
)

bf.imp <- bf(formula.imp, family="poisson")
fit.imp = brm(bf.imp, pertransecti,
              cores=4, chains = 4,
              iter = 4000,
              control = list(max_treedepth = 15),
              verbose = TRUE
)
p = plot(pp_check(fit.imp))
p + xlim(0, 50)
fullmodelR2 = bayes_R2(fit.imp)
marginalR2 = performance::r2(fit.imp, component = "fixed")


# Plot all together
models = list(
  mixtus= fit.mix,
  impatiens = fit.imp
)
params2 = c("b_idwSN_imp", "b_idwHAY_imp", "b_idwIND_imp", "b_idwURB_imp", "b_idwANN_imp", "b_idwBLU_imp", "b_idwPER_imp", "b_idwFAL_imp", "b_landscape_iji_imp_mean")
params1 = c("b_idwSN_mix", "b_idwHAY_mix", "b_idwIND_mix", "b_idwURB_mix", "b_idwANN_mix", "b_idwBLU_mix", "b_idwPER_mix", "b_idwFAL_mix", "b_landscape_iji_mix_mean")
names  = c("seminatural", "hay", "industrial", "suburban", "annual", "blueberry", "perennial", "fallow", "patch \n interspersion", "chain", "iteration", "draw", "model")

d1 = as_draws_df(fit.mix)
d1 = posterior::subset_draws(d1, variable = params1)
d1 = as.data.frame(d1) %>% 
  mutate(model = "B. mixtus")
colnames(d1) = names

d2 = as_draws_df(fit.imp)
d2 = posterior::subset_draws(d2, variable = params2)
d2 = as.data.frame(d2) %>% 
  mutate(model = "B. impatiens")
colnames(d2) = names

d = rbind(d1, d2)

df_long = d %>%
  pivot_longer(cols = 
                 all_of(c("seminatural", "hay", "industrial", "suburban", "annual", "blueberry", "perennial", "fallow", "patch \n interspersion")), 
               names_to = "parameter", values_to = "value")

# get which 95% CI overlap 0
df_summ = df_long %>%
  group_by(model, parameter) %>%
  median_qi(value, .width = 0.95) %>%
  ungroup() %>%
  mutate(less_than_zero = (.upper <= 0)) %>%
  mutate(greater_than_zero = (.lower >= 0))
df_summ$CredibleInterval = NA

for (i in 1:nrow(df_summ)){
  if (df_summ$less_than_zero[i] == TRUE){
    df_summ$CredibleInterval[i] = "less than zero"
  }else if(df_summ$greater_than_zero[i] == TRUE){
    df_summ$CredibleInterval[i] = "greater than zero"
  }else{
    df_summ$CredibleInterval[i] = "nonsignificant"
  }
}

df_long2 = df_long %>%
  left_join(df_summ[c("model", "parameter", "CredibleInterval")])
df_long2$CredibleInterval = factor(df_long2$CredibleInterval,
                                  levels = c("less than zero", "nonsignificant",
                                             "greater than zero"))
df_long2$parameter <- factor(df_long2$parameter,
                             levels = c("patch \n interspersion", "seminatural", "industrial",
                                        "suburban", "fallow", "hay", "perennial", 
                                        "blueberry", "annual"))
df_long2$model <- factor(df_long2$model,
                             levels = c("B. mixtus", "B. impatiens"))
df_long2 <- df_long2 %>%
  mutate(metric_group = ifelse(parameter == "patch \n interspersion",
                               "Configuration",
                               "Distance-weighted Composition"))


params_plot = ggplot(df_long2, aes(x = value, y = parameter, fill = CredibleInterval)) +
  scale_fill_manual(values = c(gold, "lightgrey", faded_green)) +
  stat_halfeye(point_interval = median_qi, .width = c(0.95)) +
  vline_0(linetype = "dashed") +
  #facet_wrap(~model, ncol = 1) +
  facet_grid(metric_group ~ model, scales = "free_y", space = "free_y") +
  xlab("Posterior Values") +
  ylab("Landcover Metric") +
  theme(legend.position = "none") +
  guides(colour = "none", fill = "none", alpha = "none", slab_fill = "none") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )
ggsave("figures/manuscript_figures/colony_density_compositional.jpg",
       params_plot, units = "px", height = 2500, width = 2000)


# plot correlation of raw predictors (mixtus)
landscape_corr <- cor(
  column_to_rownames(landscape_metrics[c("sample_pt", "landscape_iji_mix_mean", landscape_mix)], var = "sample_pt")
)
ggcorrplot(landscape_corr, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3, ggtheme = ggplot2::theme_minimal())

# plot correlation of posterior draws (mixtus)
posterior <- as_draws_df(fit.mix)
bayes_cor <- posterior[,1:10] %>%
  cor()
bayes_cor
ggcorrplot(bayes_cor, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3, ggtheme = ggplot2::theme_minimal())




############################################
# Plot conditional effects
###########################################

mixtus.conditional = conditional_effects(fit.mix)
impatiens.conditional = conditional_effects(fit.imp)


############################################
# mixtus ~ iji
###########################################
mixiji =
  mixtus.conditional[["landscape_iji_mix_mean"]]

# compute effort adjusted raw data
pertransectm_adjusted = pertransectm
pertransectm_adjusted$num_colonies = pertransectm$num_colonies / pertransectm$effort * mixiji$effort[1]

#ggplot
mix.iji = 
  
  #plot raw data
  ggplot(pertransectm_adjusted, aes(x = landscape_iji_mix_mean, y = num_colonies)) +
  geom_point(cex = 2, alpha = 0.4) +
  labs(x="Patch interspersion", y="*B. mixtus* colony abundance") +
  theme_bw() +
  theme(legend.position = "none") +
  #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    axis.title   = element_markdown(size = 11)
  ) +
  ylim(c(0,50)) +
  #plot model prediction with credible interval
  geom_line(data = mixiji, aes(x = landscape_iji_mix_mean, y=estimate__)) +
  geom_ribbon(data = mixiji, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = faded_green)

mix.iji


############################################
# mixtus ~ perennial
###########################################
mixper =
  mixtus.conditional[["idwPER_mix"]]

# compute effort adjusted raw data
pertransectm_adjusted = pertransectm
pertransectm_adjusted$num_colonies = pertransectm$num_colonies / pertransectm$effort * mixper$effort[1]

#ggplot
mix.per = 
  
  #plot raw data
  ggplot(pertransectm_adjusted, aes(x = idwPER_mix, y = num_colonies)) +
  geom_point(cex = 2, alpha = 0.2) +
  labs(x="Distance-weighted perennial crop area", y="") +
  theme_bw() +
  theme(legend.position = "none") +
  #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    axis.title   = element_text(size = 11)
  ) +
  ylim(c(0,50)) +
  #plot model prediction with credible interval
  geom_line(data = mixper, aes(x = idwPER_mix, y=estimate__)) +
  geom_ribbon(data = mixper, aes(ymin = lower__, ymax = upper__,
                                 alpha=0.5), fill = faded_green)

mix.per


############################################
# mixtus ~ industrial
###########################################
mixind =
  mixtus.conditional[["idwIND_mix"]]

# compute effort adjusted raw data
pertransectm_adjusted = pertransectm
pertransectm_adjusted$num_colonies = pertransectm$num_colonies / pertransectm$effort * mixind$effort[1]

#ggplot
mix.ind = 
  
  #plot raw data
  ggplot(pertransectm_adjusted, aes(x = idwIND_mix, y = num_colonies)) +
  geom_point(cex = 2, alpha = 0.2) +
  labs(x="Distance-weighted industrial area", y="") +
  theme_bw() +
  theme(legend.position = "none") +
  #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    axis.title   = element_text(size = 11)
  ) +
  ylim(c(0,50)) +
  #plot model prediction with credible interval
  geom_line(data = mixind, aes(x = idwIND_mix, y=estimate__)) +
  geom_ribbon(data = mixind, aes(ymin = lower__, ymax = upper__,
                                 alpha=0.5), fill = gold)

mix.ind



############################################
# impatiens ~ annual
###########################################
impann =
  impatiens.conditional[["idwANN_imp"]]
impann$upper__[impann$upper__ > 50] = 50 #truncate for plotting

# compute effort adjusted raw data
pertransecti_adjusted = pertransecti
pertransecti_adjusted$num_colonies = pertransecti$num_colonies / pertransecti$effort * impann$effort[1]

#ggplot
imp.ann = 
  
  #plot raw data
  ggplot(pertransecti_adjusted, aes(x = idwANN_imp, y = num_colonies)) +
  geom_point(cex = 2, alpha = 0.2) +
  labs(x="Distance-weighted annual crop area", y="*B. impatiens* colony abundance") +
  theme_bw() +
  theme(legend.position = "none") +
  #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    axis.title   = element_markdown(size = 11)
  ) +
  ylim(c(0,50)) +
  #plot model prediction with credible interval
  geom_line(data = impann, aes(x = idwANN_imp, y=estimate__)) +
  geom_ribbon(data = impann, aes(ymin = lower__, ymax = upper__,
                                 alpha=0.5), fill = faded_green)

imp.ann


############################################
# impatiens ~ urban
###########################################
impurb =
  impatiens.conditional[["idwURB_imp"]]

# compute effort adjusted raw data
pertransecti_adjusted = pertransecti
pertransecti_adjusted$num_colonies = pertransecti$num_colonies / pertransecti$effort * impurb$effort[1]

#ggplot
imp.urb = 
  
  #plot raw data
  ggplot(pertransecti_adjusted, aes(x = idwURB_imp, y = num_colonies)) +
  geom_point(cex = 2, alpha = 0.2) +
  labs(x="Distance-weighted suburban area", y="") +
  theme_bw() +
  theme(legend.position = "none") +
  #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    axis.title   = element_text(size = 11)
  ) +
  ylim(c(0,50)) +
  #plot model prediction with credible interval
  geom_line(data = impurb, aes(x = idwURB_imp, y=estimate__)) +
  geom_ribbon(data = impurb, aes(ymin = lower__, ymax = upper__,
                                 alpha=0.5), fill = faded_green)

imp.urb


############################################
# impatiens ~ industrial
###########################################
impind =
  impatiens.conditional[["idwIND_imp"]]

# compute effort adjusted raw data
pertransecti_adjusted = pertransecti
pertransecti_adjusted$num_colonies = pertransecti$num_colonies / pertransecti$effort * impind$effort[1]

#ggplot
imp.ind = 
  
  #plot raw data
  ggplot(pertransecti_adjusted, aes(x = idwIND_imp, y = num_colonies)) +
  geom_point(cex = 2, alpha = 0.2) +
  labs(x="Distance-weighted industrial area", y="") +
  theme_bw() +
  theme(legend.position = "none") +
  #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    axis.title   = element_text(size = 11)
  ) +
  ylim(c(0,50)) +
  #plot model prediction with credible interval
  geom_line(data = impind, aes(x = idwIND_imp, y=estimate__)) +
  geom_ribbon(data = impind, aes(ymin = lower__, ymax = upper__,
                                 alpha=0.5), fill = faded_green)

imp.ind



############################################
# Make grid plot
############################################

grid = grid.arrange(mix.iji, mix.per, mix.ind,
                    nullGrob(), nullGrob(), nullGrob(),
                    imp.ann, imp.urb, imp.ind,
                    ncol = 3, heights = c(10, 1, 10))
grid2 = grid.arrange(params_plot, nullGrob(), grid, ncol = 3, widths = c(10,1,15))
grid3 = ggdraw() +
  draw_plot(grid2, 0, 0, 1, 1) +
  draw_plot_label(c("A", "B"), 
                  x = c(0, 0.40), 
                  y = c(1, 1))
grid3
ggsave("figures/manuscript_figures/colony_grid.jpg",
       grid3, height = 2000, width = 4000, units = "px")
