##### Some useful Stan diagnostic code #####
# Script Initiated: May 11, 2025
# Started by: Jenna Melanson
# This code is modified from various tidbits from Michael Betancourt (https://betanalpha.github.io)

#source
source("simulate_data/mcmc_analysis_tools_rstan.R")

# extract expectands
extract_expectands <- function(stan_fit) {
  nom_params <- rstan:::extract(stan_fit, permuted=FALSE)
  N <- dim(nom_params)[3] - 1
  params <- lapply(1:N, function(n) t(nom_params[,,n]))
  names(params) <- names(stan_fit)[1:N]
  (params)
}

extract_hmc_diagnostics <- function(stan_fit) {
  diagnostic_names <- c('divergent__', 'treedepth__', 'n_leapfrog__', 
                        'stepsize__', 'energy__', 'accept_stat__')
  
  nom_params <- get_sampler_params(stan_fit, inc_warmup=FALSE)
  C <- length(nom_params)
  params <- lapply(diagnostic_names, 
                   function(name) t(sapply(1:C, function(c) 
                     nom_params[c][[1]][,name])))
  names(params) <- diagnostic_names
  (params)
}

# check distribution of divergent transitions in parameter space
x = c("mu", "beta", "theta")
y = c("mu", "beta", "theta")
samples = extract_expectands(stanFitGQ)
diagnostics = extract_hmc_diagnostics(stanFitGQ)
plot_div_pairs(x, y, samples, diagnostics)

