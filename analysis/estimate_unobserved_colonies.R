####### Fit truncated Poisson to estimate number of unobserved colonies per site
### Nov 4 2025
### J Melanson
### based on tutorial: https://freerangestats.info/blog/2018/03/20/truncated-poisson

# Load packages
library(tidyverse)
library(scales)
library(fitdistrplus)
library(rstan)
library(truncdist)


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
### Now incorporate multiple sites using sim data
###################################################

simdata = read.csv("simulate_data/colony_assignments/sim_data/true_data/mixtus_set5_sub1.csv")


colonies_by_landscape = simdata %>%
  group_by(truecolony, landscape_id) %>%
  summarize(count = n())
total_counts = colonies_by_landscape %>%
  group_by(landscape_id) %>%
  summarize(total_colonies = n())

#select stan model to fit
stanfile = paste("models/trunc_poisson.stan")
# updated so more complex than what I ran above

# prep data for stan
data = list()
data$num_sites = 6
data$num_colonies = nrow(colonies_by_landscape)
data$lower_limit = 1
data$sib_counts = colonies_by_landscape$count
data$site_ids = colonies_by_landscape$landscape_id
data$lambda_mu = mean(colonies_by_landscape$count)
data$lambda_sigma = 1


# fit stan model!
truncfit = stan(file = stanfile,
                data = data, seed = 5838299,
                chains = 4, cores = 4,
                verbose = TRUE)


true_lambda = sum(colonies_by_landscape$count)/(6*1792)

# visualize credible interval
plot(truncfit) + 
  ggtitle("Credibility interval for lambda, estimated by Stan from truncated data",
          "(Correct value is 1.3)") +
  labs(y = "Estimated parameters") +
  theme_minimal()

# these overestimate a bit the lambdas


# Perhaps try again with negative binomial, to account for overdispersion due to differing distances from colony sources
#select stan model to fit
stanfile = paste("models/trunc_negbin.stan")

# prep data for stan
data = list()
data$num_sites = 6
data$num_colonies = nrow(colonies_by_landscape)
data$lower_limit = 1
data$sib_counts = colonies_by_landscape$count
data$site_ids = colonies_by_landscape$landscape_id
data$alpha_mu = mean(colonies_by_landscape$count)
data$alpha_sigma = 1


# fit stan model!
truncfit = stan(file = stanfile,
                data = data, seed = 5838299,
                chains = 4, cores = 4,
                verbose = TRUE)
