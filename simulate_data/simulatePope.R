##### Simulate data for pope_consgenetics model #####
# Script Initiated: April 2, 2025
# By: Jenna Melanson
# Goal: reproduce simulated data from 2017 Conservation Genetics Paper

##### Load packages #####
library(matrixStats)
library(sp)
library(gstat)
library(ggplot2)
library(reshape2)
library(raster)
library(rasterVis)

##### Load source code #####
source("simulate_data/sim_src.R")

##### Required Data #####
### number of colonies (num_col)
### total number of bees (num_bees)
### number of bees per colony (col_size)
### landscape bounds (e.g., define survey space)
### colony locations (start with uniform distribution)

### number of traps
### trap locations
### floral quality at each trap (normal(0,1) for now)

### rate parameter beta (how willing are bees to travel long distances)
### rate parameter theta (how attracted are bees to flowers)

### eps (trap specific random effect) --> ignore initially
### psi (colony-specific random effect) --> ignore initially


##### Define landscape bounds, number of bees & traps #####
colony_lower = 1
colony_upper = 3500
trap_lower = 1000
trap_upper = 2500
num_col = 10
num_bees = 200
num_traps = 30

##### Define parameters #####

#what is a reasonable value for beta??
x = 1:1000
y = exp(-(1/300)*x)
plot(x, y, xlab = "distance from nest", ylab = "visitation rate")
#this is really small but gives a reasonable decay of foraging distance?
#might be better to fit 1/beta? not sure if this makes any difference....
inv_beta = 300

#what is a reasonable value for theta? (increased visitation to high floral quality)
#how many more bees would visit the *best* flower patch compared to an average patch?
y = exp(-1/inv_beta*x + (1)*-2)
plot(x, y, xlab = "distance from nest", ylab = "visitation rate")
####(realizing that this function form does not modulate floral attractiveness as a fxn
#### of distance...e.g., the ratio of bees on a good compared vs bees on a bad patch
#### stays the same regardless of distance from the nest)
#### Not sure if this is the best assumption but let's fly with it for now...
theta = 1


##### Define colony characteristics #####
col_id = 1:num_col
colsize = rep(num_bees/num_col, num_col)
col_xcoord = round(runif(num_col, colony_lower, colony_upper)) #colony locations are integers (e.g., round to nearest meter)
col_ycoord = round(runif(num_col, colony_lower, colony_upper)) #colony locations are integers (e.g., round to nearest meter)

colony_data = as.data.frame(cbind(col_id, colsize, col_xcoord, col_ycoord))

#plot colonies
landscape = plot(colony_data$col_xcoord, colony_data$col_ycoord)

##### Define trap locations #####
trap_id = 1:num_traps
trap_xcoord = round(runif(num_traps, trap_lower, trap_upper))
trap_ycoord = round(runif(num_traps, trap_lower, trap_upper))
trap_data = as.data.frame(cbind(trap_id, trap_xcoord, trap_ycoord))

#plot traps
landscape = landscape + plot(trap_data$trap_xcoord, trap_data$trap_ycoord, col = "red")

##### Define floral landscape #####
# use Brownian variogram to simulate spatial distribution of resources (as in Pope & Jha)
# this block of code (until "Simulate Data") modified from ChatGPT suggestion

# create a regular grid of points to simulate over
grid_spacing = 2
x.range <- seq(colony_lower, colony_upper, by = grid_spacing)  # X coordinates
y.range <- seq(colony_lower, colony_upper, by = grid_spacing)  # Y coordinates

grid <- expand.grid(x = x.range, y = y.range)
coordinates(grid) <- ~x + y
gridded(grid) <- TRUE

# define a Brownian variogram model
# "nugget" = 0, "sill" = total variance
vgm_model <- vgm(psill = 100, model = "Lin", nugget = 0, range = 1) #range is not actually used in a Lin model but we need it as a place holder so gstat doesn't get mad

# simulate one realization of a Gaussian random field
simulated <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, 
                   beta = 0, model = vgm_model, nmax = 20)

simulated_field <- predict(simulated, newdata = grid, nsim = 1) #change nsim to make multiple landscapes
simulated_field2 <- predict(simulated, newdata = grid, nsim = 1) #change nsim to make multiple landscapes
# this takes a while for a landscape of my size ... 

# plot the simulated resource distribution
spplot(simulated_field, zcol = "sim1", main = "Simulated Resource Distribution (Brownian variogram)")

# make a matrix of floral quality values
sim_values <- simulated_field$sim1
n_cells <- length(x.range)
coarse_mat <- matrix(sim_values, nrow = n_cells, ncol = n_cells)

# if you don't have sim values for every map unit, then expand (duplicate) each 
# coarse cell into 2x2 block in a fine grid using Kronecker product
fine_mat <- kronecker(coarse_mat, matrix(1, nrow = grid_spacing, ncol = grid_spacing))

##### Simulate data #####
# uncorrelated floral quality for now
fq = matrix(data = rnorm(colony_upper^2, 0, 1), nrow = colony_upper, ncol = colony_upper)

# efficient calculation of Pr(s = k | s in Kappa) via convolution
# with help from ChatGPT

# compute visitation rates for each colony at each grid cell using nested functions,
# compute_distances and compute_visitation_rates 
visitation_rates <- compute_visitation_rates(colonies = colony_data, 
                                             resource_quality = fq, 
                                             beta = -1/inv_beta, 
                                             theta = theta, 
                                             landscape_size = c(3500, 3500))

# quick visual check of one foraging kernel
# create a large raster object from the matrix
r <- raster(visitation_rates[[1]])

png("figures/simfigs/foragingkernel_single.png", width = 3500, height = 3500, res = 300)
levelplot(r, col.regions = viridis::viridis(100))
dev.off()

ggsave("figures/simfigs/foragingkernel_single.jpg", 
       plot = forasinglecolony, 
       width = 7000, height = 7000,
       units = "px",
       dpi = 600)


# calculate D_i for each colony (total visitation over all cells)
D_i <- sapply(visitation_rates, function(x) sum(x))

# compute probabilities of sampling for each trap k in kappa
colony_data$w_i = colony_data$colsize/num_bees


# IN PROGRESS
for (k in trap_data$trap_id){
  #Pr(s = k)
  numerator <- sum(w_i * sapply(visitation_rates, function(x) x[trap_data$trap_xcoord[trap_data$trap_id ==k],
                                                                trap_data$trap_ycoord[trap_data$trap_id == k]]))
  trap_data$
}

# Numerator for a specific trap k
numerator <- sum(w_i * sapply(visitation_rates, function(x) x[trap_data$trap_xcoord[trap_data$trap_id ==k],
                                                              trap_data$trap_ycoord[trap_data$trap_id == k]]))

# Denominator for all traps in kappa (we assume kappa = 1:3500 here for simplicity)
kappa <- 1:(landscape_size[1] * landscape_size[2])
denominator <- sum(w_i * sapply(visitation_rates, function(x) sum(x[kappa])))

# Final probability of sampling at k
probability_k <- numerator / denominator

print(probability_k)















