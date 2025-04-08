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
# pretend that each map unit = 5 meters
colony_lower = 1
colony_upper = 700
trap_lower = 200
trap_upper = 500
num_col = 100
num_bees = 2000
num_traps = 30

##### Define parameters #####

#what is a reasonable value for beta??
x = 1:200
y = exp(-(1/60)*x)
plot(x, y, xlab = "distance from nest", ylab = "visitation rate")
#this is really small but gives a reasonable decay of foraging distance?
#might be better to fit 1/beta? not sure if this makes any difference....
inv_beta = 50

#what is a reasonable value for theta? (increased visitation to high floral quality)
#how many more bees would visit the *best* flower patch compared to an average patch?
y = exp(-1/inv_beta*x + (1)*-2)
plot(x, y, xlab = "distance from nest", ylab = "visitation rate")
####(realizing that this function form does not modulate floral attractiveness as a fxn
#### of distance...e.g., the ratio of bees on a good compared vs bees on a bad patch
#### stays the same regardless of distance from the nest)
#### Not sure if this is the best assumption but let's fly with it for now...
theta = 0.5


##### Define colony characteristics #####
col_id = 1:num_col
colsize = rep(num_bees/num_col, num_col)
col_xcoord = round(runif(num_col, colony_lower, colony_upper)) #colony locations are integers (e.g., round to nearest meter)
col_ycoord = round(runif(num_col, colony_lower, colony_upper)) #colony locations are integers (e.g., round to nearest meter)
w_i = colsize/num_bees

colony_data = as.data.frame(cbind(col_id, colsize, w_i, col_xcoord, col_ycoord))

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
grid_spacing = 1
x.range <- seq(colony_lower, colony_upper, by = grid_spacing)  # X coordinates
y.range <- seq(colony_lower, colony_upper, by = grid_spacing)  # Y coordinates

grid <- expand.grid(x = x.range, y = y.range)
coordinates(grid) <- ~x + y
gridded(grid) <- TRUE

# define a Brownian variogram model
# "nugget" = 0, "sill" = total variance
vgm_model <- vgm(psill = 10, model = "Lin", nugget = 0, range = 10) #range is not actually used in a Lin model but we need it as a place holder so gstat doesn't get mad

# simulate one realization of a Gaussian random field
simulated <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, 
                   beta = 0, model = vgm_model, nmax = 20)

simulated_field <- predict(simulated, newdata = grid, nsim = 1) #change nsim to make multiple landscapes

# plot the simulated resource distribution
spplot(simulated_field, zcol = "sim1", main = "Simulated Resource Distribution (Brownian variogram)")

# make a matrix of floral quality values
sim_values <- simulated_field$sim1
n_cells <- length(x.range)
fq <- matrix(sim_values, nrow = n_cells, ncol = n_cells)

# if you don't have sim values for every map unit, then expand (duplicate) each 
# coarse cell into 2x2 block in a fine grid using Kronecker product
fq <- kronecker(coarse_mat, matrix(1, nrow = grid_spacing, ncol = grid_spacing))

# normalize / z-score the floral quality values
fq_norm = fq/sd(fq)

##### Simulate data #####
# efficient calculation of Pr(s = k | s in Kappa) via convolution
# with help from ChatGPT

# compute visitation rates for each colony at each grid cell using nested functions,
# compute_distances and compute_visitation_rates 
visitation_rates <- compute_visitation_rates(colonies = colony_data, 
                                             resource_quality = fq_norm, 
                                             beta = -1/inv_beta, 
                                             theta = theta, 
                                             landscape_size = c(700, 700))

# quick visual check of one foraging kernel
# create a large raster object from the matrix
r <- raster(visitation_rates[[1]])

png("figures/simfigs/foragingkernel_single.png", width = 3500, height = 3500, res = 300)
levelplot(r, col.regions = viridis::viridis(100))
dev.off()

# calculate D_i for each colony (total visitation over all cells)
colony_data$D_i <- sapply(visitation_rates, function(x) sum(x))

# compute probabilities of sampling from each trap k in kappa
# e.g., Pr(s = k)
# sum over C:
#### lambda_ik * w_i / D_i
lambda_ik = allocMatrix(nrow = num_col, ncol = num_traps, value = 0)

for (i in colony_data$col_id){
  for (k in trap_data$trap_id){
    lambda_ik[i,k] = visitation_rates[[i]][trap_data$trap_xcoord[trap_data$trap_id == k],
                                           trap_data$trap_ycoord[trap_data$trap_id == k]]
  }
}

trap_data$prob_s_eq_k = colSums(lambda_ik*colony_data$w_i/colony_data$D_i)

# calculate prob of sampling from any trap kappa
prob_kappa = sum(rowSums(lambda_ik)*colony_data$w_i/D_i)

# prob of sampling from a particular trap given k in kappa
trap_data$conditional_prob = trap_data$prob_s_eq_k/prob_kappa
# sanity check: conditional probs sum to 1!! yay :)

#draw a trap_id
trap_id = sample(trap_data$trap_id, size = 1, prob = trap_data$conditional_prob)

# next up: draw a value of c from Pr(c = i | s = k)
colony_data$colony_prob = ((lambda_ik[,trap_id]/colony_data$D_i)*colony_data$w_i)/trap_data$prob_s_eq_k[trap_data$trap_id == trap_id]

colony = sample(colony_data$col_id, size = 1, prob = colony_data$colony_prob)

#update ni and N
N = N-1
colony_data$colsize[colony_data$col_id == colony] = colony_data$colsize[colony_data$col_id == colony] - 1


# next: make code iterative (repeat these steps until reaching a stopping point)
