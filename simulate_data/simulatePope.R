##### Simulate data for pope_consgenetics model #####
# Script Initiated: April 2, 2025
# By: Jenna Melanson
# Goal: reproduce simulated data from 2017 Conservation Genetics Paper

##### Load packages #####
library(matrixStats)
library(sp)
library(gstat)

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
colony_lower = 0
colony_upper = 3500
trap_lower = 1000
trap_upper = 2500
num_col = 1000
num_bees = 20000
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
trap_xcoord = runif(num_traps, trap_lower, trap_upper)
trap_ycoord = runif(num_traps, trap_lower, trap_upper)
trap_data = as.data.frame(cbind(trap_id, trap_xcoord, trap_ycoord))

#plot traps
landscape = landscape + plot(trap_data$trap_xcoord, trap_data$trap_ycoord, col = "red")

##### Define floral landscape #####
# use Brownian variogram to simulate spatial distribution of resources (as in Pope & Jha)
# this block of code (until "Simulate Data") modified from ChatGPT suggestion

# create a regular grid of points to simulate over
x.range <- seq(colony_lower, colony_upper, by = 2)  # X coordinates
y.range <- seq(colony_lower, colony_upper, by = 2)  # Y coordinates

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
View(simulated_field)


##### Simulate data #####

# create lambda_full -> a landscape-wide visitation table (visitation rate between each pair of pixels based on distance and floral quality)

# start = starting location (xs, ys)
# end = ending location (xe, ye, fe) (fqe = floral quality)

lambda_full = allocMatrix(nrow = length(colony_lower:colony_upper), ncol = length(colony_lower:colony_upper), value = 0)



distance_mat = allocMatrix(nrow = num_col, ncol = num_traps, value = 0)
lambda = allocMatrix(nrow = num_col, ncol = num_traps, value = 0)
ypred = allocMatrix(nrow = num_col, ncol = num_traps, value = 0)

for (k in trap_data$trap_id){
  for (i in colony_data$col_id){
    distance_mat[i, k] = calculate_dist(colony_coords = list(colony_data$col_xcoord[i], colony_data$col_ycoord[i]),
                                  trap_coords = list(trap_data$trap_xcoord[k], trap_data$trap_ycoord[k]))
    lambda[i, k] = exp(-1/inv_beta*distance_mat[i, k] + theta*trap_data$trap_quality[k])
    
    #draw trap number using eqn 3
    
    ypred[i, k] = rpois(1, lambda[i, k])
  }
}

