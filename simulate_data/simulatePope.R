##### Simulate data for pope_consgenetics model #####
# Script Initiated: April 2, 2025
# By: Jenna Melanson
# Goal: reproduce simulated data from 2017 Conservation Genetics Paper

##### Load packages #####
library(matrixStats)

##### Load source code #####
source("/simulation/sim_src")

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
col_xcoord = runif(num_col, colony_lower, colony_upper)
col_ycoord = runif(num_col, colony_lower, colony_upper)

colony_data = as.data.frame(cbind(col_id, colsize, col_xcoord, col_ycoord))

#plot colonies
landscape = plot(colony_data$col_xcoord, colony_data$col_ycoord)

##### Define trap characteristics #####
trap_id = 1:num_traps
trap_quality = rnorm(num_traps, 0, 1)
trap_xcoord = runif(num_traps, trap_lower, trap_upper)
trap_ycoord = runif(num_traps, trap_lower, trap_upper)

trap_data = as.data.frame(cbind(trap_id, trap_quality, trap_xcoord, trap_ycoord))

#plot traps
landscape = landscape + plot(trap_data$trap_xcoord, trap_data$trap_ycoord, col = "red")


##### Simulate data #####
distance_mat = allocMatrix(nrow = num_col, ncol = num_traps, value = 0)
lambda = allocMatrix(nrow = num_col, ncol = num_traps, value = 0)

for (k in trap_data$trap_id){
  for (i in colony_data$col_id){
    distance_mat = dist
    lambda = 
      
      dis[i, k] <- distance(delta[i], trap[k]); #is distance a function in stan??
    lambda[i, k] <- -beta*dis[i, k] + theta*floral[k] + mu + zeta_scale[i] + eps_scale[k];
    
  }
}
lambda = exp(-1/inv_beta)

