##### Simulate data for pope_consgenetics model #####
# Script Initiated: April 2, 2025
# By: Jenna Melanson
# Goal: reproduce simulated data from 2017 Conservation Genetics Paper


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
landscape_lower = 0
landscape_upper = 500
num_col = 100
num_bees = 2000
num_traps = 30

##### Define parameters #####
#what is a reasonable value for beta??
x = 1:500
y = exp(-0.005*x)
plot(x, y)
#this is really small but gives a reasonab
inv_beta = 200

##### Define colony characteristics #####
col_id = 1:num_col
colsize = rep(num_bees/num_col, num_col)
col_xcoord = runif(num_col, landscape_lower, landscape_upper)
col_ycoord = runif(num_col, landscape_lower, landscape_upper)

colony_data = as.data.frame(cbind(col_id, colsize, col_xcoord, col_ycoord))

#plot colonies
landscape = plot(colony_data$col_xcoord, colony_data$col_ycoord)

##### Define trap characteristics #####
trap_id = 1:num_traps
trap_quality = rnorm(num_traps, 0, 1)
trap_xcoord = runif(num_traps, landscape_lower, landscape_upper)
trap_ycoord = runif(num_traps, landscape_lower, landscape_upper)

trap_data = as.data.frame(cbind(trap_id, trap_quality, trap_xcoord, trap_ycoord))

#plot traps
landscape = landscape + plot(trap_data$trap_xcoord, trap_data$trap_ycoord, col = "red")
