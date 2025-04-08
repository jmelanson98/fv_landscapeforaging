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

draw_N_bees = function(sample_size, # number of bees to sample
                       landscape_size, # vector of two values (xmax and ymax)
                       trapgrid_size, # will be centered within landscape
                       resource_range, # integer for variogram range parameter (e.g., 100)
                       number_traps, # integer
                       number_colonies, # a square number
                       colony_sizes, #vector of length number_colonies
                       beta, # set param value
                       theta # set param value
){
  ##### Define colony characteristics #####
  colonyid = 1:number_colonies
  numbees_start = sum(colony_sizes)
  colony_x = round(runif(number_colonies, 1, landscape_size[1])) #colony locations are integers
  colony_y = round(runif(number_colonies, 1, landscape_size[2])) #colony locations are integers
  colony_data = as.data.frame(cbind(colonyid, colony_sizes, colony_x, colony_y))
  
  ##### Define trap characteristics #####
  trapid = 1:number_traps
  grid_size = sqrt(number_traps)
  x_step = trapgrid_size[1]/(grid_size-1)
  y_step = trapgrid_size[2]/(grid_size-1)
  trap_x = (landscape_size[1] - trapgrid_size[1])/2 + x_step*(0:(grid_size-1))
  trap_y = (landscape_size[2] - trapgrid_size[2])/2 + y_step*(0:(grid_size-1))
  coords = expand.grid(trap_x = trap_x, trap_y = trap_y)
  trap_data = as.data.frame(cbind(trapid, coords)) 
  
  # optional plotting step to visualize traps and colonies
  plot(colony_data$colony_x, colony_data$colony_y, col = "black", xlab = "Long", ylab = "Lat", main = "Colonies (black) and Traps (red)")
  points(trap_data$trap_x, trap_data$trap_y, col = "red")
  
  
  ##### Define floral landscape #####
  # use Brownian variogram to simulate spatial distribution of resources (as in Pope & Jha)
  # this block of code (until "Simulate Data") modified from ChatGPT
  
  # create a regular grid of points to simulate over
  grid_spacing = 1
  x.range <- seq(1, landscape_size[1], by = grid_spacing)  # x coordinates
  y.range <- seq(1, landscape_size[2], by = grid_spacing)  # y coordinates
  
  grid <- expand.grid(x = x.range, y = y.range)
  coordinates(grid) <- ~x + y
  gridded(grid) <- TRUE
  
  # define a Brownian variogram model
  # "nugget" = 0, "sill" = total variance
  vgm_model <- vgm(psill = 1, model = "Lin", nugget = 0, resource_range)
  
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
  if (grid_spacing > 1){
    fq <- kronecker(coarse_mat, matrix(1, nrow = grid_spacing, ncol = grid_spacing))
  }
  
  
  ##### Start sampling #####
  yik = allocMatrix(nrow = number_colonies, ncol = number_traps, value = 0)
  N = numbees_start
  colony_data$ni = colony_sizes
  
  while(sum(yik) < sample_size){
    # set weights based on colony sizes
    colony_data$w_i = colony_data$ni/N
    
    # calculate Pr(s = k | s in kappa) via matrix convolution
    # code with help from chatgpt
    
    # compute visitation rates for each colony at each grid cell using nested functions,
    # compute_distances and compute_visitation_rates 
    visitation_rates <- compute_visitation_rates(colonies = colony_data, 
                                                 resource_quality = fq, 
                                                 beta = beta, 
                                                 theta = theta, 
                                                 landscape_size = landscape_size)
    
    
    # calculate D_i for each colony (total visitation over all cells)
    colony_data$D_i <- sapply(visitation_rates, function(x) sum(x))
    
    # compute probabilities of sampling from each trap k in kappa
    # e.g., Pr(s = k)
    # sum over C:
    #### lambda_ik * w_i / D_i
    lambda_ik = allocMatrix(nrow = number_colonies, ncol = number_traps, value = 0)
    
    for (i in colony_data$colonyid){
      for (k in trap_data$trapid){
        lambda_ik[i,k] = visitation_rates[[i]][trap_data$trap_x[trap_data$trapid == k],
                                               trap_data$trap_y[trap_data$trapid == k]]
      }
    }
    
    trap_data$prob_s_eq_k = colSums(lambda_ik*colony_data$w_i/colony_data$D_i)
    trap_data$prob_s_eq_k[is.na(trap_data$prob_s_eq_k)] = 0
    
    # calculate prob of sampling from any trap kappa (Pr(s in kappa))
    prob_kappa = sum(rowSums(lambda_ik)*colony_data$w_i/colony_data$D_i)
    
    # prob of sampling from a particular trap given k in kappa
    trap_data$trap_prob = trap_data$prob_s_eq_k/prob_kappa
    # sanity check: conditional probs sum to 1!! yay :)
    
    #draw a trap_id
    trap = sample(trap_data$trapid, size = 1, prob = trap_data$trap_prob)
    
    # next up: draw a value of c from Pr(c = i | s = k)
    colony_data$colony_prob = ((lambda_ik[,trap]/colony_data$D_i)*colony_data$w_i)/trap_data$prob_s_eq_k[trap_data$trapid == trap]
    # sanity check: conditional probs sum to 1!! yay :)
    
    colony = sample(colony_data$colonyid, size = 1, prob = colony_data$colony_prob)
    
    # record visitation event to yik
    yik[colony, trap] = yik[colony, trap] + 1
    
    #update ni and N
    N = N-1
    colony_data$ni[colony_data$colonyid == colony] = colony_data$ni[colony_data$colonyid == colony] - 1
    
    print(paste("Now sampling bee number", sum(yik)+1))
  }
  
  return(yik)
}



# Simulate 10 draws
obs = draw_N_bees(sample_size = 10,
                  landscape_size = c(700,700),
                  trapgrid_size = c(300,300),
                  resource_range = 10,
                  number_traps = 25,
                  number_colonies = 10,
                  colony_sizes = rep(20, 10),
                  beta = -1/50,
                  theta = 0.5)
                                        


# for plotting of a foraging kernel
# create a large raster object from the matrix
r <- raster(visitation_rates[[1]])
png("figures/simfigs/foragingkernel_single.png", width = 3500, height = 3500, res = 300)
levelplot(r, col.regions = viridis::viridis(100))
dev.off()




# next: make code iterative (repeat these steps until reaching a stopping point)
