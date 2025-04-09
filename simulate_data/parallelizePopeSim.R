##### Run parallel simulations of Pope model on server #####
# Script Initiated: April 8, 2025
# By: Jenna Melanson
# Goal: self contained script for running multiple iterations of Pope simulation on AllianceCan server

##### Load packages #####
library(matrixStats)
library(sp)
library(gstat)
library(ggplot2)
library(reshape2)
library(raster)
library(rasterVis)
library(parallel)
library(future)
library(furrr)
library(gstat)

##### Required Functions #####

##### helper function 1 #####
compute_distances <- function(colonies, landscape_size) {
  # function to compute distance matrix for each colony to each grid cell
  # via ChatGPT
  dist_matrix <- list()
  
  for (i in 1:nrow(colonies)) {
    colony_x <- colonies$colony_x[i]
    colony_y <- colonies$colony_y[i]
    
    # Create a grid of coordinates for the landscape
    x_coords <- 1:landscape_size
    y_coords <- 1:landscape_size
    
    # Compute the distance matrix for this colony (distance from colony to all grid cells)
    dist_x <- outer(x_coords, rep(colony_x, landscape_size), "-")  # Distance in x-direction
    dist_y <- outer(rep(colony_y, landscape_size), y_coords, "-")  # Distance in y-direction
    
    # Compute Euclidean distance
    dist_matrix[[i]] <- sqrt(dist_x^2 + dist_y^2)
  }
  
  return(dist_matrix)
}


##### helper function 2 #####
compute_visitation_rates <- function(colonies, resource_quality, beta, theta, landscape_size) {
  # compute visitation rates for each colony
  # via ChatGPT
  dist_matrix <- compute_distances(colonies, landscape_size)
  
  visitation_rates <- list()
  
  for (i in 1:nrow(colonies)) {
    dist <- dist_matrix[[i]]
    resource <- resource_quality
    
    # Apply the kernel function to compute visitation rates
    visitation_rates[[i]] <- exp(beta * dist + theta * resource)
  }
  
  return(visitation_rates)
}


##### Simulate floral quality landscapes #####
simulateLandscape = function(landscape_size, # integer, same for x and y
                             resource_range,
                             grid_spacing = 1) {
  # use Brownian variogram to simulate spatial distribution of resources (as in Pope & Jha)
  # this block of code (until "Simulate Data") modified from ChatGPT
  
  # create a regular grid of points to simulate over
  grid_spacing = 1
  x.range <- seq(1, landscape_size, by = grid_spacing)  # x coordinates
  y.range <- seq(1, landscape_size, by = grid_spacing)  # y coordinates
  
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
  fq <- matrix(simulated_field$sim1, nrow = length(x.range), ncol = length(y.range))
  
  # if you don't have sim values for every map unit, then expand (duplicate) each 
  # coarse cell into 2x2 block in a fine grid using Kronecker product
  if (grid_spacing > 1){
    fq <- kronecker(coarse_mat, matrix(1, nrow = grid_spacing, ncol = grid_spacing))
  }
  return(fq)
}


##### Simulate bee draws #####
draw_N_bees = function(sample_size, # number of bees to sample
                       landscape_size, # integer, size of full resource landscape
                       colonygrid_size, # integer, size of colony distribution
                       trapgrid_size, # integer, size of trap grid
                       resource_landscape, # matrix containing floral quality across landscape
                       number_traps, # integer
                       number_colonies, # a square number
                       colony_sizes, #vector of length number_colonies
                       beta, # set param value
                       theta, # set param value
                       batch_size # number of bees to draw before updating ni/N
){
  ##### Define colony characteristics #####
  colonyid = 1:number_colonies
  numbees_start = sum(colony_sizes)
  colony_x = round(runif(number_colonies, (landscape_size-colonygrid_size)/2, (landscape_size-colonygrid_size)/2 + colonygrid_size)) #colony locations are integers
  colony_y = round(runif(number_colonies, (landscape_size-colonygrid_size)/2, (landscape_size-colonygrid_size)/2 + colonygrid_size)) #colony locations are integers
  colony_data = as.data.frame(cbind(colonyid, colony_sizes, colony_x, colony_y))
  
  ##### Define trap characteristics #####
  trapid = 1:number_traps
  grid_size = sqrt(number_traps)
  step = trapgrid_size/(grid_size-1)
  trap_x = (landscape_size - trapgrid_size)/2 + step*(0:(grid_size-1))
  trap_y = (landscape_size - trapgrid_size)/2 + step*(0:(grid_size-1))
  coords = expand.grid(trap_x = trap_x, trap_y = trap_y)
  trap_data = as.data.frame(cbind(trapid, coords)) 
  
  # optional plotting step to visualize traps and colonies
  plot(colony_data$colony_x, colony_data$colony_y, col = "black", xlab = "Long", ylab = "Lat", main = "Colonies (black) and Traps (red)")
  points(trap_data$trap_x, trap_data$trap_y, col = "red")
  
  
  
  ##### Create landscape-wide visitation matrix #####
  # compute visitation rates for each colony at each grid cell using nested functions,
  # compute_distances and compute_visitation_rates 
  visitation_rates <- compute_visitation_rates(colonies = colony_data, 
                                               resource_quality = resource_landscape, 
                                               beta = beta, 
                                               theta = theta, 
                                               landscape_size = landscape_size)
  
  
  
  ##### Compute fixed values lambda_ik and D_i #####
  # calculate lambda_ik, e.g., visitation rates of each colony to specific traps
  lambda_ik = allocMatrix(nrow = number_colonies, ncol = number_traps, value = 0)
  colony_data$D_i = rep(0, number_colonies)
  
  for (i in seq_len(number_colonies)) {
    visit_mat <- visitation_rates[[i]]  # each is a matrix of size [xmax x ymax]
    
    # calculate D_i for each colony (total visitation over all cells)
    colony_data$D_i[i] <- sum(visit_mat)
    
    # pull out all trap visitation values in one vectorized step
    lambda_ik[i, ] <- mapply(function(x, y) visit_mat[x, y], trap_data$trap_x, trap_data$trap_y)
  }
  

  ##### Start sampling #####
  yik = allocMatrix(nrow = number_colonies, ncol = number_traps, value = 0)
  N = numbees_start
  colony_data$ni = colony_sizes
  
  while(sum(yik) < sample_size){
    remaining = sample_size - sum(yik)
    n_draws = min(batch_size, remaining)
    
    # set weights based on colony sizes
    colony_data$w_i = colony_data$ni/N
    
    # calculate Pr(s = k | s in kappa)
    
    #  first compute Pr(s = k)
    # to do this, sum over C:
    #### lambda_ik * w_i / D_i
    lambda_ik_scaled <- sweep(lambda_ik, 1, colony_data$w_i / colony_data$D_i, "*")  # elementwise scale by w_i / D_i
    trap_data$prob_s_eq_k <- colSums(lambda_ik_scaled)
    
    # calculate prob of sampling from any trap kappa (Pr(s in kappa))
    prob_kappa = sum(trap_data$prob_s_eq_k)
    
    # prob of sampling from a particular trap given k in kappa
    trap_data$trap_prob = trap_data$prob_s_eq_k/prob_kappa
    
    # sample traps
    traps = sample(trap_data$trapid, size = n_draws, replace = TRUE, prob = trap_data$trap_prob)
    
    # sample bees
    for (draw in seq_len(n_draws)) {
      trap = traps[draw]
      
      # calculate Pr(c = i | s = k)
      numer = (lambda_ik[, trap] / colony_data$D_i) * colony_data$w_i
      denom = trap_data$prob_s_eq_k[trap_data$trapid == trap]
      colony_probs <- numer / denom
      
      colony = sample(colony_data$colonyid, size = 1, prob = colony_probs)
      
      # record visitation event to yik
      yik[colony, trap] = yik[colony, trap] + 1
      
      #update ni and N
      N = N-1
      colony_data$ni[colony_data$colonyid == colony] = colony_data$ni[colony_data$colonyid == colony] - 1
    }
    print(paste("Now sampled", sum(yik), "of", sample_size, "bees"))
  }
  
  return(yik)
}



##### Prepare to run in parallel ! #####
# set parameter combinations
landscape_ids  = 1:10
betas <- c(-1/10, -1/25, -1/50, -1/75, -1/100)
sample_sizes <- c(250, 500, 1000)

param_grid <- expand.grid(
  landscape_id = landscape_ids,
  beta = betas,
  sample_size = sample_sizes,
  stringsAsFactors = FALSE
)

##### Set up run #####
# Get task ID from SLURM environment variable
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
params <- param_grid[task_id, ]
saveRDS(param_grid, "param_grid.rds")

# Simulate landscape based on landscape_id
fq <- simulateLandscape(landscape_size = 1100, resource_range = 10)

# Run simulation
result <- draw_N_bees(
  sample_size     = 100,
  landscape_size  = 1100,
  colonygrid_size = 700,
  trapgrid_size   = 300,
  number_traps    = 25,
  number_colonies = 100,
  colony_sizes    = rep(100, 100),
  beta            = -1/50,
  theta           = 0.5,
  resource_landscape = fq,
  batch_size = 1
)



result <- draw_N_bees(
  sample_size     = params$sample_size,
  landscape_size  = c(700, 700),
  trapgrid_size   = c(300, 300),
  resource_range  = 10,
  number_traps    = 25,
  number_colonies = 1000,
  colony_sizes    = rep(100, 1000),
  beta            = params$beta,
  theta           = 0.5,
  fq              = fq
)

# Save output
outfile <- sprintf("data/popesim2017/sim_result_%03d.rds", task_id)
saveRDS(result, outfile)


# for later
# files <- list.files("results", pattern = "sim_result_.*rds", full.names = TRUE)
# results <- lapply(files, readRDS)
# param_grid <- readRDS("param_grid.rds")
# 
# full_results <- bind_cols(param_grid, tibble(result = results))
