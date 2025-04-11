##### Pope Simulation Functions #####
##### Simulate data for pope_consgenetics model #####
# Script Initiated: April 2, 2025
# By: Jenna Melanson
# Goal: reproduce simulated data from 2017 Conservation Genetics Paper

compute_visitation_rates_and_avg_distance <- function(colonies, resource_quality, beta, theta, landscape_size) {
  x_coords <- 1:landscape_size
  y_coords <- 1:landscape_size
  
  visitation_rates <- list()
  avg_distances <- numeric(nrow(colonies))
  total_vis <- numeric(nrow(colonies))
  
  for (i in 1:nrow(colonies)) {
    colony_x <- colonies$colony_x[i]
    colony_y <- colonies$colony_y[i]
    
    # Efficiently generate coordinate grids
    dist_x <- outer(x_coords, rep(colony_x, landscape_size), "-")
    dist_y <- outer(rep(colony_y, landscape_size), y_coords, "-")
    dist <- sqrt(dist_x^2 + dist_y^2)
    
    # Compute visitation matrix
    visitation <- exp(beta * dist + theta * resource_quality)
    
    # Compute weighted average distance
    total_visitation <- sum(visitation)
    avg_distance <- if (total_visitation > 0) {
      sum(visitation * dist) / total_visitation
    } else {
      NA_real_
    }
    
    visitation_rates[[i]] <- visitation
    avg_distances[i] <- avg_distance
    total_vis[i] = total_visitation
  }
  
  return(list(visitation_rates = visitation_rates,
              avg_distances = avg_distances,
              total_vis = total_vis))
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
  trap_data$fq <- mapply(function(x, y) resource_landscape[x, y], trap_data$trap_x, trap_data$trap_y)
 
  
  # optional plotting step to visualize traps and colonies
  plot(colony_data$colony_x, colony_data$colony_y, col = "black", xlab = "Long", ylab = "Lat", main = "Colonies (black) and Traps (red)")
  points(trap_data$trap_x, trap_data$trap_y, col = "red")
  
  
  ##### Create landscape-wide visitation matrix #####
  # compute visitation rates for each colony at each grid cell, total visitation, and average foraging range
  visitation_and_foraging_range <- compute_visitation_rates_and_avg_distance(colonies = colony_data, 
                                               resource_quality = resource_landscape, 
                                               beta = beta, 
                                               theta = theta, 
                                               landscape_size = landscape_size)
  
  visitation_rates = visitation_and_foraging_range[[1]]
  colony_data$foraging_range = visitation_and_foraging_range[[2]]
  colony_data$D_i = visitation_and_foraging_range[[3]]

  ##### Compute fixed values lambda_ik #####
  # calculate lambda_ik, e.g., visitation rates of each colony to specific traps
  lambda_ik = allocMatrix(nrow = number_colonies, ncol = number_traps, value = 0)
  
  for (i in seq_len(number_colonies)) {
    visit_mat <- visitation_rates[[i]]  # each is a matrix of size [xmax x ymax]
    
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
  
  output = list(yik, colony_data, trap_data)
  return(output)
}