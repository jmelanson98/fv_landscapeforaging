##### Generalized Simulation Functions #####
# Script Initiated: May 26, 2025
# By: Jenna Melanson
# Goal: simulate data using either distance decay function (specify)
# using raster (not matrix) implementation


##### Simulate floral quality landscapes as raster #####
simulateLandscapeRaster = function(landscape_size, # integer, same for x and y
                                   resource_range) {
  # use Brownian variogram to simulate spatial distribution of resources (as in Pope & Jha
  # create a regular grid of points to simulate over
  grid = rast(nrows = landscape_size, ncols = landscape_size,
              xmin = 0, xmax = landscape_size, 
              ymin = 0, ymax = landscape_size)
  
  # get xy coordinates of center of each pixel
  xy <- as.data.frame(xyFromCell(grid, 1:ncell(grid)))
  
  # create a spatial points data frame for gstat
  sp_points <- SpatialPointsDataFrame(coords = xy, data = data.frame(id = 1:nrow(xy)),
                                      proj4string = CRS(NA_character_))
  
  # define a Brownian variogram model
  # "nugget" = 0, "sill" = total variance
  vgm_model <- vgm(psill = 1, model = "Lin", nugget = 0, resource_range)
  
  # simulate one realization of a Gaussian random field
  simulated <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, 
                     beta = 0, model = vgm_model, nmax = 20)
  simulated_field <- predict(simulated, newdata = sp_points, nsim = 1)
  
  # return as a raster
  sim_df <- as.data.frame(simulated_field)
  sim_rast <- grid
  values(sim_rast) <- sim_df$sim1
  
  # plot
  plot(sim_rast, main = "Simulated Resource Distribution (Brownian variogram)")
  
  return(sim_rast)
}



##### compute visitation rates of each colony to each raster pixel
compute_visitation_on_raster <- function(colonies, 
                                         resource_quality_rast, 
                                         rho, 
                                         theta,
                                         distance_decay) {
  
  # get coordinates for center of each pixel
  xy_coords <- crds(resource_quality_rast, df = TRUE)  # matrix of x and y
  
  # flatten raster values into a vector
  rq_vals <- values(resource_quality_rast)
  
  lambda_ik = allocMatrix(nrow = number_colonies, ncol = number_traps, value = 0)
  trapcoords = as.matrix(trap_data[, c("trap_x", "trap_y")])
  avg_distances <- numeric(nrow(colonies))
  total_vis <- numeric(nrow(colonies))
  
  for (i in 1:nrow(colonies)) {
    colony_x <- colonies$colony_x[i]
    colony_y <- colonies$colony_y[i]
    
    # compute Euclidean distance from colony to each pixel
    dx <- xy_coords[, 1] - colony_x
    dy <- xy_coords[, 2] - colony_y
    dist <- sqrt(dx^2 + dy^2)
    
    print(rho)
    print(distance_decay)
    # compute visitation rate of colony at each pixel
    if (distance_decay == "exponentiated_quadratic"){
      visitation_vals <- exp(-0.5 * (dist / rho)^2 + theta * rq_vals)
    } else if (distance_decay == "exponential"){
      visitation_vals <- exp(dist / rho + theta * rq_vals)
    } else {
      print("Sorry, not a valid decay function.")
    }
    
    total_visitation <- sum(visitation_vals, na.rm = TRUE)
    avg_distance <- if (total_visitation > 0) {
      sum(visitation_vals * dist, na.rm = TRUE) / total_visitation
    } else {
      NA_real_
    }
    
    visit_rast <- resource_quality_rast  # clone geometry
    values(visit_rast) <- visitation_vals
    
    
    # calculate lambda_ik, e.g., visitation rates of each colony to specific traps
    # pull out all trap visitation values for that colony
    trap_vals <- terra::extract(visit_rast, trapcoords)
    lambda_ik[i, ] <- vals[, 1]
    
    # record average foraging distance and total visitation
    avg_distances[i] <- avg_distance
    total_vis[i] <- total_visitation
    
    # report
    print(paste("Computed visitation for colony", i, sep = " "))
  }
  
  return(list(
    lambda_ik = lambda_ik,
    avg_distances = avg_distances,
    total_vis = total_vis
  ))
}


##### Simulate bee draws #####
draw_bees_colony_restricted = function(sample_size, # number of bees to sample
                                       landscape_size, # integer, size of full resource landscape
                                       colonygrid_size, # integer, size of colony distribution
                                       trapgrid_size, # integer, size of trap grid
                                       resource_landscape, # raster containing floral quality across landscape
                                       nesting_landscape = NULL, # matrix containing nest habitat quality across landscape
                                       number_traps, # square integer
                                       number_colonies, # positive integer
                                       colony_sizes, #vector cd ..of length number_colonies
                                       rho, # set param value
                                       theta, # set param value
                                       distance_decay # current options: "exponentiated_quadratic" and "exponential"
){
  ##### Define colony characteristics #####
  colonyid = 1:number_colonies
  numbees_start = sum(colony_sizes)
  
  colony_x = c()
  colony_y = c()
  
  ##### Simulate colony locations #####
  while (length(colony_x) < number_colonies){
    # propose a colony location
    xprop = runif(1, (landscape_size-colonygrid_size)/2, (landscape_size-colonygrid_size)/2 + colonygrid_size)
    yprop = runif(1, (landscape_size-colonygrid_size)/2, (landscape_size-colonygrid_size)/2 + colonygrid_size)
    
    # accept or reject colony location
    if (is.null(nesting_landscape)){
      nesting_val = 1
    } else {
      nesting_val = extract(nesting_landscape, cbind(xprop, yprop))
    }
    
    rand = runif(1, 0, 1)
    if(rand < nesting_val){
      colony_x = c(colony_x, xprop)
      colony_y = c(colony_y, yprop)
    }
  }
  
  colony_data = as.data.frame(cbind(colonyid, colony_sizes, colony_x, colony_y))
  
  ##### Define trap characteristics #####
  trapid = 1:number_traps
  grid_size = sqrt(number_traps)
  step = trapgrid_size/(grid_size-1)
  trap_x = (landscape_size - trapgrid_size)/2 + step*(0:(grid_size-1))
  trap_y = (landscape_size - trapgrid_size)/2 + step*(0:(grid_size-1))
  coords = expand.grid(trap_x = trap_x, trap_y = trap_y)
  trap_data = as.data.frame(cbind(trapid, coords))
  trap_data$fq <- terra::extract(resource_landscape, trap_data[, c("trap_x", "trap_y")])[, 2]
  
  
  # optional plotting step to visualize traps and colonies
  #plot(nesting_landscape, col = c("white", "black"))
  #points(colony_data$colony_x, colony_data$colony_y, col = "red", pch = 16, cex = 1.5)
  #points(trap_data$trap_x, trap_data$trap_y, col = "blue", pch = 16, cex = 1.5)
  
  
  ##### Create landscape-wide visitation matrix #####
  # compute visitation rates for each colony at each grid cell, total visitation, and average foraging range
  visitation_and_foraging_range <- compute_visitation_on_raster(colonies = colony_data,
                                                                trap_data = trap_data,
                                                                resource_quality_rast = resource_landscape,
                                                                rho = rho,
                                                                theta = theta,
                                                                distance_decay = distance_decay)
  
  
  lambda_ik = visitation_and_foraging_range[[1]]
  colony_data$foraging_range = visitation_and_foraging_range[[2]]
  colony_data$D_i = visitation_and_foraging_range[[3]]
  
  ##### Start sampling #####
  yik = allocMatrix(nrow = number_colonies, ncol = number_traps, value = 0)
  N = numbees_start
  colony_data$ni = colony_sizes
  
  while(sum(yik) < sample_size){
    
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
    
    # sample trap
    trap = sample(trap_data$trapid, size = 1, replace = TRUE, prob = trap_data$trap_prob)
    
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
    
    print(paste("Now sampled", sum(yik), "of", sample_size, "bees"))
  }
  
  output = list(yik, colony_data, trap_data, visitation_rates)
  return(output)
}
