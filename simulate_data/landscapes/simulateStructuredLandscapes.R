##### Simulate structured landscapes #####
# Script Initiated: May 22, 2025
# By: J Melanson
# Goal: Simulate landscapes with (1) grid-like configuration of +/- habitat;
# (2) simulate agricultural-esque landscapes following [[insert authors here]]

# load packages
library(matrixStats)
library(ggplot2)


# define scales
landscape_size = 1100
grid_size = 50

# create empty matrix of landscape_size x landscape_size
empty_landscape = allocMatrix(nrow = landscape_size,
                              ncol = landscape_size)

# fill in grid pattern according to grid size
is_odd = function(number){
  if (number - ( 2* (number %/% 2)) == 1){
    return(TRUE)
  } else if (number - ( 2* (number %/% 2)) == 0){
    return(FALSE)
  }
}

for (x in 1:landscape_size){
  for (y in 1:landscape_size){
    xbool = is_odd(x %/% grid_size)
    ybool = is_odd(y %/% grid_size)
    if (xbool + ybool == 1){
      empty_landscape[x,y] = 1
    } else {
      empty_landscape[x,y] = 0
    }
  }
}

filled_landscape = empty_landscape


# plot grid landscape to check
df <- data.frame(
  x = rep(1:ncol(filled_landscape), each = nrow(filled_landscape)),
  y = rep(nrow(filled_landscape):1, times = ncol(filled_landscape)),  # reverse y for image-like plot
  value = as.vector(filled_landscape)
)

ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_raster() +  # much faster than geom_tile for large matrices
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_void()
# voila!



# Put some functions here for now, might move later
compute_visitation_rates_and_avg_distance <- function(colonies, resource_quality, rho, theta, landscape_size) {
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
    visitation <- exp(-0.5*(dist/rho)^2 + theta * resource_quality)
    
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



##### Simulate bee draws #####
draw_bees_colony_restricted = function(sample_size, # number of bees to sample
                       landscape_size, # integer, size of full resource landscape
                       colonygrid_size, # integer, size of colony distribution
                       trapgrid_size, # integer, size of trap grid
                       resource_landscape, # matrix containing floral quality across landscape
                       nesting_landscape, # matrix containing nest habitat quality across landscape
                       number_traps, # square integer
                       number_colonies, # positive integer
                       colony_sizes, #vector of length number_colonies
                       rho, # set param value
                       theta # set param value
){
  ##### Define colony characteristics #####
  colonyid = 1:number_colonies
  numbees_start = sum(colony_sizes)
  
  colony_x = c()
  colony_y = c()
  
  while (length(colony_x) < number_colonies){
    # propose a colony location
    xprop = round(runif(1, (landscape_size-colonygrid_size)/2, (landscape_size-colonygrid_size)/2 + colonygrid_size)) #colony locations are integers
    yprop = round(runif(1, (landscape_size-colonygrid_size)/2, (landscape_size-colonygrid_size)/2 + colonygrid_size)) #colony locations are integers
    
    # accept or reject colony location
    nesting_val = nesting_landscape[xprop, yprop]
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
  trap_data$fq <- mapply(function(x, y) resource_landscape[x, y], trap_data$trap_x, trap_data$trap_y)
  
  
  # optional plotting step to visualize traps and colonies
  # df <- data.frame(
  #   x = rep(1:ncol(nesting_landscape), each = nrow(nesting_landscape)),
  #   y = rep(1:nrow(nesting_landscape), times = ncol(nesting_landscape)),  # reverse y for image-like plot
  #   value = as.vector(nesting_landscape)
  # )
  # 
  # a = ggplot(df) +
  #   geom_raster(aes(x = x, y = y, fill = value)) +
  #   scale_fill_viridis_c() +
  #   coord_fixed() +
  #   geom_point(data = colony_data, aes(x = colony_x, y = colony_y)) +
  #   geom_point(data = trap_data, aes(x = trap_x, y = trap_y)) +
  #   theme_void()
    
  
  
  ##### Create landscape-wide visitation matrix #####
  # compute visitation rates for each colony at each grid cell, total visitation, and average foraging range
  visitation_and_foraging_range <- compute_visitation_rates_and_avg_distance(colonies = colony_data,
                                                                             resource_quality = resource_landscape,
                                                                             rho = rho,
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





