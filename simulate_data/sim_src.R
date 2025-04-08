##### Functions for simulation #####

#calculate euclidean distance between pairs of coords
#this is dumb as shit but I'm on an airplane and I don't know the R command :'(
#replace later!
calculate_dist = function(colony_coords, trap_coords){
  xdist = abs(trap_coords[[1]] - colony_coords[[1]])
  ydist = abs(trap_coords[[2]] - colony_coords[[2]])
  hypotenuse = sqrt(xdist^2 + ydist^2)
  return(hypotenuse)
}


# function to compute distance matrix for each colony to each grid cell
# via ChatGPT
compute_distances <- function(colonies, landscape_size) {
  dist_matrix <- list()
  
  for (i in 1:nrow(colonies)) {
    colony_x <- colonies$colony_x[i]
    colony_y <- colonies$colony_y[i]
    
    # Create a grid of coordinates for the landscape
    x_coords <- 1:landscape_size[1]
    y_coords <- 1:landscape_size[2]
    
    # Compute the distance matrix for this colony (distance from colony to all grid cells)
    dist_x <- outer(x_coords, rep(colony_x, landscape_size[2]), "-")  # Distance in x-direction
    dist_y <- outer(rep(colony_y, landscape_size[1]), y_coords, "-")  # Distance in y-direction
    
    # Compute Euclidean distance
    dist_matrix[[i]] <- sqrt(dist_x^2 + dist_y^2)
  }
  
  return(dist_matrix)
}

# compute visitation rates for each colony
# via ChatGPT
compute_visitation_rates <- function(colonies, resource_quality, beta, theta, landscape_size) {
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