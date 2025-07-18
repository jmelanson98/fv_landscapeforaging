##### Simulate structured landscapes #####
# Script Initiated: May 22, 2025
# By: J Melanson
# Goal: Simulate landscapes with (1) grid-like configuration of +/- habitat;
# (2) simulate agricultural-esque landscapes following [[insert authors here]]

# load packages
library(matrixStats)
library(ggplot2)
library(terra)



#######################################################################
### Simulate a checkerboard landscape of nesting habitat quality
#######################################################################


# define scales
landscape_size = 1100
grid_size = 50

# create an empty raster
nest_raster <- rast(nrows = landscape_size, ncols = landscape_size, 
                     xmin = 0, xmax = landscape_size, 
                     ymin = 0, ymax = landscape_size)

# get xy coordinates of center of each pixel
xy <- xyFromCell(nest_raster, 1:ncell(nest_raster))

# determine which checkerboard grid they're in
col_index <- floor(xy[,1] / 50)
row_index <- floor(xy[,2] / 50)

# assign values
eps = 10^-12 # use small epsilon values so that we can later calculate log probability of habitat suitability
checker_values <- eps + (col_index + row_index) %% 2
values(nest_raster) <- checker_values

# plot
plot(nest_raster, col = c("white", "black"))
# voila!


#######################################################################
### Put some functions here for now, might move later
#######################################################################


##### simulate floral quality landscapes as raster #####
simulateLandscapeRaster = function(landscape_size, # integer, same for x and y
                             resource_range) {
  # use Brownian variogram to simulate spatial distribution of resources (as in Pope & Jha)
  
  
 
  # get xy coordinates of center of each pixel
  xy <- xyFromCell(nest_raster, 1:ncell(nest_raster))
  
  # determine which checkerboard grid they're in
  col_index <- floor(xy[,1] / 50)
  row_index <- floor(xy[,2] / 50)
  
  # assign values
  checker_values <- (col_index + row_index) %% 2
  values(nest_raster) <- checker_values
  
  
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
compute_visitation_on_raster <- function(colonies, resource_quality_rast, rho, theta) {
  
  # get coordinates for center of each pixel
  xy_coords <- crds(resource_quality_rast, df = TRUE)  # matrix of x and y
  
  # flatten raster values into a vector
  rq_vals <- values(resource_quality_rast)
  
  visitation_rates <- list() # will be a list of rasters
  avg_distances <- numeric(nrow(colonies))
  total_vis <- numeric(nrow(colonies))
  
  for (i in 1:nrow(colonies)) {
    colony_x <- colonies$colony_x[i]
    colony_y <- colonies$colony_y[i]
    
    # compute Euclidean distance from colony to each pixel
    dx <- xy_coords[, 1] - colony_x
    dy <- xy_coords[, 2] - colony_y
    dist <- sqrt(dx^2 + dy^2)
    
    # compute visitation rate of colony at each pixel
    visitation_vals <- exp(-0.5 * (dist / rho)^2 + theta * rq_vals)
    
    total_visitation <- sum(visitation_vals, na.rm = TRUE)
    avg_distance <- if (total_visitation > 0) {
      sum(visitation_vals * dist, na.rm = TRUE) / total_visitation
    } else {
      NA_real_
    }
    
    vis_rast <- resource_quality_rast  # clone geometry
    values(vis_rast) <- visitation_vals
    
    print(paste("Saving visitation for colony", i, sep = " "))
    visitation_rates[[i]] <- vis_rast
    avg_distances[i] <- avg_distance
    total_vis[i] <- total_visitation
  }
  
  return(list(
    visitation_rates = visitation_rates,
    avg_distances = avg_distances,
    total_vis = total_vis
  ))
}


##### simulate bee draws #####
draw_bees_colony_restricted = function(sample_size, # number of bees to sample
                       landscape_size, # integer, size of full resource landscape
                       colonygrid_size, # integer, size of colony distribution
                       trapgrid_size, # integer, size of trap grid
                       resource_landscape, # raster containing floral quality across landscape
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
    xprop = runif(1, (landscape_size-colonygrid_size)/2, (landscape_size-colonygrid_size)/2 + colonygrid_size)
    yprop = runif(1, (landscape_size-colonygrid_size)/2, (landscape_size-colonygrid_size)/2 + colonygrid_size)
    
    # accept or reject colony location
    nesting_val = extract(nesting_landscape, cbind(xprop, yprop))
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
  plot(nesting_landscape, col = c("white", "black"))
  points(colony_data$colony_x, colony_data$colony_y, col = "red", pch = 16, cex = 1.5)
  points(trap_data$trap_x, trap_data$trap_y, col = "blue", pch = 16, cex = 1.5)
  
  
  ##### Create landscape-wide visitation matrix #####
  # compute visitation rates for each colony at each grid cell, total visitation, and average foraging range
  visitation_and_foraging_range <- compute_visitation_on_raster(colonies = colony_data,
                                                                             resource_quality_rast = resource_landscape,
                                                                             rho = rho,
                                                                             theta = theta)
  

  visitation_rates = visitation_and_foraging_range[[1]]
  colony_data$foraging_range = visitation_and_foraging_range[[2]]
  colony_data$D_i = visitation_and_foraging_range[[3]]

  ##### Compute fixed values lambda_ik #####
  # calculate lambda_ik, e.g., visitation rates of each colony to specific traps
  lambda_ik = allocMatrix(nrow = number_colonies, ncol = number_traps, value = 0)
  trapcoords = as.matrix(trap_data[, c("trap_x", "trap_y")])

  for (i in seq_len(number_colonies)) {
    print(i)
    print(visitation_rates[[i]])
    class(visitation_rates[[i]])
    nlyr(visitation_rates[[i]])
    visit_rast <- visitation_rates[[i]]  # each is a raster of size [xmax x ymax]

    # pull out all trap visitation values for that colony
    vals <- terra::extract(visit_rast, trapcoords)
    lambda_ik[i, ] <- vals[, 1]
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




#######################################################################
### simulate one iteration of bee data for model
#######################################################################
landscape1 = simulateLandscapeRaster(landscape_size = 1100,
                                     resource_range = 10)

sample1 = draw_bees_colony_restricted(sample_size = 1000,
                                      landscape_size = 1100,
                                      colonygrid_size = 700,
                                      trapgrid_size = 300,
                                      resource_landscape = landscape1,
                                      nesting_landscape = nest_raster,
                                      number_traps = 25,
                                      number_colonies = 500,
                                      colony_sizes = rep(100,500),
                                      rho = 75,
                                      theta = 0.5)

#######################################################################
### Fit model with informative colony prior
#######################################################################

###### Prep trap data, colony data, yobs


###### Prep nest habitat matrix
# CHECK RASTER RESOLUTION AND MIN/MAX!!!
xres <- terra::xres(nest_raster)
yres <- terra::yres(nest_raster)
xmin <- terra::xmin(nest_raster)
ymin <- terra::ymin(nest_raster)

# Convert nesting landscape raster to a matrix
nest_mat = as.matrix(nest_raster, wide = TRUE)

# flip along y axis so that [row, column] indices match to (x, y) coordinates
nest_mat = nest_mat[nrow(nest_mat):1, ]

# convert to a log probability matrix (normalize and take the log)
log_prob_nest_mat = log(nest_mat/sum(nest_mat))


###### Fit GAM
# Convert raster to data.frame
df <- as.data.frame(nest_raster, xy = TRUE)
colnames(df) <- c("x", "y", "suitability")

# Scale coordinates (important for numerical stability in Stan)
df$x_scaled <- scale(df$x)[,1]
df$y_scaled <- scale(df$y)[,1]

# Fit GAM
gam_fit <- gam(suitability ~ s(x_scaled, y_scaled, bs = "tp", k = 100), data = df)


df$predicted <- predict(gam_fit, newdata = df)
ggplot(df, aes(x = x, y = y)) +
  geom_raster(aes(fill = predicted)) +
  coord_equal() +
  scale_fill_viridis_c() +
  ggtitle("GAM-Predicted Suitability")

### hahaha that looks like shit


###### Create data list for stan
data = list()
data$C = nrow(yobs)
data$K = 25
data$L = 1100
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$y = yobs
data$lowerbound = 200
data$upperbound = 900
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$rho_center = 4.5
data$rho_sd = 0.5
data$res = 1 #resolution of matrix--important!!
data$nesting_landscape = log_prob_nest_mat


###### Fit model in Stan
# not possible in stan lol
stanFit = stan(file = "models/colony_prior.stan",
               data = data, seed = 5838299,
               chains = 4, cores = 4,
               iter = 10000,
               verbose = TRUE)
