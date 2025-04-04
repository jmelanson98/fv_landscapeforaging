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
