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





