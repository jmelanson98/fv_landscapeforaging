# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
load.packages <- function(){
  
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("here",
                "readr",
                "readxl",
                "writexl", 
                "openxlsx",
                "janitor",
                "filesstrings",
                "textshaping",
                "extrafont",
                "qrcode", 
                "lubridate", 
                "naniar", 
                "Hmisc",
                "sf",
                "maps",
                "broom",
                #"rgdal",
                "raster",
                "dismo",
                "terra",
                "spData",
                "tmap",
                "leaflet",
                #"rgeos",
                #"maptools",
                "mapview",
                "stars",
                "landscapemetrics",
                "ggrepel",
                "vctrs",
                "measurements", 
                "arsenal", 
                "vegan",
                "taxize",
                "iNEXT", 
                "tidyverse",
                "stringr",
                "ggplot2", 
                "purrr",
                "magrittr",
                "dplyr",
                "tidylog",
                "data.table",
                "diffr",
                "RColorBrewer")
  
  ipak(packages)
}

