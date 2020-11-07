# Date: 11.6.2020
# Author: Luke Wilde

# Script to learn ARGOS data structure and subsample GPS to that structure

foo <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}

#load or install packages
foo( c("dplyr", "sp", "ggplot2", "geosphere", "adehabitatHR", "rgdal", "ggmap", "mapproj", "maptools", "rgeos", "raster", "leaflet", "move", "ctmm", "BBMM", "PBSmapping", "shiny", "plyr", "lubridate", "viridisLite", "rworldmap", "adehabitatLT", "adehabitatMA"))