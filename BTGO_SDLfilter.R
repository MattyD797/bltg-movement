#Author: Luke Wilde
#Date: 11.13.20

#Script to filter Argos/GPS locations for errors

#create function to load and install (missing) packages
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
foo( c("dplyr", "sp", "ggplot2", "geosphere", "adehabitatHR", "rgdal", "ggmap", "mapproj", "maptools", "rgeos", "raster", "leaflet", "move", "ctmm", "BBMM", "PBSmapping", "shiny", "plyr", "lubridate", "viridisLite", "rworldmap", "adehabitatLT", "adehabitatMA", "SDLfilter"))



#location of the working directory (input data)
setwd("C:/Users/14064/Dropbox/BTGO Movmement Study/Data/Data from Senner 2015/")

#### Load and format the movement data ####

#read in the data
movedata <- read.csv("./Continental black-tailed godwits (data from Senner et al. 2015).csv"); head(movedata)

#clean the data and put in correct format
#filter
data <- movedata[,c(3,4,5,13,26)]; head(data)
names(data)[1] <- "DateTime"
names(data)[2] <- "lon"
names(data)[3] <- "lat"
names(data)[4] <- "qi"
names(data)[5] <- "id"
#remove missing data
data <- na.omit(data)
#create time object
data$DateTime <- as.POSIXct(strptime(as.character(data$DateTime),"%Y-%m-%d %H:%M:%S", tz="GMT"))


#will want a for loop for this

#select a bird ID
bird_134755 <- data %>% filter(id == "134755")

#filter duplicates in time and in space, if the dup in space fall within the time threshold (conditional = TRUE). step parameters are in hours and km resp.
bird_134755 <- dupfilter(bird_134755, step.time = 0.01, step.dist = 0.3, conditional = TRUE)
#next, remove fixes from B,Z,0 class
bird_134755 <- ddfilter(bird_134755, qi = c(3,5,6,7))

#calculate the track parameters from the data:
  #pTime - hours from hours from a previous fix
  #sTime - are hours to a subsequent fix 
  #pDist - straight distances in kilometres from a previous
  #sDist - straight distances in kilometres to a subsequent fix 
  #pSpeed - linear speed (km/h) from a previous
  #sSpeed - to a subsequent fix 
  #inAng - the degree between the bearings of lines joining successive location points
  # meanSpeed - calculated over 'days' length of time
  #meanAngle - daily

track_param(bird_134755, param = c("time", "distance", "speed", "angle", "mean speed", "mean angle"), days = 1)

