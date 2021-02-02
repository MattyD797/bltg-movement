# Date: 2/1/20
# Author: Matt Duggan

#Script to process and implement HMM on GPS birds (argos to follow)
#based on two states

#### Load packages and data ####

#load and install (missing) packages
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
foo( c("dplyr", "sp", "ggplot2", "ggmap", "plyr", "lubridate", "SDLfilter", "moveHMM", "bsam", "parallel", "adehabitatLT", "recurse", "momentHMM"))

#clean old objects from R
rm(list=ls(all=TRUE))


#see outputs in standard (not scientific notation)
options(scipen = 999)
options(max.print=100) #if you need to see the complete output

#load Hanmeer GPS data
setwd("./CSV/")
GPS_data <- read.csv("./GPS_BTGO_Haanmeer_nestR.csv"); head(GPS_data)

#Convert to Easting and Northing
cord.dec = SpatialPoints(cbind(GPS_data$long, GPS_data$lat), proj4string=CRS("+proj=longlat"))
GPS_data[,c(2,3)] <- data.frame(spTransform(cord.dec, CRS("+proj=utm +zone=31 +datum=WGS84"))); head(GPS_data)

#rename columns
names(GPS_data)[1] <- "id"
names(GPS_data)[2] <- "x"
names(GPS_data)[3] <- "y"
names(GPS_data)[4] <- "t"

GPS_data$t <- as.POSIXct(as.character(GPS_data$t), format = "%m/%d/%Y %H:%M:%S")
GPS_data$id <- as.factor(GPS_data$id)
GPS_data <- GPS_data[,c(2,3,4,1)]

#52.91564255769956, 5.438472498396251 are the reference coordinates of Haanmeer
#filter out appropriate UTM values for nesting
GPS_data.filtered <- filter(GPS_data, y >= 5765288)
GPS_data.filtered <- filter(GPS_data.filtered, x >= 294071)

animals <- GPS_data.filtered

plot(animals$x, animals$y, col = c(1:24)[as.numeric(animals$id)], 
     pch = ".", xlab = "x", ylab = "y", asp = 1)

#A circle of radius R is drawn around each point in the trajectory. The number of revisits is calculated
#as the number of segments of the trajectory passing through that circle.
popvisit = getRecursions(animals, 7)

#plot(popvisit, animals)
#Now we have two vectors and additional revisit stats in popvisit...Time to move to momentHMM. 

#practice run through with one animal track

colnames(animals) <- c("x", "y", "time", "ID")
animals <- animals[,c(4,1,2,3)]; head(animals)

#convert back to lon and lat
cord.dec = SpatialPoints(cbind(animals$x, animals$y), proj4string=CRS("+proj=utm +zone=31 +datum=WGS84"))
animals[,c(5,6)] <- data.frame(spTransform(cord.dec, CRS("+proj=longlat +zone=31 +datum=WGS84"))); head(animals)
colnames(animals) <- c("ID", "x", "y", "time", "lon", "lat")

#take the 3rd ID
rawData <- subset(animals, ID ==unique(ID)[3]); head(rawData)

#Need to set up an even amount of time for every hou, so we predict locations
crwOut <- crawlWrap(obsData=rawData, timeStep="hour",
                    theta=c(6.855, -0.007), fixPar=c(NA,NA))

birdData <- momentuHMM::prepData(data=crwOut)

birdData$hour <- as.integer(strftime(birdData$time, format = "%H", tz="CET"))
acf(birdData$step[!is.na(birdData$step)],lag.max=300)

stateNames <- c("migrating","nesting", "foraging")

dist = list(step = "gamma", angle = "wrpcauchy")

Par0_m1 <- list(step=c(13745.1594342787, 35.4968291160, 142.1692815640, 11807.5643056707, 0.0007897712, 206.3706396916),angle=c(0.00000000000002777020918302507471,0.70603171596917158048256624169881, 0.00000000000000000000000001813604))

m1 <- momentuHMM::fitHMM(
                         data = birdData, 
                         nbStates = 3, 
                         dist = dist, 
                         Par0 = Par0_m1, 
                         estAngleMean = list(angle=FALSE), 
                         stateNames = stateNames)
plot(m1, plotCI = TRUE, covs = data.frame(hour=12))
momentuHMM::plotStates(m1)
momentuHMM::plotPR(m1)









