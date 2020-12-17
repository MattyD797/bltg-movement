# Date: 12/16/20
# Author: Luke Wilde

#Script to process and implement HMM on GPS birds (argos to follow)
#based on two states

#### Load packages and data ####

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
foo( c("dplyr", "sp", "ggplot2", "ggmap", "plyr", "lubridate", "SDLfilter", "moveHMM", "bsam", "parallel", "adehabitatLT"))

#clean old objects from R
rm(list=ls(all=TRUE))

#see outputs in standard (not scientific notation)
options(scipen = 999)
options(max.print=1600) #if you need to see the complete output
#

#read in the data and preview

setwd("C:/Users/ralph/OneDrive/Desktop/Birds/bltg-movement/CSV/")
data1 <- read.csv("./BTGO_Haanmeer_nestR.csv"); head(data1)

#these two birds had too much missing data
data2 <- dplyr::filter(data1, burst != "1009-2012")
data <- dplyr::filter(data2, burst != "2002-2014")

#### Linearly interpolate missing data ####

#convert time object to time-date format
data$date <- as.POSIXct(as.character(data$date), format = "%m/%d/%Y %H:%M:%S")

#getting error message Error in as.ltraj(xy = data[, c("long", "lat")], date = data$date, id = data$burst) : non unique dates for a given burst
#so I will remove duplicates with the dupfilter from SDLfilter

#clean the data and put in correct format
#filter
names(data)[4] <- "DateTime"
names(data)[2] <- "lon"
names(data)[3] <- "lat"
names(data)[1] <- "id"
data$qi <- 3
head(data)

#remove missing data
data <- na.omit(data)
#create time object
data$DateTime <- as.POSIXct(strptime(as.character(data$DateTime),"%Y-%m-%d %H:%M:%S", tz="GMT"))

head(data)

coordinates(data) <- c("lon", "lat")
proj4string(data) <- CRS("+proj=longlat +zone=31 +datum=WGS84")

#transform lat/long to UTM and rename the columns for ease down the road
#UTM is in meters, but is based on zones - the zone of the region in Netherlands is 31N or 32N (western Netherlands)
utm <- spTransform(data, CRSobj="+proj=utm +zone=31 +datum=WGS84")
utm_locs <- data.frame(as(utm, "SpatialPoints"))
colnames(utm_locs) <- c("UTM_W", "UTM_N")

DOY <- as.POSIXct(strptime(as.character(data$DateTime), "%Y-%m-%d", tz="GMT"))

BTGO_move_1 <- data.frame(id = data$id,date = data$DateTime, DOY = DOY, UTM_W = utm_locs$UTM_W, UTM_N = utm_locs$UTM_N); head(BTGO_move_1)

#convert to a traj object
BTGO_move <- as.ltraj(xy = BTGO_move_1[, c("UTM_W", "UTM_N")], date = BTGO_move_1$date, id = BTGO_move_1$id)



