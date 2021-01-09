# Date: 1/6/20
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
foo( c("dplyr", "sp", "ggplot2", "ggmap", "plyr", "lubridate", "SDLfilter", "moveHMM", "bsam", "parallel", "adehabitatLT"))

#clean old objects from R
rm(list=ls(all=TRUE))

#see outputs in standard (not scientific notation)
options(scipen = 999)
options(max.print=100) #if you need to see the complete output
#

#load data
setwd("./CSV/")
data1 <- read.csv("./GPS_BTGO_Haanmeer_nestR.csv"); head(data1)

#Convert to Easting and Northing
cord.dec = SpatialPoints(cbind(data1$long, data1$lat), proj4string=CRS("+proj=longlat"))
locs <- data.frame(spTransform(cord.dec, CRS("+proj=utm +zone=31 +datum=WGS84")))
head(locs)

data1$long <- locs$coords.x1
data1$lat <- locs$coords.x2

#rename columns
names(data1)[1] <- "id"
names(data1)[2] <- "x"
names(data1)[3] <- "y"
names(data1)[4] <- "t"

data1$t <- as.POSIXct(as.character(data1$t), format = "%m/%d/%Y %H:%M:%S")
data1$id <- as.factor(data1$id)

data1 <- data1[,c(2,3,4,1)]
##52.91564255769956, 5.438472498396251 are the reference coordinates of Haanmeer
#longH<-5.438472498396251
#data1$y <- latH - data1$y
#data1$x <- longH - data1$x

individual1 <- droplevels(dplyr::filter(data1, id == "1002-2012"))

#Find recursive locations
indvisit <- getRecursions(individual1, 7)

#Plot recursions
plot(indvisit, individual1, 
     xlim = c(660000, 670000),
     ylim = c(5855000, 5875000),
     legendPos = c(670000, 5860000))
drawCircle(655000, 5865000, 7)
hist(indvisit$revisits, breaks = 20, main = "", xlab = "Revisits (radius = 7)")
summary(indvisit$revisits)
head(indvisit$revisitStats)

#filter a second individual
individual2 <- droplevels(dplyr::filter(data1, id == "2004-2013"))

#Find recursive locations
indvisit <- getRecursions(individual2, 7)

#bind the two individuals
birds <- rbind(individual1, individual2)

#getRecursions
head(popvisit$revisitStats)
plot(popvisit, birds, legendPos = c(600000, 5500000))

#Use K-means clustering to cluster the (x,y) coordinates 
#of the top 20% of the locations
visitThreshold = quantile(popvisit$revisits, 0.8)
popCluster = kmeans(birds[popvisit$revisits > visitThreshold,c("x", "y")], centers = 3)
plot(birds$x, birds$y, col = c("red", "darkblue")[as.numeric(birds$id)], 
     pch = ".", xlab = "x", ylab = "y", asp = 1, cex = 2)
with(birds[popvisit$revisits > visitThreshold,],
     points(x, y, col = c(alpha("red", 0.5), alpha("darkblue", 0.5))[as.numeric(id)], 
            pch = c(15:17)[popCluster$cluster]) )
legend("topleft", pch = 15:17, legend = paste("cluster", 1:3), bty = "n")

#filter out appropriate UTM values for nesting
data1 <- filter(data1, y >= 5765288)
data1 <- filter(data1, x >= 294071)
count <- 1
x <- list()
par(mfrow=c(7, 5), mar = c(1,1,1,1))
for(val in unique(data1$id)){
  individual <- droplevels(dplyr::filter(data1, id == val))
  
  #Find recursive locations
  indvisit <- getRecursions(individual, 6)
  x[[count]] <- indvisit
  head(indvisit$revisitStats)
  plot(indvisit, 
       individual, 
       xlim = c(600000, 720000),
       ylim = c(5780000, 5880000), 
       legendPos = c(700000, 5800000))
  count <- count + 1
}
legend("bottomright", legend=levels(factor(individual$)), text.col=seq_along(levels(factor(data.c))))
count
head(indvisit)



