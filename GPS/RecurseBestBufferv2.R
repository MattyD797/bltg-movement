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
foo( c("ggspatial", "tibble", "rnaturalearth","rnaturalearthdata","dplyr", "sp", "ggplot2", "ggmap", "plyr", "lubridate", "SDLfilter", "moveHMM", "bsam", "parallel", "adehabitatLT", "recurse", "geosphere", "raster"))

#clean old objects from R
rm(list=ls(all=TRUE))

#see outputs in standard (not scientific notation)
options(scipen = 999)
options(max.print=100) #if you need to see the complete output
#

#load data
setwd("./CSV/")
data1 <- read.csv("./GPS_BTGO_Haanmeer_nestR.csv"); head(data1)

#these two birds had too much missing data
data1 <- dplyr::filter(data1, burst != "1009-2012")
data1 <- dplyr::filter(data1, burst != "2002-2014")
data1 <- dplyr::filter(data1, burst != "1002-2012")
data1 <- dplyr::filter(data1, burst != "2016-2013")
data1 <- dplyr::filter(data1, burst != "2010-2013")
data1 <- dplyr::filter(data1, burst != "2002-2015")
data1 <- dplyr::filter(data1, burst != "2041-2013")
data1 <- dplyr::filter(data1, burst != "2004-2014")

#filter out the nesting area coordinates
data1.filt <- data1[data1[,3] > 52,]

#world <- ne_countries(scale = "medium", returnclass = "sf")
NLD <- readRDS("gadm36_NLD_0_sp.rds")

#ggplot(data=world) + geom_sf() + geom_point(data=data1.filt, mapping = aes(x=long, y=lat), color="red") + coord_sf()
ggplot(NLD) +
  theme_minimal()+  # no backgroundcolor
  geom_polygon( aes(x = long, y = lat, group = group),
                color = "white",   # color is the lines of the region
                fill = "#9C9797")+ # fill is the fill of every polygon.
  geom_point(data=data1.filt, mapping = aes(x=long, y=lat, color = burst))+
  coord_map()

#Convert lat and lon to Easting and Northing
cord.dec = SpatialPoints(cbind(data1.filt$long, data1.filt$lat), proj4string=CRS("+proj=longlat"))
locs <- data.frame(spTransform(cord.dec, CRS("+proj=utm +zone=31 +datum=WGS84")))
head(locs)

data1.filt$long <- locs$coords.x1
data1.filt$lat <- locs$coords.x2

#rename columns and reformat
names(data1.filt)[1] <- "id"
names(data1.filt)[2] <- "x"
names(data1.filt)[3] <- "y"
names(data1.filt)[4] <- "t"

#change the data types for columns
data1.filt$t <- as.POSIXct(as.character(data1.filt$t), format = "%m/%d/%Y %H:%M:%S")
data1.filt$id <- as.factor(data1.filt$id)

#reorder columns
data1.utm <- data1.filt[,c(2,3,4,1)]


logbook <- data.frame()


#get the nest data for the individual
indnestloc <- read.csv("./IndNestLocation.csv")

for(j in unique(data1.utm$id)){
  individual1 <- droplevels(dplyr::filter(data1.utm, id == j))
  individual1.nest <- droplevels(dplyr::filter(indnestloc, ID == j))
 
  cord.dec = SpatialPoints(cbind(individual1.nest[1,]$lon, individual1.nest[1,]$lat), proj4string=CRS("+proj=longlat"))
  coords.known <- data.frame(spTransform(cord.dec, CRS("+proj=utm +zone=31 +datum=WGS84")))
  for (i in seq(2,50,1)){
    #This is getting the recursion at a buffer of i between 2 - 200
    indvisit <- getRecursions(individual1, i)
    coords.pred <- individual1[which.max(indvisit$revisits),c(1,2)]
    #find the euclidean distance between the predicted site and the known site
    logbook[i-1,j] <- pointDistance(coords.known, coords.pred, type="Euclidean", lonlat = FALSE)
    
  }
}

which.min(rowMeans(logbook))+1
plot(rowMeans(logbook))
rowMeans(logbook)[which.min(rowMeans(logbook))]














