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


#set the colors for each individual tag
# get domain of numeric data
#domain <- factor(movedata.breed$tag.local.identifier)

# make palette
#factpal <- colorFactor(viridis(36), domain)

#plot the data to visualize on a map
#leaflet(data=movedata.breed) %>% addTiles(group="Map")%>%
# addProviderTiles("Esri.WorldImagery",group="Sattelite") %>%
# addCircleMarkers(~location.long, ~location.lat, color = ~factpal(tag.local.identifier),label =~as.character(timestamp),radius = 3,stroke = F,fillOpacity = .7)%>%
# addLayersControl(baseGroups = c("Map", "Satellite"),options=layersControlOptions(collapsed=F))


#for Mo's data only! ----> Nate's birds may need different corrections (I don't think they do though)
#
#there are two mistakes in the ID column, this corrects these
movedata$tag.local.identifier[movedata$tag.local.identifier == "151862"] <- "151682"
movedata$tag.local.identifier[movedata$tag.local.identifier == "157539"] <- "157535"


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



#erases duplicate locations that are both less than 36 seconds apart in time AND less than 10 meters apart in space
loop_object <- vector("list", length(unique(data$id)))

for(i in unique(data$id)){

  x <- data %>% dplyr::filter(id == i)
  loop_object[[i]] <- dupfilter(x, step.time = 0.01, step.dist = 0.01, conditional = TRUE)
}

#this binds them back together to get the fixed data.
BTGO_dup.rem <- as.data.frame(do.call(rbind, loop_object)); head(BTGO_dup.rem)


#Hays et al. 2001 (doi: 10.1006/anbe.2001.1685) recommends you use these only

BTGO_dup.rem_qi <- ddfilter_loop(BTGO_dup.rem, qi = c(-1,1,2,3)); head(BTGO_dup.rem_qi) #this then is the argos data we can use further on!