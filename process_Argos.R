# Date: 10/15/20
# Author: Luke Wilde

# Script for processing the Argos data from Verhoeven et al. 2020 and calculating summary statistics on kernel area

#project: inferring breeding success from tags
#species: Black-tailed Godwit (BTGO, Limosa limosa)

#
#begin script
#

#### Load the packages ####

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
foo( c("dplyr", "sp", "ggplot2", "geosphere", "adehabitatHR", "rgdal", "ggmap", 
       "mapproj", "maptools", "rgeos", "raster", "leaflet", "move", "ctmm", "BBMM", 
       "PBSmapping", "shiny", "plyr", "lubridate", "viridisLite", "rworldmap", "sp", "adehabitatLT"))

#clean old objects from R
rm(list=ls(all=TRUE))

#see outputs in standard (not scientific notation)
options(scipen = 999)
#


#location of the working directory (input data)
setwd("bltg-movement")

#### Load and format the movement data ####

#read in the data
movedata <- read.csv("./CSV/Argos_BTGO_Haanmeer.csv"); head(movedata)

#Verhoeven et al. 2020 defines breeding area for Haanmeer to be > 52N, visually, all breeders in our study ought to have stayed E of the prime meridian (0 longitude)
#subset the data to include only the breeding location
movedata_breed <- filter(movedata, location.lat >= 52)
movedata_breed <- filter(movedata_breed, location.long >= 0)

#check
min(movedata_breed$location.lat)

#set the colors for each individual tag
# get domain of numeric data
domain <- factor(movedata_breed$tag.local.identifier)
# make palette
factpal <- colorFactor(viridis(36), domain)

#plot the data to visualize on a map
leaflet(data=movedata_breed) %>% addTiles(group="Map")%>%
  addProviderTiles("Esri.WorldImagery",group="Sattelite") %>%
  addCircleMarkers(~location.long, ~location.lat, color = ~factpal(tag.local.identifier),label =~as.character(timestamp),radius = 3,stroke = F,fillOpacity = .7)%>%
  addLayersControl(baseGroups = c("Map", "Satellite"),options=layersControlOptions(collapsed=F))

#### organize the data - formatting and projection ####
movedata <- movedata_breed[with(movedata_breed,order(movedata_breed$tag.local.identifier,movedata_breed$timestamp)),]

#define coordinates within the data - Argos is always WGS 84
coordinates(movedata) <- c("location.long", "location.lat")
proj4string(movedata) <- CRS("+proj=longlat +datum=WGS84")

#transform lat/long to UTM and rename the columns for ease down the road
#UTM is in meters, but is based on zones - the zone of the region in Netherlands is 31N or 32N (western Netherlands)
utm <- spTransform(movedata, CRS("+proj=utm +zone=31 +ellps=bessel"))
utm_locs <- data.frame(as(utm, "SpatialPoints"))
colnames(utm_locs) <- c("UTM_W", "UTM_N")

#inspect
head(movedata)

movedata$timestamp <- substr(movedata$timestamp, 1, nchar(as.character(movedata$timestamp)))

#define date and time
date_time <- as.POSIXct(strptime(as.character(movedata$timestamp),"%Y-%m-%d %H:%M:%S", tz="GMT"))
fix_date <- as.Date(format(date_time, "%Y-%m-%d"))

month <- as.numeric(format(date_time, "%m"))
year <- as.numeric(format(date_time, "%Y")) #only saves the unique digits of the year so 2015 looks like 0015
day <- as.numeric(format(date_time, "%d"))
time <- strftime(date_time, format="%H:%M:%S")

#create new data and imporved frame for HR analysis
final_move1 <- data.frame(ID=as.factor(movedata$tag.local.identifier),Height=movedata$argos.altitude,UTM_W=utm_locs$UTM_W,UTM_N=utm_locs$UTM_N,fix_date,month,year,day,time,timestamp=date_time)
final_move <- na.omit(final_move1)

#inspect
head(final_move)

#define the coordinates and projection of the new df
coordinates(final_move) <- c("UTM_W", "UTM_N")
projection(final_move) <- CRS("+proj=utm +zone=31 ellps=bessel")




#### Incorporate ID information and Subset to reproductive cycles ####
#load the individual information of tagged birds
individual <- read.csv("./CSV/Argos_Individual_Information_BTGO_Haanmeer.csv"); head(individual)
nest <- read.csv("./CSV/Argos_Nesting_Information_BTGO_Haanmeer.csv");head(nest)

#rename datasets to have consistent names for: tags, nest, and UTM coordinates
names(individual)[4] <- "ID"; head(individual); head(final_move1)
names(individual)[6] <- "NestID"; head(individual); head(nest)
names(nest)[c(3,4)] <- c("nestloc_UTM_W", "nestloc_UTM_N"); head(nest)

#now, join the nest and individual datasets first
nesters <- right_join(nest, individual, by = "NestID"); head(nesters)

#change ID to factor to be compatible with movement data
nesters$ID <- as.factor(nesters$ID)
nesters_final_move <- right_join(nesters, final_move1, by = "ID"); head(nesters_final_move)

#### Format recorded nesting dates and join to movements ####

#reformat the nesting dates to numerical
nesters_final_move$StartKnownActive <- gsub('Apr','04', nesters_final_move$StartKnownActive)
nesters_final_move$StartKnownActive <- gsub('May','05', nesters_final_move$StartKnownActive)
nesters_final_move$StartKnownActive <- gsub('Jun','06', nesters_final_move$StartKnownActive)
nesters_final_move$EndKnownActive <- gsub('Apr','04', nesters_final_move$EndKnownActive)
nesters_final_move$EndKnownActive <- gsub('May','05', nesters_final_move$EndKnownActive)
nesters_final_move$EndKnownActive <- gsub('Jun','06', nesters_final_move$EndKnownActive)

#then, reformat the dates as date objects to allow easy subsetting
nesters_final_move$nest_start_date <- as.POSIXct(strptime(as.character(nesters_final_move$StartKnownActive),"%d-%m-%Y", tz="GMT"))
nesters_final_move$nest_end_date <- as.POSIXct(strptime(as.character(nesters_final_move$EndKnownActive),"%d-%m-%Y", tz="GMT"))

nesters_final_move$nest_start_date <- as.Date(strptime(nesters_final_move$nest_start_date,"%Y-%m-%d"))
nesters_final_move$nest_end_date <- as.Date(strptime(nesters_final_move$nest_end_date,"%Y-%m-%d"))

#I havent found an easier way to do this, so I simply ask R to substitue the incorrectly formatted years with ones that agree with other formats in the code
nesters_final_move$nest_start_date <- as.Date(gsub('0016','2016',nesters_final_move$nest_start_date))
nesters_final_move$nest_start_date <- as.Date(gsub('0015','2015', nesters_final_move$nest_start_date))
nesters_final_move$nest_start_date <- as.Date(gsub('0017','2017', nesters_final_move$nest_start_date))
nesters_final_move$nest_end_date <- as.Date(gsub('0016','2016',nesters_final_move$nest_end_date))
nesters_final_move$nest_end_date <- as.Date(gsub('0015','2015', nesters_final_move$nest_end_date))
nesters_final_move$nest_end_date <- as.Date(gsub('0017','2017', nesters_final_move$nest_end_date))

#estimate the days a nest survived
nesters_final_move$nest_min_days_survived <- as.numeric(nesters_final_move$nest_end_date - nesters_final_move$nest_start_date); hist(nesters_final_move$nest_min_days_survived)

#break each down into their own columns
nesters_final_move$nest_start_month <- as.numeric(format(nesters_final_move$nest_start_date, "%m"))
nesters_final_move$nest_start_year <- as.numeric(format(nesters_final_move$nest_start_date, "%Y")) #only saves the unique digits of the year so 2015 looks like 0015
nesters_final_move$nest_start_day <- as.numeric(format(nesters_final_move$nest_start_date, "%d"))

nesters_final_move$nest_end_month <- as.numeric(format(nesters_final_move$nest_end_date, "%m"))
nesters_final_move$nest_end_year <- as.numeric(format(nesters_final_move$nest_end_date, "%Y")) #only saves the unique digits of the year so 2015 looks like 0015
nesters_final_move$nest_end_day <- as.numeric(format(nesters_final_move$nest_end_date, "%d"))

#erase second Colourcode column and rename the remaining column
nesters_final_move$Colourcode.y <- NULL
names(nesters_final_move)[1] <- "Colourcode"

#we finally have a well formatted datafile with UTM coordinates, nesting and ID information
head(nesters_final_move); str(nesters_final_move)