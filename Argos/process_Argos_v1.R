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
foo( c("dplyr", "sp", "ggplot2", "geosphere", "adehabitatHR", "rgdal", "ggmap", "mapproj", "maptools", "rgeos", "raster", "leaflet", "move", "ctmm", "BBMM", "PBSmapping", "shiny", "plyr", "lubridate", "viridisLite", "rworldmap", "adehabitatLT", "stringr"))


#clean old objects from R
rm(list=ls(all=TRUE))

#see outputs in standard (not scientific notation)
options(scipen = 999)
#options(max.print=999999) #if you need to see the complete output
#


#location of the working directory (input data)
setwd("./CSV/")

#### Load and format the movement data ####

#read in the data
movedata <- read.csv("./Haanmeer_Argos_MattLukeNathan.csv"); head(movedata)

#breakdown of fix quality:
#class 3 (error < 150m) = 24020 | 26%
#class 2 (150 m < error < 350 m) = 11876 | 13%
#class 1 (350 m < error < 1000 m) = 8424 | 9%
#class 0 (error > 1000) = 6081 | 7%
#class A (no estimate) = 13593 | 15%
#class B (no estimate) = 27688 | 30%
#class Z (invalid fix, failed) = 94 | < 1%

#movedata <- read.csv("Data from Senner 2015/Continental black-tailed godwits (data from Senner et al. 2015).csv"); head(movedata)

#Verhoeven et al. 2020 defines breeding area for Haanmeer to be > 52N, visually, all breeders in our study ought to have stayed E of the prime meridian (0 longitude)
#subset the data to include only the breeding location
#movedata.breed <- filter(movedata, location.lat >= 52)
movedata.breed <- filter(movedata, location.long >= 0)


#check
min(movedata.breed$location.lat)

unique(movedata.breed$tag.local.identifier)


#set the colors for each individual tag
# get domain of numeric data
domain <- factor(movedata.breed$tag.local.identifier)
# make palette
factpal <- colorFactor(viridis(36), domain)

#plot the data to visualize on a map
#leaflet(data=movedata.breed) %>% addTiles(group="Map")%>%
  # addProviderTiles("Esri.WorldImagery",group="Sattelite") %>%
  # addCircleMarkers(~location.long, ~location.lat, color = ~factpal(tag.local.identifier),label =~as.character(timestamp),radius = 3,stroke = F,fillOpacity = .7)%>%
  # addLayersControl(baseGroups = c("Map", "Satellite"),options=layersControlOptions(collapsed=F))

#### organize the data - formatting and projection ####
movedata <- movedata.breed[with(movedata.breed,order(movedata.breed$tag.local.identifier,movedata.breed$timestamp)),]

#define coordinates within the data - Argos is always WGS 84
coordinates(movedata) <- c("location.long", "location.lat")
proj4string(movedata) <- CRS("+proj=longlat +datum=WGS84")

#transform lat/long to UTM and rename the columns for ease down the road
#UTM is in meters, but is based on zones - the zone of the region in Netherlands is 31N or 32N (western Netherlands)
utm <- spTransform(movedata, CRSobj="+proj=utm +zone=31 +datum=WGS84")
utm_locs <- data.frame(as(utm, "SpatialPoints"))
colnames(utm_locs) <- c("UTM_W", "UTM_N")

#inspect
head(movedata)

movedata$timestamp <- substr(movedata$timestamp, 1, nchar(as.character(movedata$timestamp))-4)

#there are two mistakes in the ID column, this corrects these
movedata$tag.local.identifier[movedata$tag.local.identifier == "151862"] <- "151682"
movedata$tag.local.identifier[movedata$tag.local.identifier == "157539"] <- "157535"

#define date and time
date_time <- as.POSIXct(strptime(as.character(movedata$timestamp),"%Y-%m-%d %H:%M:%S", tz="GMT"))
fix_date <- as.Date(format(date_time, "%Y-%m-%d"))
month <- as.numeric(format(date_time, "%m"))
year <- as.numeric(format(date_time, "%Y")) #only saves the unique digits of the year so 2015 looks like 0015
day <- as.numeric(format(date_time, "%d"))
time <- strftime(date_time, format="%H:%M:%S", tz="GMT")

#create new data and improved frame for HR analysis
final_move1 <- data.frame(ID=as.factor(movedata$tag.local.identifier),Height=movedata$argos.altitude,UTM_W=utm_locs$UTM_W,UTM_N=utm_locs$UTM_N,fix_date,month,year,day,time,timestamp=date_time)
final_move <- na.omit(final_move1)

#inspect
head(final_move)

#define the coordinates and projection of the new df
coordinates(final_move) <- c("UTM_W", "UTM_N")
projection(final_move) <- CRS("+proj=utm +zone=31 +datum=WGS84")




#### Incorporate ID information and Subset to reproductive cycles ####
#load the individual information of tagged birds
individual <- read.csv("./Haanmer_Argos_IndividualInformation.csv"); head(individual)
nest <- read.csv("./Haanmer_Argos_NestInformation.csv");head(nest)

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

# need to place default value for NA's?

nesters_final_move$fix_success <- if(nesters_final_move$UTM_W == "NA") {0} else {1}

  
#####  Quick summary statistics of Argos Data ##### 
hour <- hour(nesters_final_move$timestamp) #create temporary object by hour function in lubridate
new <- as.data.frame(table(hour)) # make dataframe of frequency table
ggplot() + geom_bar(aes(hour), fill = "red", color = "red", alpha = 0.1) + theme_classic() + lims(x = c(0,24)) #visualize

count(hour)

#we see that the vast majority of fixes occur between 4-11 and 15-20, when the sattelites pass over Haanmeer

# calculate the average number of fixes per tag in each day
nesters_final_move$bird_day <- paste(nesters_final_move$ID, "-", nesters_final_move$fix_date)
bird_days <- as.data.frame(table(nesters_final_move$bird_day))

#visualize
ggplot() + geom_density(aes(bird_days$Freq), fill = "red", color = "red", alpha = 0.1) + geom_vline(xintercept = mean(bird_days$Freq), size = 1.2) + theme_classic() + lims(x = c(0,20))

#what percentage are lost with minimum of 5 fixes
quantile(bird_days$Freq, probs = seq(0,.99,0.01)) #25%

#write.csv(nesters_final_move, "C:/Users/14064/Dropbox/BTGO Movmement Study/Processed Data/Haanmeer_Argos_NestIndiv.csv")

#### convert to nestR format ####

head(nesters_final_move)


coordinates(nesters_final_move) <- c("UTM_W", "UTM_W")
proj4string(nesters_final_move) <- CRS("+proj=utm +zone=31 +datum=WGS84")

#transform lat/long to UTM and rename the columns for ease down the road
#UTM is in meters, but is based on zones - the zone of the region in Netherlands is 31N or 32N (western Netherlands)
utm <- spTransform(nesters_final_move, CRSobj="+proj=longlat +datum=WGS84")
utm_locs <- data.frame(as(utm, "SpatialPoints"))
colnames(utm_locs) <- c("long", "lat")

burst <- paste(nesters_final_move$ID, nesters_final_move$year, sep = "-")

reform_argos_nestR <- data.frame("burst" = burst, "long" = utm_locs$long, "lat" = utm_locs$lat, date = nesters_final_move$timestamp); head(reform_argos_nestR)

write.csv(reform_argos_nestR, "C:/Users/14064/Dropbox/BTGO Movmement Study/BTGO git/bltg-movement/CSV/Reformatted_Argos_Haanmeer_v1.csv")






##### Filter errors and smooth ##### 
 # see https://jmlondon.github.io/crawl-workshop/crawl-practical.html#analysis-and-coding-principles
# see https://movementecologyjournal-biomedcentral-com.pallas2.tcl.sc.edu/articles/10.1186/s40462-017-0104-2#Sec14

#this code will resample to a consistent resolution - in our case, this should be based upon either the GPS data or the fixes that Argos get during the day. 3-hr intervals suffice

