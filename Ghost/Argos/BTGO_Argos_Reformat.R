# Date: 11/3/2020
# Author: Matt Duggan

# Script for processing the Argos data from Verhoeven et al. 2020 and formatting for nestR

#project: inferring breeding success from tags
#species: Black-tailed Godwit (BTGO, Limosa limosa)

#
#begin script
#

#### Load the packages ####

#clean old objects from R
rm(list=ls(all=TRUE))

#see outputs in standard (not scientific notation)
options(scipen = 999)

#location of the working directory (input data)
setwd("bltg-movement")

#read in the data
movedata <- read.csv("./CSV/Argos_BTGO_Haanmeer.csv"); head(movedata)

#Reformat to nestR
library(dplyr)

library(nestR)

vignette("nestR")
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")

library(stats)



movedata <- filter(movedata, movedata$location.lat>= 52)
movedata <- filter(movedata, movedata$location.long >= 0)

#check
min(movedata$location.lat)


#take individuals identificatio, timestamp, long, and lat 
movedata_nestR <- select(movedata, tag.local.identifier, timestamp, location.long, location.lat)

#rename columns to standard format
movedata_nestR <- rename(movedata_nestR, c("burst" = "tag.local.identifier", "date" = "timestamp", "long" = "location.long", "lat" = "location.lat"))

#add year to the burst column
movedata_nestR$burst <- paste(movedata$tag.local.identifier, substring(movedata$timestamp, 1, 4))

#replace " " with "-"
movedata_nestR$burst <- gsub(" ", "-", movedata_nestR$burst);head(movedata_nestR)

#Put date into POSIXct format
movedata_nestR$date <- as.POSIXct(as.character(movedata_nestR$date), tz="GMT"); head(movedata_nestR)

#Run nestR
BTGO_Argos_output <- find_nests(gps_data = movedata_nestR,
                            sea_start = "03-01",
                            sea_end = "07-31",
                            nest_cycle = 27,
                            buffer = 20,
                            min_pts = 2,
                            min_d_fix = 5,
                            min_consec = 2,
                            min_top_att = 1,
                            min_days_att = 1,
                            discard_overlapping = FALSE)

#Checking first individual and nest
BTGO_Argos_output$nests %>% filter(burst == "123424-2016") %>% head()

#Take the lat and long coordinates from the nest prediction
coords_cand <- BTGO_Argos_output$nests %>% filter(burst == "123424-2016") %>% slice(1) %>% select(long, lat)

#Need to grab the 31U and UTM nest locations for each individual
individual <- read.csv("./CSV/Argos_Individual_Information_BTGO_Haanmeer.csv"); head(individual)
nest <- read.csv("./CSV/Haanmer_Argos_NestInformation.csv");head(nest)

#Change Capture_NestID to match in both documents
individual <- rename(individual, c("NestID" = "Capture_NestID")); head(individual)

#merge the data based on NestID
nest_individual <- merge(nest, individual, by = "NestID");head(nest_individual)

#filter out only the ringnumber of the individual, NestID and coordinates
nestForEachIndividual <- select(nest_individual,Ringnumber, NestID, X31U, UTM)

nestCoords <- select(nestForEachIndividual, X31U, UTM)
nestCoords$X31U <- as.numeric(nestCoords$X31U)
nestCoords$UTM <- as.numeric(nestCoords$UTM)
complete_coords <- na.omit(nestCoords)

#Make UTM and 31U into lat and long
library(rgdal)

sputm <- SpatialPoints(complete_coords, proj4string=CRS("+proj=utm +zone=31 ellps=bessel"))
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))

lat_long_locs <- data.frame(as(spgeo , "SpatialPoints"))

nestForEachIndividualLATLONG <- data.frame("Ringnumber" = 0, "NestID" = 0, "lat" = 0, "long" = 0)

copy <- na.omit(nestForEachIndividual)

lat_long_locs$RingNumber <- copy$Ringnumber
lat_long_locs$NestID <- copy$NestID

lat_long_locs <- lat_long_locs %>% select(RingNumber, NestID, everything())
lat_long_locs <- rename(lat_long_locs, c("long" = "X31U", "lat" = "UTM")); head(lat_long_locs)

coords_known <- lat_long_locs %>% filter(NestID == "1606905605") %>% select(long, lat)

geosphere::distGeo(coords_cand, coords_known)
