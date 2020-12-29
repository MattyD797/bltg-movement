#Author: Matt Duggan
#Date: 12.21.2020

#Script to filter Argos locations for errors and run summary statistics with nestR

#clean old objects from R
rm(list=ls(all=TRUE))

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
foo( c("nestR", "dplyr", "sp", "ggplot2", "geosphere", "adehabitatHR", "rgdal", "ggmap", "mapproj", "maptools", "rgeos", "raster", "leaflet", "move", "ctmm", "BBMM", "PBSmapping", "shiny", "plyr", "lubridate", "viridisLite", "rworldmap", "adehabitatLT", "adehabitatMA", "SDLfilter"))



#location of the working directory (input data)
setwd("./CSV/")

#### Load and format the movement data ####

#read in the data
movedata <- read.csv("./Argos_BTGO_Haanmeer.csv"); head(movedata)


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

#breakdown of fix quality:
#class 3 (error < 150m) = 24020 | 26%
#class 2 (150 m < error < 350 m) = 11876 | 13%
#class 1 (350 m < error < 1000 m) = 8424 | 9%
#class 0 (error > 1000) = 6081 | 7%
#class A (no estimate) = 13593 | 15%
#class B (no estimate) = 27688 | 30%
#class Z (invalid fix, failed) = 94 | < 1%
movedata.breed <- filter(movedata.breed, argos.lc == "A" | argos.lc == 1 | argos.lc == 2 | argos.lc == 3)

#clean the data and put in correct format
#filter
data <- movedata[,c(3,4,5,13,27)]; head(data)
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

#Transform into dataframe of burst, long, lat, and date. 
godwit <- BTGO_dup.rem_qi[, c(5, 2, 3, 1)]
colnames(godwit)[1] <- "burst"
colnames(godwit)[2] <- "long"
colnames(godwit)[3] <- "lat"
colnames(godwit)[4] <- "date"

#Organize data based on burst Id and time
godwit <- godwit[with(godwit,order(godwit$burst,godwit$date)),]

str(godwit)

#Verhoeven et al. 2020 defines breeding area for Haanmeer to be > 52N, visually, all breeders in our study ought to have stayed E of the prime meridian (0 longitude)
#subset the data to include only the breeding location
godwit <- filter(godwit, lat >= 52)
godwit <- filter(godwit, long >= 0)


#Prepare for nestR
long <- godwit$long
lat <-  godwit$lat
godwit$date <- substr(godwit$date, 1, nchar(as.character(godwit$date)))
date <- as.POSIXct(godwit$date, tx = "GMT", format = "%Y-%m-%d %H:%M:%S")
godwit <- na.omit(godwit)
year <- as.numeric(format(date, "%Y"))
burst <- paste(as.character(godwit$burst),as.character(year),sep="-")

godwit<- data.frame("burst" = burst,"date" =date,"long" =long,"lat" = lat)
godwit <- na.omit(godwit)

vignette("nestR")

#Parameters from Mo's 2020 paper (https://www.frontiersin.org/articles/10.3389/fevo.2019.00031/full) breeding is Mar to July, average incubation is 24 days after last egg laid (which takes 3 to 4 days)
# we will classify a complete breeding cycle as 27 d (+/- 3.4 d)

#ARGOS tags start with 400m buffer, UvA GPS with 20m buffer

#UvA  GPS trackers have error of < 7m
#data were recorded every 5 min when fully charged and every 15-30 minutes when not
#Senner et al. 2018 scrubbed data for only times where interval did not exceed 30 minutes
# since they are solar powered, fixes only happen in the day so that the max fixes were 64.

BTGO_Argos_Haanmeer_output_1 <- find_nests(gps_data = godwit,
                                           sea_start = "03-01",
                                           sea_end = "07-31",
                                           nest_cycle = 27,
                                           buffer = 400,
                                           min_pts = 2,
                                           min_d_fix = 5,
                                           min_consec = 2,
                                           min_top_att = 1,
                                           min_days_att = 1,
                                           discard_overlapping = FALSE)


#### Format known nest information ####   

#here we deterimine the best set of parameters to discriminate nests from non-nests
#specifically, we want to tune min_consec, min_top_att and min_days_att.
#

#
#to inform our choice, we compare the values obtained in output for nests and non-nests for consec_days, perc_top_vis, and perc_days_vis

#get_explodata() automates the process of selecting nests and non-nests to compare using the output above and known nest locations from the user
#choose the most visited location (#1 in output nest list)

#if true nest locations are known, you load in the known location, and then use:
#load the individual information of tagged birds
individual <- read.csv("./Haanmer_Argos_Individual_Information.csv"); head(individual)
nest <- read.csv("./Haanmer_Argos_Nest_Information.csv");head(nest)

#rename datasets to have consistent names for: tags, nest, and UTM coordinates
names(individual)[4] <- "ID"; head(individual); head(nest)
names(individual)[6] <- "NestID"; head(individual); head(nest)
names(nest)[c(3,4)] <- c("nestloc_UTM_W", "nestloc_UTM_N"); head(nest)

#now, join the nest and individual datasets first
nesters <- right_join(nest, individual, by = "NestID"); head(nesters)

#change ID to factor
nesters$ID <- as.factor(nesters$ID)


#reformat the nesting dates to numerical
nesters$StartKnownActive <- gsub('Apr','04', nesters$StartKnownActive)
nesters$StartKnownActive <- gsub('May','05', nesters$StartKnownActive)
nesters$StartKnownActive <- gsub('Jun','06', nesters$StartKnownActive)
nesters$EndKnownActive <- gsub('Apr','04', nesters$EndKnownActive)
nesters$EndKnownActive <- gsub('May','05', nesters$EndKnownActive)
nesters$EndKnownActive <- gsub('Jun','06', nesters$EndKnownActive)

#then, reformat the dates as date objects to allow easy subsetting
nesters$nest_start_date <- as.POSIXct(strptime(as.character(nesters$StartKnownActive),"%d-%m-%Y", tz="GMT"))
nesters$nest_end_date <- as.POSIXct(strptime(as.character(nesters$EndKnownActive),"%d-%m-%Y", tz="GMT"))

nesters$nest_start_date <- as.Date(strptime(nesters$nest_start_date,"%Y-%m-%d"))
nesters$nest_end_date <- as.Date(strptime(nesters$nest_end_date,"%Y-%m-%d"))

#I havent found an easier way to do this, so I simply ask R to substitue the incorrectly formatted years with ones that agree with other formats in the code
nesters$nest_start_date <- as.Date(gsub('0016','2016',nesters$nest_start_date))
nesters$nest_start_date <- as.Date(gsub('0015','2015', nesters$nest_start_date))
nesters$nest_start_date <- as.Date(gsub('0017','2017', nesters$nest_start_date))
nesters$nest_end_date <- as.Date(gsub('0016','2016',nesters$nest_end_date))
nesters$nest_end_date <- as.Date(gsub('0015','2015', nesters$nest_end_date))
nesters$nest_end_date <- as.Date(gsub('0017','2017', nesters$nest_end_date))

#estimate the days a nest survived
nesters$nest_min_days_survived <- as.numeric(nesters$nest_end_date - nesters$nest_start_date); hist(nesters$nest_min_days_survived)

#break each down into their own columns
nesters$nest_start_month <- as.numeric(format(nesters$nest_start_date, "%m"))
nesters$nest_start_year <- as.numeric(format(nesters$nest_start_date, "%Y")) #only saves the unique digits of the year so 2015 looks like 0015
nesters$nest_start_day <- as.numeric(format(nesters$nest_start_date, "%d"))

nesters$nest_end_month <- as.numeric(format(nesters$nest_end_date, "%m"))
nesters$nest_end_year <- as.numeric(format(nesters$nest_end_date, "%Y")) #only saves the unique digits of the year so 2015 looks like 0015
nesters$nest_end_day <- as.numeric(format(nesters$nest_end_date, "%d"))

#erase second Colourcode column and rename the remaining column
nesters$Colourcode.y <- NULL
names(nesters)[1] <- "Colourcode"

head(nesters)

##
##

#filter only burst with known nest location - for training pusposes
BTGO_known_nests<-nesters[!is.na(nesters$nestloc_UTM_W),]

#define coords and projection - currently in UTM
coordinates(BTGO_known_nests) <- c("nestloc_UTM_W", "nestloc_UTM_N")
proj4string(BTGO_known_nests) <- CRS("+proj=utm +zone=31 +datum=WGS84")

#transform to latlong
locs <- spTransform(BTGO_known_nests, CRSobj="+proj=longlat +zone=31 +datum=WGS84")
BTGO_known_nests_locs <- data.frame(as(locs, "SpatialPoints"))
colnames(BTGO_known_nests_locs) <- c("long", "lat")

known_nests_burst <- paste(BTGO_known_nests$ID, BTGO_known_nests$YearDeployed, sep = "-")

BTGO_known_nests1 <- data.frame("burst" = as.character(known_nests_burst), "fate" = BTGO_known_nests$NestSuccess, "start_known" = BTGO_known_nests$nest_start_date, "end_known" = BTGO_known_nests$nest_end_date, "long" = as.numeric(BTGO_known_nests_locs$long), "lat" = as.numeric(BTGO_known_nests_locs$lat)); head(BTGO_known_nests1)

#### Discriminating between nests and non-nests ####  

# first, in the case of known nests

#the output of find_nests() is a list of nests and visits
BTGO_Argos_Haanmeer_output_1$nests
BTGO_Argos_Haanmeer_output_1$visits
table(BTGO_Argos_Haanmeer_output_1$nests$burst)

# Looking at the disstance of the first coordinate from the known

output1 <- BTGO_Argos_Haanmeer_output_1$nests %>%
  dplyr::filter(burst == "144420-2015") %>%
  slice(1)

output1 <- output1[,c(3,4)]

known1 <- BTGO_known_nests1 %>%
  dplyr::filter(burst == "144420-2015") 


known1_coords <- known1[,c(5,6)]
geosphere::distGeo(output1, known1_coords)
# 

#Find nest versus not nests

(explodata <- get_explodata(candidate_nests = BTGO_Argos_Haanmeer_output_1$nests, 
                            known_coords = BTGO_known_nests1,
                            buffer = 500,
                            pick_overlapping = TRUE))


(BTGO_cart <- discriminate_nests(explodata = explodata, train_frac = .5))

BTGO_Haanmeer_output_2 <- find_nests(gps_data = godwit,
                                     sea_start = "03-01",
                                     sea_end = "07-31",
                                     nest_cycle = 27,
                                     buffer = 400,
                                     min_pts = 2,
                                     min_d_fix = 5,
                                     min_consec = 8,
                                     min_top_att = 1,
                                     min_days_att = 27,
                                     discard_overlapping = TRUE)


output2 <- BTGO_Haanmeer_output_2$nests %>%
  dplyr::filter(burst == "123424-2016") %>%
  slice(1)

output2 <- output2[,c(3,4)]; head(output2)

known2 <- BTGO_known_nests1 %>%
  dplyr::filter(burst == "123424-2016") 


known2_coords <- known2[,c(5,6)]
geosphere::distGeo(output2, known2_coords)


godwit_nests <- BTGO_Haanmeer_output_2

#Format the output of find nests into two matrices (fixes and visits)

godwit_attempts <- format_attempts(nest_info = godwit_nests, nest_cycle = 27)

#now we specify the nest survival model

#The function estimate_outcomes() fits a Bayesian hierarchical model to the histories of nest revisitation and estimates the probability that each nest is still active (i.e., "alive") on the last day of the attempt. Under the hood, estimate_outcomes() uses JAGS (Just Another Gibbs Sampler) to run the MCMC (Markov Chain Monte Carlo), via the package rjags. This gives the user the option of accessing the MCMC diagnostics from the coda package. The user can control parameters of the MCMC by passing them to the argument mcmc_params.

# the user can choose among four models depending on the literature
# null = a null model, where both survival and detection probability are constant,
# phi_time = a model where survival varies through time,
# p_time = a model where detection probability varies through time, and
# phi_time_p_time = a model where both survival and detection probability are a function of time.

godwit_outcomes <- estimate_outcomes(fixes = godwit_attempts$fixes, 
                                 visits = godwit_attempts$visits,
                                 model = "null")
str(godwit_outcomes)

#the output is a list of multiple objects of class mcarray
# p = stores the population-level value of p on each day of the attempt.
# p.b0 = is the value of the intercept of p.
# p.b1 = is the value of the slope of p.
# phi = stores the population-level value of ?? on each day of the attempt.
# phi.b0 = is the value of the intercept of ??.
# phi.b1 = is the value of the slope of ??.
# z = stores the daily survival for each attempt.
# names = stores the attempt ID, in the form burst_loc_id.
# model = is the name of the chosen model structure.

#see time series of survival and observation process
plot_survival(godwit_outcomes)
plot_detection(godwit_outcomes)


# output of Bayes survival model
summarize_outcomes(godwit_outcomes, ci = .95)


godwit_pb0_coda <- coda::as.mcmc.list(godwit_outcomes$p.b0)
godwit_phi_coda <- coda::as.mcmc.list(godwit_outcomes$phi.b0)

#inspect chains for mixing - a fuzzy caterpillar with a straight, flat line is good.
plot(godwit_pb0_coda)
plot(godwit_phi_coda)

BTGO_known_nests1
