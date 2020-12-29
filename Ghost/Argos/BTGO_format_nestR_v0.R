# Date: 10/15/20
# Author: Luke Wilde

# Script for processing the Argos data from Verhoeven et al. 2020 and calculating summary statistics on kernel area

#project: inferring breeding success from tags
#species: Black-tailed Godwit (BTGO, Limosa limosa)

#
#begin script from nestR
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
foo( c("dplyr", "sp", "ggplot2", "geosphere", "adehabitatHR", "rgdal", "ggmap", "mapproj", "maptools", "rgeos", "raster", "leaflet", "move", "ctmm", "BBMM", "PBSmapping", "shiny", "plyr", "lubridate", "viridisLite", "rworldmap", "adehabitatLT", "stringr", "nestR", "coda", "rjags"))


#### Load and format ARGOS data ####

#location of the working directory (input data)
setwd("C:/Users/ralph/OneDrive/Desktop/Birds/bltg-movement/CSV/")


#
## IMPORTANT :::: Always remove old objects from R before rerunning this script for formatting data ##
#


#clean old objects from R
rm(list=ls(all=TRUE))

#see outputs in standard (not scientific notation)
options(scipen = 999)
#options(max.print=999999) #if you need to see the complete output
#


#read in the data
movedata <- read.csv("./Argos_BTGO_Haanmeer.csv"); head(movedata)


#reference for ARGOS fix quality:
#class 3 (error < 150m) = 24020 | 26%
#class 2 (150 m < error < 350 m) = 11876 | 13%
#class 1 (350 m < error < 1000 m) = 8424 | 9%
#class 0 (error > 1000) = 6081 | 7%
#class A (no estimate) = 13593 | 15%
#class B (no estimate) = 27688 | 30%
#class Z (invalid fix, failed) = 94 | < 1%


#Verhoeven et al. 2020 defines breeding area for Haanmeer to be > 52N, visually, all breeders in our study ought to have stayed E of the prime meridian (0 longitude)
#subset the data to include only the breeding location
#movedata.breed <- filter(movedata, location.lat >= 52)
movedata.breed <- filter(movedata, location.long >= 0)
movedata.breed <- filter(movedata.breed, argos.lc == "A" | argos.lc == 1 | argos.lc == 2 | argos.lc == 3)

#check
min(movedata.breed$location.lat)
unique(movedata.breed$tag.local.identifier)


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

#### organize the data - formatting and projection ####
movedata <- movedata.breed[with(movedata.breed,order(movedata.breed$tag.local.identifier,movedata.breed$timestamp)),]

str(movedata)

#define coordinates within the data - Argos is always WGS 84
coordinates(movedata) <- c("location.long", "location.lat")
proj4string(movedata) <- CRS("+proj=longlat +datum=WGS84")

#transform lat/long to UTM and rename the columns for ease down the road
#UTM is in meters, but is based on zones - the zone of the region in Netherlands is 31N or 32N (western Netherlands)
utm <- spTransform(movedata, CRSobj="+proj=longlat +zone=31 +datum=WGS84")
utm_locs <- data.frame(as(utm, "SpatialPoints"))
colnames(utm_locs) <- c("UTM_W", "UTM_N")

#inspect
head(movedata)

movedata$timestamp <- substr(movedata$timestamp, 1, nchar(as.character(movedata$timestamp))-4)

head(movedata)


#
#for Mo's data only!
#
#there are two mistakes in the ID column, this corrects these
movedata$tag.local.identifier[movedata$tag.local.identifier == "151862"] <- "151682"
movedata$tag.local.identifier[movedata$tag.local.identifier == "157539"] <- "157535"

head(movedata)

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

#prep for nestR
burst <- paste(final_move$ID, final_move$year, sep = "-")
long <- final_move$UTM_W
lat <- final_move$UTM_N
date <- as.POSIXct(final_move$timestamp)

BTGO_Haanmeer_nestR <- data.frame("burst" = burst,"long" =long,"lat" = lat,"date" =as.POSIXct(date))

#### Begin nestR #####

#open vignette for workflow
vignette("nestR")

#locate nests

#Parameters from Mo's 2020 paper (https://www.frontiersin.org/articles/10.3389/fevo.2019.00031/full) breeding is Mar to July, average incubation is 24 days after last egg laid (which takes 3 to 4 days)
# we will classify a complete breeding cycle as 27 d (+/- 3.4 d)

#ARGOS tags start with 400m buffer, UvA GPS with 20m buffer

#UvA  GPS trackers have error of < 7m
#data were recorded every 5 min when fully charged and every 15-30 minutes when not
#Senner et al. 2018 scrubbed data for only times where interval did not exceed 30 minutes
# since they are solar powered, fixes only happen in the day so that the max fixes were 64.

#randomly select two birds to test with the find_nests function
#rand <- sample(levels(as.factor(BTGO_Haanmeer_nestR$burst)), size = 2)
#godwit.test <- BTGO_Haanmeer_nestR %>% filter(burst == rand[1] | burst == rand[2])
#head(godwit.test)

BTGO_Haanmeer_output_1 <- find_nests(gps_data = BTGO_Haanmeer_nestR,
                            sea_start = "03-01",
                            sea_end = "07-31",
                            nest_cycle = 27,
                            buffer = 150,
                            min_pts = 2,
                            min_d_fix = 5,
                            min_consec = 2,
                            min_top_att = 1,
                            min_days_att = 1,
                            discard_overlapping = FALSE)

# Prior information is not typically available for min_consec, min_top_att and min_days_att - and so the next section focuses on refining these
#all in all, the user should specify low values for the constraints so as to identify all recurrently visited locations


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
individual <- read.csv("./Haanmer_Argos_IndividualInformation.csv"); head(individual)
nest <- read.csv("./Haanmer_Argos_NestInformation.csv");head(nest)

#rename datasets to have consistent names for: tags, nest, and UTM coordinates
names(individual)[4] <- "ID"; head(individual); head(final_move1)
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
BTGO_Haanmeer_output_1$nests
BTGO_Haanmeer_output_1$visits
table(BTGO_Haanmeer_output_1$nests$burst)


# write a for loop for this - after practice!!!!

output1 <- BTGO_Haanmeer_output_1$nests %>%
  dplyr::filter(burst == "144420-2015")

known1 <- BTGO_known_nests1 %>%
  dplyr::filter(burst == "144420-2015")


# 
(explodata <- get_explodata(candidate_nests = BTGO_Haanmeer_output_1$nests, 
                                  known_coords = BTGO_known_nests1,
                                  buffer = 150,
                                  pick_overlapping = TRUE))

#then in the case of unknown nests, we could rely on field data or (when a colony is known), the GUI tool in nestR - not useful for BTGO

#explore_nests(). Here you could create an object and input that as the explodata:
#id_known <- data.frame(burst = "1134370-2013",
#                       loc_id = 2170)

#### Find set of parameter values to tell apart nests and non-nests - CART ####

#
#
# IMPORTANT: For this, would be good to figure out a way to iterate and extract the FPR and FNR for each. Then, statistically determine the best parameters
# As is, this only runs with 1 iteration, but the outcome changes with how the data are split

(BTGO_cart <- discriminate_nests(explodata = explodata, train_frac = .5))

#record the options
#min_consec = 12 days #(FPR = 0, FNR = 0.33)
#min_top_att = 52 #(FPR = 0.06, FNR = 0.17)


#### Identifying nests among revisited locations ####

#return to the find_nests function
#the only thing to change are the new parameters

BTGO_Haanmeer_output_2 <- find_nests(gps_data = BTGO_Haanmeer_nestR,
                                     sea_start = "03-01",
                                     sea_end = "07-31",
                                     nest_cycle = 27,
                                     buffer = 150,
                                     min_pts = 2,
                                     min_d_fix = 5,
                                     min_consec = 10,
                                     min_top_att = 1,
                                     min_days_att = 1,
                                     discard_overlapping = FALSE)

#review outputs
BTGO_Haanmeer_output_2$nests
BTGO_Haanmeer_output_2$visits
table(BTGO_Haanmeer_output_2$nests$burst)

# calculate the distance of each nest with a known location to the candidate nest



#store this second round for estimating the reproductive outcome
BTGO_nests <- BTGO_Haanmeer_output_2

#### Estimating reproductive outcome ####

#this stage now estimates if an attempt was successful or not depending on whether it lasted as long as the duration of a complete nesting cycle
#the underlying assumption is that the nest stops being revisited after an attempt fails, which is true for many bird species
# THIS IS A PROBLEM FOR ARGOS TAGS, WHICH COULD MISS MANY REVISTATION EVENTS. For these we need to update this to calculate Pr(Visitation) from movement

#Functions in nestR allow the estimation of reproductive outcome from GPS-derived nest revisitation histories while taking into account the probability of visit detection and allowing both detection and nest survival probability to vary in time.

BTGO_attempts <- format_attempts(nest_info = BTGO_nests, nest_cycle = 27)

#now we specify the nest survival model

#The function estimate_outcomes() fits a Bayesian hierarchical model to the histories of nest revisitation and estimates the probability that each nest is still active (i.e., "alive") on the last day of the attempt. Under the hood, estimate_outcomes() uses JAGS (Just Another Gibbs Sampler) to run the MCMC (Markov Chain Monte Carlo), via the package rjags. This gives the user the option of accessing the MCMC diagnostics from the coda package. The user can control parameters of the MCMC by passing them to the argument mcmc_params.

# the user can choose among four models depending on the literature
  # null = a null model, where both survival and detection probability are constant,
  # phi_time = a model where survival varies through time,
  # p_time = a model where detection probability varies through time, and
  # phi_time_p_time = a model where both survival and detection probability are a       function of time.

BTGO_outcomes <- estimate_outcomes(fixes = BTGO_attempts$fixes, 
                                   visits = BTGO_attempts$visits,
                                   model = "phi_time_p_time")

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
plot_survival(BTGO_outcomes)
plot_detection(BTGO_outcomes)

#see individual survival
plot_nest_surv(BTGO_outcomes, who = 15)

# output of Bayes survival model
summarize_outcomes(BTGO_outcomes, ci = .95)

#check mixing and sampling distributino of the intercept and slope parameters (pb0 and pb1 resp.)
BTGO_pb0_coda <- coda::as.mcmc.list(BTGO_outcomes$p.b0)
BTGO_pb1_coda <- coda::as.mcmc.list(BTGO_outcomes$p.b1)

#inspect chains for mixing - a fuzzy caterpillar with a straight, flat line is good.
plot(BTGO_pb0_coda)
plot(BTGO_pb1_coda)

# END #

#### Apply model to tags without known nest locations ####

setwd("C:/Users/14064/Dropbox/BTGO Movmement Study/Data/Data from Senner 2015/")


#
## IMPORTANT :::: Always remove old objects from R before rerunning this script for formatting data ##
#


#clean old objects from R
rm(list=ls(all=TRUE))

#see outputs in standard (not scientific notation)
options(scipen = 999)
#options(max.print=999999) #if you need to see the complete output
#


#read in the data
movedata <- read.csv("./Continental black-tailed godwits (data from Senner et al. 2015).csv"); head(movedata)


#reference for ARGOS fix quality:
#class 3 (error < 150m) = 24020 | 26%
#class 2 (150 m < error < 350 m) = 11876 | 13%
#class 1 (350 m < error < 1000 m) = 8424 | 9%
#class 0 (error > 1000) = 6081 | 7%
#class A (no estimate) = 13593 | 15%
#class B (no estimate) = 27688 | 30%
#class Z (invalid fix, failed) = 94 | < 1%


#Verhoeven et al. 2020 defines breeding area for Haanmeer to be > 52N, visually, all breeders in our study ought to have stayed E of the prime meridian (0 longitude)
#subset the data to include only the breeding location
#movedata.breed <- filter(movedata, location.lat >= 52)
movedata.breed <- filter(movedata, location.long >= 0)
movedata.breed <- filter(movedata.breed, argos.lc == "A" | argos.lc == 1 | argos.lc == 2 | argos.lc == 3)

#check
min(movedata.breed$location.lat)
unique(movedata.breed$tag.local.identifier)


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

#### organize the data - formatting and projection ####
movedata <- movedata.breed[with(movedata.breed,order(movedata.breed$tag.local.identifier,movedata.breed$timestamp)),]

str(movedata)

#define coordinates within the data - Argos is always WGS 84
coordinates(movedata) <- c("location.long", "location.lat")
proj4string(movedata) <- CRS("+proj=longlat +datum=WGS84")

#transform lat/long to UTM and rename the columns for ease down the road
#UTM is in meters, but is based on zones - the zone of the region in Netherlands is 31N or 32N (western Netherlands)
utm <- spTransform(movedata, CRSobj="+proj=longlat +zone=31 +datum=WGS84")
utm_locs <- data.frame(as(utm, "SpatialPoints"))
colnames(utm_locs) <- c("UTM_W", "UTM_N")

#inspect
head(movedata)

movedata$timestamp <- substr(movedata$timestamp, 1, nchar(as.character(movedata$timestamp))-4)

head(movedata)


#
#for Mo's data only!
#
#there are two mistakes in the ID column, this corrects these
movedata$tag.local.identifier[movedata$tag.local.identifier == "151862"] <- "151682"
movedata$tag.local.identifier[movedata$tag.local.identifier == "157539"] <- "157535"

head(movedata)

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

#prep for nestR
burst <- paste(final_move$ID, final_move$year, sep = "-")
long <- final_move$UTM_W
lat <- final_move$UTM_N
date <- as.POSIXct(final_move$timestamp)

#complete data
BTGO_elsewhere_nestR <- data.frame("burst" = burst,"long" =long,"lat" = lat,"date" =as.POSIXct(date))

#specify recursive model
BTGO_Haanmeer_output_2 <- find_nests(gps_data = BTGO_elsewhere_nestR,
                                     sea_start = "03-01",
                                     sea_end = "07-31",
                                     nest_cycle = 27,
                                     buffer = 150,
                                     min_pts = 2,
                                     min_d_fix = 5,
                                     min_consec = 10,
                                     min_top_att = 1,
                                     min_days_att = 1,
                                     discard_overlapping = FALSE)

#review outputs
BTGO_Haanmeer_output_2$nests
BTGO_Haanmeer_output_2$visits
table(BTGO_Haanmeer_output_2$nests$burst)

# calculate the distance of each nest with a known location to the candidate nest



#store this second round for estimating the reproductive outcome
BTGO_nests <- BTGO_Haanmeer_output_2

#### Estimating reproductive outcome ####

#this stage now estimates if an attempt was successful or not depending on whether it lasted as long as the duration of a complete nesting cycle
#the underlying assumption is that the nest stops being revisited after an attempt fails, which is true for many bird species
# THIS IS A PROBLEM FOR ARGOS TAGS, WHICH COULD MISS MANY REVISTATION EVENTS. For these we need to update this to calculate Pr(Visitation) from movement

#Functions in nestR allow the estimation of reproductive outcome from GPS-derived nest revisitation histories while taking into account the probability of visit detection and allowing both detection and nest survival probability to vary in time.

BTGO_attempts <- format_attempts(nest_info = BTGO_nests, nest_cycle = 27)

#now we specify the nest survival model

#The function estimate_outcomes() fits a Bayesian hierarchical model to the histories of nest revisitation and estimates the probability that each nest is still active (i.e., "alive") on the last day of the attempt. Under the hood, estimate_outcomes() uses JAGS (Just Another Gibbs Sampler) to run the MCMC (Markov Chain Monte Carlo), via the package rjags. This gives the user the option of accessing the MCMC diagnostics from the coda package. The user can control parameters of the MCMC by passing them to the argument mcmc_params.

# the user can choose among four models depending on the literature
# null = a null model, where both survival and detection probability are constant,
# phi_time = a model where survival varies through time,
# p_time = a model where detection probability varies through time, and
# phi_time_p_time = a model where both survival and detection probability are a       function of time.

BTGO_outcomes <- estimate_outcomes(fixes = BTGO_attempts$fixes, 
                                   visits = BTGO_attempts$visits,
                                   model = "phi_time_p_time")

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
plot_survival(BTGO_outcomes)
plot_detection(BTGO_outcomes)

#see individual survival
plot_nest_surv(BTGO_outcomes, who = 1)

# output of Bayes survival model
summarize_outcomes(BTGO_outcomes, ci = .95)

#check mixing and sampling distributino of the intercept and slope parameters (pb0 and pb1 resp.)
BTGO_pb0_coda <- coda::as.mcmc.list(BTGO_outcomes$p.b0)
BTGO_pb1_coda <- coda::as.mcmc.list(BTGO_outcomes$p.b1)

#inspect chains for mixing - a fuzzy caterpillar with a straight, flat line is good.
plot(BTGO_pb0_coda)
plot(BTGO_pb1_coda)

# END #


