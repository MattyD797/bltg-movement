# Author: Luke Wilde
# Script to filter down the UvABits_complete data

#========================================
#filters applied:
  #satellites used
  #median from straight line
  #speed
#========================================

#build function to load packages
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

#load all libraries
foo( c("rgdal","sf","sp","raster","ggplot2","maptools","viridis","tidyverse", "remotes", "leaflet", "shiny", "dplyr"))

#clean old objects from R
rm(list=ls(all=TRUE))

#see outputs in standard (not scientific notation)
options(scipen = 999)

#load the raw track
raw_tk_1009 <- read.csv("C:/Users/14064/Dropbox/BTGO Movmement Study/UvABits_complete_data/wetransfer-2ed909/1009.csv", skip = 4)
head(raw_tk_1009)

raw_tk_1009 <- filter(raw_tk_1009, latitude != "\\N")

#check to see if you should get rid of those without known satellite number?
#if > 10% do not use
tmp <- nrow(filter(raw_tk_1009, satellites_used == "\\N"))/nrow(raw_tk_1009); tmp

  if(tmp < 0.10){
    #remove all locations with unknown no. satellites. >4 is typical cut off.
    sat_tk_1009 <- filter(raw_tk_1009, satellites_used != "\\N")
  } else {
    sat_tk_1009 <- raw_tk_1009
  }


#then apply the speed and median filters as below

#speed based filter by Ingo Schiffner 2017
#calculates speed between consecutive localizations and filters out segments exceeding spd_lim
#assumes x and y coordinates are given in a projected format in meters and time(t) given as ATLASTimeStamp(ms)
#returns a data.frame containing filtered x,y and time(t)
SpdFilt <- function(x,y,t,spd_lim) 
{
  dx = diff(x,1)
  dy = diff(y,1)
  dt = diff(t,1)/1000
  de = (dx^2 + dy^2) ^ 0.5
  spd = de/dt
  spd = c(0,spd)
  xr= x[spd<=spd_lim]
  yr= y[spd<=spd_lim]
  tr= t[spd<=spd_lim]
  sdf <- data.frame(x=as.numeric(xr),y=as.numeric(yr),t=as.numeric(tr)) 
  return (sdf)
}

#before filtering, need to change time to ATLASTimeStamp(ms) and project to UTM

coordinates(sat_tk_1009) <- c("longitude", "latitude")
proj4string(sat_tk_1009) <- CRS("+proj=longlat +zone=31 +datum=WGS84")

#transform lat/long to UTM and rename the columns for ease down the road
#UTM is in meters, but is based on zones - the zone of the region in Netherlands is 31N or 32N (western Netherlands)
utm <- spTransform(sat_tk_1009, CRSobj="+proj=utm +zone=31 +datum=WGS84")
utm_locs <- data.frame(as(utm, "SpatialPoints"))
colnames(utm_locs) <- c("UTM_W", "UTM_N")

#DOY <- as.POSIXct(strptime(as.character(data$DateTime), "%Y-%m-%d"))

sat_reproj_tk_1009 <- data.frame(sat_tk_1009[1:ncol(sat_tk_1009)], utm_locs$UTM_W, utm_locs$UTM_N); head(sat_reproj_tk_1009)

sat_reproj_tk_1009$ATLAST <- as.numeric(as.POSIXct(sat_reproj_tk_1009$date_time, format = "%Y-%m-%d %H:%M:%S", tz = "GMT"))


sat_spd_tk_1009 <- SpdFilt(sat_reproj_tk_1009$utm_locs.UTM_W, sat_reproj_tk_1009$utm_locs.UTM_N, sat_reproj_tk_1009$ATLAST,50)




#basic median filter by Ingo Schiffner 2017
#calculates median x y localizations and time (t) over a predefined window (win) in ms
#assumes x and y coordinates are given in a projected format in meters and time(t) given as ATLASTimeStamp(ms)
#returns a data.frame containing filtered x,y and time(t)
MedFilt <- function(x,y,t,win) 
{
  #set bins
  tb = trunc((t - t[1])/win)+1
  mx <- aggregate(x ~ tb, FUN = median)
  my <- aggregate(y ~ tb, FUN = median)
  mt <- aggregate(t ~ tb, FUN = median)
  mtf <- data.frame(x=as.numeric(mx$x),y=as.numeric(my$y),t=as.numeric(mt$t)) 
  return (mtf)
}

clean_tk_1009 <- MedFilt(sat_spd_tk_1009$x, sat_spd_tk_1009$y,sat_spd_tk_1009$t,1500)

#get ATLASTtime back to POSIXct
clean_tk_1009$date_time <- as.POSIXct.numeric(clean_tk_1009$t, origin = "1970-01-01")

head(clean_tk_1009)




#for loop through all individuals
#
for(i in c(1001,1002,1008,2001,2002,2004,2006,2011,2012,2013,2014,2015,2017,2019,2020,2042,1009,2010,2016,2018,2023,2041)){
 
  #load the raw track
  raw_tk <- read.csv(print(paste("C:/Users/14064/Dropbox/BTGO Movmement Study/UvABits_complete_data/wetransfer-2ed909/",i,".csv", sep="")), skip = 4)
  head(raw_tk)

  #only keep rows with complete data  
 raw_tk <- filter(raw_tk, latitude != "\\N")
  
 #since there were missing rows for some, forces coords to remain numeric. They will be character otherwise
 raw_tk$latitude <- as.numeric(raw_tk$latitude)
 raw_tk$longitude <- as.numeric(raw_tk$longitude)
 
  #check to see if you should get rid of those without known satellite number?
  #if > 10% do not use
  tmp <- nrow(filter(raw_tk, satellites_used == "\\N"))/nrow(raw_tk); tmp
  
  if(tmp < 0.05){
    #remove all locations with unknown no. satellites. >4 is typical cut off.
    sat_tk <- filter(raw_tk, satellites_used != "\\N")
  } else {
    sat_tk <- raw_tk
  }
  
  
  #then apply the speed and median filters as below
  
  #speed based filter by Ingo Schiffner 2017
  #calculates speed between consecutive localizations and filters out segments exceeding spd_lim
  #assumes x and y coordinates are given in a projected format in meters and time(t) given as ATLASTimeStamp(ms)
  #returns a data.frame containing filtered x,y and time(t)
  SpdFilt <- function(x,y,t,spd_lim) 
  {
    dx = diff(x,1)
    dy = diff(y,1)
    dt = diff(t,1)/1000
    de = (dx^2 + dy^2) ^ 0.5
    spd = de/dt
    spd = c(0,spd)
    xr= x[spd<=spd_lim]
    yr= y[spd<=spd_lim]
    tr= t[spd<=spd_lim]
    sdf <- data.frame(x=as.numeric(xr),y=as.numeric(yr),t=as.numeric(tr)) 
    return (sdf)
  }
  
  #before filtering, need to change time to ATLASTimeStamp(ms) and project to UTM
  
  coordinates(sat_tk) <- c("longitude", "latitude")
  proj4string(sat_tk) <- CRS("+proj=longlat +zone=31 +datum=WGS84")
  
  #transform lat/long to UTM and rename the columns for ease down the road
  #UTM is in meters, but is based on zones - the zone of the region in Netherlands is 31N or 32N (western Netherlands)
  utm <- spTransform(sat_tk, CRSobj="+proj=utm +zone=31 +datum=WGS84")
  utm_locs <- data.frame(as(utm, "SpatialPoints"))
  colnames(utm_locs) <- c("UTM_W", "UTM_N")
  
  #DOY <- as.POSIXct(strptime(as.character(data$DateTime), "%Y-%m-%d"))
  
  sat_reproj_tk <- data.frame(sat_tk[1:ncol(sat_tk)], utm_locs$UTM_W, utm_locs$UTM_N); head(sat_reproj_tk)
  
  sat_reproj_tk$ATLAST <- as.numeric(as.POSIXct(sat_reproj_tk$date_time, format = "%Y-%m-%d %H:%M:%S", tz = "GMT"))
  
  
  sat_spd_tk <- SpdFilt(sat_reproj_tk$utm_locs.UTM_W, sat_reproj_tk$utm_locs.UTM_N, sat_reproj_tk$ATLAST,150)
  
  
  
  
  #basic median filter by Ingo Schiffner 2017
  #calculates median x y localizations and time (t) over a predefined window (win) in ms
  #assumes x and y coordinates are given in a projected format in meters and time(t) given as ATLASTimeStamp(ms)
  #returns a data.frame containing filtered x,y and time(t)
 MedFilt <- function(x,y,t,win) 
  {
    #set bins
    tb = trunc((t - t[1])/win)+1
    mx <- aggregate(x ~ tb, FUN = median)
    my <- aggregate(y ~ tb, FUN = median)
    mt <- aggregate(t ~ tb, FUN = median)
    mtf <- data.frame(x=as.numeric(mx$x),y=as.numeric(my$y),t=as.numeric(mt$t)) 
    return (mtf)
  }
  
  clean_tk <- sat_spd_tk #MedFilt(sat_spd_tk$x, sat_spd_tk$y,sat_spd_tk$t,1500)
  
  #get ATLASTtime back to POSIXct
  clean_tk$date_time <- as.POSIXct.numeric(clean_tk$t, origin = "1970-01-01")
  
  head(clean_tk) 
  
  write.csv(clean_tk, print(paste("C:/Users/14064/Dropbox/BTGO Movmement Study/UvABits_complete_data/cleaned_tracks/",i,"_clean.csv", sep="")))
  
  
}


1009
2010
2016
2018
2023
2041














## comments below ##










#running median filter by Ingo Schiffner 2017
#calculates median x y localizations and time (t) over a predefined window (win) in ms in predefined time steps (stp)
#assumes x and y coordinates are given in a projected format in meters and time(t) given as ATLASTimeStamp(ms)
#returns a data.frame containing filtered x,y and time(t) and n number of samples used to calculate for each time step
# #RunMedFilt <- function(x,y,t,win,stp) 
# {
#   rmt <- NULL
#   rmx <- NULL
#   rmy <- NULL
#   rmn <- NULL
#   i <- 1
#   while (i < length(t))
#   {
#     #get window
#     nt <- which(t >= (t[i]+med_win))
#     if (!is.na(nt[1])) {
#       nt <- nt[1]
#       #calculate median
#       rmt<-c(rmt,median(t[i:nt]))
#       rmx<-c(rmx,median(x[i:nt]))
#       rmy<-c(rmy,median(y[i:nt]))
#       rmn<-c(rmn,nt-i)
#       #get next step
#       ni <- which(t >= (t[i]+med_stp))
#       i <- ni[1]
#     } else {
#       break
#     }
#   }
#   rmdf <- data.frame(x=as.numeric(rmx),y=as.numeric(rmy),t=as.numeric(rmt),n=as.numeric(rmn)) 
# }




# 
# 
# source('../SQLITE/ATLAS_SQLite_GetLoc.R')
# source('RunMedFilt.R')
# source('MedFilt.R')
# source('SpdFilt.R')
# 
# #parameters
# max_spd <- 50 #maximum speed in m/s (should be adjusted to the maximum speed of your animal)
# med_win <- 30000 #median window in ms (should be adjusted to your sampling rate and the behavioral timeframe)
# med_stp <- 10000 #median step width in ms for running median filter
# 
# #select database file
# sqlfile = choose.files(default = "", caption = "Select SQLite Database", multi = FALSE, filters = c(".sqlite", "*.*"))
# SQLDATA = ATLAS_SQLite_GetLoc(sqlfile)
# TagList = unique(SQLDATA$TAG)
# 
# for(tn in 1:length(TagList))
# {
#   #get data
#   t_dat <- subset(SQLDATA,SQLDATA$TAG==TagList[tn])
#   plot(t_dat$X,t_dat$Y, col='red')
#   
#   #speed filter
#   s_dat = SpdFilt(t_dat$X, t_dat$Y, t_dat$TIME, max_spd)
#   points(s_dat$x,s_dat$y, col='green')
#   
#   #basic median filter
#   m_dat = MedFilt(s_dat$x, s_dat$y, s_dat$t, med_win)
#   lines(m_dat$x,m_dat$y)
#   
#   #runing median filter
#   rm_dat = RunMedFilt(s_dat$x, s_dat$y, s_dat$t, med_win, med_stp)
#   lines(rm_dat$x,rm_dat$y, col='blue')
# }
