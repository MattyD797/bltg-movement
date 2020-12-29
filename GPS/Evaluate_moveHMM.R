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
options(max.print=100) #if you need to see the complete output
#

#read in the data and preview

setwd("./CSV/")
data1 <- read.csv("./GPS_BTGO_Haanmeer_nestR.csv"); head(data1)

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

DOY <- as.POSIXct(strptime(as.character(data$DateTime), "%Y-%m-%d"))

BTGO_move_1 <- data.frame(id = data$id,date = data$DateTime, DOY = DOY, UTM_W = utm_locs$UTM_W, UTM_N = utm_locs$UTM_N); head(BTGO_move_1)

#convert to a traj object
BTGO_move <- as.ltraj(xy = BTGO_move_1[, c("UTM_W", "UTM_N")], date = BTGO_move_1$date, id = BTGO_move_1$id)

loop_object <- vector("list", length(1:25))
for(i in 1:25){loop_object[[i]] <- print(median(BTGO_move[[i]]$dt, na.rm = TRUE))} 
lags <- as.data.frame(do.call(rbind, loop_object))

mean(lags$V1) #this is the typical lag of all tags in seconds: 1490.14

#interpolate between known locations
# BTGO_move2 <- redisltraj(na.omit(BTGO_move),3600,type="time"); BTGO_move2
# 
# #convert to data frame
# data <- ld(BTGO_move2)
# head(data)

for(i in 1:length(BTGO_move)){
  ref <- round(BTGO_move[[i]]$date[1], "hours")
  BTGO_move[i] %>% 
    setNA(ltraj = ., date.ref = ref, dt = 1, units = "hour") %>%
    sett0(ltraj = ., date.ref = ref, dt = 1, units = "hour") -> BTGO_move[i]
}

is.regular(BTGO_move)

data <- ld(BTGO_move)

#keep only columns with id, date, x, y
data <- data[,c(1,2,3,11)]

#convert from m to km
data$x <- data$x/1000 
data$y <- data$y/1000

#id needs to be capitalized so that each animal has its own track
colnames(data)[4] <- "ID"

#create hmm data object
hmmdata <- prepData(data, type="UTM", coordNames = c("x", "y"))

head(hmmdata) #step lengths are in km

#### Visualize parameters ####
#we imagine three states: adult foraging, tending, resting

ggplot(hmmdata, aes(step)) + geom_density(color = "red", fill = "red", alpha = 0.3) + theme_classic() + coord_cartesian(xlim=c(0,10)) 

ggplot(hmmdata, aes(angle)) + geom_density(color = "red", fill = "red", alpha = 0.3) + coord_cartesian(xlim=c(-pi,pi)) + theme_classic() # TA near 0 correspond to persitant, directed movement

#check how many zeros
whichzero <- which(hmmdata$step == 0)

#proportion of ssteps of length 0 in data
length(whichzero)/nrow(hmmdata)


#### Choosing starting values in parellel ####
# good resource here is 

#running in parallel is faster as it uses different computer cores 
ncores <- detectCores() - 1
cl <- makeCluster(getOption("cl.cores", ncores))
# Export objects needed in parallelised function to cluster
clusterExport(cl, list("hmmdata", "fitHMM"))

#Number of tries with different starting values
niter <- 25
# Create list of starting values
allPar0 <- lapply(as.list(1:niter), function(x) {
  # Step length mean
  stepMean0 <- runif(3,
                     min = c(0.001, 0.3, 1),
                     max = c(0.3, 1,10))
  # Step length standard deviation - when gamma dist (default), the sd is on the same magnitude as the mean
  stepSD0 <- runif(3,
                   min = c(0.001, 0.3, 1),
                   max = c(0.3, 1,10))
  # Turning angle mean
  angleMean0 <- c(0, 0, 0)
  # Turning angle concentration
  angleCon0 <- runif(3,
                     min = c(0.1, 1,5),
                     max = c(1, 5,10))
  # Return vectors of starting values
  stepPar0 <- c(stepMean0, stepSD0)
  anglePar0 <- c(angleMean0, angleCon0)
  return(list(step = stepPar0, angle = anglePar0))
})
# Fit the niter models in parallel
allm_parallel <- parLapply(cl = cl, X = allPar0, fun = function(par0) {
  m <- fitHMM(data = hmmdata, nbStates = 3, stepPar0 = par0$step,
              anglePar0 = par0$angle)
  return(m)
})


# Extract likelihoods of fitted models
allnllk <- unlist(lapply(allm_parallel, function(m) m$mod$minimum))
allnllk


#Index of best fitting model (smallest negative log-likelihood)
whichbest <- which.min(allnllk)
# Best fitting model
mbest <- allm_parallel[[whichbest]]
mbest


#### Visualize in plots ####

#here we see the distribution of states
plot(mbest)


#state 1 is most like chick tending (or foraging) behavior, and would be difficult to seperate those two
par(mfrow=c(10,4))
plotStates(mbest)
dev.off()

#checking the global fit - this looks pretty good
par(mar=  c(1,1,1,1))
plotPR(mbest)
dev.off()

#### Printing outputs ####

#print the movement data
#write.csv(mbest$data, file = "C:/Users/14064/Dropbox/BTGO Movmement Study/Output/HMM_GPS.csv")

#print the state probabilities - For a given model, computes the probability of the process being in the different states at each time point.
#write.csv(stateProbs(mbest), file = "C:/Users/14064/Dropbox/BTGO Movmement Study/Output/HMM_GPS_states.csv")

#most probable states for a given sequence - from the viterbi algorithm 
#write.csv(viterbi(mbest), file = "C:/Users/14064/Dropbox/BTGO Movmement Study/Output/HMM_GPS_viterbi_mostlikely_states.csv")

#
moveHMM::CI(mbest, 0.95)



