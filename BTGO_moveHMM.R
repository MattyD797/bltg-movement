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
options(max.print=1600) #if you need to see the complete output
#

#read in the data and preview

setwd("C:/Users/14064/Dropbox/BTGO Movmement Study/BTGO git/bltg-movement/CSV/")
data1 <- read.csv("./BTGO_Haanmeer_nestR.csv"); head(data1)

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
BTGO_move2 <- redisltraj(na.omit(BTGO_move),3600,type="time"); BTGO_move2

#convert to data frame
data <- ld(BTGO_move2)
head(data)

#keep only columns with id, date, x, y
data <- data[,c(1,2,3,11)]

#convert from m to km
data$x <- data$x/1000 
data$y <- data$y/1000

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













#### test fitHMM ####

## initial parameters for gamma and von Mises distributions

mu0 <- c(0.01,2) # step mean (two parameters: one for each state) in km. state 1 involving relatively short steps and many turnings, state 2 involving loger steps and fewer turnings
sigma0 <- c(0.01,2) # step SD

#zeromass0 <- c(0.1,0.05,0.01,0.001) # step zero-mass
angleMean0 <- c(pi,0) # angle mean
kappa0 <- c(1,1) # angle concentration - corresponds to how concentrated the data are around the mean, large is concentrated, small is not concentrated

stepPar0 <- c(mu0,sigma0)
anglePar0 <- c(angleMean0,kappa0)


m <- fitHMM(data=hmmdata, nbStates = 2, stepPar0=stepPar0,anglePar0 = anglePar0,formula=~1)

m


#### choosing starting values following ####

#try many starting values: Here were fit the model with many different sets of starting values, and select the best model fit among those
#One possible way to do this is to generate starting values at random, from a distribution of plausible values.
# We can use a uniform distribution fitted over the range of values observed in the histograms before

# For reproducibility
set.seed(12345)
# Number of tries with different starting values
niter <- 25
# Save list of fitted models
allm <- list()
for(i in 1:niter) {
  # Step length mean
  stepMean0 <- runif(2,
                     min = c(0.001, 0.1),
                     max = c(0.1, 0.5))
  # Step length standard deviation
  stepSD0 <- runif(2,
                   min = c(0.001, 0.1),
                   max = c(0.1, 0.5))
  # Turning angle mean
  angleMean0 <- c(pi, 0)
  # Turning angle concentration
  angleCon0 <- runif(2,
                     min = c(0.1, 1),
                     max = c(5, 10))
  # Fit model
  stepPar0 <- c(stepMean0, stepSD0)
  anglePar0 <- c(angleMean0, angleCon0)
  allm[[i]] <- fitHMM(data = hmmdata, nbStates = 2, stepPar0 = stepPar0,
                      anglePar0 = anglePar0)
}


# Extract likelihoods of fitted models
allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
allnllk

#these clearly have issues with convergence!!!!!!!!! MLE should not be in the billions of trillions!

#Index of best fitting model (smallest negative log-likelihood)
whichbest <- which.min(allnllk)
# Best fitting model
mbest <- allm[[whichbest]]
mbest
