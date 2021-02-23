# Date: 02/10/21
# Author: Luke Wilde

#Script to process and implement momentuHMM on GPS birds 

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
foo( c("dplyr", "sp", "ggplot2", "ggmap", "plyr", "lubridate", "SDLfilter", "moveHMM", "momentuHMM", "bsam", "parallel", "adehabitatLT", "recurse", "fitdistrplus", "logspline"))

#clean old objects from R
rm(list=ls(all=TRUE))

#see outputs in standard (not scientific notation)
options(scipen = 999)
options(max.print=100) #if you need to see the complete output
#

#read in the data and preview


data1 <- read.csv("./CSV/GPS_BTGO_Haanmeer_nestR.csv"); head(data1)

#these two birds had too much missing data
data <- dplyr::filter(data1, burst != "1009-2012" | burst != "2002-2014" | burst != "1002-2012"| burst != "2002-2015" | burst != "2010-2013" | burst != "2004-2014" | burst != "2041-2013")
#data <- dplyr::filter(data2, burst != "2002-2014")

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

#DOY <- as.POSIXct(strptime(as.character(data$DateTime), "%Y-%m-%d"))

BTGO_move_1 <- data.frame(id = data$id,date = data$DateTime, UTM_W = utm_locs$UTM_W, UTM_N = utm_locs$UTM_N); head(BTGO_move_1)

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

test <- filter(data, ID == "1009-2013" | ID ==  "2004-2013" | ID == "2010-2013")

head(test)


#### calculate recursion and residence time ####

Num <- length(unique(as.vector(test$ID)))

#loop_object <- vector("list", length(1:Num)) 
# mat <- matrix(ncol=1, nrow=length(test)) # create empty data frame to populate with for-loop
# for(i in 1:Num){mat[i,] <- getRecursions(subset(test, ID == unique(ID)[i]), 7)} 
# recursivity <- as.data.frame(do.call(rbind, loop_object))
# head(recursivity)

# for(i in length(unique(animals$id))){
#   recursions_all = getRecursions(subset(animals, id == unique(id)[i]), 7)
# }

recursions_all = getRecursions(data, radius = 3, timeunits = c("hours"))

data <- cbind(data, recursions_all$revisits)
data <- cbind(data, recursions_all$residenceTime)

#rename
names(data)[5] <- "revisits"; names(data)[6] <- "residence_time" 
head(data)

#### create HMM object for momentuHMM ####
#create hmm data object
hmmdata <- momentuHMM::prepData(data, type="UTM", coordNames = c("x", "y"))

head(hmmdata) #step lengths are in km

#### Visualize parameters ####
#we imagine four states: adult foraging, chick-tending, nesting, migrating
par(mfrow = c(1,1,1,1))
hist(hmmdata$step);quantile(hmmdata$step, probs = seq(0,1,.05), na.rm = T)

hist(hmmdata$angle);quantile(hmmdata$angle, probs = seq(0,1,.05), na.rm = T) # TA near 0 correspond to persistent, directed movement

hist(hmmdata$revisits); quantile(hmmdata$revisits, probs = seq(0,1,.05))

hist(hmmdata$residence_time); quantile(hmmdata$residence_time, probs = seq(0,1,.05))


#check how many zeros
whichzero <- which(hmmdata$step == 0)

#proportion of ssteps of length 0 in data
length(whichzero)/nrow(hmmdata)

#### choosing the correct distribution for revisits and residence_time ####

#these packages allow you to visualize how close your observed distribution fits the theoretical distributions
library(fitdistrplus)
library(logspline) #install.packages("logspline")

#plot - the blue dot is your data, and depending on how close it is to markers (which can be shaded regions, points, or lines), you get an estimate of the proper distribution
descdist(hmmdata$revisits, discrete = F)
descdist(hmmdata$residence_time, discrete = F)


#it would seem that revisits is closest to a normal or uniform distribution and residence time is closest to either a beta or exponential
#we can data the fit of either and compare quantitatively to choose

revisit.fit.normal <- fitdist(hmmdata$revisits, "norm")
revisit.fit.uniform <- fitdist(hmmdata$revisits, "unif")
revisit.fit.gamma <- fitdist(hmmdata$revisits, "gamma")

residence_time.fit.lnorm <- fitdist(hmmdata$residence_time, "lnorm")
residence_time.fit.gamma<- fitdist(hmmdata$residence_time, "gamma")

plot(revisit.fit.normal)
plot(revisit.fit.uniform)
plot(revisit.fit.gamma)
plot(residence_time.fit.lnorm)
plot(residence_time.fit.weibull)

revisit.fit.normal$AIC ; revisit.fit.uniform$AIC ; revisit.fit.gamma$AIC 

residence_time.fit.lnorm$AIC ;residence_time.fit.gamma$AIC 


# revisits = gamma distribution, residence_time = weibull
#

#### attempt with recurse data ####
nb <- 2
stateNames <- c("nonnesting", "nesting")
dist = list(step = "gamma", angle = "wrpcauchy", revisits = "gamma", residence_time = "lnorm")


#these step and angle came from the initials in model 7 without recurse data, the revisits and residence time were chosen from the quantile vector
#angle concentration needs to be between 0 - 1
Par0_m1 <- list(step=c(0.0274685066089, 0.4155683122342, 0.0156864326750, 0.9640296170488),angle=c(0.202341507835, 0.9370594799590), revisits=c(50,13,50,13), residence_time =c(2612,1400,2612,1400))

m1 <- momentuHMM::fitHMM( data = hmmdata,   nbStates = nb,
  dist = dist, 
  Par0 = Par0_m1, 
  estAngleMean = list(angle=FALSE), 
  stateNames = stateNames)

m1

plot(m1)



#### Choosing starting values in parallel - without recurse (give sense of step and angle) ####
# good resource here is 

#create hmm data object
hmmdata_moveHMM <- moveHMM::prepData(data[,c(1,2,4)], type="UTM", coordNames = c("x", "y"))

head(hmmdata_moveHMM); class(hmmdata_moveHMM)

# Create cluster of size ncores
ncores <- detectCores() - 1
cl <- makeCluster(getOption("cl.cores", ncores))
# Export objects needed in parallelised function to cluster
clusterExport(cl, list("hmmdata_moveHMM", "fitHMM"))
# Number of tries with different starting values
niter <- 25
# Create list of starting values
allPar0 <- lapply(as.list(1:niter), function(x) {
  
  ##
  # stateNames: nesting, chick tending, foraging, migrating
  ##
  
  # Step length mean
  stepMean0 <- runif(4,
                     min = c(0.005, 0.03, 0.3,1),
                     max = c(0.03,0.3, 1,5))
  # Step length standard deviation
  stepSD0 <- runif(4,
                   min = c(0.005, 0.03, 0.3,1),
                   max = c(0.03,0.3, 1,5))
  # Turning angle mean
  angleMean0 <- c(0, 0, 0, 0)
  # Turning angle concentration - a large concentration indicates strong persistence in direction if the mean is 0, closer to 0 means undirected movement
  angleCon0 <- runif(4,
                     min = c(0.01,0.03,1, 3),
                     max = c(0.03,1, 3, 5))
  #revisits
  revisitsMean0 = runif(4, 
                        min = c(64,39,15,3),
                        max = c(150, 64,39,15))
  revisitsSD0 = runif(4, 
                        min = c(64,39,15,3),
                        max = c(150, 64,39,15))
  
  #residence time
  residence_timeMean0 = runif(4, 
                              min = c(2612,1800,1200,830),
                              max = c(3000,2600,1800,1200))
  residence_timeSD0 = runif(4, 
                              min = c(2612,1800,1200,830),
                              max = c(3000,2600,1800,1200))      
  
  # Return vectors of starting values
  stepPar0 <- c(stepMean0, stepSD0)
  anglePar0 <- c(angleMean0, angleCon0)
  revisitsPar0 <- c(revisitsMean0, revisitsSD0)
  residence_timePar0 <- c(residence_timeMean0, residence_timeSD0)
  return(list(step = stepPar0, angle = anglePar0, revisits = revisitsPar0, residence_time = residence_timePar0))
})

# Fit the niter models in parallel
allm_parallel <- parLapply(cl = cl, X = allPar0, fun = function(par0) {
  m <- moveHMM::fitHMM(hmmdata_moveHMM, nbStates = 4, stepPar0 = par0$step,  anglePar0 = par0$angle)
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


#### Choosing starting values in parallel - with recurse ####
# good resource here is 

# Create cluster of size ncores
ncores <- detectCores() - 1
cl <- makeCluster(getOption("cl.cores", ncores))
# Export objects needed in parallelised function to cluster
clusterExport(cl, list("hmmdata", "momentuHMM::fitHMM"))
# Number of tries with different starting values
niter <- 25
# Create list of starting values
allPar0 <- lapply(as.list(1:niter), function(x) {
  
  ##
  # stateNames: nesting, chick tending, foraging, migrating
  ##
  
  # Step length mean
  stepMean0 <- runif(2,
                     min = c(0.005, 0.3),
                     max = c(0.2,0.9))
  # Step length standard deviation
  stepSD0 <- runif(2,
                   min = c(0.005, 0.3),
                   max = c(0.2,0.9))
  # Turning angle mean
  angleMean0 <- c(0, 0)
  # Turning angle concentration - a large concentration indicates strong persistance in direction if the mean is 0, closer to 0 means undirected movement
  angleCon0 <- runif(2,
                     min = c(0.01,1),
                     max = c(0.03,3))
  # Return vectors of starting values
  stepPar0 <- c(stepMean0, stepSD0)
  anglePar0 <- c(angleMean0, angleCon0)
  return(list(step = stepPar0, angle = anglePar0))
})

# Fit the niter models in parallel
allm_parallel <- parLapply(cl = cl, X = allPar0, fun = function(par0) {
  m <- moveHMM::fitHMM(hmmdata_moveHMM, nbStates = 2, stepPar0 = par0$step,  anglePar0 = par0$angle)
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



