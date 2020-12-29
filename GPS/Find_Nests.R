# date: 10.27.20
# Author: Luke Wilde

# Script for estimating nesting location and fate from GPS data of BTGO using nestR

#install the package nestR from online source:
remotes::install_github("picardis/nestR", build_vignettes = TRUE)

#load the packages
x <- c("nestR", "ggplot2", "dplyr", "devtools", "lubridate"); lapply(x, require, character.only= TRUE)

#open vignette for workflow
vignette("nestR")

#### Process GPS data into format of nestR data ####
#clean old objects from R
rm(list=ls(all=TRUE))

#see outputs in standard (not scientific notation)
options(scipen = 999)
#

#set working directory
setwd("./CSV/")

#load the data
godwit <- read.csv("./GPS_BTGO_Haanmeer_nestR.csv"); head(godwit)

#the package has the date column as class POSIXct and the lat, long as numerics
godwit$date <- substr(godwit$date, 1, nchar(as.character(godwit$date)))

options(digits=10)
godwit$lat <- as.numeric(godwit$lat) #need the coords in numeric
godwit$long <- as.numeric(godwit$long)
godwit$date <- as.POSIXct(as.character(godwit$date), tz="GMT", format = "%m/%d/%Y %H:%M:%S") #converted date to a POSIX object

#now we have the data in the same format as the example datasets
#can compare to example data in the nestR package using: data(woodstorks); str(woodstorks); head(woodstorks)
str(godwit); head(godwit)

#randomly select two birds to test with the find_nests function
rand <- sample(levels(as.factor(godwit$burst)), size = 2)
godwit.test <- godwit %>% filter(burst == rand[1] | burst == rand[2])
head(godwit.test)

#from Mo's 2020 paper (https://www.frontiersin.org/articles/10.3389/fevo.2019.00031/full) breeding is Mar to July, average incubation is 24 days after last egg laid (which takes 3 to 4 days)
# we will classify a complete breeding cycle as 27 d (+/- 3.4 d)

#UvA-bits trackers have error of < 7m
#data were recorded every 5 min when fully charged and every 15-30 minutes when not
#Senner et al. 2018 scrubbed data for only times where interval did not exceed 30 minutes

BTGO_output_1 <- find_nests(gps_data = godwit,
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

# Prior information is not typically available for min_consec, min_top_att and min_days_att - and so the next section focuses on refining these
#all in all, the user should specify low values for the constraints so as to identify all recurrently visited locations


#the output of find_nests() is a list of nests and visits
BTGO_output_1$nests
BTGO_output_1$visits
           
