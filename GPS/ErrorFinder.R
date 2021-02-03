#Read in the movement data
movement_data <- read.csv("CSV/GPS_BTGO_Haanmeer_nestR.csv", header = TRUE)
head(movement_data)

#Read in the known nest locations
nest_loc <- read.csv("CSV/GPS Data/UvA_nest_locs.csv", header = TRUE)
head(nest_loc)

#Reformat nest locations into lat and long
library(proj4)
proj4string <- "+proj=utm +zone=31 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
xy <- nest_loc[,c(2,3)]
head(xy)
#Transform
pj <- project(xy, proj4string, inverse=TRUE)
latlon <- data.frame(lat=pj$y, lon=pj$x)
head(latlon)

#See more digits in lat and long
options("digits"=12)

#Combine NestID back to find the individual associatiated
nest_loc[,c(2,3)] <- latlon
head(nest_loc)
#rename to unit lat and long
library(dplyr)
nest_loc <- nest_loc %>% 
  dplyr::rename(
    lat = Easting,
    lon = Northing
  )
#remove noise
nest_loc <- nest_loc[,c(1:3)]
head(nest_loc)

#Read in transmitter information
transmitter_info <- read.csv("CSV/GPS Data/UvA_Transmitter_Info.csv", header = T)
head(transmitter_info)

#Reformat transmitter info to only include ID of individual, ID of nest and Year
transmitter_info.INY <- transmitter_info[,c(4,5,7)]
transmitter_info.INY$ID <- paste0(transmitter_info.INY$TransmitterID, "-",transmitter_info.INY$YearDeployed)
indAndNest <- transmitter_info.INY[,c(4,3)]
head(indAndNest)

#Add the known Hatch Dates
indAndNest.loc <- cbind(indAndNest, nest_loc[,c(2,3)], transmitter_info$HatchDate)
names(indAndNest.loc)[names(indAndNest.loc) == "transmitter_info$HatchDate"] <- "Hatch Date"
head(indAndNest.loc)
head(movement_data)


#reformat the date for the movement data
movement_data$date <- as.Date(as.character(movement_data$date), format = "%m/%d/%Y")

head(movement_data)


#reformat the date data
indAndNest.loc$`Hatch Date` <- gsub('Apr','04', indAndNest.loc$`Hatch Date`)
indAndNest.loc$`Hatch Date` <- gsub('May','05', indAndNest.loc$`Hatch Date`)
indAndNest.loc$`Hatch Date` <- gsub('Jun','06', indAndNest.loc$`Hatch Date`)
indAndNest.loc$`Hatch Date` <- gsub('Apr','04', indAndNest.loc$`Hatch Date`)
indAndNest.loc$`Hatch Date` <- gsub('May','05', indAndNest.loc$`Hatch Date`)
indAndNest.loc$`Hatch Date` <- gsub('Jun','06', indAndNest.loc$`Hatch Date`)


indAndNest.loc$`Hatch Date` <- as.Date(indAndNest.loc$`Hatch Date`, "%d-%m-%Y")
indAndNest.loc$`Hatch Date` <- as.Date(gsub('00','20', indAndNest.loc$`Hatch Date`))
head(indAndNest.loc)



#find smallest distance

results <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(results) <- c("ID", "dist", "knownlat", "knownlon")

#Run a for loop through each year of each individual
library(geosphere)
for (val in unique(movement_data$burst))
{
  if(val %in% indAndNest.loc$ID){
    
    output2 <- indAndNest.loc %>%
      dplyr::filter(ID == val) 
    
    known_coords <- output2[,c(4,3)]
    
    output1 <- movement_data %>%
      dplyr::filter(burst == val)
    
    output1.d <- output1 %>%
      dplyr::filter(date <= output2$`Hatch Date`)
    
    coords <- output1[,c(2,3)]
    
    if(nrow(output1.d) > 0)
    {
      de <- data.frame(val, min(distm(known_coords, coords)), output2$lat, output2$lon)
      results <- rbind(results, de)
    } else {
      de <- data.frame(val, -1, -1, -1)
      de <- de %>% 
        dplyr::rename(
          min.distm.known_coords..coords.. = X.1,
          output2.lat = X.1.1,
          output2.lon = X.1.2
        )
      results <- rbind(results, de)
    }
  }
}

results <- results %>% 
  dplyr::rename(
    "Distance (m)" = min.distm.known_coords..coords..,
    "Nest Latitude" = output2.lat,
    "Nest Longitude" = output2.lon
  )
str(results)


write.csv(results, "CSV/ErrorRateFiltered.csv")









