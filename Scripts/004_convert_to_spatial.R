##############################################################
# Authors: 
# Ana Benítez-López (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 25th of September 2019
# Modified July 2020

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script is to merge all data used in the analyses and convert it to spatial data for Figure 2
# in the paper: 

# Benítez-López et al. The island rule explains consistent patterns of 
# body size evolution across terrestrial vertebrates. 

##############################################################
# Packages needed                                         ####
##############################################################

#load libraries
library(ggplot2)
library(dplyr)
library(sf)
library(sp)
library(mapview)

#clean memory
 rm(list=ls())

##############################################################
# Importing datasets                                      ####
##############################################################

#load data
mamdata<-read.csv("Data/Final data/mamdata_ph.csv", header = TRUE) 
birddata<-read.csv("Data/Final data/birddata_ph.csv", header = TRUE)
reptdata<-read.csv("Data/Final data/reptdata_ph.csv", header = TRUE)
amphdata<-read.csv("Data/Final data/amphdata_ph.csv", header = TRUE)

mamdata<-mamdata[,c("Class", "Island", "Long_i", "Lat_i")]
birddata<-birddata[,c("Class", "Island", "Long_i", "Lat_i")]
reptdata<-reptdata[,c("Class", "Island", "Long_i", "Lat_i")]
amphdata<-amphdata[,c("Class", "Island", "Long_i", "Lat_i")]

# merge datasets
data<-rbind(mamdata,birddata, reptdata,amphdata)

# save full dataset
# write.csv(data, file=" Data/Final data/islanddata_ms.csv")

##############################################################
# Converting data to SpatialPointsDataFrame               ####
##############################################################
xy <- data[,c("Long_i", "Lat_i")]
spdata <- SpatialPointsDataFrame(coords = xy, data = data,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

spdata <- st_as_sf(data, coords = c('Long_i', 'Lat_i'),crs = "+init=epsg:4326")

#import countries
# country<-st_read("Spatial/country1mVMAP_un2008.shp")
# mapview(spdata, zcol = "Class", legend = TRUE)+ mapview(country)

#save spatialpoints dataframe --> island locations, give full path
st_write(spdata, "Spatial/spdata.shp", update=TRUE)

#############################################################################################
# Make points of different size depending on the number of insular populations included ####
############################################################################################
### analysis per group

mamdata_group<- mamdata %>% 
        group_by(Island) %>%
        summarize(Class=first(Class),
                  Lat_i= first(Lat_i),
                  Long_i = first(Long_i),
                  count =n())

birddata_group<- birddata %>% 
  group_by(Island) %>%
  summarize(Class=first(Class),
            Lat_i= first(Lat_i),
            Long_i = first(Long_i),
            count =n())

reptdata_group<- reptdata %>% 
  group_by(Island) %>%
  summarize(Class=first(Class),
            Lat_i= first(Lat_i),
            Long_i = first(Long_i),
            count =n())

amphdata_group<- amphdata %>% 
  group_by(Island) %>%
  summarize(Class=first(Class),
            Lat_i= first(Lat_i),
            Long_i = first(Long_i),
            count =n())

data_group<-rbind(mamdata_group,birddata_group, reptdata_group,amphdata_group)

# convert to spatial
xy <- data_group[,c("Long_i", "Lat_i")]
spdata_group <- SpatialPointsDataFrame(coords = xy, data = data_group,
                                 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

#import countries
# country<-st_read("Spatial/country1mVMAP_un2008.shp")
# mapview(spdata_group, zcol = "Class", legend = TRUE)+ mapview(country)

spdata_group <- st_as_sf(data_group, coords = c('Long_i', 'Lat_i'),crs = "+init=epsg:4326")

#save spatialpoints dataframe --> island locations, give full path
st_write(spdata_group, "Spatial/spdata_group.shp")

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/convert_to_spatial_data_R_session.txt")
sessionInfo()
sink()

### End of script ####