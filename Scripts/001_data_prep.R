##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estacion Biologica de Donana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 28th of August 2019
# Modified in July 2020

##############################################################
# Description of script and instructions
##############################################################

# This script is to load, extract environmental variables (NDVI - resource availability; 
# temperature and precipitation), calculate population-level body size indices, 
# and subset the island rule database for the four 
# taxonomic groups: Mammals, Birds, Amphibians, Reptiles

# Benitez-Lopez et al. The island rule explains consistent patterns of 
# body size evolution

##############################################################
# Packages needed
##############################################################
#load libraries
library(sf)
library(ggplot2)
library(raster)

# Clear memory 
# rm(list=ls())

##############################################################
# Importing datasets
##############################################################

#load data
data<-read.csv("Data/islandrule_clean.csv", stringsAsFactors = FALSE)

#load environmental variables
env.var_0.1<-stack("Spatial/env.var_0.1.tif")
env.var_0.5<-stack("Spatial/env.var_0.5.tif")
env.var_1<- stack("Spatial/env.var_1.tif")
names(env.var_0.1)<-c("NDVI_0.1", "SDNDVI_0.1","tmean_0.1", "tseas_0.1","prec_0.1")
names(env.var_0.5)<-c("NDVI_0.5", "SDNDVI_0.5","tmean_0.5", "tseas_0.5","prec_0.5")
names(env.var_1)<-c("NDVI_1", "SDNDVI_1", "tmean_1", "tseas_1","prec_1")

##############################################################
# Extracting environmental variables
##############################################################
# extract coordinates first
data<-data[!is.na(data$Lat_i),] #remove NAs in spatial coordinates
data<-data[!is.na(data$Mean_m),] #remove datapoints that do not have mainland size
data<-data[!is.na(data$Mean_i),] #remove datapoints that do not have island size

coord<-data[,c("Long_i", "Lat_i")]

env.var.df_0.1<- raster::extract(env.var_0.1, coord) 
env.var.df_0.5<-raster::extract(env.var_0.5, coord)
env.var.df_1<-raster::extract(env.var_1, coord) 

data_temp<-data.frame(data,env.var.df_0.1,env.var.df_0.5,env.var.df_1)

# We need to multiply NDVI, SDNDVI, EVI and SDEVI by 0.0001 to get the original values
data_temp[,c("NDVI_0.1", "SDNDVI_0.1","NDVI_0.5", "SDNDVI_0.5","NDVI_1", "SDNDVI_1")]<-data_temp[,c("NDVI_0.1", "SDNDVI_0.1", "NDVI_0.5", "SDNDVI_0.5","NDVI_1", "SDNDVI_1")]*0.0001

##############################################################
# Assigning environmental variables based on island size  ####
##############################################################
data_temp$NDVI  <-ifelse(data_temp$Island_km2 <= 1000 & !is.na(data_temp$NDVI_0.1), data_temp$NDVI_0.1, 
                         ifelse(data_temp$Island_km2 > 1000 & data_temp$Island_km2 <= 10000 & !is.na(data_temp$NDVI_0.5), data_temp$NDVI_0.5,data_temp$NDVI_1))

data_temp$SDNDVI<-ifelse(data_temp$Island_km2 <= 1000 & !is.na(data_temp$SDNDVI_0.1), data_temp$SDNDVI_0.1, 
                         ifelse(data_temp$Island_km2 > 1000 & data_temp$Island_km2 <= 10000 & !is.na(data_temp$SDNDVI_0.5), data_temp$SDNDVI_0.5,data_temp$SDNDVI_1))

data_temp$tmean <-ifelse(data_temp$Island_km2 <= 1000 & !is.na(data_temp$tmean_0.1), data_temp$tmean_0.1, 
                         ifelse(data_temp$Island_km2 > 1000 & data_temp$Island_km2 <= 10000 & !is.na(data_temp$tmean_0.5), data_temp$tmean_0.5,data_temp$tmean_1))

data_temp$tseas <-ifelse(data_temp$Island_km2 <= 1000 & !is.na(data_temp$tseas_0.1), data_temp$tseas_0.1, 
                         ifelse(data_temp$Island_km2 > 1000 & data_temp$Island_km2 <= 10000 & !is.na(data_temp$tseas_0.5), data_temp$tseas_0.5,data_temp$tseas_1))

data_temp$prec  <-ifelse(data_temp$Island_km2 <= 1000 & !is.na(data_temp$prec_0.1), data_temp$prec_0.1, 
                         ifelse(data_temp$Island_km2 > 1000 & data_temp$Island_km2 <= 10000 & !is.na(data_temp$prec_0.5), data_temp$prec_0.5,data_temp$prec_1))

# There are some islands without environmental information. Assign based on lower resolution rasters.
data_temp[is.na(data_temp$tmean),'Island']


#there's not climatic raster data for the island St. Helena. I assign it from values reporte on climatic reports: 
# Feistel, R., Hagen, E., & Grant, K. (2003). Climatic changes in the subtropical Southeast Atlantic: the St. Helena Island climate index (1893–1999). 
# Progress in Oceanography, 59(2-3), 321-337.https://www.sciencedirect.com/science/article/pii/S0079661103001678?via%3Dihub

data_temp[data_temp$Island == "St Helena","tmean"]<-mean(c(19.48,20.53,20.67,20.13,18.92,17.45,16.25,15.7,15.72,16.09,16.89,18.05)) #Table 5
data_temp[data_temp$Island == "St Helena","tseas"]<-sd(c(19.48,20.53,20.67,20.13,18.92,17.45,16.25,15.7,15.72,16.09,16.89,18.05))*100 #Table 5
data_temp[data_temp$Island == "St Helena","prec"]<-sum(33.63,50.31,64.93,48.63,50.03,58.43,61.86,48.40,34.86,24.34,18.23,22.21) #Table 5
data_temp[data_temp$Island == "St Helena","pseas"]<-sd(c(33.63,50.31,64.93,48.63,50.03,58.43,61.86,48.40,34.86,24.34,18.23,22.21))/mean(c(33.63,50.31,64.93,48.63,50.03,58.43,61.86,48.40,34.86,24.34,18.23,22.21))*100 #Table 5

#2935 rows

##############################################################
# Making sure the data is normally distributed            ####
##############################################################
data<-data_temp

data$Island_km2<-log10(data$Island_km2)
data$Dist_near_mainland<-log10(data$Dist_near_mainland +0.001)
data$SDNDVI<-log10(data$SDNDVI+0.001)
data$prec<-log10(data$prec)
data$tseas<-log10(data$tseas)

#########################################################################
# Calculating population-level body size based on male and female data ##
#########################################################################
data$Mean_m_pooled <- ifelse(data$Sex == "Both" & c(!is.na(data$Mean_m_male) & !is.na(data$Mean_m_female)), 
                             (data$Mean_m_male*data$N_m_male + data$Mean_m_female*data$N_m_female)/(data$N_m_male+data$N_m_female),
                             data$Mean_m)

data[is.na(data$Mean_m_pooled), "Mean_m_pooled"] <- data[is.na(data$Mean_m_pooled), "Mean_m"] #fix 2 for which sample size was not given by sex

ggplot(data) + geom_point(aes(log10(Mean_m), log10(Mean_m_pooled))) + geom_smooth(aes(log10(Mean_m), log10(Mean_m_pooled)), method="lm") + 
  geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed",  size = 1.2)

summary(lm(log10(Mean_m_pooled) ~log10(Mean_m), data=data)) #similar to normal average

data$Mean_i_pooled <- ifelse(data$Sex == "Both" & c(!is.na(data$Mean_i_male) & !is.na(data$Mean_i_female)),
                             (data$Mean_i_male*data$N_i_male + data$Mean_i_female*data$N_i_female)/(data$N_i_male + data$N_i_female),
                             data$Mean_i)

data[is.na(data$Mean_i_pooled), "Mean_i_pooled"] <- data[is.na(data$Mean_i_pooled), "Mean_i"] #fix 2 for which sample size was not given by sex

ggplot(data) + geom_point(aes(log10(Mean_i), log10(Mean_i_pooled))) + geom_smooth(aes(log10(Mean_i), log10(Mean_i_pooled)), method="lm") + 
  geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed",  size = 1.2)

summary(lm(log10(Mean_i_pooled) ~log10(Mean_i), data=data)) #similar to normal average

# Calculate pooled estimates as in Table 6.5.a (https://training.cochrane.org/handbook/current)
# usual pooled SD provides a within-subgroup SD rather than an SD for the combined group, so provides an underestimate of the desired SD.
# We use the formula for SD for the combined group

# mainland data
# data$sd_m_pooled <- ifelse(data$Sex == "Both" & c(!is.na(data$sd_m_male) & !is.na(data$sd_m_female)), 
#                            sqrt(((data$N_m_male-1)*data$sd_m_male^2 + 
#                                    (data$N_m_female-1)*data$sd_m_female^2)/  
#                                   (data$N_m_male+data$N_m_female-1)),
#                            data$sd_m)

data$sd_m_combined <- ifelse(data$Sex == "Both" & c(!is.na(data$sd_m_male) & !is.na(data$sd_m_female)), 
                           sqrt(((data$N_m_male-1)*data$sd_m_male^2 + 
                                   (data$N_m_female-1)*data$sd_m_female^2 +
                                   ((data$N_m_male*data$N_m_female)/(data$N_m_male+data$N_m_female))*
                                   (data$Mean_m_male^2 + data$Mean_m_female^2 - 2*(data$Mean_m_male*data$Mean_m_female)))/  
                                  (data$N_m_male+data$N_m_female-1)),
                           data$sd_m)

# ggplot(data) + geom_point(aes(log10(sd_m_pooled), log10(sd_m_combined))) + geom_smooth(aes(log10(sd_m_pooled), log10(sd_m_combined)), method="lm") + 
#   geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed",  size = 1.2) #using normal pooled formula are systematically underestimated
# 
# ggplot(data) + geom_point(aes(sd_m_pooled, sd_m_combined)) + geom_smooth(aes(sd_m_pooled, sd_m_combined), method="lm") + 
#   geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed",  size = 1.2) +xlim(0,1000) + ylim(0,1000) #using normal pooled formula are systematically underestimated

# island data
# data$sd_i_pooled <- ifelse(data$Sex == "Both" & c(!is.na(data$sd_i_male) & !is.na(data$sd_i_female)), 
#                            sqrt(((data$N_i_male-1)*data$sd_i_male^2 + 
#                                    (data$N_i_female-1)*data$sd_i_female^2)/  
#                                   (data$N_i_male+data$N_i_female-1)),
#                            data$sd_i)
# 

data$sd_i_combined <- ifelse(data$Sex == "Both" & c(!is.na(data$sd_i_male) & !is.na(data$sd_i_female)), 
                             sqrt(((data$N_i_male-1)*data$sd_i_male^2 + 
                                     (data$N_i_female-1)*data$sd_i_female^2 +
                                     ((data$N_i_male*data$N_i_female)/(data$N_i_male+data$N_i_female))*
                                     (data$Mean_i_male^2 + data$Mean_i_female^2 - 2*(data$Mean_i_male*data$Mean_i_female)))/  
                                    (data$N_i_male+data$N_i_female-1)),
                             data$sd_i)
# 
# ggplot(data) + geom_point(aes(log10(sd_i_pooled), log10(sd_i_combined))) + geom_smooth(aes(log10(sd_i_pooled), log10(sd_i_combined)), method="lm") + 
#   geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed",  size = 1.2) #using normal pooled formula are systematically underestimated
# 
# i<-ggplot(data) + geom_point(aes(sd_i, sd_i_combined)) + geom_smooth(aes(sd_i, sd_i_combined), method="lm") + 
#   geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed",  size = 1.2) +xlim(0,1000) + ylim(0,1000) #using normal pooled formula are systematically underestimated
# 
# library(plotly)
# ggplotly(i)

# rename variables for simplicity
data$Mean_m <- data$Mean_m_pooled
data$Mean_i <- data$Mean_i_pooled
data$sd_m <- data$sd_m_combined
data$sd_i <- data$sd_i_combined

##############################################################
# Subsetting databases by class for specific analyses     ####
##############################################################

mamdata<-subset(data, Class == "Mammals") 
birddata<-subset(data, Class == "Birds")
amphdata<-subset(data, Class == "Amphibians")
reptdata<-subset(data, Class == "Reptiles")

# save data
write.csv(mamdata,file= "Data/mamdata.csv", row.names = FALSE)
write.csv(birddata,file= "Data/birddata.csv",  row.names = FALSE)
write.csv(amphdata,file= "Data/amphdata.csv",  row.names = FALSE)
write.csv(reptdata,file= "Data/reptdata.csv",  row.names = FALSE)

# saving session information with all packages versions for reproducibility purposes
sink("~/New projects/Island rule/Data/Final data/data_prep_R_session.txt")
sessionInfo()
sink()

### End of script ####
