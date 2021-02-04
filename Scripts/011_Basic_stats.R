##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 21st of August 2020


##############################################################
# Description of script and instructions                  ####
##############################################################

# This script calculates some basic descriptive stats for the paper: 

# Benitez-Lopez et al.The island rule explains consistent patterns of 
# body size evolution across terrestrial vertebrates. 


##############################################################
# Packages needed                                         ####
##############################################################

library(dplyr)
library(stringr)

#clean memory
rm(list=ls())

##############################################################
# Importing datasets and functions                        ####
##############################################################

#Load data
mamdata<-read.csv("Data/Final data/mamdata_ph.csv", header = TRUE, stringsAsFactors = FALSE)
birddata<-read.csv("Data/Final data/birddata_ph.csv", header = TRUE, stringsAsFactors = FALSE)
reptdata<-read.csv("Data/Final data/reptdata_ph.csv", header = TRUE, stringsAsFactors = FALSE)
amphdata<-read.csv("Data/Final data/amphdata_ph.csv", header = TRUE, stringsAsFactors = FALSE)

###############################################################
# Calculate sample sizes (comparisons, species, specimens) ####
###############################################################

# Total island-mainland comparisons
sum(nrow(mamdata),nrow(birddata), nrow(reptdata), nrow(amphdata))

# Total number species
sum(length(unique(mamdata$Binomial)),
    length(unique(birddata$Binomial)),
    length(unique(reptdata$Binomial)),
    length(unique(amphdata$Binomial)))

sum(length(unique(mamdata$Species_island)),
    length(unique(birddata$Species_island)),
    length(unique(reptdata$Species_island)),
    length(unique(amphdata$Species_island)))

# Total number specimens
sum(mamdata$N_m,birddata$N_m, reptdata$N_m, amphdata$N_m)
sum(mamdata$N_i,birddata$N_i, reptdata$N_i, amphdata$N_i)


###################################################################
# Calculate number of conspecific vs interspecific comparisons ####
###################################################################

#Mammals ####
# Total number of conspecific vs interspecific comparisons
nrow(mamdata[mamdata$Phylogeny == "Conspecific",]) #999

#Birds ####
# Total number of conspecific vs interspecific comparisons
nrow(birddata[birddata$Phylogeny == "Conspecific",]) #429

#Reptiles ####
# Total number of conspecific vs interspecific comparisons
nrow(reptdata[reptdata$Phylogeny == "Conspecific",]) #483

#Amphibians ####
# Total number of conspecific vs interspecific comparisons
nrow(amphdata[amphdata$Phylogeny == "Conspecific",]) #157; 12.2%

### Calculate for all datasets

(nrow(mamdata[mamdata$Phylogeny == "Conspecific",]) + 
  nrow(birddata[birddata$Phylogeny == "Conspecific",]) +
  nrow(reptdata[reptdata$Phylogeny == "Conspecific",]) +
  nrow(amphdata[amphdata$Phylogeny == "Conspecific",]))/
  (nrow(mamdata) + nrow(birddata) + nrow(reptdata) + nrow(amphdata))*100 #83.4%


###################################################################
# Calculate number of archipelagos included                    ####
###################################################################

(nrow(mamdata[mamdata$Archipielago == "Yes",]) + 
    nrow(birddata[birddata$Archipielago == "Yes",]) +
    nrow(reptdata[reptdata$Archipielago == "Yes",]) +
    nrow(amphdata[amphdata$Archipielago == "Yes",]))/
  (nrow(mamdata) + nrow(birddata) + nrow(reptdata) + nrow(amphdata))*100 #3.22%

###############################################################################################################
# Calculate number of museum and live specimens and number of species from Pigot et al 2020 and Borja Milá ####
###############################################################################################################

# Number specimens
sum(birddata[birddata$Reference == "Pigot et al. 2020",]$N_m) + sum(birddata[birddata$Reference == "Pigot et al. 2020",]$N_i) 
+ sum(birddata[birddata$Reference == "Milá unpublished",]$N_m) +sum(birddata[birddata$Reference == "Milá unpublished",]$N_i)

# Number species
length(unique(birddata[birddata$Reference == "Pigot et al. 2020" | birddata$Reference == "Milá unpublished",]$Species_main))
length(unique(birddata[birddata$Reference == "Pigot et al. 2020" | birddata$Reference == "Milá unpublished",]$Species_island))

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/basic_stats_R_session.txt")
sessionInfo()
sink()

### End of script ####