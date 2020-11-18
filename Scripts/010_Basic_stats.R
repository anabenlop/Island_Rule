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
all.equal(mamdata$Species_main, mamdata$Species_island) # 187 mismatches

mam_comp <- mamdata %>% 
            mutate(Conspecific = if_else(Species_main == Species_island, "yes", "no"))

# View(mam_comp[mam_comp$Conspecific == "no",]) # some are just subspecies

# Upgrade subspecies to species level for comparison

# mam_subsp <- mam_comp[mam_comp$Conspecific == "no",]

mam_subsp_main <- str_split(mam_comp$Species_main, " ")
mam_subsp_island<- str_split(mam_comp$Species_island, " ")

mam_subsp_all <- data.frame(Genus_Main = NA, Species_Main = NA, Subspecies_Main = NA,
                                 Genus_Island = NA, Species_Island = NA, Subspecies_Island = NA)

 # i <- 1
for (i in 1:length(mam_subsp_main)) {
    
    mam_subsp_all[i, "Genus_Main"] <- mam_subsp_main[[i]][1]
    mam_subsp_all[i, "Species_Main"] <- mam_subsp_main[[i]][2] 
    mam_subsp_all[i, "Subspecies_Main"] <- mam_subsp_main[[i]][3]
    mam_subsp_all[i, "Genus_Island"] <- mam_subsp_island[[i]][1]
    mam_subsp_all[i, "Species_Island"] <- mam_subsp_island[[i]][2] 
    mam_subsp_all[i, "Subspecies_Island"] <- mam_subsp_island[[i]][3]
}     

# Compare species names, detect which are subspecies and which interspecific comparisons

mam_subsp_all <- mam_subsp_all %>% 
    mutate(Conspecific = if_else(Species_Main == Species_Island, "yes", "no"))

nrow(mam_subsp_all[mam_subsp_all$Conspecific == "no",]) # 59 comparisons (5.6%)

length(unique(mam_subsp_all[mam_subsp_all$Conspecific == "no",]$Species_Main)) #42 mainland species
length(unique(mam_subsp_all[mam_subsp_all$Conspecific == "no",]$Species_Island)) #53 insular endemic species (16.16%)

#Birds ####
# Total number of conspecific vs interspecific comparisons
all.equal(birddata$Species_main, birddata$Species_island) # 382 mismatches

bird_comp <- birddata %>% 
    mutate(Conspecific = if_else(Species_main == Species_island, "yes", "no"))

# View(bird_comp[bird_comp$Conspecific == "no",]) # some are just subspecies

# Upgrade subspecies to species level for comparison

bird_subsp_main <- str_split(bird_comp$Species_main, " ")
bird_subsp_island<- str_split(bird_comp$Species_island, " ")

bird_subsp_all <- data.frame(Genus_Main = NA, Species_Main = NA, Subspecies_Main = NA,
                            Genus_Island = NA, Species_Island = NA, Subspecies_Island = NA)

# i <- 1
for (i in 1:length(bird_subsp_main)) {
    
    bird_subsp_all[i, "Genus_Main"] <- bird_subsp_main[[i]][1]
    bird_subsp_all[i, "Species_Main"] <- bird_subsp_main[[i]][2] 
    bird_subsp_all[i, "Subspecies_Main"] <- bird_subsp_main[[i]][3]
    bird_subsp_all[i, "Genus_Island"] <- bird_subsp_island[[i]][1]
    bird_subsp_all[i, "Species_Island"] <- bird_subsp_island[[i]][2] 
    bird_subsp_all[i, "Subspecies_Island"] <- bird_subsp_island[[i]][3]
    
}     


# Compare species names, detect which are subspecies and which interspecific comparisons

bird_subsp_all <- bird_subsp_all %>% 
    mutate(Conspecific = if_else(Species_Main == Species_Island, "yes", "no"))

nrow(bird_subsp_all[bird_subsp_all$Conspecific == "no",]) # 282 comparisons (39.49%)

length(unique(bird_subsp_all[bird_subsp_all$Conspecific == "no",]$Species_Main)) #216 mainland species
length(unique(bird_subsp_all[bird_subsp_all$Conspecific == "no",]$Species_Island)) #242 insular endemic species (39.80%)

#Reptiles ####
# Total number of conspecific vs interspecific comparisons
all.equal(reptdata$Species_main, reptdata$Species_island) # 160 mismatches

rept_comp <- reptdata %>% 
    mutate(Conspecific = if_else(Species_main == Species_island, "yes", "no"))

# View(rept_comp[rept_comp$Conspecific == "no",]) # some are just subspecies

# Upgrade subspecies to species level for comparison

rept_subsp_main <- str_split(rept_comp$Species_main, " ")
rept_subsp_island<- str_split(rept_comp$Species_island, " ")

rept_subsp_all <- data.frame(Genus_Main = NA, Species_Main = NA, Subspecies_Main = NA,
                             Genus_Island = NA, Species_Island = NA, Subspecies_Island = NA)

# i <- 1
for (i in 1:length(rept_subsp_main)) {
    
    rept_subsp_all[i, "Genus_Main"] <- rept_subsp_main[[i]][1]
    rept_subsp_all[i, "Species_Main"] <- rept_subsp_main[[i]][2] 
    rept_subsp_all[i, "Subspecies_Main"] <- rept_subsp_main[[i]][3]
    rept_subsp_all[i, "Genus_Island"] <- rept_subsp_island[[i]][1]
    rept_subsp_all[i, "Species_Island"] <- rept_subsp_island[[i]][2] 
    rept_subsp_all[i, "Subspecies_Island"] <- rept_subsp_island[[i]][3]
    
}     


# Compare species names, detect which are subspecies and which interspecific comparisons

rept_subsp_all <- rept_subsp_all %>% 
    mutate(Conspecific = if_else(Species_Main == Species_Island, "yes", "no"))

nrow(rept_subsp_all[rept_subsp_all$Conspecific == "no",]) # 79 comparisons (14.42%)

length(unique(rept_subsp_all[rept_subsp_all$Conspecific == "no",]$Species_Main)) #35 mainland species
length(unique(rept_subsp_all[rept_subsp_all$Conspecific == "no",]$Species_Island)) #46 insular endemic species (26.59%)

#Amphibians ####
# Total number of conspecific vs interspecific comparisons
all.equal(amphdata$Species_main, amphdata$Species_island) # 22 mismatches

amph_comp <- amphdata %>% 
    mutate(Conspecific = if_else(Species_main == Species_island, "yes", "no"))

# View(amph_comp[amph_comp$Conspecific == "no",]) # some are just subspecies

# Upgrade subspecies to species level for comparison

amph_subsp_main <- str_split(amph_comp$Species_main, " ")
amph_subsp_island<- str_split(amph_comp$Species_island, " ")

amph_subsp_all <- data.frame(Genus_Main = NA, Species_Main = NA, Subspecies_Main = NA,
                             Genus_Island = NA, Species_Island = NA, Subspecies_Island = NA)

# i <- 1
for (i in 1:length(amph_subsp_main)) {
    
    amph_subsp_all[i, "Genus_Main"] <- amph_subsp_main[[i]][1]
    amph_subsp_all[i, "Species_Main"] <- amph_subsp_main[[i]][2] 
    amph_subsp_all[i, "Subspecies_Main"] <- amph_subsp_main[[i]][3]
    amph_subsp_all[i, "Genus_Island"] <- amph_subsp_island[[i]][1]
    amph_subsp_all[i, "Species_Island"] <- amph_subsp_island[[i]][2] 
    amph_subsp_all[i, "Subspecies_Island"] <- amph_subsp_island[[i]][3]
    
}     


# Compare species names, detect which are subspecies and which interspecific comparisons

amph_subsp_all <- amph_subsp_all %>% 
    mutate(Conspecific = if_else(Species_Main == Species_Island, "yes", "no"))

nrow(amph_subsp_all[amph_subsp_all$Conspecific == "no",]) # 21 comparisons (11.73%)

length(unique(amph_subsp_all[amph_subsp_all$Conspecific == "no",]$Species_Main)) #15 mainland species
length(unique(amph_subsp_all[amph_subsp_all$Conspecific == "no",]$Species_Island)) #20 insular endemic species (34.48%)


### Calculate for all datasets

(nrow(mam_subsp_all[mam_subsp_all$Conspecific == "no",]) + 
  nrow(bird_subsp_all[bird_subsp_all$Conspecific == "no",]) +
  nrow(rept_subsp_all[rept_subsp_all$Conspecific == "no",]) +
  nrow(amph_subsp_all[amph_subsp_all$Conspecific == "no",]) )/
  (nrow(mamdata) + nrow(birddata) + nrow(reptdata) + nrow(amphdata))*100

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/basic_stats_R_session.txt")
sessionInfo()
sink()

### End of script ####