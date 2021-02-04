##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 11th of December 2020

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script runs the phylogenetic meta-regressions for all complete cases (no imputed SD) to test the island rule
# for the paper: 

# Benitez-Lopez et al.The island rule explains consistent patterns of 
# body size evolution across terrestrial vertebrates. 


##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)

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

# loading phylogenetic matrixes 
load("Data/Final data/mam_phylo_cor.Rdata") #mam_phylo_cor
load("Data/Final data/bird_phylo_cor.Rdata") #bird_phylo_cor
load("Data/Final data/rept_phylo_cor.Rdata") #rept_phylo_cor
load("Data/Final data/amph_phylo_cor.Rdata") #amph_phylo_cor

# load necessary functions
source("Scripts/000_Functions.R")

# Calculate varcovar matrix
Vmam<- bldiag(lapply(split(mamdata, mamdata$CommonControl), calc.v))
is.positive.definite(Vmam) # FALSE
Vmam<-PDfunc(Vmam)
is.positive.definite(Vmam) # TRUE

Vbird<- bldiag(lapply(split(birddata, birddata$CommonControl), calc.v))
is.positive.definite(Vbird) # FALSE
Vbird<-PDfunc(Vbird)
is.positive.definite(Vbird) # TRUE

Vrept<- bldiag(lapply(split(reptdata, reptdata$CommonControl), calc.v))
is.positive.definite(Vrept) # TRUE

Vamph<- bldiag(lapply(split(amphdata, amphdata$CommonControl), calc.v))
is.positive.definite(Vamph) # TRUE

################################################
# Testing insular size shifts: Island rule  ####
################################################
RE = list(~ 1 | Reference,~1|ID, ~1|SPID, ~1| Binomial)

#mammals####
phylocor<-list(Binomial= mam_phylo_cor)
metamam_bias<-rma.mv(RR~logmass,V=Vmam, subset = imputed == "No", data=mamdata, random= RE,  R = phylocor)
summary(metamam_bias)
mR2.func(metamam_bias)
cR2.func(metamam_bias)

#birds#### 
phylocor<-list(Binomial= bird_phylo_cor)
metabird_bias<-rma.mv(RR~logmass,V=Vbird, subset = imputed == "No",  data=birddata, random= RE, R = phylocor)
summary(metabird_bias)
mR2.func(metabird_bias)
cR2.func(metabird_bias)

#reptiles####
phylocor<-list(Binomial= rept_phylo_cor)
metarept_bias<-rma.mv(RR~logmass,V=Vrept, subset = imputed == "No", data=reptdata,  random= RE, R=phylocor) 
summary(metarept_bias)
mR2.func(metarept_bias)
cR2.func(metarept_bias)

#amphibians####
phylocor<-list(Binomial= amph_phylo_cor)
metaamph_bias<-rma.mv(RR~logmass ,V=Vamph, subset = imputed == "No", data=amphdata, random= RE, R = phylocor) 
summary(metaamph_bias)
mR2.func(metaamph_bias)
cR2.func(metaamph_bias)


# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/checkimputationbias_R_session.txt")
sessionInfo()
sink()

### End of script ###