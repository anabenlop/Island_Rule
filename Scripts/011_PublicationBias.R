##############################################################
# Authors: 
# Ana Benítez-López (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 08th of September 2020

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script is to test the influence of data source on the island rule analyses. Here we see test if island-mainland comparisons 
# from studies testing the island rule are biased towards a significant effect compared to comparisons from studies not testing the 
# island rule. We test the effect of data source on the intercept, and on the intercept and slope of the lnRR ~ mainland mass regression. 
# Analyses done for the paper: 

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
# Vrept<-PDfunc(Vrept)
# is.positive.definite(Vrept) # TRUE

Vamph<- bldiag(lapply(split(amphdata, amphdata$CommonControl), calc.v))
is.positive.definite(Vamph) # TRUE

####################################################
# FIT MODELS WITH DATA SOURCE TYPE AS MODERATOR ####
####################################################
# mammals####
RE = list(~ 1 | Reference,~1|ID, ~1|SPID, ~1| Binomial)

#mammals####
phylocor<-list(Binomial= mam_phylo_cor)

# rma.mv
metamam_bias <- rma.mv(RR ~ logmass + Data_source_type,V=Vmam,  data=mamdata,
                  random= RE, R= phylocor) 
summary(metamam_bias)
anova(metamam_bias, btt = 3)

saveRDS(metamam_bias, file = "Data/Final data/metamam_bias.Rdata")

metamam_bias2 <- rma.mv(RR ~ logmass*Data_source_type,V=Vmam,  data=mamdata,
                       random= RE, R= phylocor) 
summary(metamam_bias2)
anova(metamam_bias2, btt = 4)

saveRDS(metamam_bias2, file = "Data/Final data/metamam_bias2.Rdata")


#birds####
phylocor<-list(Binomial= bird_phylo_cor)
metabird_bias<-rma.mv(RR ~ logmass + Data_source_type, V=Vbird,  data=birddata,
                  random= RE, R= phylocor) 
summary(metabird_bias)
anova(metabird_bias, btt = 3)

saveRDS(metabird_bias, file = "Data/Final data/metabird_bias.Rdata")

metabird_bias2 <- rma.mv(RR ~ logmass*Data_source_type,V=Vbird,  data=birddata,
                        random= RE, R= phylocor) 
summary(metabird_bias2)
anova(metabird_bias2, btt = 4)

saveRDS(metabird_bias2, file = "Data/Final data/metabird_bias2.Rdata")


#reptiles####
phylocor<-list(Binomial= rept_phylo_cor)
metarept_bias<-rma.mv(RR~logmass+Data_source_type,V=Vrept,  data=reptdata,  random= RE, R= phylocor) 
summary(metarept_bias)

anova(metarept_bias, btt = 3)

saveRDS(metarept_bias, file = "Data/Final data/metarept_bias.Rdata")

metarept_bias2<-rma.mv(RR~logmass*Data_source_type,V=Vrept,  data=reptdata,  random= RE, R= phylocor) 
summary(metarept_bias2)

anova(metarept_bias2, btt = 4)

saveRDS(metarept_bias2, file = "Data/Final data/metarept_bias2.Rdata")

#amphibians####
phylocor<-list(Binomial=amph_phylo_cor)

metaamph_bias<-rma.mv(RR~logmass +Data_source_type ,V = Vamph,  data=amphdata, random=RE, R = phylocor)
summary(metaamph_bias)

anova(metaamph_bias, btt = 3)

saveRDS(metaamph_bias, file = "Data/Final data/metaamph_bias.Rdata")

metaamph_bias2<-rma.mv(RR~logmass*Data_source_type,V = Vamph,  data=amphdata,  random= RE, R= phylocor) 
summary(metaamph_bias2)

anova(metaamph_bias2, btt = 4)

saveRDS(metaamph_bias2, file = "Data/Final data/metaamph_bias2.Rdata")

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/publication_bias_R_session.txt")
sessionInfo()
sink()

### End of script ####


