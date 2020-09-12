##############################################################
# Authors: 
# Ana Benitez-López (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 09 of February 2020

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script is to produce the Extended Data Figure 3 
# in the paper: 

# Benitez-Lopez et al.The island rule explains consistent patterns of 
# body size evolution across terrestrial vertebrates

##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(RColorBrewer)

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

# Load models
metamam <- readRDS(file = "Data/Final data/metamam.Rdata")
metabird <- readRDS(file = "Data/Final data/metabird.Rdata")
metarept <- readRDS(file = "Data/Final data/metarept.Rdata")
metaamph <- readRDS(file = "Data/Final data/metaamph.Rdata")

# load necessary functions
source("Scripts/000_Functions.R")

# Mammals ####
# var exp by moderators
mR2.func(metamam)
cR2.func(metamam)

#variance components
sigma_s<-metamam$sigma2[1]/sum(metamam$sigma2)
sigma_s*100 # 16.56% variance due to between-study effects
sigma_o<-metamam$sigma2[2]/sum(metamam$sigma2)
sigma_o*100 # 35.03.48% variance due within-study effects 
sigma_sp<-metamam$sigma2[3]/sum(metamam$sigma2)
sigma_sp*100 # 24.36% variance due to between-species variation (non-phylogenetic variaton at species level)
H2 <-metamam$sigma2[4]/sum(metamam$sigma2) 
H2 # lamdba or phylogenetic heritability - 0.14
H2*100

mam_variance <-data.frame(class = "Mammals", study = sigma_s, species = sigma_sp, phylogeny = H2, observ = sigma_o)

#birds#### 
# var exp by moderators
mR2.func(metabird)
cR2.func(metabird)

#variance components
sigma_s<-metabird$sigma2[1]/sum(metabird$sigma2)
sigma_s*100 # 1.42% variance due to between-study effects
sigma_o<-metabird$sigma2[2]/sum(metabird$sigma2)
sigma_o*100 # 54.1 % variance due within-study effects(residual)
sigma_sp<-metabird$sigma2[3]/sum(metabird$sigma2)
sigma_sp*100 # 38.9% variance due to between-species variation (non-phylogenetic variaton at species level)
H2 <-metabird$sigma2[4]/sum(metabird$sigma2) 
H2 # lamdba or phylogenetic heritability - 0.06
H2*100

bird_variance <-data.frame(class = "Birds", study = sigma_s, species = sigma_sp, phylogeny = H2, observ = sigma_o)

#reptiles####
# var exp by moderators
mR2.func(metarept)
cR2.func(metarept)

#variance components
sigma_s<-metarept$sigma2[1]/sum(metarept$sigma2)
sigma_s*100 # 23.5% variance due to between-study effects
sigma_o<-metarept$sigma2[2]/sum(metarept$sigma2)
sigma_o*100 # 19.4 % variance due within-study effects(residual)
sigma_sp<-metarept$sigma2[3]/sum(metarept$sigma2)
sigma_sp*100 # 29.6% variance due to between-species variation (non-phylogenetic variaton at species level)
H2 <-metarept$sigma2[4]/sum(metarept$sigma2) 
H2 # lamdba or phylogenetic heritability - 0.28
H2*100

rept_variance <-data.frame(class = "Reptiles", study = sigma_s, species = sigma_sp, phylogeny = H2, observ = sigma_o)

#amphibians####
# var exp by moderators
mR2.func(metaamph)
cR2.func(metaamph)

#variance components
sigma_s<-metaamph$sigma2[1]/sum(metaamph$sigma2)
sigma_s*100 # 10.1% variance due to between-study effects
sigma_o<-metaamph$sigma2[2]/sum(metaamph$sigma2)
sigma_o*100 # 31.7 % variance due within-study effects(residual)
sigma_sp<-metaamph$sigma2[3]/sum(metaamph$sigma2)
sigma_sp*100 # 58.2% variance due to between-species variation (non-phylogenetic variaton at species level)
H2 <-metaamph$sigma2[4]/sum(metaamph$sigma2) 
H2 # lamdba or phylogenetic heritability - < 0.01
H2*100

amph_variance<-data.frame(class = "Amphibians", study = sigma_s, species = sigma_sp, phylogeny = H2, observ = sigma_o)

#merge database
vardata<-rbind(mam_variance, bird_variance, rept_variance,amph_variance)
colnames(vardata)<-c("Class", "Source", "Species", "Phylogeny", "Residual")

# adjust dataframe
vardata_t<-melt(vardata, id.vars = c(colnames(vardata)[1]))
vardata_t$value<-round(vardata_t$value*100,2)

#figure
display.brewer.all()

varexp_island<-ggplot(vardata_t) + geom_bar(aes(Class, value, fill = variable), stat ="identity") +
  scale_fill_brewer(palette = "Spectral") + theme_bw(base_size=24)+
  theme(element_blank(), axis.text=element_text(size=24)) +
  guides(fill=guide_legend(title=" "))+ ylab("Variation accounted by random factors (%)")+
  xlab(" ") 

tiff('Results/Figures/Fig S3_var_RE.tif', res=300, width=3100, height=3000)
varexp_island
dev.off()

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/variance_expl_RE_R_session.txt")
sessionInfo()
sink()

# End of script ####
