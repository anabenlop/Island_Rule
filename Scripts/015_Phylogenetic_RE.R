##############################################################
# Authors: 
# Ana Benitez-López (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 15 of December 2020

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script is extracts and plots the phylogenetic RE
# for the paper: 

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
# var exp by phylogeny
mR2.func(metamam)
cR2.func(metamam)

RE_mam <- ranef(metamam)
RE_mam_sp <- RE_mam$Binomial
RE_mam_sp$Binomial <- row.names(RE_mam_sp)

RE_mam_temp <- left_join(RE_mam_sp, mamdata[,c("Binomial", "Order", "Family")], by = "Binomial")
RE_mam_sp <- RE_mam_temp[!duplicated(RE_mam_temp),]
RE_mam_sp$size <- ifelse(RE_mam_sp$intrcpt < 0, "Dwarf", "Giant")

RE_mam_sp$Order <- as.factor(RE_mam_sp$Order)

mamplot <- ggplot(RE_mam_sp) + geom_boxplot(aes(Order, intrcpt, colour = size)) +
  scale_colour_viridis_d() + theme_bw() + coord_flip() + ylab("lnRR") + 
  scale_x_discrete(limits = rev(levels(RE_mam_sp$Order))) + ylim(-0.2, 0.2) + labs(tag = "a")

# tiff('Results/Figures/Fig_mam_clade.tif', res=300, width=1500, height=2000)
# mamplot
# dev.off()

#birds#### 
# var exp by phylogeny
mR2.func(metabird)
cR2.func(metabird)

RE_bird <- ranef(metabird)
RE_bird_sp <- RE_bird$Binomial
RE_bird_sp$Binomial <- row.names(RE_bird_sp)

RE_bird_temp <- left_join(RE_bird_sp, birddata[,c("Binomial", "Order", "Family")], by = "Binomial")
RE_bird_sp <- RE_bird_temp[!duplicated(RE_bird_temp),]
RE_bird_sp$size <- ifelse(RE_bird_sp$intrcpt < 0, "Dwarf", "Giant")

RE_bird_sp$Order <- as.factor(RE_bird_sp$Order)

birdplot <- ggplot(RE_bird_sp) + geom_boxplot(aes(Order, intrcpt, colour = size)) +
  scale_colour_viridis_d() + theme_bw() + coord_flip() + ylab("lnRR") + 
  scale_x_discrete(limits = rev(levels(RE_bird_sp$Order))) + ylim(-0.2, 0.2) + labs(tag = "b")

# tiff('Results/Figures/Fig_bird_clade.tif', res=300, width=1500, height=2000)
# birdplot
# dev.off()

#reptiles####
# var exp by phylogeny
mR2.func(metarept)
cR2.func(metarept)

RE_rept <- ranef(metarept)
RE_rept_sp <- RE_rept$Binomial
RE_rept_sp$Binomial <- row.names(RE_rept_sp)

RE_rept_temp <- left_join(RE_rept_sp, reptdata[,c("Binomial", "Order", "Family")], by = "Binomial")
RE_rept_sp <- RE_rept_temp[!duplicated(RE_rept_temp),]
RE_rept_sp$size <- ifelse(RE_rept_sp$intrcpt < 0, "Dwarf", "Giant")

RE_rept_sp$Family <- as.factor(RE_rept_sp$Family)

reptplot <- ggplot(RE_rept_sp) + geom_boxplot(aes(Family, intrcpt, colour = size)) +
  scale_colour_viridis_d() + theme_bw() + coord_flip() + ylab("lnRR") + 
  scale_x_discrete(limits = rev(levels(RE_rept_sp$Family))) + ylim(-0.5, 0.5) + labs(tag = "c")

# tiff('Results/Figures/Fig_rept_clade.tif', res=300, width=1500, height=2000)
# reptplot
# dev.off()

#amphibians####
# var exp by phylogeny
mR2.func(metaamph)
cR2.func(metaamph)

RE_amph <- ranef(metaamph)
RE_amph_sp <- RE_amph$Binomial
RE_amph_sp$Binomial <- row.names(RE_amph_sp)

RE_amph_temp <- left_join(RE_amph_sp, amphdata[,c("Binomial", "Order", "Family")], by = "Binomial")
RE_amph_sp <- RE_amph_temp[!duplicated(RE_amph_temp),]
RE_amph_sp$size <- ifelse(RE_amph_sp$intrcpt < 0, "Dwarf", "Giant")

RE_amph_sp$Family <- as.factor(RE_amph_sp$Family)

amphplot <- ggplot(RE_amph_sp) + geom_boxplot(aes(Family, intrcpt, colour = size)) +
  scale_colour_viridis_d() + theme_bw() + coord_flip() + ylab("lnRR") + 
  scale_x_discrete(limits = rev(levels(RE_amph_sp$Family))) + ylim(-0.2, 0.2) + labs(tag = "d")

# tiff('Results/Figures/Fig_amph_clade.tif', res=300, width=1500, height=2000)
# amphplot
# dev.off()

#### FIGURE 3 ####
multiplot<-ggarrange(mamplot,birdplot,reptplot,amphplot, ncol = 2, nrow = 2, align = "v")
tiff('Results/Figures/ED_Fig_4.tif', res=300, width=3100, height=3000)
multiplot
dev.off()

pdf('Results/Figures/ED_Fig_4.pdf', width=10, height=9.5)
multiplot
dev.off()

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/phylogenetic_effects_R_session.txt")
sessionInfo()
sink()

# End of script ####
