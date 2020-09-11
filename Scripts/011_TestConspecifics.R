##############################################################
# Authors: 
# Ana Benítez-López (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 20th of July 2020

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script is to run alternative models to test the island rule
# for the paper: 

# Benítez-López et al. The island rule explains consistent patterns of 
# body size evolution across terrestrial vertebrates. 


##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)
library(ggplot2)
library(ggpubr)
library(lme4)
library(jtools)
library(effects)
# library(phylo_gmm)
library(phytools)

#clean memory
# rm(list=ls())

##############################################################
# Importing datasets and functions                        ####
##############################################################

#Load data
mamdata<-read.csv("~/New projects/Island rule/Data/Final data/mamdata_ph.csv", header = TRUE, stringsAsFactors = FALSE)
birddata<-read.csv("~/New projects/Island rule/Data/Final data/birddata_ph.csv", header = TRUE, stringsAsFactors = FALSE)
reptdata<-read.csv("~/New projects/Island rule/Data/Final data/reptdata_ph.csv", header = TRUE, stringsAsFactors = FALSE)
amphdata<-read.csv("~/New projects/Island rule/Data/Final data/amphdata_ph.csv", header = TRUE, stringsAsFactors = FALSE)

# loading phylogenetic matrixes 
load("~/New projects/Island rule/Data/Final data/mam_phylo_cor.Rdata") #mam_phylo_cor
load("~/New projects/Island rule/Data/Final data/bird_phylo_cor.Rdata") #bird_phylo_cor
load("~/New projects/Island rule/Data/Final data/rept_phylo_cor.Rdata") #rept_phylo_cor
load("~/New projects/Island rule/Data/Final data/amph_phylo_cor.Rdata") #amph_phylo_cor

# load functions
# source(file= "~/New projects/Island rule/Scripts/phyloglmm-v1.0.0/wzmli-phyloglmm-e7aee53/lme4_phylo_setup.R")
source("~/New projects/Island rule/Scripts/000_Functions.R")



### Run models with permutations