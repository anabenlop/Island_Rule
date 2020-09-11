##############################################################
# Authors: 
# Ana Benítez-López (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 1th of October 2019
# Last modification 4th of March 2020

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script is to run naive models that do not account for intrapopulation variability and/or phylogeny 
# for the paper: 

# Benitez-Lopez et al.The island rule explains consistent patterns of 
# body size evolution across terrestrial vertebrates. 

##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)
library(ggplot2)

#clean memory
rm(list=ls())

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

# load necessary functions
source("~/New projects/Island rule/Scripts/000_Functions.R")

####FIT MODELS####
# metaregressions
RE = list(~ 1 | Reference,~1|ID, ~1|SPID, ~1| Binomial)

#mammals####
phylocor<-list(Binomial= mam_phylo_cor)

# rma.mv
metamam2<-rma.mv(RR~logmass,V=1,  data=mamdata,
                  random= RE, R= phylocor, control=list(optimizer="bobyqa")) 
summary(metamam2)
saveRDS(metamam2, file = "~/New projects/Island rule/Data/Final data/metamam2.Rdata")

logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

#predict for vector of mass
df_m<-predict(metamam2, newmods = cbind(logmass), addx=TRUE)
df_m<-data.frame(df_m)
df_m$logmass<-df_m$X.logmass

M2<-ggplot(mamdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RR), colour= "#0072B2",size = 4,shape=20, alpha=I(.3)) +scale_shape_identity()+
  geom_line(data=df_m,aes(logmass,pred),color="#0072B2", size = 1.2)+
  geom_ribbon(data=df_m, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#0072B2", alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18))+xlab("log10(mass mainland (g))")+
  scale_x_continuous(breaks=seq(0,6,1)) + 
  scale_y_continuous(breaks=seq(-1.6,1.6,0.4), limits= c(-1.6,1.6))+
  labs(tag = "a")

#birds####
phylocor<-list(Binomial= bird_phylo_cor)
metabird2<-rma.mv(RR~logmass,V=1,  data=birddata,
                  random= RE, R= phylocor) 
summary(metabird2)
saveRDS(metabird2, file = "~/New projects/Island rule/Data/Final data/metabird2.Rdata")

logmass <- seq(from = min(birddata$logmass), to =  max(birddata$logmass) , length.out = 1000)

#predict for vector of mass
df_b<-predict(metabird2, newmods = cbind(logmass), addx=TRUE)
df_b<-data.frame(df_b)
df_b$logmass<-df_b$X.logmass

B2<-ggplot(birddata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RR), size = 4,shape=20, color="#CC0000",alpha=I(.3)) +scale_shape_identity()+
  geom_line(data=df_b,aes(logmass,pred),color="#CC0000", size = 1.2)+
  geom_ribbon(data=df_b, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#CC0000", alpha=I(.4))+
  theme(element_blank(),axis.text=element_text(size=18))+ xlab("log10(mass mainland (g))") + 
  scale_x_continuous(breaks=seq(0.6,3.7,0.8), limits = c(0.6,3.7))+ 
  scale_y_continuous(breaks=seq(-1.2,1.2,0.3), limits= c(-1.2,1.2))+
  labs(tag = "b")

#reptiles####
phylocor<-list(Binomial= rept_phylo_cor)
metarept2<-rma.mv(RR~logmass,V=1,  data=reptdata,  random= RE, R= phylocor) 
summary(metarept2)
saveRDS(metarept2, file = "~/New projects/Island rule/Data/Final data/metarept2.Rdata")

logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

#predict for vector of mass
df_r<-predict(metarept2, newmods = cbind(logmass), addx=TRUE)
df_r<-data.frame(df_r)
df_r$logmass<-df_r$X.logmass

R2<-ggplot(reptdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RR),colour="#E69F00", size = 4, shape = 20, alpha=I(.3)) + 
  geom_line(data=df_r,aes(logmass,pred),color="#E69F00", size = 1.2)+
  geom_ribbon(data=df_r, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#E69F00", alpha=I(.4))+
  theme(element_blank(),axis.text=element_text(size=18))+ xlab("log10(mass mainland (g))") +
  #scale_x_continuous(breaks=seq(-1,4,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, 0.6), limits = c(-2.4,2.4)) +
  labs(tag = "c")

#amphibians####
phylocor<-list(Binomial=amph_phylo_cor)

metaamph2<-rma.mv(RR~logmass ,V=1,  data=amphdata, random=RE, R = phylocor)
summary(metaamph2)

logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

#predict for vector of mass
df_a<-predict(metaamph2, newmods = cbind(logmass), addx=TRUE)
df_a<-data.frame(df_a)
df_a$logmass<-df_a$X.logmass

A2<-ggplot(amphdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RR),colour="#009E73", size = 4,shape = 20, alpha=I(.3)) + 
  geom_line(data=df_a,aes(logmass,pred),color="#009E73", size = 1.2)+
  geom_ribbon(data=df_a, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#009E73", alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18))+ xlab("log10(mass mainland (g))") + 
  scale_x_continuous(breaks=seq(-0.5,1.7,0.5), limits = c(-0.5,1.7))+ 
  scale_y_continuous(breaks=seq(-1.6,1.6,0.4), limits =c(-1.6,1.6)) + 
  labs(tag = "d")

tiff('~/New projects/Island rule/Results/Figures/Figure_nosamplingvar.tif', res=300, width=3100, height=3000)
multiplot(M2,R2,B2,A2, cols=2)
dev.off()
