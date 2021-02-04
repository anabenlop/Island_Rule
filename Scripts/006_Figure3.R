##############################################################
# Authors: 
# Ana Benitez-López (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 27th of January 2021

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script is used to produce the Figure 3 
# in the paper: 

# Benitez-Lopez et al.The island rule explains consistent patterns of 
# body size evolution in terrestrial vertebrates

##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(png)
library(ggimage)
library(grid)
library(rsvg)
library(grImport2)

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

# Mammals ####
logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

#predict for vector of mass
df_m<-predict(metamam, newmods = cbind(logmass), addx=TRUE)
df_m<-data.frame(df_m)
df_m$logmass<-df_m$X.logmass

#calculate size of points based on sampling variance only
wi    <- 1/sqrt(mamdata$var)
size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

# import silhouette
# raster format
sil_M <- readPNG("Silhouettes/PhyloPic.72f2f997.Steven-Traver.Cervus-elaphus.png")

M<-ggplot(mamdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  annotation_custom(rasterGrob(sil_M,
                               x = unit(0.14, "npc"),
                               y = unit(0.15, "npc"),
                               width = unit(0.22,"npc"),
                               height = unit(0.27,"npc")),
                              -Inf, Inf, -Inf, Inf) +
  # annotation_custom(pictureGrob(svg_M,
  #                               x = unit(0.15, "npc"), 
  #                               y = unit(0.15, "npc"),
  #                               width = unit(0.18,"npc"),
  #                               height = unit(0.2,"npc")), 
  #                               -Inf, Inf, -Inf, Inf) +
  geom_point(aes(logmass,RR), colour= "#0072B2",size = size,shape=20, alpha=I(.3)) +scale_shape_identity()+
  geom_line(data=df_m,aes(logmass,pred),color="#0072B2", size = 1.2)+
  geom_ribbon(data=df_m, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#0072B2", alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+xlab("log10(mass mainland (g))")+ 
  ylab("lnRR")+ 
  scale_x_continuous(breaks=seq(0,6,1)) + 
  scale_y_continuous(breaks=seq(-1.6,1.6,0.4), limits= c(-1.6,1.6))+ 
  labs(tag = "a")

#birds#### 
logmass <- seq(from = min(birddata$logmass), to =  max(birddata$logmass) , length.out = 1000)

#predict for vector of mass
df_b<-predict(metabird, newmods = cbind(logmass), addx=TRUE)
df_b<-data.frame(df_b)
df_b$logmass<-df_b$X.logmass

#calculate size of points
wi    <- 1/sqrt(birddata$var)
size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

# import silhouette
# raster format
sil_B<- readPNG("Silhouettes/PhyloPic.67a9ecfd.Sylviidae_Sylvioidea_PublicDom1.0_flipped.png")

B<-ggplot(birddata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  annotation_custom(rasterGrob(sil_B,
                               x = unit(0.14, "npc"),
                               y = unit(0.13, "npc"),
                               width = unit(0.17,"npc"),
                               height = unit(0.16,"npc")),
                               -Inf, Inf, -Inf, Inf) +
  geom_point(aes(logmass,RR), size = size,shape=20, color="#FF5412",alpha=I(.3)) +scale_shape_identity()+
  geom_line(data=df_b,aes(logmass,pred),color="#FF5412", size = 1.2)+
  geom_ribbon(data=df_b, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#FF5412", alpha=I(.4))+
  theme(element_blank(),axis.text=element_text(size=18,colour ="black"))+ xlab("log10(mass mainland (g))") + ylab("lnRR")+ 
  scale_x_continuous(breaks=seq(0.6,3.7,0.8), limits = c(0.6,3.7))+ 
  scale_y_continuous(breaks=seq(-1.4,1.4,0.4), limits= c(-1.4,1.4))+
  #ylim(-1.3,1.3)+
  labs(tag = "b")

#reptiles####
logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

#predict for vector of mass
df_r<-predict(metarept, newmods = cbind(logmass), addx=TRUE)
df_r<-data.frame(df_r)
df_r$logmass<-df_r$X.logmass

#calculate size of points
wi    <- 1/sqrt(reptdata$var)
size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

# import silhouette
# raster format
sil_R<- readPNG("Silhouettes/Steven-Traver.Varanus_Varanus.png")

R<-ggplot(reptdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  annotation_custom(rasterGrob(sil_R,
                               x = unit(0.16, "npc"),
                               y = unit(0.14, "npc"),
                               width = unit(0.26,"npc"),
                               height = unit(0.16,"npc")),
                    -Inf, Inf, -Inf, Inf) +
  geom_point(aes(logmass,RR),colour="#E69F00", size = size, shape = 20, alpha=I(.3)) + 
  geom_line(data=df_r,aes(logmass,pred),color="#E69F00", size = 1.2)+
  geom_ribbon(data=df_r, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#E69F00", alpha=I(.4))+
  theme(element_blank(),axis.text=element_text(size=18, colour ="black"))+ xlab("log10(mass mainland (g))") + ylab("lnRR")+ 
  #scale_x_continuous(breaks=seq(-1,4,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, 0.6), limits = c(-2.4,2.4)) +
  #ylim(-2.3,2.3)+
  labs(tag = "c")

#amphibians####
logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

#predict for vector of mass
df_a<-predict(metaamph, newmods = cbind(logmass), addx=TRUE)
df_a<-data.frame(df_a)
df_a$logmass<-df_a$X.logmass

#calculate size of points
wi    <- 1/sqrt(amphdata$var)
size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

# import silhouette
# raster format
sil_A<- readPNG("Silhouettes/Will-Booker.Hyla-versicolor_CC0.1.0.png")

A<-ggplot(amphdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  annotation_custom(rasterGrob(sil_A,
                               x = unit(0.14, "npc"),
                               y = unit(0.15, "npc"),
                               width = unit(0.16,"npc"),
                               height = unit(0.18,"npc")),
                               -Inf, Inf, -Inf, Inf) +
  geom_point(aes(logmass,RR),colour="#009E73", size = size,shape = 20, alpha=I(.3)) + 
  geom_line(data=df_a,aes(logmass,pred),color="#009E73", size = 1.2)+
  geom_ribbon(data=df_a, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#009E73", alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+ xlab("log10(mass mainland (g))") + 
  ylab("lnRR")+
  # scale_x_continuous(breaks=seq(-0.6,1.7,0.5), limits = c(-0.6,1.7))+
  scale_y_continuous(breaks=seq(-1.6,1.6,0.4), limits =c(-1.6,1.6)) +
  labs(tag = "d")

#### FIGURE 3 ####
fig3<-ggarrange(M,B,R,A, ncol = 2, nrow = 2, align = "v")
# tiff('Results/Figures/Figure3.tif', res=300, width=3100, height=3000)
# fig3
# dev.off()

pdf('Results/Figures/Figure3.pdf', width=10, height=9.5)
fig3
dev.off()


# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/Figure3_R_session.txt")
sessionInfo()
sink()

# End of script ####
