##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 1th of October 2019
# Modified July 2020

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script is to produce the Figure 4 
# for the paper: 

# Benitez-Lopez et al.The island rule explains consistent patterns of 
# body size evolution across terrestrial vertebrates. 


##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)
library(ggplot2)
library(ggpubr)


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
metamam4 <- readRDS(file = "Data/Final data/metamam4.Rdata")
metabird4 <- readRDS(file = "Data/Final data/metabird4.Rdata")
metarept4 <- readRDS(file = "Data/Final data/metarept4.Rdata")
metaamph4 <- readRDS(file = "Data/Final data/metaamph4.Rdata")

# load necessary functions
source("Scripts/000_Functions.R")

#Island area and remoteness##### 
#MAMMALS####
logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

# island isolation
small_island<-quantile(mamdata$Island_km2, prob = 0.1, names =FALSE) #4 km2
large_distance<-quantile(mamdata$Dist_near_mainland, prob = 0.9,names =FALSE) #150 km
large_island<-quantile(mamdata$Island_km2, prob = 0.9, names =FALSE) #32900
small_distance<-quantile(mamdata$Dist_near_mainland, prob = 0.1, names =FALSE) #1.5km

l<-predict(metamam4, newmods = cbind(logmass, large_distance, small_island,
                                  logmass*large_distance, logmass*small_island), addx=TRUE) 
h<-predict(metamam4, newmods = cbind(logmass,small_distance, large_island,
                                  logmass*small_distance, logmass*large_island), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Small and remote islands", nrow(l)) 
h<-data.frame(h)
h$Islandtype<-rep("Close and large islands", nrow(h))

df<-rbind(l,h)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Close and large islands" = "#9966FF","Small and remote islands"="#E69F00") #  orange palette

# import silhouette
# raster format
# sil_M <- readPNG("Silhouettes/PhyloPic.72f2f997.Steven-Traver.Cervus-elaphus.png")

Mc<-ggplot(mamdata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  # annotation_custom(rasterGrob(sil_M,
  #                              x = unit(0.14, "npc"),
  #                              y = unit(0.15, "npc"),
  #                              width = unit(0.22,"npc"),
  #                              height = unit(0.27,"npc")),
  #                             -Inf, Inf, -Inf, Inf) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ 
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), axis.text=element_text(size=18, colour = "black"), legend.position = c(0.65,0.82), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+ 
   guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+
  ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ 
  scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  labs(tag = "a")
Mc

# BIRDS ####
logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

# island isolation
small_island<-quantile(birddata$Island_km2, prob = 0.1, names =FALSE) #4 km2
large_distance<-quantile(birddata$Dist_near_mainland, prob = 0.9,names =FALSE) #150 km
large_island<-quantile(birddata$Island_km2, prob = 0.9, names =FALSE) #32900
small_distance<-quantile(birddata$Dist_near_mainland, prob = 0.1, names =FALSE) #1.5km

l<-predict(metabird4, newmods = cbind(logmass, large_distance, small_island,
                                      logmass*large_distance, logmass*small_island), addx=TRUE) 
h<-predict(metabird4, newmods = cbind(logmass,small_distance, large_island,
                                      logmass*small_distance, logmass*large_island), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Small and remote islands", nrow(l)) 
h<-data.frame(h)
h$Islandtype<-rep("Close and large islands", nrow(h))

df<-rbind(l,h)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Close and large islands" = "#9966FF","Small and remote islands"="#E69F00") #  orange palette

# import silhouette
# raster format
# sil_B<- readPNG("Silhouettes/PhyloPic.67a9ecfd.Sylviidae_Sylvioidea_PublicDom1.0_flipped.png")

Bc<-ggplot(birddata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  # annotation_custom(rasterGrob(sil_B,
  #                              x = unit(0.14, "npc"),
  #                              y = unit(0.13, "npc"),
  #                              width = unit(0.17,"npc"),
  #                              height = unit(0.16,"npc")),
  #                             -Inf, Inf, -Inf, Inf) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ 
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), axis.text=element_text(size=18, colour = "black"), legend.position = "none")+ 
  # guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+
  ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ 
  scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  labs(tag = "b") #+Qmtext
Bc

# REPTILES####
logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

# create variables
small_island<-quantile(reptdata$Island_km2, prob = 0.1, names =FALSE) #4 km2
large_distance<-quantile(reptdata$Dist_near_mainland, prob = 0.9,names =FALSE) #150 km
large_island<-quantile(reptdata$Island_km2, prob = 0.9, names =FALSE) #32900
small_distance<-quantile(reptdata$Dist_near_mainland, prob = 0.1, names =FALSE) #1.5km

l<-predict(metarept4, newmods = cbind(logmass, large_distance, small_island,
                                      logmass*large_distance, logmass*small_island), addx=TRUE) 
h<-predict(metarept4, newmods = cbind(logmass,small_distance, large_island,
                                      logmass*small_distance, logmass*large_island), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Small and remote islands", nrow(l)) 
h<-data.frame(h)
h$Islandtype<-rep("Close and large islands", nrow(h))

df<-rbind(l,h)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Close and large islands" = "#9966FF","Small and remote islands"="#E69F00") #  orange palette

# import silhouette
# raster format
# sil_R<- readPNG("Silhouettes/Steven-Traver.Varanus_Varanus.png")

Rc<-ggplot(reptdata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  # annotation_custom(rasterGrob(sil_R,
  #                              x = unit(0.16, "npc"),
  #                              y = unit(0.14, "npc"),
  #                              width = unit(0.26,"npc"),
  #                              height = unit(0.16,"npc")),
  #                             -Inf, Inf, -Inf, Inf) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ 
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), axis.text=element_text(size=18, colour = "black"), legend.position = "none")+ 
  # guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+
  ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ 
  scale_y_continuous(breaks=seq(-2.4,2.4, by=0.8), limits= c(-2.4,2.4), labels = scales::number_format(accuracy = 0.1))+
  labs(tag = "c") #+ Qmtext
Rc

# AMPHIBIANS ####
logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

# create variables
small_island<-quantile(amphdata$Island_km2, prob = 0.1, names =FALSE) #4 km2
large_distance<-quantile(amphdata$Dist_near_mainland, prob = 0.9,names =FALSE) #150 km
large_island<-quantile(amphdata$Island_km2, prob = 0.9, names =FALSE) #32900
small_distance<-quantile(amphdata$Dist_near_mainland, prob = 0.1, names =FALSE) #1.5km

l<-predict(metaamph4, newmods = cbind(logmass, large_distance, small_island,
                                      logmass*large_distance, logmass*small_island), addx=TRUE) 
h<-predict(metaamph4, newmods = cbind(logmass,small_distance, large_island,
                                      logmass*small_distance, logmass*large_island), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Small and remote islands", nrow(l)) 
h<-data.frame(h)
h$Islandtype<-rep("Close and large islands", nrow(h))

df<-rbind(l,h)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Close and large islands" = "#9966FF","Small and remote islands"="#E69F00") #  orange palette

# import silhouette
# raster format
# sil_A<- readPNG("Silhouettes/Will-Booker.Hyla-versicolor_CC0.1.0.png")

Ac<-ggplot(amphdata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  # annotation_custom(rasterGrob(sil_A,
  #                              x = unit(0.14, "npc"),
  #                              y = unit(0.15, "npc"),
  #                              width = unit(0.16,"npc"),
  #                              height = unit(0.18,"npc")),
  #                               -Inf, Inf, -Inf, Inf) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ 
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), axis.text=element_text(size=18, colour = "black"), 
        legend.position = "none" )+ 
  # guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+
  ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ 
  #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  labs(tag = "d") #+ Qmtext
Ac

####FIGURE 4 ####
fig4<-ggarrange(Mc,Bc,Rc,Ac, ncol = 2, nrow = 2, align = "v")
# tiff('Results/Figures/Figure4.tif', res=300, width=3100, height=3000)
# fig4
# dev.off()

pdf('Results/Figures/Figure4.pdf', width=10, height=9.5)
fig4
dev.off()

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/figure_4_R_session.txt")
sessionInfo()
sink()


### End of script ####
