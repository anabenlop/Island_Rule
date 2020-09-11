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
# rm(list=ls())

##############################################################
# Importing datasets and functions                        ####
##############################################################

#Load data
mamdata<-read.csv("~/New projects/Island rule/Data/Final data/mamdata_ph.csv", header = TRUE, stringsAsFactors = FALSE)
birddata<-read.csv("~/New projects/Island rule/Data/Final data/birddata_ph.csv", header = TRUE, stringsAsFactors = FALSE)
reptdata<-read.csv("~/New projects/Island rule/Data/Final data/reptdata_ph.csv", header = TRUE, stringsAsFactors = FALSE)
amphdata<-read.csv("~/New projects/Island rule/Data/Final data/amphdata_ph.csv", header = TRUE, stringsAsFactors = FALSE)

# Load models
meta4 <- readRDS(file = "~/New projects/Island rule/Data/Final data/meta4.Rdata")
metabird4 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metabird4.Rdata")
metarept4 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metarept4.Rdata")
metaamph4 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metaamph4.Rdata")

# load necessary functions
source("~/New projects/Island rule/Scripts/000_Functions.R")

#Island area and remoteness##### 
#MAMMALS####
logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

#extract Qm
prednames<-row.names(meta4$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(meta4, btt = 2) 
test_2pred<-anova(meta4, btt = 3) 
test_3pred<-anova(meta4, btt = 4) 
test_int<-anova(meta4, btt = 5) 
test_int2<-anova(meta4, btt = 6) 
test_int3<-anova(meta4, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

# island isolation
small_island<-quantile(mamdata$Island_km2, prob = 0.1, names =FALSE) #4 km2
large_distance<-quantile(mamdata$Dist_near_mainland, prob = 0.9,names =FALSE) #150 km
large_island<-quantile(mamdata$Island_km2, prob = 0.9, names =FALSE) #32900
small_distance<-quantile(mamdata$Dist_near_mainland, prob = 0.1, names =FALSE) #1.5km

l<-predict(meta4, newmods = cbind(logmass, large_distance, small_island,
                                  logmass*large_distance, logmass*small_island), addx=TRUE) 
h<-predict(meta4, newmods = cbind(logmass,small_distance, large_island,
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

Qmtext<-annotate(geom="text", x=1.6, y= -0.8, label= Qm, size = 5)

Mc<-ggplot(mamdata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ 
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), axis.text=element_text(size=18, colour = "black"), legend.position = c(0.68,0.82), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+ 
   guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+
  ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ 
  scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  labs(tag = "a")#+Qmtext
Mc

# BIRDS ####
logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

#extract Qm
prednames<-row.names(metabird4$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(metabird4, btt = 2) 
test_2pred<-anova(metabird4, btt = 3) 
test_3pred<-anova(metabird4, btt = 4) 
test_int<-anova(metabird4, btt = 5) 
test_int2<-anova(metabird4, btt = 6) 
test_int3<-anova(metabird4, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1], metabird4$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2], metabird4$QM))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

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

Qmtext<-annotate(geom="text", x=1.6, y= -0.8, label= Qm, size = 5)

Bc<-ggplot(birddata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ 
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), axis.text=element_text(size=18, colour = "black"), legend.position = c(0.68,0.82), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+ 
  guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+
  ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ 
  scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  labs(tag = "b") #+Qmtext
Bc

# REPTILES####
logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

#extract Qm
prednames<-row.names(metarept4$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(metarept4, btt = 2) 
test_2pred<-anova(metarept4, btt = 3) 
test_3pred<-anova(metarept4, btt = 4) 
test_int<-anova(metarept4, btt = 5) 
test_int2<-anova(metarept4, btt = 6) 
test_int3<-anova(metarept4, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

# write.csv(Qm_tot, "~/New projects/Island rule/Results/Reptiles/Coef/anova_dist_area.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 5)


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

Rc<-ggplot(reptdata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ 
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), axis.text=element_text(size=18, colour = "black"), legend.position = c(0.68,0.82), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+ 
  guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+
  ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ 
  scale_y_continuous(breaks=seq(-2.4,2.4, by=0.8), limits= c(-2.4,2.4), labels = scales::number_format(accuracy = 0.1))+
  labs(tag = "c") #+ Qmtext
Rc

# AMPHIBIANS ####
logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

#extract Qm
prednames<-row.names(metaamph4$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(metaamph4, btt = 2) 
test_2pred<-anova(metaamph4, btt = 3) 
test_3pred<-anova(metaamph4, btt = 4) 
test_int<-anova(metaamph4, btt = 5) 
test_int2<-anova(metaamph4, btt = 6) 
test_int3<-anova(metaamph4, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1], metaamph4$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2], metaamph4$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

# write.csv(Qm_tot, "~/New projects/Island rule/Results/Amphibians/Coef/anova_dist_area.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.6, y= -1.6, label= Qm, size = 5)

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

Ac<-ggplot(amphdata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ 
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), axis.text=element_text(size=18, colour = "black"), legend.position = c(0.68,0.82), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+ 
  guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+
  ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ 
  #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  labs(tag = "d") #+ Qmtext
Ac

####FIGURE 4 ####
fig4<-ggarrange(Mc,Bc,Rc,Ac, ncol = 2, nrow = 2, align = "v")
tiff('~/New projects/Island rule/Results/Figures/Figure4.tif', res=300, width=3100, height=3000)
fig4
dev.off()

### End of script ####
