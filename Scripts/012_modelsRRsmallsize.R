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

# This script runs the phylogenetic meta-regressions to produce the results 
# in the paper: 

# Benitez-Lopez et al.The island rule explains consistent patterns of 
# body size evolution across terrestrial vertebrates. 


##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)
library(ggplot2)
library(ggpubr)
library(tictoc)

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


#check if small sample size could be a problem
hist(mamdata$N_i + mamdata$N_m, breaks = 30)
summary(mamdata$N_i + mamdata$N_m, breaks = 30)
nrow(mamdata[(mamdata$N_i + mamdata$N_m<= 4),])

summary(mamdata[(mamdata$N_i + mamdata$N_m<= 4),]$Mean_m) #min 21.79
summary(mamdata[(mamdata$N_i + mamdata$N_m<= 4),]$Mean_i) #min 14.45

summary(birddata[(birddata$N_i + birddata$N_m<= 4),]$Mean_m) #min 6.648
summary(birddata[(birddata$N_i + birddata$N_m<= 4),]$Mean_i) #min 3.898

summary(reptdata[(reptdata$N_i + reptdata$N_m<= 4),]$Mean_m) #min 92.63
summary(reptdata[(reptdata$N_i + reptdata$N_m<= 4),]$Mean_i) #min 46.96

summary(amphdata[(amphdata$N_i + amphdata$N_m<= 16),]$Mean_m) #min 0.7276
summary(amphdata[(amphdata$N_i + amphdata$N_m<= 16),]$Mean_i) #min 0.7690

#probably just problematic for birds and maybe amphibians

# Calculate RR small and sampling var
mamdata$RRsm <- mamdata$RR + 0.5*(mamdata$sd_i^2/(mamdata$N_i*mamdata$Mean_i^2) - mamdata$sd_m^2/(mamdata$N_m*mamdata$Mean_m^2))
mamdata$var_sm <- mamdata$var + 0.5*(mamdata$sd_i^4/(mamdata$N_i*mamdata$Mean_i^4) - mamdata$sd_m^4/(mamdata$N_m*mamdata$Mean_m^4))

birddata$RRsm <- birddata$RR + 0.5*(birddata$sd_i^2/(birddata$N_i*birddata$Mean_i^2) - birddata$sd_m^2/(birddata$N_m*birddata$Mean_m^2))
birddata$var_sm <- birddata$var + 0.5*(birddata$sd_i^4/(birddata$N_i*birddata$Mean_i^4) - birddata$sd_m^4/(birddata$N_m*birddata$Mean_m^4))

reptdata$RRsm <- reptdata$RR + 0.5*(reptdata$sd_i^2/(reptdata$N_i*reptdata$Mean_i^2) - reptdata$sd_m^2/(reptdata$N_m*reptdata$Mean_m^2))
reptdata$var_sm <- reptdata$var + 0.5*(reptdata$sd_i^4/(reptdata$N_i*reptdata$Mean_i^4) - reptdata$sd_m^4/(reptdata$N_m*reptdata$Mean_m^4))

amphdata$RRsm <- amphdata$RR + 0.5*(amphdata$sd_i^2/(amphdata$N_i*amphdata$Mean_i^2) - amphdata$sd_m^2/(amphdata$N_m*amphdata$Mean_m^2))
amphdata$var_sm <- amphdata$var + 0.5*(amphdata$sd_i^4/(amphdata$N_i*amphdata$Mean_i^4) - amphdata$sd_m^4/(amphdata$N_m*amphdata$Mean_m^4))

# Check how RR and RRsm align
ggplot(mamdata) + geom_point(aes(RR,RRsm)) + geom_smooth(aes(RR,RRsm), method = "lm")
ggplot(birddata) + geom_point(aes(RR,RRsm)) + geom_smooth(aes(RR,RRsm), method = "lm")
ggplot(reptdata) + geom_point(aes(RR,RRsm)) + geom_smooth(aes(RR,RRsm), method = "lm")
ggplot(amphdata) + geom_point(aes(RR,RRsm)) + geom_smooth(aes(RR,RRsm), method = "lm")

# all of them pretty similar

# ggplot(mamdata) + geom_point(aes(var,var_sm)) + geom_smooth(aes(var, var_sm), method = "lm")
# ggplot(birddata) + geom_point(aes(var,var_sm)) + geom_smooth(aes(var, var_sm), method = "lm")
# ggplot(reptdata) + geom_point(aes(var,var_sm)) + geom_smooth(aes(var,var_sm), method = "lm")
# ggplot(amphdata) + geom_point(aes(var,var_sm)) + geom_smooth(aes(var,var_sm), method = "lm") 

#some variation for reptiles and amphibians

#rename varsm so that the common control covariance formula works
mamdata$var <- mamdata$var_sm
birddata$var <- birddata$var_sm
reptdata$var <- reptdata$var_sm
amphdata$var <- amphdata$var_sm

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
is.positive.definite(Vrept) # FALSE
Vrept<-PDfunc(Vrept)
is.positive.definite(Vrept) # TRUE

Vamph<- bldiag(lapply(split(amphdata, amphdata$CommonControl), calc.v))
is.positive.definite(Vamph) # TRUE

#################################################################################
# Testing insular size shifts: Island rule, correction for small sample size ####
#################################################################################
RE = list(~ 1 | Reference,~1|ID, ~1|SPID, ~1| Binomial)

#mammals####
phylocor<-list(Binomial= mam_phylo_cor)
metamam<-rma.mv(RRsm~logmass,V=Vmam, data=mamdata, random= RE,  R = phylocor)
summary(metamam)
mR2.func(metamam)
cR2.func(metamam)

logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

#predict for vector of mass
df_m<-predict(metamam, newmods = cbind(logmass), addx=TRUE)
df_m<-data.frame(df_m)
df_m$logmass<-df_m$X.logmass

#calculate size of points
wi    <- 1/sqrt(mamdata$var)
size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

M<-ggplot(mamdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RRsm), colour= "#0072B2",size = size,shape=20, alpha=I(.3)) +scale_shape_identity()+
  geom_line(data=df_m,aes(logmass,pred),color="#0072B2", size = 1.2)+
  geom_ribbon(data=df_m, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#0072B2", alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+xlab("log10(mass mainland (g))")+ ylab("lnRR(delta)")+ 
  scale_x_continuous(breaks=seq(0,6,1)) + 
  scale_y_continuous(breaks=seq(-1.6,1.6,0.4), limits= c(-1.6,1.6))+ 
  labs(tag = "a")

#birds#### 
phylocor<-list(Binomial= bird_phylo_cor)
metabird<-rma.mv(RRsm~logmass,V=Vbird,  data=birddata, random= RE, R = phylocor)
summary(metabird)
mR2.func(metabird)
cR2.func(metabird)

logmass <- seq(from = min(birddata$logmass), to =  max(birddata$logmass) , length.out = 1000)

#predict for vector of mass
df_b<-predict(metabird, newmods = cbind(logmass), addx=TRUE)
df_b<-data.frame(df_b)
df_b$logmass<-df_b$X.logmass

#calculate size of points
wi    <- 1/sqrt(birddata$var)
size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

B<-ggplot(birddata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RRsm), size = size,shape=20, color="#CC0000",alpha=I(.3)) +scale_shape_identity()+
  geom_line(data=df_b,aes(logmass,pred),color="#CC0000", size = 1.2)+
  geom_ribbon(data=df_b, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#CC0000", alpha=I(.4))+
  theme(element_blank(),axis.text=element_text(size=18,colour ="black"))+ xlab("log10(mass mainland (g))") + ylab("lnRR(delta)")+ 
  scale_x_continuous(breaks=seq(0.6,3.7,0.8), limits = c(0.6,3.7))+ 
  scale_y_continuous(breaks=seq(-1.4,1.4,0.4), limits= c(-1.4,1.4))+
  #ylim(-1.3,1.3)+
  labs(tag = "b")
B

#reptiles####
phylocor<-list(Binomial= rept_phylo_cor)
metarept<-rma.mv(RRsm~logmass,V=Vrept,  data=reptdata,  random= RE, R=phylocor) 
summary(metarept)
mR2.func(metarept)
cR2.func(metarept)

logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

#predict for vector of mass
df_r<-predict(metarept, newmods = cbind(logmass), addx=TRUE)
df_r<-data.frame(df_r)
df_r$logmass<-df_r$X.logmass

#calculate size of points
wi    <- 1/sqrt(reptdata$var)
size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

R<-ggplot(reptdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RRsm),colour="#E69F00", size = size, shape = 20, alpha=I(.3)) + 
  geom_line(data=df_r,aes(logmass,pred),color="#E69F00", size = 1.2)+
  geom_ribbon(data=df_r, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#E69F00", alpha=I(.4))+
  theme(element_blank(),axis.text=element_text(size=18, colour ="black"))+ xlab("log10(mass mainland (g))") + ylab("lnRR(delta)")+ 
  #scale_x_continuous(breaks=seq(-1,4,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, 0.6), limits = c(-2.4,2.4)) +
  #ylim(-2.3,2.3)+
  labs(tag = "c")

#amphibians####
phylocor<-list(Binomial= amph_phylo_cor)
metaamph<-rma.mv(RRsm~logmass ,V=Vamph,  data=amphdata, random= RE, R = phylocor) 
summary(metaamph)
mR2.func(metaamph)
cR2.func(metaamph)

logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

#predict for vector of mass
df_a<-predict(metaamph, newmods = cbind(logmass), addx=TRUE)
df_a<-data.frame(df_a)
df_a$logmass<-df_a$X.logmass

#calculate size of points
wi    <- 1/sqrt(amphdata$var)
size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

A<-ggplot(amphdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RRsm),colour="#009E73", size = size,shape = 20, alpha=I(.3)) + 
  geom_line(data=df_a,aes(logmass,pred),color="#009E73", size = 1.2)+
  geom_ribbon(data=df_a, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#009E73", alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+ xlab("log10(mass mainland (g))") + ylab("lnRR(delta)") +
  # scale_x_continuous(breaks=seq(-0.6,1.7,0.5), limits = c(-0.6,1.7))+
  scale_y_continuous(breaks=seq(-1.6,1.6,0.4), limits =c(-1.6,1.6)) +
  labs(tag = "d")

saveRDS(metamam, file = "Data/Final data/metamam_sm.Rdata")
saveRDS(metabird, file = "Data/Final data/metabird_sm.Rdata")
saveRDS(metarept, file = "Data/Final data/metarept_sm.Rdata")
saveRDS(metaamph, file = "Data/Final data/metaamph_sm.Rdata")

#### FIGURE 3 ####
multiplot<-ggarrange(M,B,R,A, ncol = 2, nrow = 2, align = "v")
tiff('Results/Figures/Figure_RRsmall.tif', res=300, width=3100, height=3000)
multiplot
dev.off()