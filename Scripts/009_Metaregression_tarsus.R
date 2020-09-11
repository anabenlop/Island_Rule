##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 21st of April 2020

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script runs the phylogenetic meta-regression for birds using
# tarsus length as morphometric measure for the paper: 

# Benitez-Lopez et al.The island rule explains consistent patterns of 
# body size evolution across terrestrial vertebrates. 


##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)
library(ggplot2)

#clean memory
# rm(list=ls())

##############################################################
# Importing datasets and functions                        ####
##############################################################

#Load data
birddata<-read.csv("~/New projects/Island rule/Data/Final data/birddata_ph2.csv", header = TRUE)

# loading phylogenetic matrixes 
load("~/New projects/Island rule/Data/Final data/bird_phylo_cor2.Rdata") #bird_phylo_cor

# load necessary functions
source("~/New projects/Island rule/Scripts/000_Functions.R")

# Calculate varcovar matrix
Vbird<- bldiag(lapply(split(birddata, birddata$CommonControl), calc.v))
is.positive.definite(Vbird) # FALSE
Vbird<-PDfunc(Vbird)
is.positive.definite(Vbird) # TRUE

#ISLAND RULE PLOT####
RE = list(~ 1 | Reference,~1|ID, ~1|SPID, ~1| Binomial)

#birds#### 
phylocor<-list(Binomial= bird_phylo_cor)
metabird<-rma.mv(RR~logmass,V=Vbird,  data=birddata, random= RE, 
                 R = phylocor#,  control=list(optimizer= "bobyqa")
) 

summary(metabird)
mR2.func(metabird)

logmass <- seq(from = min(birddata$logmass), to =  max(birddata$logmass) , length.out = 1000)

#predict for vector of mass
df_b<-predict(metabird, newmods = cbind(logmass), addx=TRUE)
df_b<-data.frame(df_b)
df_b$logmass<-df_b$X.logmass

#calculate size of points
wi    <- 1/sqrt(birddata$var)
size  <- 1 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

B<-ggplot(birddata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RR), size = size,shape=20, color="#CC0000",alpha=I(.3)) +scale_shape_identity()+
  geom_line(data=df_b,aes(logmass,pred),color="#CC0000", size = 1.2)+
  geom_ribbon(data=df_b, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#CC0000", alpha=I(.4))+
  theme(element_blank(),axis.text=element_text(size=18))+ xlab("log10(mass mainland (g))") + 
  scale_x_continuous(breaks=seq(0.6,3.7,0.8), limits = c(0.6,3.7))+ 
  scale_y_continuous(breaks=seq(-1.2,1.2,0.3), limits= c(-1.2,1.2))

tiff('~/New projects/Island rule/Results/Figures/Figure_Birds.tif', res=300, width=1550, height=1500)
B
dev.off()

