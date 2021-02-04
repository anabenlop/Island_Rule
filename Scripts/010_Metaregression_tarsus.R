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
 rm(list=ls())

##############################################################
# Importing datasets and functions                        ####
##############################################################

#Load data
birddata<-read.csv("Data/Final data/birddata_ph_tarsus.csv", header = TRUE)

# loading phylogenetic matrixes 
load("Data/Final data/bird_phylo_cor_tarsus.Rdata") #bird_phylo_cor

# load necessary functions
source("Scripts/000_Functions.R")

# Calculate varcovar matrix
Vbird<- bldiag(lapply(split(birddata, birddata$CommonControl), calc.v))
is.positive.definite(Vbird) # FALSE
Vbird<-PDfunc(Vbird)
is.positive.definite(Vbird) # TRUE

#ISLAND RULE PLOT####
RE = list(~ 1 | Reference,~1|ID, ~1|SPID, ~1| Binomial)

#birds#### 
phylocor<-list(Binomial= bird_phylo_cor)
metabird_tarsus<-rma.mv(RR~logmass,V=Vbird,  data=birddata, random= RE, 
                 R = phylocor) 
summary(metabird_tarsus)
mR2.func(metabird_tarsus)

logmass <- seq(from = min(birddata$logmass), to =  max(birddata$logmass) , length.out = 1000)

#predict for vector of mass
df_b<-predict(metabird_tarsus, newmods = cbind(logmass), addx=TRUE)
df_b<-data.frame(df_b)
df_b$logmass<-df_b$X.logmass

#calculate size of points
wi    <- 1/sqrt(birddata$var)
size  <- 1 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

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
  theme(element_blank(),axis.text=element_text(size=18))+ xlab("log10(mass mainland (g))") + ylab("lnRR")+
  scale_x_continuous(breaks=seq(0.6,3.7,0.8), limits = c(0.6,3.7))+ 
  scale_y_continuous(breaks=seq(-1.2,1.2,0.3), limits= c(-1.2,1.2))

tiff('Results/Figures/Figure S2.tif', res=300, width=1550, height=1500)
B
dev.off()

pdf('Results/Figures/Figure S2.pdf', width=5, height=4.75)
B
dev.off()

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/bird_meta-regression_tarsus_R_session.txt")
sessionInfo()
sink()

# End of script ####

