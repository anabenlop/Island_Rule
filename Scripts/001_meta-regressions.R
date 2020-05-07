##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estacion Biologica de Doñana (EBD-CSIC, Spain) 
# Email: abenitez81@gmail.com

# Script first created on the 11th of November 2019
# Last modification April 2020

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


#clean memory
# rm(list=ls())

##############################################################
# Importing datasets and functions                        ####
##############################################################

#Load data
mamdata<-read.csv("~/GitHub/Island_Rule/Data/mamdata_ph.csv", header = TRUE)
birddata<-read.csv("~/GitHub/Island_Rule/Data//birddata_ph.csv", header = TRUE)
reptdata<-read.csv("~/GitHub/Island_Rule/Data/reptdata_ph.csv", header = TRUE)
amphdata<-read.csv("~/GitHub/Island_Rule/Data/amphdata_ph.csv", header = TRUE)

#Convert diet to Carn -Non Carn
mamdata$guild2<- ifelse(mamdata$guild == "Carn", "Carn", "No Carn")
mamdata$guild2<- factor(mamdata$guild2)


# loading phylogenetic matrixes 
load("~/GitHub/Island_Rule/Data/mam_phylo_cor.Rdata") #mam_phylo_cor
load("~/GitHub/Island_Rule/Data/bird_phylo_cor.Rdata") #bird_phylo_cor
load("~/GitHub/Island_Rule/Data/rept_phylo_cor.Rdata") #rept_phylo_cor
load("~/GitHub/Island_Rule/Data/amph_phylo_cor.Rdata") #amph_phylo_cor

# load necessary functions
source("~/GitHub/Island_Rule/Scripts/000_Functions.R")

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
is.positive.definite(Vrept) # TRUE

Vamph<- bldiag(lapply(split(amphdata, amphdata$CommonControl), calc.v))
is.positive.definite(Vamph) # TRUE

#ISLAND RULE ANALYSIS####
RE = list(~ 1 | Reference,~1|ID, ~1|SPID, ~1| Binomial)

#mammals####
phylocor<-list(Binomial= mam_phylo_cor)
metamam<-rma.mv(RR~logmass,V=Vmam,  data=mamdata, random= RE,  R = phylocor)
summary(metamam)
mR2.func(metamam)

logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

#predict for vector of mass
df_m<-predict(metamam, newmods = cbind(logmass), addx=TRUE)
df_m<-data.frame(df_m)
df_m$logmass<-df_m$X.logmass

#calculate size of points
wi    <- 1/sqrt(mamdata$var)
size  <- 3 + 30.0 * (wi - min(wi))/(max(wi) - min(wi))

M<-ggplot(mamdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RR), colour= "#0072B2",size = size,shape=20, alpha=I(.3)) +scale_shape_identity()+
  geom_line(data=df_m,aes(logmass,pred),color="#0072B2", size = 1.2)+
  geom_ribbon(data=df_m, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#0072B2", alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18))+xlab("log10(mass mainland (g))")+
  scale_x_continuous(breaks=seq(0,6,1)) + 
  scale_y_continuous(breaks=seq(-1.6,1.6,0.4), limits= c(-1.6,1.6))+
  labs(tag = "a")

#birds#### 
phylocor<-list(Binomial= bird_phylo_cor)
metabird<-rma.mv(RR~logmass,V=Vbird,  data=birddata, random= RE, 
                   R = phylocor) 
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
  scale_y_continuous(breaks=seq(-1.2,1.2,0.3), limits= c(-1.2,1.2))+
  labs(tag = "b")

#reptiles####
phylocor<-list(Binomial= rept_phylo_cor)
metarept<-rma.mv(RR~logmass,V=Vrept,  data=reptdata,  random= RE, R=phylocor) 
summary(metarept)
mR2.func(metarept)

logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

#predict for vector of mass
df_r<-predict(metarept, newmods = cbind(logmass), addx=TRUE)
df_r<-data.frame(df_r)
df_r$logmass<-df_r$X.logmass

#calculate size of points
wi    <- 1/sqrt(reptdata$var)
size  <- 3 + 30.0 * (wi - min(wi))/(max(wi) - min(wi))

R<-ggplot(reptdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RR),colour="#E69F00", size = size, shape = 20, alpha=I(.3)) + 
  geom_line(data=df_r,aes(logmass,pred),color="#E69F00", size = 1.2)+
  geom_ribbon(data=df_r, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#E69F00", alpha=I(.4))+
  theme(element_blank(),axis.text=element_text(size=18))+ xlab("log10(mass mainland (g))") +
  scale_x_continuous(breaks=seq(-1,4,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, 0.6), limits = c(-2.4,2.4)) +
  labs(tag = "c")

#amphibians####
phylocor<-list(Binomial= amph_phylo_cor)
metaamph<-rma.mv(RR~logmass ,V=var,  data=amphdata, random= RE, R = phylocor) 
summary(metaamph)
mR2.func(metaamph)

logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

#predict for vector of mass
df_a<-predict(metaamph, newmods = cbind(logmass), addx=TRUE)
df_a<-data.frame(df_a)
df_a$logmass<-df_a$X.logmass

#calculate size of points
wi    <- 1/sqrt(amphdata$var)
size  <- 0.5 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

A<-ggplot(amphdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RR),colour="#009E73", size = size,shape = 20, alpha=I(.3)) + 
  geom_line(data=df_a,aes(logmass,pred),color="#009E73", size = 1.2)+
  geom_ribbon(data=df_a, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#009E73", alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18))+ xlab("log10(mass mainland (g))") + 
  scale_x_continuous(breaks=seq(-0.5,1.7,0.5), limits = c(-0.5,1.7))+
  scale_y_continuous(breaks=seq(-1.6,1.6,0.4), limits =c(-1.6,1.6)) +
  labs(tag = "d")

#### FIGURE 3 ####
multiplot<-ggarrange(M,B,R,A, ncol = 2, nrow = 2, align = "v")
tiff('~/GitHub/Island_Rule/Figure3.tif', res=300, width=3100, height=3000)
multiplot
dev.off()

###ECOLOGICAL HYPOTHESIS ANALYSES#####

# MAMMALS####
phylocor<-list(Binomial=mam_phylo_cor)
logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

#island area ####
metamam2<-rma.mv(RR~logmass*Island_km2,V = Vmam, data=mamdata, random= RE,
              R = phylocor, control=list(optimizer= "bobyqa"),method = "REML") 
summary(metamam2) 

#extract Qm
prednames<-row.names(metamam2$beta)[-1]

test_1pred<-anova(metamam2, btt = 2) 
test_2pred<-anova(metamam2, btt = 3) 
test_int<-anova(metamam2, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))
Qmtext<-annotate(geom="text", x=1.6, y= -0.8, label= Qm, size = 6)

# Calculation R2
mR2.func(metamam2) 

# island area 
logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

Island_km2<-quantile(mamdata$Island_km2, prob = 0.1) #63 km2
s<-predict(metamam2, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

Island_km2<-quantile(mamdata$Island_km2, prob = 0.9) #147910.8
l<-predict(metamam2, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Small", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Large", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Large"="#9966FF","Small" = "#E69F00") #red yellow palette "Warm no pred"= "#FFCC33", "Cold no pred"= "#0072B2

Ma<-ggplot(mamdata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island type"),colour=guide_legend(title="Island type"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "a")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Mammals/Hyp/area.tif")
Ma
dev.off()

##distance####
metamam3<-rma.mv(RR~logmass*Dist_near_mainland, V=Vmam,  data=mamdata,  random= RE,
              R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metamam3) 

coef<-data.frame(b =metamam3$b, lci = metamam3$ci.lb, uci =  metamam3$ci.ub)

#extract Qm
prednames<-row.names(metamam3$beta)[-1]

test_1pred<-anova(metamam3, btt = 2) 
test_2pred<-anova(metamam3, btt = 3) 
test_int<-anova(metamam3, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))
Qmtext<-annotate(geom="text", x=1.6, y= -0.8, label= Qm, size = 6)

# Calculation R2
mR2.func(metamam3) 

logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

Dist_near_mainland<-quantile(mamdata$Dist_near_mainland, prob = 0.1) #63 km2
s<-predict(metamam3, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

Dist_near_mainland<-quantile(mamdata$Dist_near_mainland, prob = 0.9) #147910.8
l<-predict(metamam3, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Close", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Remote", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Remote"="#9966FF","Close" = "#E69F00") 

Mb<-ggplot(mamdata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island type"),colour=guide_legend(title="Island type"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "b")+ Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Mammals/Hyp/distance.tif")
Mb
dev.off()

##Island area and remoteness####
metamam4<-rma.mv(RR~logmass*Dist_near_mainland +logmass*Island_km2,V=Vmam,  data=mamdata, random= RE,
              R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metamam4)

coef<-data.frame(b =metamam4$b, lci = metamam4$ci.lb, uci =  metamam4$ci.ub)

#extract Qm
prednames<-row.names(metamam4$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(metamam4, btt = 2) 
test_2pred<-anova(metamam4, btt = 3) 
test_3pred<-anova(metamam4, btt = 4) 
test_int<-anova(metamam4, btt = 5) 
test_int2<-anova(metamam4, btt = 6) 
test_int3<-anova(metamam4, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 4))

# Calculation R2
mR2.func(metamam4) #12.3

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
l$Islandtype<-rep("Small remote", nrow(l)) 
h<-data.frame(h)
h$Islandtype<-rep("Close large", nrow(h))

df<-rbind(l,h)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Close large" = "#9966FF","Small remote"="#E69F00") #  orange palette

Qmtext<-annotate(geom="text", x=1.6, y= -0.8, label= Qm, size = 6)

Mc<-ggplot(mamdata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island type"),colour=guide_legend(title="Island type"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ 
  scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "c")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Mammals/Hyp/distance_area.tif")
Mc
dev.off()

# diet ####
metamam5<-rma.mv(RR~logmass*guild2,V=Vmam,  data=mamdata,  random= RE,  
              R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metamam5)

coef<-data.frame(b =metamam5$b, lci = metamam5$ci.lb, uci =  metamam5$ci.ub)

#extract Qm
prednames<-row.names(metamam5$beta)[-1]

test_1pred<-anova(metamam5, btt = 2) 
test_2pred<-anova(metamam5, btt = 3) 
test_int<-anova(metamam5, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))
Qmtext<-annotate(geom="text", x=1.6, y= -0.8, label= Qm, size = 6)

# Calculation R2
mR2.func(metamam5)

c<-predict(metamam5, newmods = cbind(logmass, 0, 0), addx=TRUE)
nc<-predict(metamam5, newmods = cbind(logmass, 1, logmass), addx=TRUE)

### merge data frames and plot all together
c<-data.frame(c)
c$Islandtype<-rep("Carn", nrow(c))
nc<-data.frame(nc)
nc$Islandtype<-rep("Non Carn", nrow(nc))

df<-rbind(c,nc) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Carn" = "#E69F00", "Non Carn" = "#9966FF")

Md<-ggplot(mamdata)+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values=colpalette) + scale_fill_manual(values=colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Diet"),colour=guide_legend(title="Diet"))+
  ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("")+ labs(tag = "d")+ Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Mammals/Hyp/diet.tif")
Md
dev.off()

# temperature####
metamam6<-rma.mv(RR~logmass*tmean,V=Vmam,  data=mamdata,random= RE,  
               R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metamam6) 

coef<-data.frame(b =metamam6$b, lci = metamam6$ci.lb, uci =  metamam6$ci.ub)

#extract Qm
prednames<-row.names(metamam6$beta)[-1]

test_1pred<-anova(metamam6, btt = 2) 
test_2pred<-anova(metamam6, btt = 3) 
test_int<-anova(metamam6, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))
Qmtext<-annotate(geom="text", x=1.6, y= -0.8, label= Qm, size = 6)

# Calculation R2
mR2.func(metamam6) 

#temperature
tmean<-quantile(mamdata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metamam6, newmods = cbind(logmass, tmean, logmass*tmean), addx=TRUE) 

tmean<-quantile(mamdata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metamam6, newmods = cbind(logmass, tmean,logmass*tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Me<-ggplot(mamdata)+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "e")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Mammals/Hyp/tmean.tif") #interactive effect
Me
dev.off()

# temperature intercept####
metamam6b<-rma.mv(RR~logmass + tmean,V=Vmam,  data=mamdata,random= RE,  
              R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metamam6b)

coef<-data.frame(b =metamam6b$b, lci = metamam6b$ci.lb, uci =  metamam6b$ci.ub)

#extract Qm
prednames<-row.names(metamam6b$beta)[-1]

test_1pred<-anova(metamam6b, btt = 2) 
test_2pred<-anova(metamam6b, btt = 3) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_2pred[1]), digits = 3), ", p-val = ", round(as.numeric(test_2pred[2]), digits = 3))
Qmtext<-annotate(geom="text", x=1.6, y= -0.8, label= Qm, size = 6)

# Calculation R2
mR2.func(metamam6b) 

#temperature
tmean<-quantile(mamdata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metamam6b, newmods = cbind(logmass, tmean), addx=TRUE) 

tmean<-quantile(mamdata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metamam6b, newmods = cbind(logmass, tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Me2<-ggplot(mamdata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ 
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "g")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Mammals/Hyp/tmean_int.tif") #interactive effect
Me2
dev.off()

# temp seas####
meta6h<-rma.mv(RR~logmass*tseas,V=Vmam,  data=mamdata,  random= RE,  
               R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(meta6h) 

coef<-data.frame(b =meta6h$b, lci = meta6h$ci.lb, uci =  meta6h$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Mammals/Coef/coef_tseas.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(meta6h$beta)[-1]

test_1pred<-anova(meta6h, btt = 2) 
test_2pred<-anova(meta6h, btt = 3) 
test_int<-anova(meta6h, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Mammals/Coef/anova_tseas.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))


# Calculation R2
mR2.func(meta6h) #12.72

tseas<-quantile(mamdata$tseas, prob = 0.1, names=FALSE) #1.7
l<-predict(meta6h, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

tseas<-quantile(mamdata$tseas, prob = 0.9, names=FALSE) #81.79
h<-predict(meta6h, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Qmtext<-annotate(geom="text", x=1.6, y= -0.8, label= Qm, size = 6)

Mh<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.76,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Temperature season."),colour=guide_legend(title="Temperature season."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  # ylim(-0.7,0.7)+
  ggtitle("")+ labs(tag = "h")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Mammals/Hyp/tseas.tif")
Mh
dev.off()

# resource availability####
meta5e<-rma.mv(RR~logmass*NDVI,V=Vmam,  data=mamdata,  random= RE,  
               R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(meta5e)

coef<-data.frame(b =meta5e$b, lci = meta5e$ci.lb, uci =  meta5e$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Mammals/Coef/coef_ndvi.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(meta5e$beta)[-1]

test_1pred<-anova(meta5e, btt = 2) 
test_2pred<-anova(meta5e, btt = 3) 
test_int<-anova(meta5e, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Mammals/Coef/anova_ndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(meta5e) #12.34

NDVI<-quantile(mamdata$NDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(meta5e, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

NDVI<-quantile(mamdata$NDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(meta5e, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Qmtext<-annotate(geom="text", x=1.6, y= -0.8, label= Qm, size = 6)

Me<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource availab."),colour=guide_legend(title="Resource availab."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("")+ labs(tag = "e")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Mammals/Hyp/ndvi.tif")
Me
dev.off()

# seasonality in resources####
meta5g<-rma.mv(RR~logmass*SDNDVI,V=Vmam,  data=mamdata,  random= RE,  
               R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(meta5g)

coef<-data.frame(b =meta5g$b, lci = meta5g$ci.lb, uci =  meta5g$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Mammals/Coef/coef_sdndvi.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(meta5g$beta)[-1]

test_1pred<-anova(meta5g, btt = 2) 
test_2pred<-anova(meta5g, btt = 3) 
test_int<-anova(meta5g, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Mammals/Coef/anova_sdndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))


# Calculation R2
mR2.func(meta5g) #11.11

SDNDVI<-quantile(mamdata$SDNDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(meta5g, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

SDNDVI<-quantile(mamdata$SDNDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(meta5g, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Qmtext<-annotate(geom="text", x=1.6, y= -0.8, label= Qm, size = 6)

Mf<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource season."),colour=guide_legend(title="Resource season."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("")+ labs(tag = "f")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Mammals/Hyp/sdndvi.tif")
Mf
dev.off()



# multiplot####
multimam<-ggarrange(Ma,Mb,Mc,Md1,Me,Mf,Mg2,Mh, Mi, ncol = 3, nrow = 3, align = "v")
tiff('~/GitHub/Island_Rule/Results/Figures/Figure_hyp_mam_upd1.tif', res=300, width=6000, height=6000, compression = "lzw")
multimam
dev.off()

# BIRDS ####
phylocor<-list(Binomial=bird_phylo_cor)
logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

#island area ####
metabird2<-rma.mv(RR~logmass*Island_km2,V = Vbird, data=birddata, random= RE,
              R = phylocor, control=list(optimizer= "bobyqa"),method = "REML") 
summary(metabird2) #this is the actual hypothesis, QM(df = 3) = 30.6821, p-val < .0001

coef<-data.frame(b =metabird2$b, lci = metabird2$ci.lb, uci =  metabird2$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Birds/Coef/coef_area.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metabird2$beta)[-1]

test_1pred<-anova(metabird2, btt = 2) 
test_2pred<-anova(metabird2, btt = 3) 
test_int<-anova(metabird2, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Birds/Coef/anova_area.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x=1.6, y= -1, label= Qm, size = 6)

# Calculation R2
mR2.func(metabird2) #7.91

# island area combined with predators

Island_km2<-quantile(birddata$Island_km2, prob = 0.1) #63 km2
s<-predict(metabird2, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

Island_km2<-quantile(birddata$Island_km2, prob = 0.9) #147910.8
l<-predict(metabird2, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Small", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Large", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

#colpalette<-c("#FFCC33","#FF6600","#CC0000")
colpalette<-c("Large"="#9966FF","Small" = "#E69F00") #red yellow palette "Warm no pred"= "#FFCC33", "Cold no pred"= "#0072B2

Ba<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island area"),colour=guide_legend(title="Island area"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1,1, by=0.2), limits= c(-1,1))+
  ggtitle("") + labs(tag = "a") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Birds/Hyp/area.tif")
Ba
dev.off()

##distance####
metabird3<-rma.mv(RR~logmass*Dist_near_mainland, V=Vbird,  data=birddata,  random= RE,
              R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metabird3) #QM(df = 3) = 22.2331, p-val < .0001

coef<-data.frame(b =metabird3$b, lci = metabird3$ci.lb, uci =  metabird3$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Birds/Coef/coef_distance.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metabird3$beta)[-1]

test_1pred<-anova(metabird3, btt = 2) 
test_2pred<-anova(metabird3, btt = 3) 
test_int<-anova(metabird3, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Birds/Coef/anova_dist.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x=1.6, y= -1, label= Qm, size = 6)

# Calculation R2
mR2.func(metabird3) #7.09

logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

Dist_near_mainland<-quantile(birddata$Dist_near_mainland, prob = 0.1) #63 km2
s<-predict(metabird3, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

Dist_near_mainland<-quantile(birddata$Dist_near_mainland, prob = 0.9) #147910.8
l<-predict(metabird3, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Near", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Far", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Far"="#9966FF","Near" = "#E69F00") # orange purple

Bb<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island isolation"),colour=guide_legend(title="Island isolation"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1,1, by=0.2), limits= c(-1,1))+
  ggtitle("") + labs(tag = "b") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Birds/Hyp/distance.tif")
Bb
dev.off()

##island area and distance####
metabird4<-rma.mv(RR~logmass*Dist_near_mainland +logmass*Island_km2,V=Vbird,  data=birddata, random= RE,
              R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metabird4)

coef<-data.frame(b =metabird4$b, lci = metabird4$ci.lb, uci =  metabird4$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Birds/Coef/coef_distance_area.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metabird4$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(metabird4, btt = 2) 
test_2pred<-anova(metabird4, btt = 3) 
test_3pred<-anova(metabird4, btt = 4) 
test_int<-anova(metabird4, btt = 5) 
test_int2<-anova(metabird4, btt = 6) 
test_int3<-anova(metabird4, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Birds/Coef/anova_dist_area.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 4))

Qmtext<-annotate(geom="text", x=1.6, y= -1, label= Qm, size = 6)

# Calculation R2
mR2.func(metabird4) #9.41

logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

# create variables
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
l$Islandtype<-rep("Small remote", nrow(l)) 
h<-data.frame(h)
h$Islandtype<-rep("Close large", nrow(h))

#df<-rbind(s,m,l)
df<-rbind(l,h)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

#colpalette<-c("#FFCC33","#FF6600","#CC0000")
colpalette<-c("Close large" = "#9966FF","Small remote"="#E69F00") #  orange palette

Bc<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ 
  scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1,1, by=0.2), limits= c(-1,1))+
  ggtitle("") + labs(tag = "c") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Birds/Hyp/distance_area.tif")
Bc
dev.off()

##predation pressure####
# richness
metabird7d<-rma.mv(RR~logmass*TotPredrichEndo,V=Vbird,  data=birddata,  random= RE,  
               R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metabird7d)

coef<-data.frame(b =metabird7d$b, lci = metabird7d$ci.lb, uci =  metabird7d$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Birds/Coef/coef_predrich.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metabird7d$beta)[-1]

test_1pred<-anova(metabird7d, btt = 2) 
test_2pred<-anova(metabird7d, btt = 3) 
test_int<-anova(metabird7d, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Birds/Coef/anova_predrich.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x=1.6, y= -1, label= Qm, size = 6)

# Calculation R2
mR2.func(meta7d) #11.93

# richness of predators
logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

TotPredrichEndo<-quantile(birddata$TotPredrichEndo, prob = 0.1) 
s<-predict(metabird7d, newmods = cbind(logmass,TotPredrichEndo, logmass*TotPredrichEndo), addx=TRUE)

TotPredrichEndo<-quantile(birddata$TotPredrichEndo, prob = 0.9) #147910.8
l<-predict(metabird7d, newmods = cbind(logmass,TotPredrichEndo, logmass*TotPredrichEndo), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Low", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("High", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("High"="#9966FF","Low" = "#E69F00") 

Bd1<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Predator richness"),colour=guide_legend(title="Predator richness"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1,1, by=0.2), limits= c(-1,1))+
  ggtitle("") + labs(tag = "d")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Birds/Hyp/predrich.tif")
Bd1
dev.off()

###mammalian richness
metabird7d2<-rma.mv(RR~logmass*MPredRichEndo,V=Vbird,  data=birddata,  random= RE,  
                R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metabird7d2)

coef<-data.frame(b =metabird7d2$b, lci = metabird7d2$ci.lb, uci =  metabird7d2$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Birds/Coef/coef_predrich_mam.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metabird7d2$beta)[-1]

test_1pred<-anova(metabird7d2, btt = 2) 
test_2pred<-anova(metabird7d2, btt = 3) 
test_int<-anova(metabird7d2, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Birds/Coef/anova_predrich_mam.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metabird7d2) #11.93

# richness of predators
logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

MPredRichEndo<-quantile(birddata$MPredRichEndo, prob = 0.1) 
s<-predict(metabird7d2, newmods = cbind(logmass,MPredRichEndo, logmass*MPredRichEndo), addx=TRUE)

MPredRichEndo<-quantile(birddata$MPredRichEndo, prob = 0.9) #147910.8
l<-predict(metabird7d2, newmods = cbind(logmass,MPredRichEndo, logmass*MPredRichEndo), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Low", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("High", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("High"="#9966FF","Low" = "#E69F00") 

Qmtext<-annotate(geom="text", x=1.6, y= -0.8, label= Qm, size = 6)

Bd2<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Predator richness"),colour=guide_legend(title="Predator richness"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "d")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Birds/Hyp/predrich_mam.tif")
Bd2
dev.off()

# presence
metabird7e<-rma.mv(RR~logmass*PredEndo,V=Vbird,  data=birddata,  random= RE,  
               R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metabird7e)

coef<-data.frame(b =metabird7e$b, lci = metabird7e$ci.lb, uci =  metabird7e$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Birds/Coef/coef_pred.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metabird7e$beta)[-1]

test_1pred<-anova(metabird7e, btt = 2) 
test_2pred<-anova(metabird7e, btt = 3) 
test_int<-anova(metabird7e, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Birds/Coef/anova_pred.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metabird7e) #6.42

PredEndo<-0
nv<-predict(metabird7e, newmods = cbind(logmass,PredEndo, logmass*PredEndo), addx=TRUE) 
PredEndo<-1
v<-predict(metabird7e, newmods = cbind(logmass,PredEndo,logmass*PredEndo), addx=TRUE) 

### merge data frames and plot all together
nv<-data.frame(nv)
nv$Islandtype<-rep("No predators", nrow(l)) 
v<-data.frame(v)
v$Islandtype<-rep("Predators", nrow(h))

df<-rbind(v,nv)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Predators"="#9966FF","No predators" = "#E69F00") 

Qmtext<-annotate(geom="text", x=1.6, y= -1, label= Qm, size = 6)

Bd3<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island type"),colour=guide_legend(title="Island type"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1,1, by=0.2), limits= c(-1,1))+
  ggtitle("") + labs(tag = "d")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Birds/Hyp/pred_pres.tif")
Bd3
dev.off() #makes no sense --> 21 zeroes, 685 ones

# resource availability####
metabird5e<-rma.mv(RR~logmass*NDVI,V=Vbird,  data=birddata,  random= RE,  
               R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metabird5e)

coef<-data.frame(b =metabird5e$b, lci = metabird5e$ci.lb, uci =  metabird5e$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Birds/Coef/coef_ndvi.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metabird5e$beta)[-1]

test_1pred<-anova(metabird5e, btt = 2) 
test_2pred<-anova(metabird5e, btt = 3) 
test_int<-anova(metabird5e, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Birds/Coef/anova_ndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x=1.6, y= -1, label= Qm, size = 6)

# Calculation R2
mR2.func(metabird5e) #6.3

NDVI<-quantile(birddata$NDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(metabird5e, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

NDVI<-quantile(birddata$NDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(metabird5e, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Be<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource availab."),colour=guide_legend(title="Resource availab."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1,1, by=0.2), limits= c(-1,1))+
  ggtitle("")+ labs(tag = "e") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Birds/Hyp/ndvi.tif")
Be
dev.off()

# seasonality in resources####
metabird5g<-rma.mv(RR~logmass*SDNDVI,V=Vbird,  data= birddata,  random= RE,  
               R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metabird5g)

coef<-data.frame(b =metabird5g$b, lci = metabird5g$ci.lb, uci =  metabird5g$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Birds/Coef/coef_sdndvi.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metabird5g$beta)[-1]

test_1pred<-anova(metabird5g, btt = 2) 
test_2pred<-anova(metabird5g, btt = 3) 
test_int<-anova(metabird5g, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Birds/Coef/anova_sdndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x=1.6, y= -1, label= Qm, size = 6)

# Calculation R2
mR2.func(metabird5g) #6.56

SDNDVI<-quantile(birddata$SDNDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(metabird5g, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

SDNDVI<-quantile(birddata$SDNDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(metabird5g, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Bf<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource season."),colour=guide_legend(title="Resource season."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1,1, by=0.2), limits= c(-1,1))+
  ggtitle("")+ labs(tag = "f") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Birds/Hyp/sdndvi.tif")
Bf
dev.off()

# island temperature ####
metabird6f<-rma.mv(RR~logmass*tmean,V=Vbird,  data=birddata,random= RE,  
               R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metabird6f) #QM(df = 3) = 34.1169, p-val < .0001

coef<-data.frame(b =metabird6f$b, lci = metabird6f$ci.lb, uci =  metabird6f$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Birds/Coef/coef_tmean.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metabird6f$beta)[-1]

test_1pred<-anova(metabird6f, btt = 2) 
test_2pred<-anova(metabird6f, btt = 3) 
test_int<-anova(metabird6f, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Birds/Coef/anova_tmean.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x=1.6, y= -1, label= Qm, size = 6)

# Calculation R2
mR2.func(metabird6f) #9.46

#temperature
tmean<-quantile(birddata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metabird6f, newmods = cbind(logmass, tmean, logmass*tmean), addx=TRUE) 

tmean<-quantile(birddata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metabird6f, newmods = cbind(logmass, tmean,logmass*tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Bg<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1,1, by=0.2), limits= c(-1,1))+
  ggtitle("") + labs(tag = "g") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Birds/Hyp/tmean.tif")
Bg
dev.off()

# temperature intercept####
metabird6<-rma.mv(RR~logmass + tmean,V=Vbird,  data=birddata,random= RE,  
              R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metabird6)

coef<-data.frame(b =metabird6$b, lci = metabird6$ci.lb, uci =  metabird6$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Birds/Coef/coef_tmean_int.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metabird6$beta)[-1]

test_1pred<-anova(metabird6, btt = 2) 
test_2pred<-anova(metabird6, btt = 3) 
# test_int<-anova(meta6, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Birds/Coef/anova_tmean_int.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_2pred[1]), digits = 3), ", p-val = ", round(as.numeric(test_2pred[2]), digits = 3))

Qmtext<-annotate(geom="text", x=1.6, y= -1, label= Qm, size = 6)

# Calculation R2
mR2.func(metabird6) #14.06

#temperature
tmean<-quantile(birddata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metabird6, newmods = cbind(logmass, tmean), addx=TRUE) 

tmean<-quantile(birddata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metabird6, newmods = cbind(logmass, tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Bg2<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1,1, by=0.2), limits= c(-1,1))+
  ggtitle("") + labs(tag = "g") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Birds/Hyp/tmean_intercept.tif") #interactive effect
Bg2
dev.off()

# temp seas####
metabird6h<-rma.mv(RR~logmass*tseas,V=Vbird,  data=birddata,  random= RE,  
               R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metabird6h) 

coef<-data.frame(b =metabird6h$b, lci = metabird6h$ci.lb, uci =  metabird6h$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Birds/Coef/coef_tseas.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metabird6h$beta)[-1]

test_1pred<-anova(metabird6h, btt = 2) 
test_2pred<-anova(metabird6h, btt = 3) 
test_int<-anova(metabird6h, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Birds/Coef/anova_tseas.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x=1.6, y= -1, label= Qm, size = 6)

# Calculation R2
mR2.func(metabird6h) #9.1

tseas<-quantile(birddata$tseas, prob = 0.1, names=FALSE) #1.68
l<-predict(metabird6h, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

tseas<-quantile(birddata$tseas, prob = 0.9, names=FALSE) #2.53
h<-predict(metabird6h, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Bh<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.76,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Temperature season."),colour=guide_legend(title="Temperature season."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1,1, by=0.2), limits= c(-1,1))+
  ggtitle("")+ labs(tag = "h") +Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Birds/Hyp/tseas.tif")
Bh
dev.off()

# diet ####
birddata[birddata$Species_main == "Passer italiae", "guild"] <- "PlantSeed"
birddata$guild2<- ifelse(birddata$guild == "VertFishScav", "Carn", "No Carn")
birddata$guild2<- factor(birddata$guild2)

metabird8<-rma.mv(RR~logmass*guild2 ,V=Vbird,  data=birddata,  random= RE,  
              R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metabird8)

coef<-data.frame(b =metabird8$b, lci = metabird8$ci.lb, uci =  metabird8$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Birds/Coef/coef_diet.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metabird8$beta)[-1]

test_1pred<-anova(metabird8, btt = 2) 
test_2pred<-anova(metabird8, btt = 3) 
test_int<-anova(metabird8, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Birds/Coef/anova_diet.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x=1.6, y= -1, label= Qm, size = 6)

# Calculation R2
mR2.func(metabird8) #6.36

c<-predict(metabird8, newmods = cbind(logmass, 0, 0), addx=TRUE)
nc<-predict(metabird8, newmods = cbind(logmass, 1, logmass), addx=TRUE)

### merge data frames and plot all together
c<-data.frame(c)
c$Islandtype<-rep("Carn", nrow(c))
nc<-data.frame(nc)
nc$Islandtype<-rep("Non Carn", nrow(nc))

df<-rbind(c,nc) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Carn" = "#E69F00", "Non Carn" = "#9966FF")

Bi<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values=colpalette) + scale_fill_manual(values=colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Diet"),colour=guide_legend(title="Diet"))+
  ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1,1, by=0.2), limits= c(-1,1))+
  ggtitle("")+ labs(tag = "i") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Birds/Hyp/diet.tif")
Bi
dev.off()

# multiplot####
multibirds<-ggarrange(Ba,Bb,Bc,Bd1,Be,Bf,Bg,Bh, Bi, ncol = 3, nrow = 3, align = "v")
tiff('~/GitHub/Island_Rule/Results/Figures/Figure_hyp_birds_upd1.tif', res=300, width=6000, height=6000, compression = "lzw")
multibirds
dev.off()

multibirds<-ggarrange(Ba,Bb,Bc,Bd2,Be,Bf,Bg,Bh, Bi, ncol = 3, nrow = 3, align = "v")
tiff('~/GitHub/Island_Rule/Results/Figures/Figure_hyp_birds_upd2.tif', res=300, width=6000, height=6000, compression = "lzw")
multibirds
dev.off()


#REPTILES####
phylocor<-list(Binomial=rept_phylo_cor)
logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

#island area ####
metarept2<-rma.mv(RR~logmass*Island_km2,V = Vrept, data=reptdata, random= RE,
                  R = phylocor, control=list(optimizer= "bobyqa"),method = "REML") 
summary(metarept2) #this is the actual hypothesis, QM(df = 3) = 30.6821, p-val < .0001

coef<-data.frame(b =metarept2$b, lci = metarept2$ci.lb, uci =  metarept2$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Reptiles/Coef/coef_area.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept2$beta)[-1]

test_1pred<-anova(metarept2, btt = 2) 
test_2pred<-anova(metarept2, btt = 3) 
test_int<-anova(metarept2, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Reptiles/Coef/anova_area.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept2) #24.25

# island area
Island_km2<-quantile(reptdata$Island_km2, prob = 0.1) #63 km2
s<-predict(metarept2, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

Island_km2<-quantile(reptdata$Island_km2, prob = 0.9) #147910.8
l<-predict(metarept2, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Small", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Large", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Large"="#9966FF","Small" = "#E69F00") #red yellow palette "Warm no pred"= "#FFCC33", "Cold no pred"= "#0072B2

Ra<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island area"),colour=guide_legend(title="Island area"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "a") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Reptiles/Hyp/area.tif")
Ra
dev.off()

##distance####
metarept3<-rma.mv(RR~logmass*Dist_near_mainland, V=Vrept,  data=reptdata,  random= RE,
                  R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metarept3) #QM(df = 3) = 22.2331, p-val < .0001

coef<-data.frame(b =metarept3$b, lci = metarept3$ci.lb, uci =  metarept3$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Reptiles/Coef/coef_distance.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept3$beta)[-1]

test_1pred<-anova(metarept3, btt = 2) 
test_2pred<-anova(metarept3, btt = 3) 
test_int<-anova(metarept3, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Reptiles/Coef/anova_dist.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept3) #20.3

logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

Dist_near_mainland<-quantile(reptdata$Dist_near_mainland, prob = 0.1) #63 km2
s<-predict(metarept3, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

Dist_near_mainland<-quantile(reptdata$Dist_near_mainland, prob = 0.9) #147910.8
l<-predict(metarept3, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)


### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Near", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Far", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Far"="#9966FF","Near" = "#E69F00") # orange purple

Rb<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island isolation"),colour=guide_legend(title="Island isolation"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "b") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Reptiles/Hyp/distance.tif")
Rb
dev.off()

##island area and distance####
metarept4<-rma.mv(RR~logmass*Dist_near_mainland +logmass*Island_km2,V=Vrept,  data=reptdata, random= RE,
                  R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metarept4)

coef<-data.frame(b =metarept4$b, lci = metarept4$ci.lb, uci =  metarept4$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Reptiles/Coef/coef_distance_area.csv", row.names = FALSE)

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

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Reptiles/Coef/anova_dist_area.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept4) #23.85

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
l$Islandtype<-rep("Small remote", nrow(l)) 
h<-data.frame(h)
h$Islandtype<-rep("Close large", nrow(h))

df<-rbind(l,h)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Close large" = "#9966FF","Small remote"="#E69F00") #  orange palette

# point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE) # turn off sci notation

Rc<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ 
  #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "c") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Reptiles/Hyp/distance_area.tif")
Rc
dev.off()

##predation pressure####
# richness
metarept7d1<-rma.mv(RR~logmass*TotPredrichEcto,V=Vrept,  data=reptdata,  random= RE,  
                   R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metarept7d1)

coef<-data.frame(b =metarept7d1$b, lci = metarept7d1$ci.lb, uci =  metarept7d1$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Reptiles/Coef/coef_predrich.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept7d1$beta)[-1]

test_1pred<-anova(metarept7d1, btt = 2) 
test_2pred<-anova(metarept7d1, btt = 3) 
test_int<-anova(metarept7d1, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Reptiles/Coef/anova_predrich.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept7d1) #11.93

# richness of predators
logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

TotPredrichEcto<-quantile(reptdata$TotPredrichEcto, prob = 0.1) 
s<-predict(metarept7d1, newmods = cbind(logmass,TotPredrichEcto, logmass*TotPredrichEcto), addx=TRUE)

TotPredrichEcto<-quantile(reptdata$TotPredrichEcto, prob = 0.9) 
l<-predict(metarept7d1, newmods = cbind(logmass,TotPredrichEcto, logmass*TotPredrichEcto), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Low", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("High", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("High"="#9966FF","Low" = "#E69F00") 

Rd1<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Predator richness"),colour=guide_legend(title="Predator richness"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "d")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Reptiles/Hyp/predrich.tif")
Rd1
dev.off()

###mammalian richness
metarept7d2<-rma.mv(RR~logmass*MPredRichEcto,V=Vrept,  data=reptdata,  random= RE,  
                    R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metarept7d2)

coef<-data.frame(b =metarept7d2$b, lci = metarept7d2$ci.lb, uci =  metarept7d2$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Reptiles/Coef/coef_predrich_mam.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept7d2$beta)[-1]

test_1pred<-anova(metarept7d2, btt = 2) 
test_2pred<-anova(metarept7d2, btt = 3) 
test_int<-anova(metarept7d2, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Reptiles/Coef/anova_predrich_mam.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept7d2) #11.93

# richness of predators
logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

MPredRichEcto<-quantile(reptdata$MPredRichEndo, prob = 0.1) 
s<-predict(metarept7d2, newmods = cbind(logmass,MPredRichEcto, logmass*MPredRichEcto), addx=TRUE)

MPredRichEcto<-quantile(reptdata$MPredRichEcto, prob = 0.9) #147910.8
l<-predict(metarept7d2, newmods = cbind(logmass,MPredRichEcto, logmass*MPredRichEcto), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Low", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("High", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("High"="#9966FF","Low" = "#E69F00") 

Rd2<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Predator richness"),colour=guide_legend(title="Predator richness"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "d")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Reptiles/Hyp/predrich_mam.tif")
Rd2
dev.off()


# pred pres
metarept7d<-rma.mv(RR~logmass*PredEcto,V=Vrept,  data=reptdata,  random= RE,  
                   R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metarept7d) #QM(df = 3) = 28.7268, p-val < .0001, 221 0 and 234 1

coef<-data.frame(b =metarept7d$b, lci = metarept7d$ci.lb, uci =  metarept7d$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Reptiles/Coef/coef_pred.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept7d$beta)[-1]

test_1pred<-anova(metarept7d, btt = 2) 
test_2pred<-anova(metarept7d, btt = 3) 
test_int<-anova(metarept7d, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Reptiles/Coef/anova_pred.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)


# Calculation R2
mR2.func(metarept7d) #20.59

PredEcto<-0
nv<-predict(metarept7d, newmods = cbind(logmass,PredEcto, logmass*PredEcto), addx=TRUE) 
PredEcto<-1
v<-predict(metarept7d, newmods = cbind(logmass,PredEcto,logmass*PredEcto), addx=TRUE) 

### merge data frames and plot all together
nv<-data.frame(nv)
nv$Islandtype<-rep("No predators", nrow(l)) 
v<-data.frame(v)
v$Islandtype<-rep("Predators", nrow(h))

df<-rbind(v,nv)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Predators"="#9966FF","No predators" = "#E69F00") 

Rd3<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island type"),colour=guide_legend(title="Island type"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "d") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Reptiles/Hyp/pred_pres.tif")
Rd3
dev.off()

# resource availability####
metarept5e<-rma.mv(RR~logmass*NDVI,V=Vrept,  data=reptdata,  random= RE,  
                   R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metarept5e)

coef<-data.frame(b =metarept5e$b, lci = metarept5e$ci.lb, uci =  metarept5e$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Reptiles/Coef/coef_ndvi.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept5e$beta)[-1]

test_1pred<-anova(metarept5e, btt = 2) 
test_2pred<-anova(metarept5e, btt = 3) 
test_int<-anova(metarept5e, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Reptiles/Coef/anova_ndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept5e) #22.06

NDVI<-quantile(reptdata$NDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(metarept5e, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

NDVI<-quantile(reptdata$NDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(metarept5e, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Re<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource availab."),colour=guide_legend(title="Resource availab."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "e") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Reptiles/Hyp/ndvi.tif")
Re
dev.off()

# seasonality in resources####
metarept5g<-rma.mv(RR~logmass*SDNDVI,V=Vrept,  data= reptdata,  random= RE,  
                   R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metarept5g)

coef<-data.frame(b =metarept5g$b, lci = metarept5g$ci.lb, uci =  metarept5g$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Reptiles/Coef/coef_sdndvi.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept5g$beta)[-1]

test_1pred<-anova(metarept5g, btt = 2) 
test_2pred<-anova(metarept5g, btt = 3) 
test_int<-anova(metarept5g, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Reptiles/Coef/anova_sdndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept5g) #20.65

SDNDVI<-quantile(reptdata$SDNDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(metarept5g, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

SDNDVI<-quantile(reptdata$SDNDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(metarept5g, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Rf<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource season."),colour=guide_legend(title="Resource season."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "f") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Reptiles/Hyp/sdndvi.tif")
Rf
dev.off()

# island temperature ####
metarept6f<-rma.mv(RR~logmass*tmean,V=Vrept,  data=reptdata,random= RE,  
                   R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metarept6f) #QM(df = 3) = 34.1169, p-val < .0001

coef<-data.frame(b =metarept6f$b, lci = metarept6f$ci.lb, uci =  metarept6f$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Reptiles/Coef/coef_tmean.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept6f$beta)[-1]

test_1pred<-anova(metarept6f, btt = 2) 
test_2pred<-anova(metarept6f, btt = 3) 
test_int<-anova(metarept6f, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Reptiles/Coef/anova_tmean.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept6f) #25.71

#temperature
tmean<-quantile(reptdata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metarept6f, newmods = cbind(logmass, tmean, logmass*tmean), addx=TRUE) 

tmean<-quantile(reptdata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metarept6f, newmods = cbind(logmass, tmean,logmass*tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Rg<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "g") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Reptiles/Hyp/tmean.tif")
Rg
dev.off()

# temperature intercept####
metarept6<-rma.mv(RR~logmass + tmean,V=Vrept,  data=reptdata,random= RE,  
                  R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metarept6)

coef<-data.frame(b =metarept6$b, lci = metarept6$ci.lb, uci =  metarept6$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Reptiles/Coef/coef_tmean_int.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept6$beta)[-1]

test_1pred<-anova(metarept6, btt = 2) 
test_2pred<-anova(metarept6, btt = 3) 
# test_int<-anova(meta6, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Reptiles/Coef/anova_tmean_int.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_2pred[1]), digits = 3), ", p-val = ", round(as.numeric(test_2pred[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept6) #22.33

#temperature
tmean<-quantile(reptdata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metarept6, newmods = cbind(logmass, tmean), addx=TRUE) 

tmean<-quantile(reptdata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metarept6, newmods = cbind(logmass, tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Rg2<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #cale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "g") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Reptiles/Hyp/tmean_intercept.tif") #interactive effect
Rg2
dev.off()

# temp seas####
metarept6h<-rma.mv(RR~logmass*tseas,V=Vrept,  data=reptdata,  random= RE,  
                   R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metarept6h) 

coef<-data.frame(b =metarept6h$b, lci = metarept6h$ci.lb, uci =  metarept6h$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Reptiles/Coef/coef_tseas.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept6f$beta)[-1]

test_1pred<-anova(metarept6h, btt = 2) 
test_2pred<-anova(metarept6h, btt = 3) 
test_int<-anova(metarept6h, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Reptiles/Coef/anova_tseas.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept6h) #21.92

tseas<-quantile(reptdata$tseas, prob = 0.1, names=FALSE) #1.68
l<-predict(metarept6h, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

tseas<-quantile(reptdata$tseas, prob = 0.9, names=FALSE) #2.53
h<-predict(metarept6h, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Rh<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.76,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Temperature season."),colour=guide_legend(title="Temperature season."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "h") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Reptiles/Hyp/tseas.tif")
Rh
dev.off()

# diet ####
reptdata$guild2<- ifelse(reptdata$guild == "Carnivorous", "Carn", "No Carn")
reptdata$guild2<- factor(reptdata$guild2)

metarept8<-rma.mv(RR~logmass*guild2 ,V=Vrept,  data=reptdata,  random= RE,  
                  R = phylocor, control=list(optimizer= "bobyqa"),method = "REML")  
summary(metarept8)

coef<-data.frame(b =metarept8$b, lci = metarept8$ci.lb, uci =  metarept8$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Reptiles/Coef/coef_diet.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept8$beta)[-1]

test_1pred<-anova(metarept8, btt = 2) 
test_2pred<-anova(metarept8, btt = 3) 
test_int<-anova(metarept8, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Reptiles/Coef/anova_diet.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept8) #6.36

c<-predict(metarept8, newmods = cbind(logmass, 0, 0), addx=TRUE)
nc<-predict(metarept8, newmods = cbind(logmass, 1, logmass), addx=TRUE)

### merge data frames and plot all together
c<-data.frame(c)
c$Islandtype<-rep("Carn", nrow(c))
nc<-data.frame(nc)
nc$Islandtype<-rep("Non Carn", nrow(nc))

df<-rbind(c,nc) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Carn" = "#E69F00", "Non Carn" = "#9966FF")

Ri<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values=colpalette) + scale_fill_manual(values=colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Diet"),colour=guide_legend(title="Diet"))+
  ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "i") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Reptiles/Hyp/diet.tif")
Ri
dev.off()

# multiplot####
multirept<-ggarrange(Ra,Rb,Rc,Rd1,Re,Rf,Rg2,Rh,Ri, ncol = 3, nrow = 3, align = "v")
tiff('~/GitHub/Island_Rule/Results/Figures/Figure_hyp_rept_upd1.tif', res=300, width=6000, height=6000, compression = "lzw")
multirept
dev.off()

#AMPHIBIANS####
phylocor<-list(Binomial=amph_phylo_cor)
logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

#island area ####
metaamph2<-rma.mv(RR~logmass*Island_km2,V = Vamph, data=amphdata, random= RE,
                  R = phylocor,method = "REML") 
summary(metaamph2) #this is the actual hypothesis, QM(df = 3) = 30.6821, p-val < .0001

coef<-data.frame(b =metaamph2$b, lci = metaamph2$ci.lb, uci =  metaamph2$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_area.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph2$beta)[-1]

test_1pred<-anova(metaamph2, btt = 2) 
test_2pred<-anova(metaamph2, btt = 3) 
test_int<-anova(metaamph2, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_area.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph2) #4.1

# island area
Island_km2<-quantile(amphdata$Island_km2, prob = 0.1) #63 km2
s<-predict(metaamph2, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

Island_km2<-quantile(amphdata$Island_km2, prob = 0.9) #147910.8
l<-predict(metaamph2, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Small", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Large", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Large"="#9966FF","Small" = "#E69F00") #red yellow palette "Warm no pred"= "#FFCC33", "Cold no pred"= "#0072B2

Aa<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island area"),colour=guide_legend(title="Island area"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "a") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/area.tif")
Aa
dev.off()

##distance####
metaamph3<-rma.mv(RR~logmass*Dist_near_mainland, V=Vamph,  data=amphdata,  random= RE,
                  R = phylocor,method = "REML")  
summary(metaamph3) #QM(df = 3) = 22.2331, p-val < .0001

coef<-data.frame(b =metaamph3$b, lci = metaamph3$ci.lb, uci =  metaamph3$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_distance.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph3$beta)[-1]

test_1pred<-anova(metaamph3, btt = 2) 
test_2pred<-anova(metaamph3, btt = 3) 
test_int<-anova(metaamph3, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_dist.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph3) #20.3

logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

Dist_near_mainland<-quantile(amphdata$Dist_near_mainland, prob = 0.1) #63 km2
s<-predict(metaamph3, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

Dist_near_mainland<-quantile(amphdata$Dist_near_mainland, prob = 0.9) #147910.8
l<-predict(metaamph3, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Near", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Far", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Far"="#9966FF","Near" = "#E69F00") # orange purple

Ab<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island isolation"),colour=guide_legend(title="Island isolation"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "b") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/distance.tif")
Ab
dev.off()

##island area and distance####
metaamph4<-rma.mv(RR~logmass*Dist_near_mainland +logmass*Island_km2,V=Vamph,  data=amphdata, random= RE,
                  R = phylocor,method = "REML")  
summary(metaamph4)

coef<-data.frame(b =metaamph4$b, lci = metaamph4$ci.lb, uci =  metaamph4$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_distance_area.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph4$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(metaamph4, btt = 2) 
test_2pred<-anova(metaamph4, btt = 3) 
test_3pred<-anova(metaamph4, btt = 4) 
test_int<-anova(metaamph4, btt = 5) 
test_int2<-anova(metaamph4, btt = 6) 
test_int3<-anova(metaamph4, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_dist_area.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph4) #23.85

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
l$Islandtype<-rep("Small remote", nrow(l)) 
h<-data.frame(h)
h$Islandtype<-rep("Close large", nrow(h))

df<-rbind(l,h)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Close large" = "#9966FF","Small remote"="#E69F00") #  orange palette

Ac<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ 
  #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "c") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/distance_area.tif")
Ac
dev.off()

##predation pressure####
# richness
metaamph7d1<-rma.mv(RR~logmass*TotPredrichEcto,V=Vamph,  data=amphdata,  random= RE,  
                    R = phylocor,method = "REML")  
summary(metaamph7d1)

coef<-data.frame(b =metaamph7d1$b, lci = metaamph7d1$ci.lb, uci =  metaamph7d1$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_predrich.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph7d1$beta)[-1]

test_1pred<-anova(metaamph7d1, btt = 2) 
test_2pred<-anova(metaamph7d1, btt = 3) 
test_int<-anova(metaamph7d1, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_predrich.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph7d1) #11.93

# richness of predators
logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

TotPredrichEcto<-quantile(amphdata$TotPredrichEcto, prob = 0.1) 
s<-predict(metaamph7d1, newmods = cbind(logmass,TotPredrichEcto, logmass*TotPredrichEcto), addx=TRUE)

TotPredrichEcto<-quantile(amphdata$TotPredrichEcto, prob = 0.9) 
l<-predict(metaamph7d1, newmods = cbind(logmass,TotPredrichEcto, logmass*TotPredrichEcto), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Low", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("High", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("High"="#9966FF","Low" = "#E69F00") 

Ad1<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Predator richness"),colour=guide_legend(title="Predator richness"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "d")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/predrich.tif")
Ad1
dev.off()

###mammalian richness
metaamph7d2<-rma.mv(RR~logmass*MPredRichEcto,V=Vamph,  data=amphdata,  random= RE,  
                    R = phylocor, method = "REML")  
summary(metaamph7d2)

coef<-data.frame(b =metaamph7d2$b, lci = metaamph7d2$ci.lb, uci =  metaamph7d2$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_predrich_mam.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph7d2$beta)[-1]

test_1pred<-anova(metaamph7d2, btt = 2) 
test_2pred<-anova(metaamph7d2, btt = 3) 
test_int<-anova(metaamph7d2, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_predrich_mam.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.9, y= -2.3, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph7d2) #11.93

# richness of predators
logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

MPredRichEcto<-quantile(amphdata$MPredRichEndo, prob = 0.1) 
s<-predict(metaamph7d2, newmods = cbind(logmass,MPredRichEcto, logmass*MPredRichEcto), addx=TRUE)

MPredRichEcto<-quantile(amphdata$MPredRichEcto, prob = 0.9) #147910.8
l<-predict(metaamph7d2, newmods = cbind(logmass,MPredRichEcto, logmass*MPredRichEcto), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Low", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("High", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("High"="#9966FF","Low" = "#E69F00") 

Ad2<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Predator richness"),colour=guide_legend(title="Predator richness"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.3,2.3, by=0.4), limits= c(-2.3,2.3))+
  ggtitle("") + labs(tag = "d")+Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/predrich_mam.tif")
Ad2
dev.off()

# pred pres
metaamph7d<-rma.mv(RR~logmass*PredEcto,V=Vamph,  data=amphdata,  random= RE,  
                   R = phylocor,method = "REML")  
summary(metaamph7d) #QM(df = 3) = 28.7268, p-val < .0001, 221 0 and 234 1

coef<-data.frame(b =metaamph7d$b, lci = metaamph7d$ci.lb, uci =  metaamph7d$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_pred.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph7d$beta)[-1]

test_1pred<-anova(metaamph7d, btt = 2) 
test_2pred<-anova(metaamph7d, btt = 3) 
test_int<-anova(metaamph7d, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_pred.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)


# Calculation R2
mR2.func(metaamph7d) #20.59

PredEcto<-0
nv<-predict(metaamph7d, newmods = cbind(logmass,PredEcto, logmass*PredEcto), addx=TRUE) 
PredEcto<-1
v<-predict(metaamph7d, newmods = cbind(logmass,PredEcto,logmass*PredEcto), addx=TRUE) 

### merge data frames and plot all together
nv<-data.frame(nv)
nv$Islandtype<-rep("No predators", nrow(l)) 
v<-data.frame(v)
v$Islandtype<-rep("Predators", nrow(h))

df<-rbind(v,nv)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Predators"="#9966FF","No predators" = "#E69F00") 

Ad3<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island type"),colour=guide_legend(title="Island type"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "d") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/pred_pres.tif")
Ad3
dev.off()

# resource availability####
metaamph5e<-rma.mv(RR~logmass*NDVI,V=Vamph,  data=amphdata,  random= RE,  
                   R = phylocor, method = "REML")  
summary(metaamph5e)

coef<-data.frame(b =metaamph5e$b, lci = metaamph5e$ci.lb, uci =  metaamph5e$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_ndvi.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph5e$beta)[-1]

test_1pred<-anova(metaamph5e, btt = 2) 
test_2pred<-anova(metaamph5e, btt = 3) 
test_int<-anova(metaamph5e, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_ndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph5e) #22.06

NDVI<-quantile(amphdata$NDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(metaamph5e, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

NDVI<-quantile(amphdata$NDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(metaamph5e, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Ae<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource availab."),colour=guide_legend(title="Resource availab."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "e") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/ndvi.tif")
Ae
dev.off()

# seasonality in resources####
metaamph5g<-rma.mv(RR~logmass*SDNDVI,V=Vamph,  data= amphdata,  random= RE,  
                   R = phylocor, method = "REML")  
summary(metaamph5g)

coef<-data.frame(b =metaamph5g$b, lci = metaamph5g$ci.lb, uci =  metaamph5g$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_sdndvi.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph5g$beta)[-1]

test_1pred<-anova(metaamph5g, btt = 2) 
test_2pred<-anova(metaamph5g, btt = 3) 
test_int<-anova(metaamph5g, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_sdndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph5g) #20.65

SDNDVI<-quantile(amphdata$SDNDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(metaamph5g, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

SDNDVI<-quantile(amphdata$SDNDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(metaamph5g, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Af<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource season."),colour=guide_legend(title="Resource season."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "f") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/sdndvi.tif")
Af
dev.off()

#### seasonality in resources intercept ####
metaamph5h<-rma.mv(RR~logmass + SDNDVI,V=Vamph,  data= amphdata,  random= RE,  
                   R = phylocor, method = "REML")  
summary(metaamph5h)

coef<-data.frame(b =metaamph5h$b, lci = metaamph5h$ci.lb, uci =  metaamph5h$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_sdndvi_int.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph5h$beta)[-1]

test_1pred<-anova(metaamph5h, btt = 2) 
test_2pred<-anova(metaamph5h, btt = 3) 
# test_int<-anova(meta6, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_sdndvi_int.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_2pred[1]), digits = 3), ", p-val = ", round(as.numeric(test_2pred[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph5h) #6.12

SDNDVI<-quantile(amphdata$SDNDVI, prob = 0.1, names=FALSE) 
l<-predict(metaamph5h, newmods = cbind(logmass, SDNDVI), addx=TRUE)

SDNDVI<-quantile(amphdata$SDNDVI, prob = 0.9, names=FALSE) 
h<-predict(metaamph5h, newmods = cbind(logmass, SDNDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Af2<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource season."),colour=guide_legend(title="Resource season."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "f") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/sdndvi_int.tif")
Af2
dev.off()

# island temperature ####
metaamph6f<-rma.mv(RR~logmass*tmean,V=Vamph,  data=amphdata,random= RE,  
                   R = phylocor,method = "REML")  
summary(metaamph6f)

coef<-data.frame(b =metaamph6f$b, lci = metaamph6f$ci.lb, uci =  metaamph6f$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_tmean.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph6f$beta)[-1]

test_1pred<-anova(metaamph6f, btt = 2) 
test_2pred<-anova(metaamph6f, btt = 3) 
test_int<-anova(metaamph6f, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_tmean.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph6f) #25.71

#temperature
tmean<-quantile(amphdata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metaamph6f, newmods = cbind(logmass, tmean, logmass*tmean), addx=TRUE) 

tmean<-quantile(amphdata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metaamph6f, newmods = cbind(logmass, tmean,logmass*tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Ag<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "g") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/tmean.tif")
Ag
dev.off()

# temperature intercept####
metaamph6<-rma.mv(RR~logmass + tmean,V=Vamph,  data=amphdata,random= RE,  
                  R = phylocor,method = "REML")  
summary(metaamph6)

coef<-data.frame(b =metaamph6$b, lci = metaamph6$ci.lb, uci =  metaamph6$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_tmean_int.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph6$beta)[-1]

test_1pred<-anova(metaamph6, btt = 2) 
test_2pred<-anova(metaamph6, btt = 3) 
# test_int<-anova(meta6, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_tmean_int.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_2pred[1]), digits = 3), ", p-val = ", round(as.numeric(test_2pred[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph6) #22.33

#temperature
tmean<-quantile(amphdata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metaamph6, newmods = cbind(logmass, tmean), addx=TRUE) 

tmean<-quantile(amphdata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metaamph6, newmods = cbind(logmass, tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Ag2<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #cale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "g") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/tmean_intercept.tif") #interactive effect
Ag2
dev.off()

# temp seas####
metaamph6h<-rma.mv(RR~logmass*tseas,V=Vamph,  data=amphdata,  random= RE,  
                   R = phylocor,method = "REML")  
summary(metaamph6h) 

coef<-data.frame(b =metaamph6h$b, lci = metaamph6h$ci.lb, uci =  metaamph6h$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_tseas.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph6f$beta)[-1]

test_1pred<-anova(metaamph6h, btt = 2) 
test_2pred<-anova(metaamph6h, btt = 3) 
test_int<-anova(metaamph6h, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_tseas.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph6h) #21.92

tseas<-quantile(amphdata$tseas, prob = 0.1, names=FALSE) #1.68
l<-predict(metaamph6h, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

tseas<-quantile(amphdata$tseas, prob = 0.9, names=FALSE) #2.53
h<-predict(metaamph6h, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Ah<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.76,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Temperature season."),colour=guide_legend(title="Temperature season."))+ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "h") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/tseas.tif")
Ah
dev.off()

# prec ####
metaamph8<-rma.mv(RR~logmass + prec ,V=Vamph,  data=amphdata,  random= RE,  
                  R = phylocor,method = "REML")  
summary(metaamph8)

coef<-data.frame(b =metaamph8$b, lci = metaamph8$ci.lb, uci =  metaamph8$ci.ub)
write.csv(coef, "~/GitHub/Island_Rule/Results/Amphibians/Coef/coef_prec.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph8$beta)[-1]

test_1pred<-anova(metaamph8, btt = 2) 
test_2pred<-anova(metaamph8, btt = 3) 
#test_int<-anova(metaamph8, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1]))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2]))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-prednames

write.csv(Qm_tot, "~/GitHub/Island_Rule/Results/Amphibians/Coef/anova_prec_int.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_2pred[1]), digits = 3), ", p-val = ", round(as.numeric(test_2pred[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph8) #2.12

prec<-quantile(amphdata$prec, prob = 0.1, names=FALSE) 
l<-predict(metaamph8, newmods = cbind(logmass, prec), addx=TRUE)

prec<-quantile(amphdata$prec, prob = 0.9, names=FALSE) 
h<-predict(metaamph8, newmods = cbind(logmass, prec), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Ai<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values=colpalette) + scale_fill_manual(values=colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Precipitation"),colour=guide_legend(title="Precipitation"))+
  ylab("RR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "i") + Qmtext

tiff(filename = "~/GitHub/Island_Rule/Results/Amphibians/Hyp/prec.tif")
Ai
dev.off()

# multiplot####
multiamph<-ggarrange(Aa,Ab,Ac,Ad1,Ae,Af2,Ag2,Ah,Ai, ncol = 3, nrow = 3, align = "v")
tiff('~/GitHub/Island_Rule/Results/Figures/Figure_hyp_amph_upd1.tif', res=300, width=6000, height=6000, compression = "lzw")
multiamph
dev.off()


