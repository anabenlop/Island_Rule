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

# load necessary functions
source("~/New projects/Island rule/Scripts/000_Functions.R")

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
# Vrept<-PDfunc(Vrept)
# is.positive.definite(Vrept) # TRUE

Vamph<- bldiag(lapply(split(amphdata, amphdata$CommonControl), calc.v))
is.positive.definite(Vamph) # TRUE

#ISLAND RULE PLOTS####
RE = list(~ 1 | Reference,~1|ID, ~1|SPID, ~1| Binomial)

#mammals####
phylocor<-list(Binomial= mam_phylo_cor)
metamam<-rma.mv(RR~logmass,V=Vmam, data=mamdata, random= RE,  R = phylocor)
summary(metamam)
mR2.func(metamam)
cR2.func(metamam)

# ggplot(mamdata) + geom_point(aes(logmass, RR)) + geom_hline(yintercept = 0, col = "red")+
#   geom_smooth(aes(logmass, RR), method ="lm")

logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

#predict for vector of mass
df_m<-predict(metamam, newmods = cbind(logmass), addx=TRUE)
df_m<-data.frame(df_m)
df_m$logmass<-df_m$X.logmass

#calculate size of points
wi    <- 1/sqrt(mamdata$var)
size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))

M<-ggplot(mamdata)+ geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_point(aes(logmass,RR), colour= "#0072B2",size = size,shape=20, alpha=I(.3)) +scale_shape_identity()+
  geom_line(data=df_m,aes(logmass,pred),color="#0072B2", size = 1.2)+
  geom_ribbon(data=df_m, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#0072B2", alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+xlab("log10(mass mainland (g))")+ ylab("lnRR")+ 
  scale_x_continuous(breaks=seq(0,6,1)) + 
  scale_y_continuous(breaks=seq(-1.6,1.6,0.4), limits= c(-1.6,1.6))+ 
  labs(tag = "a")

#birds#### 
phylocor<-list(Binomial= bird_phylo_cor)
metabird<-rma.mv(RR~logmass,V=Vbird,  data=birddata, random= RE, R = phylocor)
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
  geom_point(aes(logmass,RR), size = size,shape=20, color="#CC0000",alpha=I(.3)) +scale_shape_identity()+
  geom_line(data=df_b,aes(logmass,pred),color="#CC0000", size = 1.2)+
  geom_ribbon(data=df_b, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#CC0000", alpha=I(.4))+
  theme(element_blank(),axis.text=element_text(size=18,colour ="black"))+ xlab("log10(mass mainland (g))") + ylab("lnRR")+ 
  scale_x_continuous(breaks=seq(0.6,3.7,0.8), limits = c(0.6,3.7))+ 
  scale_y_continuous(breaks=seq(-1.4,1.4,0.4), limits= c(-1.4,1.4))+
  #ylim(-1.3,1.3)+
  labs(tag = "b")
B

#reptiles####
phylocor<-list(Binomial= rept_phylo_cor)
metarept<-rma.mv(RR~logmass,V=Vrept,  data=reptdata,  random= RE, R=phylocor) 
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
  geom_point(aes(logmass,RR),colour="#E69F00", size = size, shape = 20, alpha=I(.3)) + 
  geom_line(data=df_r,aes(logmass,pred),color="#E69F00", size = 1.2)+
  geom_ribbon(data=df_r, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#E69F00", alpha=I(.4))+
  theme(element_blank(),axis.text=element_text(size=18, colour ="black"))+ xlab("log10(mass mainland (g))") + ylab("lnRR")+ 
  #scale_x_continuous(breaks=seq(-1,4,1)) +
   scale_y_continuous(breaks=seq(-2.4,2.4, 0.6), limits = c(-2.4,2.4)) +
  #ylim(-2.3,2.3)+
  labs(tag = "c")

#amphibians####
phylocor<-list(Binomial= amph_phylo_cor)
metaamph<-rma.mv(RR~logmass ,V=var,  data=amphdata, random= RE, R = phylocor) 
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
  geom_point(aes(logmass,RR),colour="#009E73", size = size,shape = 20, alpha=I(.3)) + 
  geom_line(data=df_a,aes(logmass,pred),color="#009E73", size = 1.2)+
  geom_ribbon(data=df_a, aes(x=logmass, ymin=ci.lb,ymax=ci.ub), fill = "#009E73", alpha=I(.4))+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+ xlab("log10(mass mainland (g))") + ylab("lnRR") +
  # scale_x_continuous(breaks=seq(-0.6,1.7,0.5), limits = c(-0.6,1.7))+
  scale_y_continuous(breaks=seq(-1.6,1.6,0.4), limits =c(-1.6,1.6)) +
  labs(tag = "d")

saveRDS(metamam, file = "~/New projects/Island rule/Data/Final data/metamam.Rdata")
saveRDS(metabird, file = "~/New projects/Island rule/Data/Final data/metabird.Rdata")
saveRDS(metarept, file = "~/New projects/Island rule/Data/Final data/metarept.Rdata")
saveRDS(metaamph, file = "~/New projects/Island rule/Data/Final data/metaamph.Rdata")

#### FIGURE 3 ####
multiplot<-ggarrange(M,B,R,A, ncol = 2, nrow = 2, align = "v")
tiff('~/New projects/Island rule/Results/Figures/Figure3.tif', res=300, width=3100, height=3000)
multiplot
dev.off()

###ECOLOGICAL HYPOTHESIS ANALYSES#####
# MAMMALS####
tic("Run all models for mammals")
phylocor<-list(Binomial=mam_phylo_cor)
logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

#island area ####
# meta2<-rma.mv(RR~logmass*Island_km2,V = Vmam, data=mamdata, random= RE,
#               R = phylocor,method = "REML") 
# summary(meta2) 
# 
# saveRDS(meta2, file = "~/New projects/Island rule/Data/Final data/meta2.Rdata")
# 

meta2 <- readRDS(file = "~/New projects/Island rule/Data/Final data/meta2.Rdata")

coef<-data.frame(b =meta2$b, lci = meta2$ci.lb, uci =  meta2$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Mammals/Coef/coef_area.csv")

#extract Qm
prednames<-row.names(meta2$beta)[-1]

test_1pred<-anova(meta2, btt = 2) 
test_2pred<-anova(meta2, btt = 3) 
test_int<-anova(meta2, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], meta2$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], meta2$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Mammals/Coef/anova_area.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x= 2.5, y= -0.8, label= Qm, size = 6)

# Calculation R2
mR2.func(meta2) #11.75

# island area 
logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

Island_km2<-quantile(mamdata$Island_km2, prob = 0.1) #63 km2
s<-predict(meta2, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

Island_km2<-quantile(mamdata$Island_km2, prob = 0.9) #147910.8
l<-predict(meta2, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Small", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Large", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Large"="#9966FF","Small" = "#E69F00") #red yellow palette "Warm no pred"= "#FFCC33", "Cold no pred"= "#0072B2

Ma<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island area"),colour=guide_legend(title="Island area"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "a")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Mammals/Hyp/area.tif")
Ma
dev.off()

##distance####
# meta3<-rma.mv(RR~logmass*Dist_near_mainland, V=Vmam,  data=mamdata,  random= RE,
#               R = phylocor,method = "REML")  
# summary(meta3) #QM(df = 3) = 22.9852, p-val < .0001
# 
# saveRDS(meta3, file = "~/New projects/Island rule/Data/Final data/meta3.Rdata")

meta3 <- readRDS(file = "~/New projects/Island rule/Data/Final data/meta3.Rdata")

coef<-data.frame(b =meta3$b, lci = meta3$ci.lb, uci =  meta3$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Mammals/Coef/coef_dist.csv")

#extract Qm
prednames<-row.names(meta3$beta)[-1]

test_1pred<-anova(meta3, btt = 2) 
test_2pred<-anova(meta3, btt = 3) 
test_int<-anova(meta3, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], meta3$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], meta3$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Mammals/Coef/anova_dist.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(meta3) 

logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

Dist_near_mainland<-quantile(mamdata$Dist_near_mainland, prob = 0.1) #63 km2
s<-predict(meta2, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

Dist_near_mainland<-quantile(mamdata$Dist_near_mainland, prob = 0.9) #147910.8
l<-predict(meta2, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Close", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Remote", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Remote"="#9966FF","Close" = "#E69F00") # orange purple

Qmtext<-annotate(geom="text", x=2.5, y= -0.8, label= Qm, size = 6)

Mb<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island isolation"),colour=guide_legend(title="Island isolation"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "b")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Mammals/Hyp/distance.tif")
Mb
dev.off()

##Island area and remoteness####
# meta4<-rma.mv(RR~logmass*Dist_near_mainland +logmass*Island_km2,V=Vmam,  data=mamdata, random= RE,
#               R = phylocor, method = "REML")  
# summary(meta4)
# 
# saveRDS(meta4, file = "~/New projects/Island rule/Data/Final data/meta4.Rdata")

coef<-data.frame(b =meta4$b, lci = meta4$ci.lb, uci =  meta4$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Mammals/Coef/coef_dist_area.csv")

#extract Qm
prednames<-row.names(meta4$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(meta4, btt = 2) 
test_2pred<-anova(meta4, btt = 3) 
test_3pred<-anova(meta4, btt = 4) 
test_int<-anova(meta4, btt = 5) 
test_int2<-anova(meta4, btt = 6) 
test_int3<-anova(meta4, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1], meta4$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2], meta4$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Mammals/Coef/anova_dist_area.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

# Calculation R2
mR2.func(meta4) #12.3

logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

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
l$Islandtype<-rep("Small remote", nrow(l)) 
h<-data.frame(h)
h$Islandtype<-rep("Close large", nrow(h))

df<-rbind(l,h)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Close large" = "#9966FF","Small remote"="#E69F00") #  orange palette

Qmtext<-annotate(geom="text", x=2.6, y= -0.8, label= Qm, size = 6)

Mc<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.74,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ 
  scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "c")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Mammals/Hyp/distance_area.tif")
Mc
dev.off()

# diet ####
# mamdata$guild2<- ifelse(mamdata$guild == "Carn", "Carn", "No Carn")
# mamdata$guild2<- factor(mamdata$guild2)
# meta5<-rma.mv(RR~logmass*guild2,V=Vmam,  data=mamdata,  random= RE,  
#               R = phylocor, method = "REML")  
# summary(meta5)
# 
# saveRDS(meta5, file = "~/New projects/Island rule/Data/Final data/meta5.Rdata")

coef<-data.frame(b =meta5$b, lci = meta5$ci.lb, uci =  meta5$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Mammals/Coef/coef_diet.csv")

#extract Qm
prednames<-row.names(meta5$beta)[-1]

test_1pred<-anova(meta5, btt = 2) 
test_2pred<-anova(meta5, btt = 3) 
test_int<-anova(meta5, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], meta5$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], meta5$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Mammals/Coef/anova_diet.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(meta5) #9.93

c<-predict(meta5, newmods = cbind(logmass, 0, 0), addx=TRUE)
nc<-predict(meta5, newmods = cbind(logmass, 1, logmass), addx=TRUE)

### merge data frames and plot all together
c<-data.frame(c)
c$Islandtype<-rep("Carn", nrow(c))
nc<-data.frame(nc)
nc$Islandtype<-rep("Non Carn", nrow(nc))

df<-rbind(c,nc) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Carn" = "#E69F00", "Non Carn" = "#9966FF")

Qmtext<-annotate(geom="text", x=2.5, y= -0.8, label= Qm, size = 6)

Md<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values=colpalette) + scale_fill_manual(values=colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Diet"),colour=guide_legend(title="Diet"))+
  ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("")+ labs(tag = "d")+ Qmtext

tiff(filename = "~/New projects/Island rule/Results/Mammals/Hyp/diet.tif")
Md
dev.off()

# temperature####
# meta6<-rma.mv(RR~logmass*tmean,V=Vmam,  data=mamdata,random= RE,  
#                R = phylocor,method = "REML")  
# summary(meta6) 
# 
# saveRDS(meta6, file = "~/New projects/Island rule/Data/Final data/meta6.Rdata")

coef<-data.frame(b =meta6$b, lci = meta6$ci.lb, uci =  meta6$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Mammals/Coef/coef_tmean.csv")

#extract Qm
prednames<-row.names(meta6$beta)[-1]

test_1pred<-anova(meta6, btt = 2) 
test_2pred<-anova(meta6, btt = 3) 
test_int<-anova(meta6, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], meta6$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], meta6$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Mammals/Coef/anova_tmean.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(meta6) 

#temperature
tmean<-quantile(mamdata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(meta6, newmods = cbind(logmass, tmean, logmass*tmean), addx=TRUE) 

tmean<-quantile(mamdata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(meta6, newmods = cbind(logmass, tmean,logmass*tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Qmtext<-annotate(geom="text", x=2.5, y= -0.8, label= Qm, size = 6)

Me<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.76,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  # ylim(-0.7,0.7)+
  ggtitle("") + labs(tag = "e")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Mammals/Hyp/tmean.tif") #interactive effect
Me
dev.off()

# temperature intercept####
# meta6b<-rma.mv(RR~logmass + tmean,V=Vmam,  data=mamdata,random= RE,  
#               R = phylocor,method = "REML")  
# summary(meta6b)
# 
# saveRDS(meta6b, file = "~/New projects/Island rule/Data/Final data/meta6b.Rdata")

coef<-data.frame(b =meta6b$b, lci = meta6b$ci.lb, uci =  meta6b$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Mammals/Coef/coef_tmean_int.csv")

#extract Qm
prednames<-row.names(meta6b$beta)[-1]

test_1pred<-anova(meta6b, btt = 2) 
test_2pred<-anova(meta6b, btt = 3) 
# test_int<-anova(meta6b, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], meta6b$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], meta6b$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Mammals/Coef/anova_tmean_int.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_2pred[1]), digits = 3), ", p-val = ", round(as.numeric(test_2pred[2]), digits = 3))

# Calculation R2
mR2.func(meta6b) 

#temperature
tmean<-quantile(mamdata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(meta6b, newmods = cbind(logmass, tmean), addx=TRUE) 

tmean<-quantile(mamdata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(meta6b, newmods = cbind(logmass, tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Qmtext<-annotate(geom="text", x=2.5, y= -0.8, label= Qm, size = 6)

Me2<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.74,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  # ylim(-0.7,0.7)+
  ggtitle("") + labs(tag = "e")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Mammals/Hyp/tmean_int.tif") #interactive effect
Me2
dev.off()

# temp seas####
# meta7<-rma.mv(RR~logmass*tseas,V=Vmam,  data=mamdata,  random= RE,  
#                R = phylocor,method = "REML")  
# summary(meta7) 
# 
# saveRDS(meta7, file = "~/New projects/Island rule/Data/Final data/meta7.Rdata")

coef<-data.frame(b =meta7$b, lci = meta7$ci.lb, uci =  meta7$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Mammals/Coef/coef_tseas.csv")

#extract Qm
prednames<-row.names(meta7$beta)[-1]

test_1pred<-anova(meta7, btt = 2) 
test_2pred<-anova(meta7, btt = 3) 
test_int<-anova(meta7, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], meta7$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], meta7$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Mammals/Coef/anova_tseas.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))


# Calculation R2
mR2.func(meta7) 

tseas<-quantile(mamdata$tseas, prob = 0.1, names=FALSE) #1.7
l<-predict(meta7, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

tseas<-quantile(mamdata$tseas, prob = 0.9, names=FALSE) #81.79
h<-predict(meta7, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Qmtext<-annotate(geom="text", x=2.5, y= -0.8, label= Qm, size = 6)

Mf<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.66,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Temperature seasonality"),colour=guide_legend(title="Temperature seasonality"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  # ylim(-0.7,0.7)+
  ggtitle("")+ labs(tag = "f")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Mammals/Hyp/tseas.tif")
Mf
dev.off()

# resource availability####
# meta8<-rma.mv(RR~logmass*NDVI,V=Vmam,  data=mamdata,  random= RE,  
#               R = phylocor,method = "REML")  
# summary(meta8)
# 
# saveRDS(meta8, file = "~/New projects/Island rule/Data/Final data/meta8.Rdata")
# 
coef<-data.frame(b =meta8$b, lci = meta8$ci.lb, uci =  meta8$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Mammals/Coef/coef_ndvi.csv")

#extract Qm
prednames<-row.names(meta8$beta)[-1]

test_1pred<-anova(meta8, btt = 2) 
test_2pred<-anova(meta8, btt = 3) 
test_int<-anova(meta8, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], meta8$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], meta8$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Mammals/Coef/anova_ndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(meta8)

NDVI<-quantile(mamdata$NDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(meta8, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

NDVI<-quantile(mamdata$NDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(meta8, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Qmtext<-annotate(geom="text", x=2.5, y= -0.8, label= Qm, size = 6)

Mg<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.70,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource availability"),colour=guide_legend(title="Resource availability"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("")+ labs(tag = "g")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Mammals/Hyp/ndvi.tif")
Mg
dev.off()

# seasonality in resources####
# meta9<-rma.mv(RR~logmass*SDNDVI,V=Vmam,  data=mamdata,  random= RE,  
#                R = phylocor,method = "REML")  
# summary(meta9)
# 
# saveRDS(meta9, file = "~/New projects/Island rule/Data/Final data/meta9.Rdata")

coef<-data.frame(b =meta9$b, lci = meta9$ci.lb, uci =  meta9$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Mammals/Coef/coef_sdndvi.csv")

#extract Qm
prednames<-row.names(meta9$beta)[-1]

test_1pred<-anova(meta9, btt = 2) 
test_2pred<-anova(meta9, btt = 3) 
test_int<-anova(meta9, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], meta9$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], meta9$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Mammals/Coef/anova_sdndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))


# Calculation R2
mR2.func(meta9)

SDNDVI<-quantile(mamdata$SDNDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(meta9, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

SDNDVI<-quantile(mamdata$SDNDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(meta9, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Qmtext<-annotate(geom="text", x=2.5, y= -0.8, label= Qm, size = 6)

Mh<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.70,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource seasonality"),colour=guide_legend(title="Resource seasonality"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("")+ labs(tag = "h")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Mammals/Hyp/sdndvi.tif")
Mh
dev.off()

#### seasonality in resources intercept ####
# metamam10<-rma.mv(RR~logmass + SDNDVI,V=Vmam,  data= mamdata,  random= RE,  
#                    R = phylocor, method = "REML")  
# summary(metamam10)
# 
# saveRDS(metamam10, file = "~/New projects/Island rule/Data/Final data/metamam10.Rdata")

coef<-data.frame(b =metamam10$b, lci = metamam10$ci.lb, uci =  metamam10$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Mammals/Coef/coef_sdndvi_int.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metamam10$beta)[-1]

test_1pred<-anova(metamam10, btt = 2) 
test_2pred<-anova(metamam10, btt = 3) 
# test_int<-anova(meta6, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], metamam10$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], metamam10$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Mammals/Coef/anova_sdndvi_int.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_2pred[1]), digits = 3), ", p-val = ", round(as.numeric(test_2pred[2]), digits = 3))

Qmtext<- annotate(geom="text", x= 0.2, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metamam10) 

SDNDVI<-quantile(mamdata$SDNDVI, prob = 0.1, names=FALSE) 
l<-predict(metamam10, newmods = cbind(logmass, SDNDVI), addx=TRUE)

SDNDVI<-quantile(mamdata$SDNDVI, prob = 0.9, names=FALSE) 
h<-predict(metamam10, newmods = cbind(logmass, SDNDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Mh2<-ggplot(mamdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.70,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource seasonality"),colour=guide_legend(title="Resource seasonality"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "h") + Qmtext

tiff(filename = "~/New projects/Island rule/Results/Mammals/Hyp/sdndvi_int.tif")
Mh2
dev.off()

# multiplot####
multimam<-ggarrange(Ma,Mb,Mc,Md,Me2,Mf,Mg,Mh, ncol = 2, nrow = 4, align = "v")
tiff('~/New projects/Island rule/Results/Figures/Figure_hyp_mam_upd1.tif', res=300, width=3500, height=7000, compression = "lzw")
multimam
dev.off()

toc() #1792.58 s, 30 min

# BIRDS ####
tic("Run all models for birds")
phylocor<-list(Binomial=bird_phylo_cor)
logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

#island area ####
# metabird2<-rma.mv(RR~logmass*Island_km2,V = Vbird, data=birddata, random= RE,
#               R = phylocor,method = "REML") 
# summary(metabird2) 
# 
# saveRDS(metabird2, file = "~/New projects/Island rule/Data/Final data/metabird2.Rdata")

coef<-data.frame(b =metabird2$b, lci = metabird2$ci.lb, uci =  metabird2$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Birds/Coef/coef_area.csv")

#extract Qm
prednames<-row.names(metabird2$beta)[-1]

test_1pred<-anova(metabird2, btt = 2) 
test_2pred<-anova(metabird2, btt = 3) 
test_int<-anova(metabird2, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metabird2$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metabird2$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Birds/Coef/anova_area.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x= 1.8, y= -0.8, label= Qm, size = 6)

# Calculation R2
mR2.func(metabird2) 

# island area 
logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

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

colpalette<-c("Large"="#9966FF","Small" = "#E69F00") #red yellow palette "Warm no pred"= "#FFCC33", "Cold no pred"= "#0072B2

Ba<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.82,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island area"),colour=guide_legend(title="Island area"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "a")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Birds/Hyp/area.tif")
Ba
dev.off()

##distance####
# metabird3<-rma.mv(RR~logmass*Dist_near_mainland, V=Vbird,  data=birddata,  random= RE,
#               R = phylocor,method = "REML")  
# summary(metabird3) 
# 
# saveRDS(metabird3, file = "~/New projects/Island rule/Data/Final data/metabird3.Rdata")

metabird3 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metabird3.Rdata")

coef<-data.frame(b =metabird3$b, lci = metabird3$ci.lb, uci =  metabird3$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Birds/Coef/coef_dist.csv")

#extract Qm
prednames<-row.names(metabird3$beta)[-1]

test_1pred<-anova(metabird3, btt = 2) 
test_2pred<-anova(metabird3, btt = 3) 
test_int<-anova(metabird3, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metabird3$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metabird3$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Birds/Coef/anova_dist.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metabird3) 

logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

Dist_near_mainland<-quantile(birddata$Dist_near_mainland, prob = 0.1) #63 km2
s<-predict(metabird2, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

Dist_near_mainland<-quantile(birddata$Dist_near_mainland, prob = 0.9) #147910.8
l<-predict(metabird2, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Close", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Remote", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Remote"="#9966FF","Close" = "#E69F00") # orange purple

Qmtext<-annotate(geom="text", x=1.8, y= -0.8, label= Qm, size = 6)

Bb<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island isolation"),colour=guide_legend(title="Island isolation"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "b")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Birds/Hyp/distance.tif")
Bb
dev.off()

##Island area and remoteness####
# metabird4<-rma.mv(RR~logmass*Dist_near_mainland +logmass*Island_km2,V=Vbird,  data=birddata, random= RE,
#               R = phylocor, method = "REML")  
# summary(metabird4)
# 
# saveRDS(metabird4, file = "~/New projects/Island rule/Data/Final data/metabird4.Rdata")

metabird4 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metabird4.Rdata")

coef<-data.frame(b =metabird4$b, lci = metabird4$ci.lb, uci =  metabird4$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Birds/Coef/coef_dist_area.csv")

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

write.csv(Qm_tot, "~/New projects/Island rule/Results/Birds/Coef/anova_dist_area.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

# Calculation R2
mR2.func(metabird4) 

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
l$Islandtype<-rep("Small remote", nrow(l)) 
h<-data.frame(h)
h$Islandtype<-rep("Close large", nrow(h))

df<-rbind(l,h)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Close large" = "#9966FF","Small remote"="#E69F00") #  orange palette

Qmtext<-annotate(geom="text", x=1.8, y= -0.8, label= Qm, size = 6)

Bc<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.74,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ 
  scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("") + labs(tag = "c")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Birds/Hyp/distance_area.tif")
Bc
dev.off()

# diet ####
# birddata$guild2<- ifelse(birddata$guild == "VertFishScav", "Carn", "No Carn")
# birddata$guild2<- factor(birddata$guild2)
# metabird5<-rma.mv(RR~logmass*guild2,V=Vbird,  data=birddata,  random= RE,  
#               R = phylocor, method = "REML")  
# summary(metabird5)
# 
# saveRDS(metabird5, file = "~/New projects/Island rule/Data/Final data/metabird5.Rdata")

coef<-data.frame(b =metabird5$b, lci = metabird5$ci.lb, uci =  metabird5$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Birds/Coef/coef_diet.csv")

#extract Qm
prednames<-row.names(metabird5$beta)[-1]

test_1pred<-anova(metabird5, btt = 2) 
test_2pred<-anova(metabird5, btt = 3) 
test_int<-anova(metabird5, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metabird5$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metabird5$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Birds/Coef/anova_diet.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metabird5) 

c<-predict(metabird5, newmods = cbind(logmass, 0, 0), addx=TRUE)
nc<-predict(metabird5, newmods = cbind(logmass, 1, logmass), addx=TRUE)

### merge data frames and plot all together
c<-data.frame(c)
c$Islandtype<-rep("Carn", nrow(c))
nc<-data.frame(nc)
nc$Islandtype<-rep("Non Carn", nrow(nc))

df<-rbind(c,nc) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Carn" = "#E69F00", "Non Carn" = "#9966FF")

Qmtext<-annotate(geom="text", x=1.8, y= -0.8, label= Qm, size = 6)

Bd<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values=colpalette) + scale_fill_manual(values=colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Diet"),colour=guide_legend(title="Diet"))+
  ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("")+ labs(tag = "d")+ Qmtext

tiff(filename = "~/New projects/Island rule/Results/Birds/Hyp/diet.tif")
Bd
dev.off()

# temperature####
# metabird6<-rma.mv(RR~logmass*tmean,V=Vbird,  data=birddata,random= RE,  
#               R = phylocor,method = "REML")  
# summary(metabird6) 
# 
# saveRDS(metabird6, file = "~/New projects/Island rule/Data/Final data/metabird6.Rdata")

coef<-data.frame(b =metabird6$b, lci = metabird6$ci.lb, uci =  metabird6$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Birds/Coef/coef_tmean.csv")

#extract Qm
prednames<-row.names(metabird6$beta)[-1]

test_1pred<-anova(metabird6, btt = 2) 
test_2pred<-anova(metabird6, btt = 3) 
test_int<-anova(metabird6, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metabird6$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metabird6$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Birds/Coef/anova_tmean.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metabird6) 

#temperature
tmean<-quantile(birddata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metabird6, newmods = cbind(logmass, tmean, logmass*tmean), addx=TRUE) 

tmean<-quantile(birddata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metabird6, newmods = cbind(logmass, tmean,logmass*tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Qmtext<-annotate(geom="text", x=1.8, y= -0.8, label= Qm, size = 6)

Be<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  # ylim(-0.7,0.7)+
  ggtitle("") + labs(tag = "e")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Birds/Hyp/tmean.tif") #interactive effect
Be
dev.off()

# temperature intercept####
# metabird6b<-rma.mv(RR~logmass + tmean,V=Vbird,  data=birddata,random= RE,  
#                R = phylocor,method = "REML")  
# summary(metabird6b)
# 
# saveRDS(metabird6b, file = "~/New projects/Island rule/Data/Final data/metabird6b.Rdata")

metabird6b <- readRDS(file = "~/New projects/Island rule/Data/Final data/metabird6b.Rdata")

coef<-data.frame(b =metabird6b$b, lci = metabird6b$ci.lb, uci =  metabird6b$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Birds/Coef/coef_tmean_int.csv")

#extract Qm
prednames<-row.names(metabird6b$beta)[-1]

test_1pred<-anova(metabird6b, btt = 2) 
test_2pred<-anova(metabird6b, btt = 3) 
# test_int<-anova(metabird6b, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], metabird6b$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], metabird6b$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Birds/Coef/anova_tmean_int.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_2pred[1]), digits = 3), ", p-val = ", round(as.numeric(test_2pred[2]), digits = 3))

# Calculation R2
mR2.func(metabird6b) 

#temperature
tmean<-quantile(birddata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metabird6b, newmods = cbind(logmass, tmean), addx=TRUE) 

tmean<-quantile(birddata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metabird6b, newmods = cbind(logmass, tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Qmtext<-annotate(geom="text", x=1.7, y= -0.8, label= Qm, size = 6)

Be2<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.74,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  # ylim(-0.7,0.7)+
  ggtitle("") + labs(tag = "e")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Birds/Hyp/tmean_int.tif") #interactive effect
Be2
dev.off()

# temp seas####
# metabird7<-rma.mv(RR~logmass*tseas,V=Vbird,  data=birddata,  random= RE,  
#               R = phylocor,method = "REML")  
# summary(metabird7) 
# 
# saveRDS(metabird7, file = "~/New projects/Island rule/Data/Final data/metabird7.Rdata")

coef<-data.frame(b =metabird7$b, lci = metabird7$ci.lb, uci =  metabird7$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Birds/Coef/coef_tseas.csv")

#extract Qm
prednames<-row.names(metabird7$beta)[-1]

test_1pred<-anova(metabird7, btt = 2) 
test_2pred<-anova(metabird7, btt = 3) 
test_int<-anova(metabird7, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metabird7$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metabird7$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Birds/Coef/anova_tseas.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metabird7) 

tseas<-quantile(birddata$tseas, prob = 0.1, names=FALSE) #1.7
l<-predict(metabird7, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

tseas<-quantile(birddata$tseas, prob = 0.9, names=FALSE) #81.79
h<-predict(metabird7, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Qmtext<-annotate(geom="text", x=1.9, y= -0.8, label= Qm, size = 6)

Bf<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.66,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Temperature seasonality"),colour=guide_legend(title="Temperature seasonality"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  # ylim(-0.7,0.7)+
  ggtitle("")+ labs(tag = "f")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Birds/Hyp/tseas.tif")
Bf
dev.off()

# resource availability####
# metabird8<-rma.mv(RR~logmass*NDVI,V=Vbird,  data=birddata,  random= RE,  
#               R = phylocor,method = "REML")  
# summary(metabird8)
# 
# saveRDS(metabird8, file = "~/New projects/Island rule/Data/Final data/metabird8.Rdata")

coef<-data.frame(b =metabird8$b, lci = metabird8$ci.lb, uci =  metabird8$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Birds/Coef/coef_ndvi.csv")

#extract Qm
prednames<-row.names(metabird8$beta)[-1]

test_1pred<-anova(metabird8, btt = 2) 
test_2pred<-anova(metabird8, btt = 3) 
test_int<-anova(metabird8, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metabird8$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metabird8$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Birds/Coef/anova_ndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metabird8) 

NDVI<-quantile(birddata$NDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(metabird8, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

NDVI<-quantile(birddata$NDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(metabird8, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Qmtext<-annotate(geom="text", x=1.8, y= -0.8, label= Qm, size = 6)

Bg<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.70,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource availability"),colour=guide_legend(title="Resource availability"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("")+ labs(tag = "g")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Birds/Hyp/ndvi.tif")
Bg
dev.off()

# seasonality in resources####
# metabird9<-rma.mv(RR~logmass*SDNDVI,V=Vbird,  data=birddata,  random= RE,  
#               R = phylocor,method = "REML")  
# summary(metabird9)
# 
# saveRDS(metabird9, file = "~/New projects/Island rule/Data/Final data/metabird9.Rdata")

metabird9 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metabird9.Rdata")

coef<-data.frame(b =metabird9$b, lci = metabird9$ci.lb, uci =  metabird9$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Birds/Coef/coef_sdndvi.csv")

#extract Qm
prednames<-row.names(metabird9$beta)[-1]

test_1pred<-anova(metabird9, btt = 2) 
test_2pred<-anova(metabird9, btt = 3) 
test_int<-anova(metabird9, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metabird9$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metabird9$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Birds/Coef/anova_sdndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metabird9) #11.11

SDNDVI<-quantile(birddata$SDNDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(metabird9, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

SDNDVI<-quantile(birddata$SDNDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(metabird9, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Qmtext<-annotate(geom="text", x=1.8, y= -0.8, label= Qm, size = 6)

Bh<-ggplot(birddata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.70,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource seasonality"),colour=guide_legend(title="Resource seasonality"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-0.8,0.8, by=0.4), limits= c(-0.8,0.8))+
  ggtitle("")+ labs(tag = "h")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Birds/Hyp/sdndvi.tif")
Bh
dev.off()


# multiplot####
multibird<-ggarrange(Ba,Bb,Bc,Bd,Be2,Bf,Bg,Bh, ncol = 2, nrow = 4, align = "v")
tiff('~/New projects/Island rule/Results/Figures/Figure_hyp_birds_upd.tif', res=300, width=3500, height=7000, compression = "lzw")
multibird
dev.off()
toc() #

#REPTILES####
tic("Run all models for reptiles")
phylocor<-list(Binomial=rept_phylo_cor)
logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

#island area ####
# metarept2<-rma.mv(RR~logmass*Island_km2,V = Vrept, data=reptdata, random= RE,
#                   R = phylocor,method = "REML") 
# summary(metarept2) 
# 
# saveRDS(metarept2, file = "~/New projects/Island rule/Data/Final data/metarept2.Rdata")

metarept2 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metarept2.Rdata")

coef<-data.frame(b =metarept2$b, lci = metarept2$ci.lb, uci =  metarept2$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Reptiles/Coef/coef_area.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept2$beta)[-1]

test_1pred<-anova(metarept2, btt = 2) 
test_2pred<-anova(metarept2, btt = 3) 
test_int<-anova(metarept2, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metarept2$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metarept2$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Reptiles/Coef/anova_area.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x= 1.5, y= -2.4, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept2) #22.7

# island area 
logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

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
  guides(fill=guide_legend(title="Island area"),colour=guide_legend(title="Island area"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, by=0.6), limits= c(-2.4,2.4))+
  ggtitle("") + labs(tag = "a")+ Qmtext

tiff(filename = "~/New projects/Island rule/Results/Reptiles/Hyp/area.tif")
Ra
dev.off()

##distance####
# metarept3<-rma.mv(RR~logmass*Dist_near_mainland, V=Vrept,  data=reptdata,  random= RE,
#                   R = phylocor,method = "REML")  
# summary(metarept3) #QM(df = 3) = 22.9852, p-val < .0001
# 
# saveRDS(metarept3, file = "~/New projects/Island rule/Data/Final data/metarept3.Rdata")

metarept3 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metarept3.Rdata")

coef<-data.frame(b =metarept3$b, lci = metarept3$ci.lb, uci =  metarept3$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Reptiles/Coef/coef_dist.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept3$beta)[-1]

test_1pred<-anova(metarept3, btt = 2) 
test_2pred<-anova(metarept3, btt = 3) 
test_int<-anova(metarept3, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metarept3$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metarept3$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Reptiles/Coef/anova_dist.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metarept3) 

logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

Dist_near_mainland<-quantile(reptdata$Dist_near_mainland, prob = 0.1) #63 km2
s<-predict(metarept2, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

Dist_near_mainland<-quantile(reptdata$Dist_near_mainland, prob = 0.9) #147910.8
l<-predict(metarept2, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Close", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Remote", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Remote"="#9966FF","Close" = "#E69F00") # orange purple

Qmtext<-annotate(geom="text", x=1.5, y= -2.4, label= Qm, size = 6)

Rb<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island isolation"),colour=guide_legend(title="Island isolation"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, by=0.6), limits= c(-2.4,2.4))+
  ggtitle("") + labs(tag = "b")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Reptiles/Hyp/distance.tif")
Rb
dev.off()

##Island area and remoteness####
# metarept4<-rma.mv(RR~logmass*Dist_near_mainland +logmass*Island_km2,V=Vrept,  data=reptdata, random= RE,
#                   R = phylocor, method = "REML")  
# summary(metarept4)
# 
# saveRDS(metarept4, file = "~/New projects/Island rule/Data/Final data/metarept4.Rdata")

metarept4 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metarept4.Rdata")

coef<-data.frame(b =metarept4$b, lci = metarept4$ci.lb, uci =  metarept4$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Reptiles/Coef/coef_dist_area.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept4$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(metarept4, btt = 2) 
test_2pred<-anova(metarept4, btt = 3) 
test_3pred<-anova(metarept4, btt = 4) 
test_int<-anova(metarept4, btt = 5) 
test_int2<-anova(metarept4, btt = 6) 
test_int3<-anova(metarept4, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1], metarept4$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2],  metarept4$QM))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Reptiles/Coef/anova_dist_area.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

# Calculation R2
mR2.func(metarept4) #12.3

logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

# island isolation
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

Qmtext<-annotate(geom="text", x=1.5, y= -2.4, label= Qm, size = 6)

Rc<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.74,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ 
  scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, by=0.6), limits= c(-2.4,2.4))+
  ggtitle("") + labs(tag = "c")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Reptiles/Hyp/distance_area.tif")
Rc
dev.off()

# diet ####
# reptdata$guild2<- ifelse(reptdata$guild == "Carnivorous", "Carn", "No Carn")
# reptdata$guild2<- factor(reptdata$guild2)
# metarept5<-rma.mv(RR~logmass*guild2,V=Vrept,  data=reptdata,  random= RE,  
#                   R = phylocor, method = "REML")  
# summary(metarept5)
# 
# saveRDS(metarept5, file = "~/New projects/Island rule/Data/Final data/metarept5.Rdata")

metarept5 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metarept5.Rdata")

coef<-data.frame(b =metarept5$b, lci = metarept5$ci.lb, uci =  metarept5$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Reptiles/Coef/coef_diet.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept5$beta)[-1]

test_1pred<-anova(metarept5, btt = 2) 
test_2pred<-anova(metarept5, btt = 3) 
test_int<-anova(metarept5, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1],metarept5$QM ))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metarept5$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Reptiles/Coef/anova_diet.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metarept5) 

c<-predict(metarept5, newmods = cbind(logmass, 0, 0), addx=TRUE)
nc<-predict(metarept5, newmods = cbind(logmass, 1, logmass), addx=TRUE)

### merge data frames and plot all together
c<-data.frame(c)
c$Islandtype<-rep("Carn", nrow(c))
nc<-data.frame(nc)
nc$Islandtype<-rep("Non Carn", nrow(nc))

df<-rbind(c,nc) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Carn" = "#E69F00", "Non Carn" = "#9966FF")

Qmtext<-annotate(geom="text", x=1.5, y= -2.4, label= Qm, size = 6)

Rd<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values=colpalette) + scale_fill_manual(values=colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Diet"),colour=guide_legend(title="Diet"))+
  ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, by=0.6), limits= c(-2.4,2.4))+
  ggtitle("")+ labs(tag = "d")+ Qmtext

tiff(filename = "~/New projects/Island rule/Results/Reptiles/Hyp/diet.tif")
Rd
dev.off()

# temperature####
# metarept6<-rma.mv(RR~logmass*tmean,V=Vrept,  data=reptdata,random= RE,  
#                   R = phylocor,method = "REML")  
# summary(metarept6) 
# 
# saveRDS(metarept6, file = "~/New projects/Island rule/Data/Final data/metarept6.Rdata")

metarept6 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metarept6.Rdata")

coef<-data.frame(b =metarept6$b, lci = metarept6$ci.lb, uci =  metarept6$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Reptiles/Coef/coef_tmean.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept6$beta)[-1]

test_1pred<-anova(metarept6, btt = 2) 
test_2pred<-anova(metarept6, btt = 3) 
test_int<-anova(metarept6, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1],metarept6$QM ))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metarept6$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Reptiles/Coef/anova_tmean.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metarept6) #13.41

#temperature
tmean<-quantile(reptdata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metarept6, newmods = cbind(logmass, tmean, logmass*tmean), addx=TRUE) 

tmean<-quantile(reptdata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metarept6, newmods = cbind(logmass, tmean,logmass*tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Qmtext<-annotate(geom="text", x=1.5, y= -2.4, label= Qm, size = 6)

Re<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.74,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, by=0.6), limits= c(-2.4,2.4))+
  # ylim(-0.7,0.7)+
  ggtitle("") + labs(tag = "e")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Reptiles/Hyp/tmean.tif") #interactive effect
Re
dev.off()

# temperature intercept####
# metarept6b<-rma.mv(RR~logmass + tmean,V=Vrept,  data=reptdata,random= RE,  
#                    R = phylocor,method = "REML")  
# summary(metarept6b)
# 
# saveRDS(metarept6b, file = "~/New projects/Island rule/Data/Final data/metarept6b.Rdata")

metarept6b <- readRDS(file = "~/New projects/Island rule/Data/Final data/metarept6b.Rdata")

coef<-data.frame(b =metarept6b$b, lci = metarept6b$ci.lb, uci =  metarept6b$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Reptiles/Coef/coef_tmean_int.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept6b$beta)[-1]

test_1pred<-anova(metarept6b, btt = 2) 
test_2pred<-anova(metarept6b, btt = 3) 
# test_int<-anova(metarept6b, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], metarept6b$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], metarept6b$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Reptiles/Coef/anova_tmean_int.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_2pred[1]), digits = 3), ", p-val = ", round(as.numeric(test_2pred[2]), digits = 3))

# Calculation R2
mR2.func(metarept6b) 

#temperature
tmean<-quantile(reptdata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metarept6b, newmods = cbind(logmass, tmean), addx=TRUE) 

tmean<-quantile(reptdata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metarept6b, newmods = cbind(logmass, tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Qmtext<-annotate(geom="text", x=1.5, y= -2.4, label= Qm, size = 6)

Re2<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.74,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, by=0.6), limits= c(-2.4,2.4))+
  # ylim(-0.7,0.7)+
  ggtitle("") + labs(tag = "e")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Reptiles/Hyp/tmean_int.tif") #interactive effect
Re2
dev.off()

# temp seas####
# metarept7<-rma.mv(RR~logmass*tseas,V=Vrept,  data=reptdata,  random= RE,  
#                   R = phylocor,method = "REML")  
# summary(metarept7) 
# 
# saveRDS(metarept7, file = "~/New projects/Island rule/Data/Final data/metarept7.Rdata")

metarept7 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metarept7.Rdata")

coef<-data.frame(b =metarept7$b, lci = metarept7$ci.lb, uci =  metarept7$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Reptiles/Coef/coef_tseas.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept7$beta)[-1]

test_1pred<-anova(metarept7, btt = 2) 
test_2pred<-anova(metarept7, btt = 3) 
test_int<-anova(metarept7, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metarept7$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metarept7$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Reptiles/Coef/anova_tseas.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metarept7) 

tseas<-quantile(reptdata$tseas, prob = 0.1, names=FALSE) #1.7
l<-predict(metarept7, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

tseas<-quantile(reptdata$tseas, prob = 0.9, names=FALSE) #81.79
h<-predict(metarept7, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Qmtext<-annotate(geom="text", x=1.5, y= -2.4, label= Qm, size = 6)

Rf<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.66,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Temperature seasonality"),colour=guide_legend(title="Temperature seasonality"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, by=0.6), limits= c(-2.4,2.4))+
  # ylim(-0.7,0.7)+
  ggtitle("")+ labs(tag = "f")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Reptiles/Hyp/tseas.tif")
Rf
dev.off()

# resource availability####
# metarept8<-rma.mv(RR~logmass*NDVI,V=Vrept,  data=reptdata,  random= RE,  
#                   R = phylocor,method = "REML")  
# summary(metarept8)
# 
# saveRDS(metarept8, file = "~/New projects/Island rule/Data/Final data/metarept8.Rdata")
# 

metarept8 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metarept8.Rdata")

coef<-data.frame(b =metarept8$b, lci = metarept8$ci.lb, uci =  metarept8$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Reptiles/Coef/coef_ndvi.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept8$beta)[-1]

test_1pred<-anova(metarept8, btt = 2) 
test_2pred<-anova(metarept8, btt = 3) 
test_int<-anova(metarept8, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metarept8$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metarept8$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Reptiles/Coef/anova_ndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metarept8) 

NDVI<-quantile(reptdata$NDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(metarept8, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

NDVI<-quantile(reptdata$NDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(metarept8, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Qmtext<-annotate(geom="text", x=1.5, y= -2.4, label= Qm, size = 6)

Rg<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.70,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource availability"),colour=guide_legend(title="Resource availability"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, by=0.6), limits= c(-2.4,2.4))+
  ggtitle("")+ labs(tag = "g")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Reptiles/Hyp/ndvi.tif")
Rg
dev.off()

# seasonality in resources####
# metarept9<-rma.mv(RR~logmass*SDNDVI,V=Vrept,  data=reptdata,  random= RE,  
#                   R = phylocor,method = "REML")  
# summary(metarept9)
# 
# saveRDS(metarept9, file = "~/New projects/Island rule/Data/Final data/metarept9.Rdata")

metarept9 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metarept9.Rdata")

coef<-data.frame(b =metarept9$b, lci = metarept9$ci.lb, uci =  metarept9$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Reptiles/Coef/coef_sdndvi.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept9$beta)[-1]

test_1pred<-anova(metarept9, btt = 2) 
test_2pred<-anova(metarept9, btt = 3) 
test_int<-anova(metarept9, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metarept9$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metarept9$QM))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Reptiles/Coef/anova_sdndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))


# Calculation R2
mR2.func(metarept9) 

SDNDVI<-quantile(reptdata$SDNDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(metarept9, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

SDNDVI<-quantile(reptdata$SDNDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(metarept9, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Qmtext<-annotate(geom="text", x=1.5, y= -2.4, label= Qm, size = 6)

Rh<-ggplot(reptdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.70,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource seasonality"),colour=guide_legend(title="Resource seasonality"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-2.4,2.4, by=0.6), limits= c(-2.4,2.4))+
  ggtitle("")+ labs(tag = "h")+Qmtext

tiff(filename = "~/New projects/Island rule/Results/Reptiles/Hyp/sdndvi.tif")
Rh
dev.off()

# multiplot####
multirept<-ggarrange(Ra,Rb,Rc,Rd,Re,Rf,Rg,Rh, ncol = 2, nrow = 4, align = "v")
tiff('~/New projects/Island rule/Results/Figures/Figure_hyp_rept_upd.tif', res=300, width=3500, height=7000, compression = "lzw")
multirept
dev.off()
toc()

#AMPHIBIANS####
tic("Run all models for amphibians")
phylocor<-list(Binomial=amph_phylo_cor)
logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

#island area ####
# metaamph2<-rma.mv(RR~logmass*Island_km2,V = Vamph, data=amphdata, random= RE,
#                   R = phylocor,method = "REML") 
# summary(metaamph2) 
# 
# saveRDS(metaamph2, file = "~/New projects/Island rule/Data/Final data/metaamph2.Rdata")

metaamph2 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metaamph2.Rdata")

coef<-data.frame(b =metaamph2$b, lci = metaamph2$ci.lb, uci =  metaamph2$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Amphibians/Coef/coef_area.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph2$beta)[-1]

test_1pred<-anova(metaamph2, btt = 2) 
test_2pred<-anova(metaamph2, btt = 3) 
test_int<-anova(metaamph2, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metaamph2$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metaamph2$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Amphibians/Coef/anova_area.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph2)

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
  guides(fill=guide_legend(title="Island area"),colour=guide_legend(title="Island area"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "a") + Qmtext

tiff(filename = "~/New projects/Island rule/Results/Amphibians/Hyp/area.tif")
Aa
dev.off()

##distance####
# metaamph3<-rma.mv(RR~logmass*Dist_near_mainland, V=Vamph,  data=amphdata,  random= RE,
#                   R = phylocor,method = "REML")  
# summary(metaamph3) 
# 
# saveRDS(metaamph3, file = "~/New projects/Island rule/Data/Final data/metaamph3.Rdata")

metaamph3 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metaamph3.Rdata")

coef<-data.frame(b =metaamph3$b, lci = metaamph3$ci.lb, uci =  metaamph3$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Amphibians/Coef/coef_distance.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph3$beta)[-1]

test_1pred<-anova(metaamph3, btt = 2) 
test_2pred<-anova(metaamph3, btt = 3) 
test_int<-anova(metaamph3, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metaamph3$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metaamph3$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Amphibians/Coef/anova_dist.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph3) 

logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

Dist_near_mainland<-quantile(amphdata$Dist_near_mainland, prob = 0.1) #63 km2
s<-predict(metaamph3, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

Dist_near_mainland<-quantile(amphdata$Dist_near_mainland, prob = 0.9) #147910.8
l<-predict(metaamph3, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

### merge data frames and plot all together
s<-data.frame(s)
s$Islandtype<-rep("Close", nrow(s)) 
l<-data.frame(l)
l$Islandtype<-rep("Remote", nrow(l))

df<-rbind(s,l)
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Remote"="#9966FF","Close" = "#E69F00") # orange purple

Ab<-ggplot(amphdata)+ 
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Island isolation"),colour=guide_legend(title="Island isolation"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "b") + Qmtext

tiff(filename = "~/New projects/Island rule/Results/Amphibians/Hyp/distance.tif")
Ab
dev.off()

##island area and remoteness####
# metaamph4<-rma.mv(RR~logmass*Dist_near_mainland +logmass*Island_km2,V=Vamph,  data=amphdata, random= RE,
#                   R = phylocor,method = "REML")  
# summary(metaamph4)
# 
# saveRDS(metaamph4, file = "~/New projects/Island rule/Data/Final data/metaamph4.Rdata")
# 

metaamph4 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metaamph4.Rdata")

coef<-data.frame(b =metaamph4$b, lci = metaamph4$ci.lb, uci =  metaamph4$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Amphibians/Coef/coef_distance_area.csv", row.names = FALSE)

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

write.csv(Qm_tot, "~/New projects/Island rule/Results/Amphibians/Coef/anova_dist_area.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph4) 

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
  theme(element_blank(), legend.position = c(0.74,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Area and isolation"),colour=guide_legend(title="Area and isolation"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ 
  #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "c") + Qmtext

tiff(filename = "~/New projects/Island rule/Results/Amphibians/Hyp/distance_area.tif")
Ac
dev.off()

# prec ####
# metaamph5<-rma.mv(RR~logmass*prec ,V=Vamph,  data=amphdata,  random= RE,  
#                   R = phylocor,method = "REML")  
# summary(metaamph5)

# saveRDS(metaamph5, file = "~/New projects/Island rule/Data/Final data/metaamph5.Rdata")

metaamph5 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metaamph5.Rdata")

coef<-data.frame(b =metaamph5$b, lci = metaamph5$ci.lb, uci =  metaamph5$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Amphibians/Coef/coef_prec.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph5$beta)[-1]

test_1pred<-anova(metaamph5, btt = 2) 
test_2pred<-anova(metaamph5, btt = 3) 
test_int<-anova(metaamph5, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metaamph5$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metaamph5$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Amphibians/Coef/anova_prec.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph5) 

prec<-quantile(amphdata$prec, prob = 0.1, names=FALSE) 
l<-predict(metaamph5, newmods = cbind(logmass, prec, logmass*prec), addx=TRUE)

prec<-quantile(amphdata$prec, prob = 0.9, names=FALSE) 
h<-predict(metaamph5, newmods = cbind(logmass, prec, logmass*prec), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Ad<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values=colpalette) + scale_fill_manual(values=colpalette)+
  theme(element_blank(), legend.position = c(0.78,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Precipitation"),colour=guide_legend(title="Precipitation"))+
  ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "d") + Qmtext

tiff(filename = "~/New projects/Island rule/Results/Amphibians/Hyp/prec.tif")
Ad
dev.off()

# island temperature ####
# metaamph6<-rma.mv(RR~logmass*tmean,V=Vamph,  data=amphdata,random= RE,  
#                    R = phylocor,method = "REML")  
# summary(metaamph6)
# 
# saveRDS(metaamph6, file = "~/New projects/Island rule/Data/Final data/metaamph6.Rdata")

metaamph6 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metaamph6.Rdata")

coef<-data.frame(b =metaamph6$b, lci = metaamph6$ci.lb, uci =  metaamph6$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Amphibians/Coef/coef_tmean.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph6$beta)[-1]

test_1pred<-anova(metaamph6, btt = 2) 
test_2pred<-anova(metaamph6, btt = 3) 
test_int<-anova(metaamph6, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metaamph6$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metaamph6$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Amphibians/Coef/anova_tmean.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph6) 

#temperature
tmean<-quantile(amphdata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metaamph6, newmods = cbind(logmass, tmean, logmass*tmean), addx=TRUE) 

tmean<-quantile(amphdata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metaamph6, newmods = cbind(logmass, tmean,logmass*tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands
df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Ae<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.74,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "g") + Qmtext

tiff(filename = "~/New projects/Island rule/Results/Amphibians/Hyp/tmean.tif")
Ae
dev.off()

# temperature intercept####
# metaamph6b<-rma.mv(RR~logmass + tmean,V=Vamph,  data=amphdata,random= RE,  
#                   R = phylocor,method = "REML")  
# summary(metaamph6b)
# 
# saveRDS(metaamph6b, file = "~/New projects/Island rule/Data/Final data/metaamph6b.Rdata")

metaamph6b <- readRDS(file = "~/New projects/Island rule/Data/Final data/metaamph6b.Rdata")


coef<-data.frame(b =metaamph6b$b, lci = metaamph6b$ci.lb, uci =  metaamph6b$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Amphibians/Coef/coef_tmean_int.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph6b$beta)[-1]

test_1pred<-anova(metaamph6b, btt = 2) 
test_2pred<-anova(metaamph6b, btt = 3) 
# test_int<-anova(meta6, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], metaamph6b$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], metaamph6b$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Amphibians/Coef/anova_tmean_int.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_2pred[1]), digits = 3), ", p-val = ", round(as.numeric(test_2pred[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph6b) 

#temperature
tmean<-quantile(amphdata$tmean, prob = 0.9, names = FALSE) #27 degrees
w<-predict(metaamph6b, newmods = cbind(logmass, tmean), addx=TRUE) 

tmean<-quantile(amphdata$tmean, prob = 0.1, names =FALSE) #6.1 degrees
c<-predict(metaamph6b, newmods = cbind(logmass, tmean), addx=TRUE) 

### merge data frames and plot all together
w<-data.frame(w)
w$Islandtype<-rep("Warm", nrow(w))
c<-data.frame(c)
c$Islandtype<-rep("Cold", nrow(c))

df<-rbind(w,c) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Warm" = "#E69F00", "Cold" = "#9966FF")

Ae2<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+ #small remote, red
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.74,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Mean temperature"),colour=guide_legend(title="Mean temperature"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ #cale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "g") + Qmtext

tiff(filename = "~/New projects/Island rule/Results/Amphibians/Hyp/tmean_intercept.tif") #interactive effect
Ae2
dev.off()

# temp seas####
# metaamph7<-rma.mv(RR~logmass*tseas,V=Vamph,  data=amphdata,  random= RE,  
#                    R = phylocor,method = "REML")  
# summary(metaamph7) 
# 
# saveRDS(metaamph7, file = "~/New projects/Island rule/Data/Final data/metaamph7.Rdata")

metaamph7 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metaamph7.Rdata")

coef<-data.frame(b =metaamph7$b, lci = metaamph7$ci.lb, uci =  metaamph7$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Amphibians/Coef/coef_tseas.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph7$beta)[-1]

test_1pred<-anova(metaamph7, btt = 2) 
test_2pred<-anova(metaamph7, btt = 3) 
test_int<-anova(metaamph7, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metaamph7$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metaamph7$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Amphibians/Coef/anova_tseas.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph7) 

tseas<-quantile(amphdata$tseas, prob = 0.1, names=FALSE) #1.68
l<-predict(metaamph7, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

tseas<-quantile(amphdata$tseas, prob = 0.9, names=FALSE) #2.53
h<-predict(metaamph7, newmods = cbind(logmass, tseas, logmass*tseas), addx=TRUE)

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
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.66,0.86),legend.background = element_blank(),
        axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Temperature seasonality"),colour=guide_legend(title="Temperature seasonality"))+ylab("lnRR")+
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
    xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "h") + Qmtext

tiff(filename = "~/New projects/Island rule/Results/Amphibians/Hyp/tseas.tif")
Af
dev.off()

# resource availability####
# metaamph8<-rma.mv(RR~logmass*NDVI,V=Vamph,  data=amphdata,  random= RE,  
#                    R = phylocor, method = "REML")  
# summary(metaamph8)
# 
# saveRDS(metaamph8, file = "~/New projects/Island rule/Data/Final data/metaamph8.Rdata")

metaamph8 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metaamph8.Rdata")

coef<-data.frame(b =metaamph8$b, lci = metaamph8$ci.lb, uci =  metaamph8$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Amphibians/Coef/coef_ndvi.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph8$beta)[-1]

test_1pred<-anova(metaamph8, btt = 2) 
test_2pred<-anova(metaamph8, btt = 3) 
test_int<-anova(metaamph8, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metaamph8$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metaamph8$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Amphibians/Coef/anova_ndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph8) #22.06

NDVI<-quantile(amphdata$NDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(metaamph8, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

NDVI<-quantile(amphdata$NDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(metaamph8, newmods = cbind(logmass, NDVI, logmass*NDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Ag<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.70,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource availability"),colour=guide_legend(title="Resource availability"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "e") + Qmtext

tiff(filename = "~/New projects/Island rule/Results/Amphibians/Hyp/ndvi.tif")
Ag
dev.off()

# seasonality in resources####
# metaamph9<-rma.mv(RR~logmass*SDNDVI,V=Vamph,  data= amphdata,  random= RE,  
#                    R = phylocor, method = "REML")  
# summary(metaamph9)
# 
# saveRDS(metaamph9, file = "~/New projects/Island rule/Data/Final data/metaamph9.Rdata")

metaamph9 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metaamph9.Rdata")

coef<-data.frame(b =metaamph9$b, lci = metaamph9$ci.lb, uci =  metaamph9$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Amphibians/Coef/coef_sdndvi.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph9$beta)[-1]

test_1pred<-anova(metaamph9, btt = 2) 
test_2pred<-anova(metaamph9, btt = 3) 
test_int<-anova(metaamph9, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metaamph9$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metaamph9$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Amphibians/Coef/anova_sdndvi.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph9) 

SDNDVI<-quantile(amphdata$SDNDVI, prob = 0.1, names=FALSE) #22.38
l<-predict(metaamph9, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

SDNDVI<-quantile(amphdata$SDNDVI, prob = 0.9, names=FALSE) #87.09
h<-predict(metaamph9, newmods = cbind(logmass, SDNDVI, logmass*SDNDVI), addx=TRUE)

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
  theme(element_blank(), legend.position = c(0.70,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource seasonality"),colour=guide_legend(title="Resource seasonality"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "f") + Qmtext

tiff(filename = "~/New projects/Island rule/Results/Amphibians/Hyp/sdndvi.tif")
Ah
dev.off()

#### seasonality in resources intercept ####
# metaamph10<-rma.mv(RR~logmass + SDNDVI,V=Vamph,  data= amphdata,  random= RE,  
#                    R = phylocor, method = "REML")  
# summary(metaamph10)
# 
# saveRDS(metaamph10, file = "~/New projects/Island rule/Data/Final data/metaamph10.Rdata")

metaamph10 <- readRDS(file = "~/New projects/Island rule/Data/Final data/metaamph10.Rdata")

coef<-data.frame(b =metaamph10$b, lci = metaamph10$ci.lb, uci =  metaamph10$ci.ub)
write.csv(coef, "~/New projects/Island rule/Results/Amphibians/Coef/coef_sdndvi_int.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph10$beta)[-1]

test_1pred<-anova(metaamph10, btt = 2) 
test_2pred<-anova(metaamph10, btt = 3) 
# test_int<-anova(meta6, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], metaamph10$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], metaamph10$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "~/New projects/Island rule/Results/Amphibians/Coef/anova_sdndvi_int.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_2pred[1]), digits = 3), ", p-val = ", round(as.numeric(test_2pred[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph10) 

SDNDVI<-quantile(amphdata$SDNDVI, prob = 0.1, names=FALSE) 
l<-predict(metaamph10, newmods = cbind(logmass, SDNDVI), addx=TRUE)

SDNDVI<-quantile(amphdata$SDNDVI, prob = 0.9, names=FALSE) 
h<-predict(metaamph10, newmods = cbind(logmass, SDNDVI), addx=TRUE)

### merge data frames and plot all together
l<-data.frame(l)
l$Islandtype<-rep("Low", nrow(l))
h<-data.frame(h)
h$Islandtype<-rep("High", nrow(h))

df<-rbind(l,h) #all islands

df$Islandtype<-as.factor(df$Islandtype)
df$logmass<-df$X.logmass

colpalette<-c("Low" = "#9966FF", "High" = "#E69F00")

Ah2<-ggplot(amphdata)+ #geom_point(aes(logmass,RR), color = "grey",
  #     size = size,shape=16, 
  #     alpha=I(.3)) +scale_shape_identity()+
  geom_hline(yintercept= 0, size=1.2, linetype="dashed", color="dark gray")+theme_bw(base_size=18) +
  geom_line(data=df,aes(logmass,pred,colour=Islandtype, group = Islandtype), size = 1)+
  geom_ribbon(data=df,aes(logmass,ymin= ci.lb ,ymax=ci.ub ,fill=Islandtype, group = Islandtype), alpha = I(.5))+
  scale_colour_manual(values = colpalette) + scale_fill_manual(values = colpalette)+
  theme(element_blank(), legend.position = c(0.70,0.86), axis.text=element_text(size=18))+ 
  guides(fill=guide_legend(title="Resource seasonality"),colour=guide_legend(title="Resource seasonality"))+ylab("lnRR")+
  xlab("log10(mass mainland (g))")+ #scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(-1.6,1.6, by=0.4), limits= c(-1.6,1.6))+
  ggtitle("") + labs(tag = "f") + Qmtext

tiff(filename = "~/New projects/Island rule/Results/Amphibians/Hyp/sdndvi_int.tif")
Ah2
dev.off()

# multiplot####
multiamph<-ggarrange(Aa,Ab,Ac,Ad,Ae,Af,Ag,Ah2, ncol = 2, nrow = 4, align = "v")
tiff('~/New projects/Island rule/Results/Figures/Figure_hyp_amph_upd.tif', res=300, width=3500, height=7000, compression = "lzw")
multiamph
dev.off()
toc() # 171.12 s

### End of script ###