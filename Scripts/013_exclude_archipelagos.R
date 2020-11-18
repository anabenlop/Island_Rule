##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 9th of November 2020

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script runs the phylogenetic meta-regressions excluding archipelagos to test the effects of island
# area and spatial isolation for the paper: 

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

############################################################################
# Testing effects of distance and area after excluding archipelagos    #####
############################################################################
RE = list(~ 1 | Reference,~1|ID, ~1|SPID, ~1| Binomial)

# MAMMALs####
phylocor<-list(Binomial=mam_phylo_cor)
logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

#island area ####
metamam2_arch<-rma.mv(RR~logmass*Island_km2, subset = Archipielago == "No", V = Vmam, data=mamdata, random= RE,
                 R = phylocor,method = "REML")
summary(metamam2_arch)

saveRDS(metamam2_arch, file = "Data/Final data/metamam2_arch.Rdata")

metamam2_arch <- readRDS(file = "Data/Final data/metamam2_arch.Rdata")

coef<-data.frame(b =metamam2_arch$b, lci = metamam2_arch$ci.lb, uci =  metamam2_arch$ci.ub)
write.csv(coef, "Results/Mammals/Coef/coef_area_arch.csv")

#extract Qm
prednames<-row.names(metamam2_arch$beta)[-1]

test_1pred<-anova(metamam2_arch, btt = 2) 
test_2pred<-anova(metamam2_arch, btt = 3) 
test_int<-anova(metamam2_arch, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metamam2_arch$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metamam2_arch$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "Results/Mammals/Coef/anova_area_arch.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x= 2.5, y= -0.8, label= Qm, size = 6)

# Calculation R2
mR2.func(metamam2_arch) #11.75

# island area 
logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

Island_km2<-quantile(mamdata$Island_km2, prob = 0.1) #63 km2
s<-predict(metamam2_arch, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

Island_km2<-quantile(mamdata$Island_km2, prob = 0.9) #147910.8
l<-predict(metamam2_arch, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

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

tiff(filename = "Results/Mammals/Hyp/area_arch.tif")
Ma
dev.off()

##distance####
metamam3_arch<-rma.mv(RR~logmass*Dist_near_mainland, subset = Archipielago == "No", V=Vmam,  data=mamdata,  random= RE,
                 R = phylocor,method = "REML")
summary(metamam3_arch) #QM(df = 3) = 22.9852, p-val < .0001

saveRDS(metamam3_arch, file = "Data/Final data/metamam3_arch.Rdata")

metamam3_arch <- readRDS(file = "Data/Final data/metamam3_arch.Rdata")

coef<-data.frame(b =metamam3_arch$b, lci = metamam3_arch$ci.lb, uci =  metamam3_arch$ci.ub)
write.csv(coef, "Results/Mammals/Coef/coef_dist_arch.csv")

#extract Qm
prednames<-row.names(metamam3_arch$beta)[-1]

test_1pred<-anova(metamam3_arch, btt = 2) 
test_2pred<-anova(metamam3_arch, btt = 3) 
test_int<-anova(metamam3_arch, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metamam3_arch$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metamam3_arch$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "Results/Mammals/Coef/anova_dist_arch.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metamam3_arch) 

logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

Dist_near_mainland<-quantile(mamdata$Dist_near_mainland, prob = 0.1) #63 km2
s<-predict(metamam2_arch, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

Dist_near_mainland<-quantile(mamdata$Dist_near_mainland, prob = 0.9) #147910.8
l<-predict(metamam2_arch, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

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

tiff(filename = "Results/Mammals/Hyp/distance_arch.tif")
Mb
dev.off()

##Island area and remoteness####
metamam4_arch<-rma.mv(RR~logmass*Dist_near_mainland +logmass*Island_km2,
                      subset = Archipielago == "No", V=Vmam,  data=mamdata, random= RE,
                 R = phylocor, method = "REML")
summary(metamam4_arch)

saveRDS(metamam4_arch, file = "Data/Final data/metamam4_arch.Rdata")

metamam4_arch <- readRDS(file = "Data/Final data/metamam4_arch.Rdata")

coef<-data.frame(b =metamam4_arch$b, lci = metamam4_arch$ci.lb, uci =  metamam4_arch$ci.ub)
write.csv(coef, "Results/Mammals/Coef/coef_dist_area_arch.csv")

#extract Qm
prednames<-row.names(metamam4_arch$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(metamam4_arch, btt = 2) 
test_2pred<-anova(metamam4_arch, btt = 3) 
test_3pred<-anova(metamam4_arch, btt = 4) 
test_int<-anova(metamam4_arch, btt = 5) 
test_int2<-anova(metamam4_arch, btt = 6) 
test_int3<-anova(metamam4_arch, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1], metamam4_arch$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2], metamam4_arch$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "Results/Mammals/Coef/anova_dist_area_arch.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

# Calculation R2
mR2.func(metamam4_arch) #12.3

logmass <- seq(from = min(mamdata$logmass), to = max(mamdata$logmass), length.out = 1000)

# island isolation
small_island<-quantile(mamdata$Island_km2, prob = 0.1, names =FALSE) #4 km2
large_distance<-quantile(mamdata$Dist_near_mainland, prob = 0.9,names =FALSE) #150 km
large_island<-quantile(mamdata$Island_km2, prob = 0.9, names =FALSE) #32900
small_distance<-quantile(mamdata$Dist_near_mainland, prob = 0.1, names =FALSE) #1.5km

l<-predict(metamam4_arch, newmods = cbind(logmass, large_distance, small_island,
                                     logmass*large_distance, logmass*small_island), addx=TRUE) 
h<-predict(metamam4_arch, newmods = cbind(logmass,small_distance, large_island,
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

tiff(filename = "Results/Mammals/Hyp/distance_area_arch.tif")
Mc
dev.off()

# BIRDS ####
phylocor<-list(Binomial=bird_phylo_cor)
logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

#island area ####
metabird2_arch<-rma.mv(RR~logmass*Island_km2,V = Vbird, subset = Archipielago == "No", data=birddata, random= RE,
                  R = phylocor,method = "REML")
summary(metabird2_arch)

saveRDS(metabird2_arch, file = "Data/Final data/metabird2_arch.Rdata")

metabird2_arch <- readRDS(file = "Data/Final data/metabird2_arch.Rdata")

coef<-data.frame(b =metabird2_arch$b, lci = metabird2_arch$ci.lb, uci =  metabird2_arch$ci.ub)
write.csv(coef, "Results/Birds/Coef/coef_area_arch.csv")

#extract Qm
prednames<-row.names(metabird2_arch$beta)[-1]

test_1pred<-anova(metabird2_arch, btt = 2) 
test_2pred<-anova(metabird2_arch, btt = 3) 
test_int<-anova(metabird2_arch, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metabird2_arch$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metabird2_arch$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "Results/Birds/Coef/anova_area_arch.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x= 1.8, y= -0.8, label= Qm, size = 6)

# Calculation R2
mR2.func(metabird2_arch) 

# island area 
logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

Island_km2<-quantile(birddata$Island_km2, prob = 0.1) #63 km2
s<-predict(metabird2_arch, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

Island_km2<-quantile(birddata$Island_km2, prob = 0.9) #147910.8
l<-predict(metabird2_arch, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

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

tiff(filename = "Results/Birds/Hyp/area_arch.tif")
Ba
dev.off()

##distance####
metabird3_arch<-rma.mv(RR~logmass*Dist_near_mainland, V=Vbird, subset = Archipielago == "No", data=birddata,  random= RE,
                  R = phylocor,method = "REML")
summary(metabird3_arch)

saveRDS(metabird3_arch, file = "Data/Final data/metabird3_arch.Rdata")

metabird3_arch <- readRDS(file = "Data/Final data/metabird3_arch.Rdata")

coef<-data.frame(b =metabird3_arch$b, lci = metabird3_arch$ci.lb, uci =  metabird3_arch$ci.ub)
write.csv(coef, "Results/Birds/Coef/coef_dist_arch.csv")

#extract Qm
prednames<-row.names(metabird3_arch$beta)[-1]

test_1pred<-anova(metabird3_arch, btt = 2) 
test_2pred<-anova(metabird3_arch, btt = 3) 
test_int<-anova(metabird3_arch, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metabird3_arch$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metabird3_arch$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "Results/Birds/Coef/anova_dist_arch.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metabird3_arch) 

logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

Dist_near_mainland<-quantile(birddata$Dist_near_mainland, prob = 0.1) #63 km2
s<-predict(metabird2_arch, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

Dist_near_mainland<-quantile(birddata$Dist_near_mainland, prob = 0.9) #147910.8
l<-predict(metabird2_arch, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

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

tiff(filename = "Results/Birds/Hyp/distance_arch.tif")
Bb
dev.off()

##Island area and remoteness####
metabird4_arch<-rma.mv(RR~logmass*Dist_near_mainland +logmass*Island_km2,V=Vbird, subset = Archipielago == "No",
                       data=birddata, random= RE, R = phylocor, method = "REML")
summary(metabird4_arch)

saveRDS(metabird4_arch, file = "Data/Final data/metabird4_arch.Rdata")

metabird4_arch <- readRDS(file = "Data/Final data/metabird4_arch.Rdata")

coef<-data.frame(b =metabird4_arch$b, lci = metabird4_arch$ci.lb, uci =  metabird4_arch$ci.ub)
write.csv(coef, "Results/Birds/Coef/coef_dist_area_arch.csv")

#extract Qm
prednames<-row.names(metabird4_arch$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(metabird4_arch, btt = 2) 
test_2pred<-anova(metabird4_arch, btt = 3) 
test_3pred<-anova(metabird4_arch, btt = 4) 
test_int<-anova(metabird4_arch, btt = 5) 
test_int2<-anova(metabird4_arch, btt = 6) 
test_int3<-anova(metabird4_arch, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1], metabird4_arch$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2], metabird4_arch$QM))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "Results/Birds/Coef/anova_dist_area_arch.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

# Calculation R2
mR2.func(metabird4_arch) 

logmass <- seq(from = min(birddata$logmass), to = max(birddata$logmass), length.out = 1000)

# island isolation
small_island<-quantile(birddata$Island_km2, prob = 0.1, names =FALSE) #4 km2
large_distance<-quantile(birddata$Dist_near_mainland, prob = 0.9,names =FALSE) #150 km
large_island<-quantile(birddata$Island_km2, prob = 0.9, names =FALSE) #32900
small_distance<-quantile(birddata$Dist_near_mainland, prob = 0.1, names =FALSE) #1.5km

l<-predict(metabird4_arch, newmods = cbind(logmass, large_distance, small_island,
                                      logmass*large_distance, logmass*small_island), addx=TRUE) 
h<-predict(metabird4_arch, newmods = cbind(logmass,small_distance, large_island,
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

tiff(filename = "Results/Birds/Hyp/distance_area_arch.tif")
Bc
dev.off()

#REPTILES####
phylocor<-list(Binomial=rept_phylo_cor)
logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

#island area ####
metarept2_arch<-rma.mv(RR~logmass*Island_km2,V = Vrept,subset = Archipielago == "No", data=reptdata, random= RE,
                  R = phylocor,method = "REML")
summary(metarept2_arch)

saveRDS(metarept2_arch, file = "Data/Final data/metarept2_arch.Rdata")

metarept2_arch <- readRDS(file = "Data/Final data/metarept2_arch.Rdata")

coef<-data.frame(b =metarept2_arch$b, lci = metarept2_arch$ci.lb, uci =  metarept2_arch$ci.ub)
write.csv(coef, "Results/Reptiles/Coef/coef_area_arch.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept2_arch$beta)[-1]

test_1pred<-anova(metarept2_arch, btt = 2) 
test_2pred<-anova(metarept2_arch, btt = 3) 
test_int<-anova(metarept2_arch, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metarept2_arch$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metarept2_arch$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "Results/Reptiles/Coef/anova_area_arch.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<-annotate(geom="text", x= 1.5, y= -2.4, label= Qm, size = 6)

# Calculation R2
mR2.func(metarept2_arch) #22.7

# island area 
logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

Island_km2<-quantile(reptdata$Island_km2, prob = 0.1) #63 km2
s<-predict(metarept2_arch, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

Island_km2<-quantile(reptdata$Island_km2, prob = 0.9) #147910.8
l<-predict(metarept2_arch, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

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

tiff(filename = "Results/Reptiles/Hyp/area_arch.tif")
Ra
dev.off()

##distance####
metarept3_arch<-rma.mv(RR~logmass*Dist_near_mainland, V=Vrept,subset = Archipielago == "No", data=reptdata,  random= RE,
                  R = phylocor,method = "REML")
summary(metarept3_arch) #QM(df = 3) = 22.9852, p-val < .0001

saveRDS(metarept3_arch, file = "Data/Final data/metarept3_arch.Rdata")

metarept3_arch <- readRDS(file = "Data/Final data/metarept3_arch.Rdata")

coef<-data.frame(b =metarept3_arch$b, lci = metarept3_arch$ci.lb, uci =  metarept3_arch$ci.ub)
write.csv(coef, "Results/Reptiles/Coef/coef_dist_arch.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept3_arch$beta)[-1]

test_1pred<-anova(metarept3_arch, btt = 2) 
test_2pred<-anova(metarept3_arch, btt = 3) 
test_int<-anova(metarept3_arch, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metarept3_arch$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metarept3_arch$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "Results/Reptiles/Coef/anova_dist_arch.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

# Calculation R2
mR2.func(metarept3_arch) 

logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

Dist_near_mainland<-quantile(reptdata$Dist_near_mainland, prob = 0.1) #63 km2
s<-predict(metarept2_arch, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

Dist_near_mainland<-quantile(reptdata$Dist_near_mainland, prob = 0.9) #147910.8
l<-predict(metarept2_arch, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

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

tiff(filename = "Results/Reptiles/Hyp/distance_arch.tif")
Rb
dev.off()

##Island area and remoteness####
metarept4_arch<-rma.mv(RR~logmass*Dist_near_mainland +logmass*Island_km2,V=Vrept,subset = Archipielago == "No",
                  data=reptdata, random= RE,
                  R = phylocor, method = "REML")
summary(metarept4_arch)

saveRDS(metarept4_arch, file = "Data/Final data/metarept4_arch.Rdata")

metarept4_arch <- readRDS(file = "Data/Final data/metarept4_arch.Rdata")

coef<-data.frame(b =metarept4_arch$b, lci = metarept4_arch$ci.lb, uci =  metarept4_arch$ci.ub)
write.csv(coef, "Results/Reptiles/Coef/coef_dist_area_arch.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metarept4_arch$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(metarept4_arch, btt = 2) 
test_2pred<-anova(metarept4_arch, btt = 3) 
test_3pred<-anova(metarept4_arch, btt = 4) 
test_int<-anova(metarept4_arch, btt = 5) 
test_int2<-anova(metarept4_arch, btt = 6) 
test_int3<-anova(metarept4_arch, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1], metarept4_arch$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2],  metarept4_arch$QM))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "Results/Reptiles/Coef/anova_dist_area_arch.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

# Calculation R2
mR2.func(metarept4_arch) #12.3

logmass <- seq(from = min(reptdata$logmass), to = max(reptdata$logmass), length.out = 1000)

# island isolation
small_island<-quantile(reptdata$Island_km2, prob = 0.1, names =FALSE) #4 km2
large_distance<-quantile(reptdata$Dist_near_mainland, prob = 0.9,names =FALSE) #150 km
large_island<-quantile(reptdata$Island_km2, prob = 0.9, names =FALSE) #32900
small_distance<-quantile(reptdata$Dist_near_mainland, prob = 0.1, names =FALSE) #1.5km

l<-predict(metarept4_arch, newmods = cbind(logmass, large_distance, small_island,
                                      logmass*large_distance, logmass*small_island), addx=TRUE) 
h<-predict(metarept4_arch, newmods = cbind(logmass,small_distance, large_island,
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

tiff(filename = "Results/Reptiles/Hyp/distance_area_arch.tif")
Rc
dev.off()

#AMPHIBIANS####
phylocor<-list(Binomial=amph_phylo_cor)
logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

#island area ####
metaamph2_arch<-rma.mv(RR~logmass*Island_km2,V = Vamph, data=amphdata, random= RE,
                  R = phylocor,method = "REML")
summary(metaamph2_arch)

saveRDS(metaamph2_arch, file = "Data/Final data/metaamph2_arch.Rdata")

metaamph2_arch <- readRDS(file = "Data/Final data/metaamph2_arch.Rdata")

coef<-data.frame(b =metaamph2_arch$b, lci = metaamph2_arch$ci.lb, uci =  metaamph2_arch$ci.ub)
write.csv(coef, "Results/Amphibians/Coef/coef_area_arch.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph2_arch$beta)[-1]

test_1pred<-anova(metaamph2_arch, btt = 2) 
test_2pred<-anova(metaamph2_arch, btt = 3) 
test_int<-anova(metaamph2_arch, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metaamph2_arch$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metaamph2_arch$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "Results/Amphibians/Coef/anova_area_arch.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph2_arch)

# island area
Island_km2<-quantile(amphdata$Island_km2, prob = 0.1) #63 km2
s<-predict(metaamph2_arch, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

Island_km2<-quantile(amphdata$Island_km2, prob = 0.9) #147910.8
l<-predict(metaamph2_arch, newmods = cbind(logmass,Island_km2, logmass*Island_km2), addx=TRUE)

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

tiff(filename = "Results/Amphibians/Hyp/area_arch.tif")
Aa
dev.off()

##distance####
metaamph3_arch<-rma.mv(RR~logmass*Dist_near_mainland, V=Vamph,  data=amphdata,  random= RE,
                  R = phylocor,method = "REML")
summary(metaamph3_arch)

saveRDS(metaamph3_arch, file = "Data/Final data/metaamph3_arch.Rdata")

metaamph3_arch <- readRDS(file = "Data/Final data/metaamph3_arch.Rdata")

coef<-data.frame(b =metaamph3_arch$b, lci = metaamph3_arch$ci.lb, uci =  metaamph3_arch$ci.ub)
write.csv(coef, "Results/Amphibians/Coef/coef_distance_arch.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph3_arch$beta)[-1]

test_1pred<-anova(metaamph3_arch, btt = 2) 
test_2pred<-anova(metaamph3_arch, btt = 3) 
test_int<-anova(metaamph3_arch, btt = 4) 

Qm_df<-t(data.frame(test_1pred[1],test_2pred[1], test_int[1], metaamph3_arch$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2], test_int[2], metaamph3_arch$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "Results/Amphibians/Coef/anova_dist_arch.csv")

Qm <- paste0("QM(df = 1) = ", round(as.numeric(test_int[1]), digits = 3), ", p-val = ", round(as.numeric(test_int[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph3_arch) 

logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

Dist_near_mainland<-quantile(amphdata$Dist_near_mainland, prob = 0.1) #63 km2
s<-predict(metaamph3_arch, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

Dist_near_mainland<-quantile(amphdata$Dist_near_mainland, prob = 0.9) #147910.8
l<-predict(metaamph3_arch, newmods = cbind(logmass,Dist_near_mainland, logmass*Dist_near_mainland), addx=TRUE)

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

tiff(filename = "Results/Amphibians/Hyp/distance_arch.tif")
Ab
dev.off()

##island area and remoteness####
metaamph4_arch<-rma.mv(RR~logmass*Dist_near_mainland +logmass*Island_km2,V=Vamph,  data=amphdata, random= RE,
                  R = phylocor,method = "REML")
summary(metaamph4_arch)

saveRDS(metaamph4_arch, file = "Data/Final data/metaamph4_arch.Rdata")
 
metaamph4_arch <- readRDS(file = "Data/Final data/metaamph4_arch.Rdata")

coef<-data.frame(b =metaamph4_arch$b, lci = metaamph4_arch$ci.lb, uci =  metaamph4_arch$ci.ub)
write.csv(coef, "Results/Amphibians/Coef/coef_distance_area_arch.csv", row.names = FALSE)

#extract Qm
prednames<-row.names(metaamph4_arch$beta)[-1]
prednames<-c(prednames,"logmass:dist & logmass:area")

test_1pred<-anova(metaamph4_arch, btt = 2) 
test_2pred<-anova(metaamph4_arch, btt = 3) 
test_3pred<-anova(metaamph4_arch, btt = 4) 
test_int<-anova(metaamph4_arch, btt = 5) 
test_int2<-anova(metaamph4_arch, btt = 6) 
test_int3<-anova(metaamph4_arch, btt = c(5,6)) 

Qm_df<- t(data.frame(test_1pred[1],test_2pred[1],test_3pred[1],test_int[1], test_int2[1], test_int3[1], metaamph4_arch$QM))
Qmp_df<-t(data.frame(test_1pred[2],test_2pred[2],test_3pred[2],test_int[2], test_int2[2], test_int3[2], metaamph4_arch$QMp))

Qm_tot<-cbind(Qm_df,Qmp_df)

colnames(Qm_tot)<-c("Qm", "p")
rownames(Qm_tot)<-c(prednames, "fullmodel")

write.csv(Qm_tot, "Results/Amphibians/Coef/anova_dist_area_arch.csv")

Qm <- paste0("QM(df = 2) = ", round(as.numeric(test_int3[1]), digits = 3), ", p-val = ", round(as.numeric(test_int3[2]), digits = 3))

Qmtext<- annotate(geom="text", x=0.35, y= -1.6, label= Qm, size = 6)

# Calculation R2
mR2.func(metaamph4_arch) 

logmass <- seq(from = min(amphdata$logmass), to = max(amphdata$logmass), length.out = 1000)

# create variables
small_island<-quantile(amphdata$Island_km2, prob = 0.1, names =FALSE) #4 km2
large_distance<-quantile(amphdata$Dist_near_mainland, prob = 0.9,names =FALSE) #150 km
large_island<-quantile(amphdata$Island_km2, prob = 0.9, names =FALSE) #32900
small_distance<-quantile(amphdata$Dist_near_mainland, prob = 0.1, names =FALSE) #1.5km

l<-predict(metaamph4_arch, newmods = cbind(logmass, large_distance, small_island,
                                      logmass*large_distance, logmass*small_island), addx=TRUE) 
h<-predict(metaamph4_arch, newmods = cbind(logmass,small_distance, large_island,
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

tiff(filename = "Results/Amphibians/Hyp/distance_area_arch.tif")
Ac
dev.off()

# End of script ####