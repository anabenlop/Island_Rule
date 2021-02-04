##############################################################
# Authors: 
# Ana Benítez-López (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 20th of July 2020

##############################################################
# Description of script and instructions                  ####
##############################################################

# This script is to run alternative models to test the island rule
# for the paper: 

# Benítez-López et al. The island rule explains consistent patterns of 
# body size evolution across terrestrial vertebrates. 


##############################################################
# Packages needed                                         ####
##############################################################
library(metafor)
library(ggplot2)
library(ggpubr)
library(car)
library(plotly)
library(png)
library(ggimage)
library(grid)
library(rsvg)
library(grImport2)

#clean memory
# rm(list=ls())

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

# load functions
source("Scripts/000_Functions.R")

# load models
# if you already ran the models and just want to modify the figures
metamam_alt <- readRDS(file = "Data/Final data/metamam_alt.Rdata")
metabird_alt <- readRDS(file = "Data/Final data/metabird_alt.Rdata")
metarept_alt <- readRDS(file = "Data/Final data/metarept_alt.Rdata")
metaamph_alt <- readRDS(file = "Data/Final data/metaamph_alt.Rdata")

# Transform data and calculate sampling variances
mamdata$logMean_i <- log(mamdata$Mean_i)
mamdata$logMean_m <- log(mamdata$Mean_m)
mamdata$var_i <- mamdata$sd_i^2/(mamdata$Mean_i^2*mamdata$N_i)
mamdata$var_m <- mamdata$sd_m^2/(mamdata$Mean_m^2*mamdata$N_m)

birddata$logMean_i <- log(birddata$Mean_i)
birddata$logMean_m <- log(birddata$Mean_m)
birddata$var_i <- birddata$sd_i^2/(birddata$Mean_i^2*birddata$N_i)
birddata$var_m <- birddata$sd_m^2/(birddata$Mean_m^2*birddata$N_m)

reptdata$logMean_i <- log(reptdata$Mean_i)
reptdata$logMean_m <- log(reptdata$Mean_m)
reptdata$var_i <- reptdata$sd_i^2/(reptdata$Mean_i^2*reptdata$N_i)
reptdata$var_m <- reptdata$sd_m^2/(reptdata$Mean_m^2*reptdata$N_m)

amphdata$logMean_i <- log(amphdata$Mean_i)
amphdata$logMean_m <- log(amphdata$Mean_m)
amphdata$var_i <- amphdata$sd_i^2/(amphdata$Mean_i^2*amphdata$N_i)
amphdata$var_m <- amphdata$sd_m^2/(amphdata$Mean_m^2*amphdata$N_m)

#######################################################################
# Phylogenetic meta-regression models island mass vs mainland mass ####
#######################################################################
RE = list(~ 1 | Reference,~1|ID, ~1|SPID, ~1| Binomial)

#mammals####
phylocor<-list(Binomial= mam_phylo_cor)
# metamam_alt <- rma.mv(logMean_i ~ logMean_m, V= var_i, subset = sd_i != 0, data=mamdata, random= RE,
#                       R = phylocor, method = "REML")
# summary(metamam_alt)

# Compute t-student H0: slope=1.Seen in:  https://stackoverflow.com/questions/33060601/test-if-the-slope-in-simple-linear-regression-equals-to-a-given-constant-in-r  
tstats <- (1-metamam_alt$beta[2,1])/metamam_alt$se[2]
# Calculates two tailed probability
pval<- 2 * pt(abs(tstats), df = df.residual(metamam_alt), lower.tail = FALSE)
print(pval) # t-test: 5.19 , p: 2.478329e-07

# Alternatively use LRT, or Wald-type test from car package
linearHypothesis(metamam_alt, "logMean_m = 1") # Chisq: 28.23   1.077e-07 ***

coef<-round(data.frame(b =metamam_alt$b, lci = metamam_alt$ci.lb, uci =  metamam_alt$ci.ub), digits = 3)
int <- paste0("a = ", coef$b[1], " (", coef$lci[1], " - ", coef$uci[1], ")")
slo <- paste0("b = ", coef$b[2], " (", coef$lci[2], " - ", coef$uci[2], ")")

int_text<-annotate(geom="text", x= 10.3, y= 1.8, label= int, size = 4)
slo_text<-annotate(geom="text", x= 10.3, y= 1, label= slo, size = 4)

logMean_m <- seq(from = min(mamdata$logMean_m), to = max(mamdata$logMean_m), length.out = 1000)

#predict for vector of mass
df_m<-predict(metamam_alt, newmods = cbind(logMean_m), addx=TRUE)
df_m<-data.frame(df_m)
df_m$logMean_m<-df_m$X.logMean_m
df_m$Mean_m<-exp(df_m$logMean_m)/1000
df_m$Mean_i_pred<-exp(df_m$pred)/1000

# import silhouette
# raster format
sil_M <- readPNG("Silhouettes/PhyloPic.72f2f997.Steven-Traver.Cervus-elaphus.png")

Malt<-ggplot(mamdata)+ geom_abline(intercept = 0, slope = 1, col = "dark gray", linetype = "dashed",  size = 0.8)+ 
  annotation_custom(rasterGrob(sil_M,
                               x = unit(0.14, "npc"),
                               y = unit(0.85, "npc"),
                               width = unit(0.22,"npc"),
                               height = unit(0.27,"npc")),
                               -Inf, Inf, -Inf, Inf) +theme_bw(base_size=18) +
  geom_line(data=df_m,aes(logMean_m, pred),color="#0072B2", size = 1)+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+xlab("ln(mass mainland (g))")+ ylab("ln(mass island (g))")+ 
  scale_x_continuous(breaks=seq(1,13,2), limits= c(1,13)) +
  scale_y_continuous(breaks=seq(1,13,2), limits= c(1,13)) +
  int_text + slo_text +
  labs(tag = "a")
# Malt

# back-transformed (in raw scale)
int_text<-annotate(geom="text", x= 200, y= 12, label= int, size = 4)
slo_text<-annotate(geom="text", x= 200, y= 1, label= slo, size = 4)

# for species between 1-100 kg
# Malt_bt<-ggplot(mamdata)+ geom_abline(intercept = 0, slope = 1, col = "dark gray", linetype = "dashed",  size = 0.8)+
#   theme_bw(base_size=18) +
#   geom_line(data=df_m,aes(Mean_m, Mean_i_pred),color="#0072B2", size = 1)+
#   theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+xlab("Mass mainland (kg))")+ ylab("Mass island (kg))")+
#   scale_x_continuous(breaks=seq(1,100,10), limits= c(1,100)) +
#   scale_y_continuous(breaks=seq(1,100,10), limits= c(1,100)) +
#   #int_text + slo_text +
#   labs(tag = "a")
# Malt_bt
# 
# ggplotly(Malt_bt)
# 
# 
# # for species between 10-200 g
# Malt_bt_sm<-ggplot(mamdata)+ geom_abline(intercept = 0, slope = 1, col = "dark gray", linetype = "dashed",  size = 0.8)+
#   theme_bw(base_size=18) +
#   geom_line(data=df_m,aes(Mean_m*1000, Mean_i_pred*1000),color="#0072B2", size = 1)+
#   theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+xlab("Mass mainland (g))")+ ylab("Mass island (g))")+
#   scale_x_continuous(breaks=seq(10,200,10), limits= c(10,200)) +
#   scale_y_continuous(breaks=seq(10,200,10), limits= c(10,200)) +
#   #↨int_text + slo_text +
#   labs(tag = "a")
# Malt_bt_sm
# 
# ggplotly(Malt_bt_sm)

#birds ####
phylocor<-list(Binomial= bird_phylo_cor)
# metabird_alt <- rma.mv(logMean_i ~ logMean_m, V= var_i, subset = sd_i != 0,  data=birddata, random= RE,
#                        R = phylocor, method = "REML")
# summary(metabird_alt) 

# Compute t-student H0: slope=1.Seen in:  https://stackoverflow.com/questions/33060601/test-if-the-slope-in-simple-linear-regression-equals-to-a-given-constant-in-r  
tstats <- (1-metabird_alt$beta[2,1])/metabird_alt$se[2]
# Calculates two tailed probability
pval<- 2 * pt(abs(tstats), df = df.residual(metabird_alt), lower.tail = FALSE)
print(pval) # t-test:  5.010759 , p: 6.898717e-07 

# Use LRT, or Wald-type test from car package
linearHypothesis(metabird_alt, "logMean_m = 1") # Chisq: 24.137  8.974e-07 ***

coef<-round(data.frame(b =metabird_alt$b, lci = metabird_alt$ci.lb, uci =  metabird_alt$ci.ub), digits = 3)
int <- paste0("a = ", coef$b[1], " (", coef$lci[1], " - ", coef$uci[1], ")")
slo <- paste0("b = ", coef$b[2], " (", coef$lci[2], " - ", coef$uci[2], ")")

int_text<-annotate(geom="text", x= 7.2, y= 1.5, label= int, size = 4)
slo_text<-annotate(geom="text", x= 7.2, y= 1, label= slo, size = 4)

logMean_m <- seq(from = min(birddata$logMean_m), to = max(birddata$logMean_m), length.out = 1000)

#predict for vector of mass
df_b<-predict(metabird_alt, newmods = cbind(logMean_m), addx=TRUE)
df_b<-data.frame(df_b)
df_b$logMean_m<-df_b$X.logMean_m
df_b$Mean_m<-exp(df_b$logMean_m)/1000
df_b$Mean_i_pred<-exp(df_b$pred)/1000

# import silhouette
# raster format
sil_B<- readPNG("Silhouettes/PhyloPic.67a9ecfd.Sylviidae_Sylvioidea_PublicDom1.0_flipped.png")

Balt<-ggplot(birddata)+ geom_abline(intercept = 0, slope = 1, col = "dark gray", linetype = "dashed",  size = 0.8)+ 
  annotation_custom(rasterGrob(sil_B,
                               x = unit(0.14, "npc"),
                               y = unit(0.85, "npc"),
                               width = unit(0.17,"npc"),
                               height = unit(0.16,"npc")),
                                -Inf, Inf, -Inf, Inf) +theme_bw(base_size=18) +
  geom_line(data=df_b,aes(logMean_m, pred),color="#CC0000", size = 1.2)+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+xlab("ln(mass mainland (g))")+ ylab("ln(mass island (g))")+ 
  scale_x_continuous(breaks=seq(1,9,2),limits= c(1,9)) + 
  scale_y_continuous(breaks=seq(1,9,2),limits= c(1,9)) + 
  int_text + slo_text + 
  labs(tag = "b")
# Balt

# back-transformed
# int_text<-annotate(geom="text", x= 3, y= 0.5, label= int, size = 4)
# slo_text<-annotate(geom="text", x= 3, y= 0.2, label= slo, size = 4)
# 
# Balt<-ggplot(birddata)+ geom_abline(intercept = 0, slope = 1, col = "dark gray", linetype = "dashed",  size = 0.8)+ 
#   theme_bw(base_size=18) +
#   geom_line(data=df_b,aes(Mean_m, Mean_i_pred),color="#CC0000", size = 1)+
#   theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+xlab("Mass mainland (kg))")+ ylab("Mass island (kg))")+ 
#   scale_x_continuous(breaks=seq(0,4,1), limits= c(0,4.1)) +
#   scale_y_continuous(breaks=seq(0,4,1), limits= c(0,4.1)) +
#   int_text + slo_text +
#   labs(tag = "b")
# Balt
# 
# ggplotly(Balt)

# reptiles ####
phylocor<-list(Binomial= rept_phylo_cor)
# metarept_alt <- rma.mv(logMean_i ~ logMean_m, V= var_i,  data=reptdata, random= RE,
#                        R = phylocor, method = "REML")
# summary(metarept_alt) 

# Compute t-student H0: slope=1.Seen in:  https://stackoverflow.com/questions/33060601/test-if-the-slope-in-simple-linear-regression-equals-to-a-given-constant-in-r  
tstats <- (1-metarept_alt$beta[2,1])/metarept_alt$se[2]
# Calculates two tailed probability
pval<- 2 * pt(abs(tstats), df = df.residual(metarept_alt), lower.tail = FALSE)
print(pval) # t-test:  5.188 , p: 3.002474e-07 

# Use LRT, or Wald-type test from car package
linearHypothesis(metarept_alt, "logMean_m = 1") # Chisq: 28.316  1.031e-07 ***

coef<-round(data.frame(b =metarept_alt$b, lci = metarept_alt$ci.lb, uci =  metarept_alt$ci.ub), digits = 3)
int <- paste0("a = ", coef$b[1],"0", " (", coef$lci[1], " - ", coef$uci[1], ")") #had to manually add a zero because the rounding does not allow decimals to end in zero
slo <- paste0("b = ", coef$b[2], " (", coef$lci[2], " - ", coef$uci[2], ")")

int_text<-annotate(geom="text", x= 9, y= -1, label= int, size = 4)
slo_text<-annotate(geom="text", x= 9, y= -2, label= slo, size = 4)

logMean_m <- seq(from = min(reptdata$logMean_m), to = max(reptdata$logMean_m), length.out = 1000)

#predict for vector of mass
df_r<-predict(metarept_alt, newmods = cbind(logMean_m), addx=TRUE)
df_r<-data.frame(df_r)
df_r$logMean_m<-df_r$X.logMean_m
df_r$Mean_m<-exp(df_r$logMean_m)/1000
df_r$Mean_i_pred<-exp(df_r$pred)/1000

# import silhouette
# raster format
sil_R<- readPNG("Silhouettes/Steven-Traver.Varanus_Varanus.png")

Ralt<-ggplot(reptdata)+ geom_abline(intercept = 0, slope = 1, col = "dark gray", linetype = "dashed",  size = 0.8)+ 
  annotation_custom(rasterGrob(sil_R,
                               x = unit(0.16, "npc"),
                               y = unit(0.85, "npc"),
                               width = unit(0.26,"npc"),
                               height = unit(0.16,"npc")),
                              -Inf, Inf, -Inf, Inf) +theme_bw(base_size=18) +
  geom_line(data=df_r,aes(logMean_m, pred),color="#E69F00", size = 1.2)+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+xlab("ln(mass mainland (g))")+ ylab("ln(mass island (g))")+ 
  scale_x_continuous(breaks=seq(-2,12,2),limits= c(-2,12)) +
  scale_y_continuous(breaks=seq(-2,12,2),limits= c(-2,12)) +
  int_text + slo_text +
  labs(tag = "c")
# Ralt

# back-transformed
# int_text<-annotate(geom="text", x= 20, y= 3, label= int, size = 4)
# slo_text<-annotate(geom="text", x= 20, y= 1, label= slo, size = 4)
# 
# #calculate size of points
# wi    <- 1/sqrt(reptdata$var_i)
# size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))
# 
# Ralt_bt<-ggplot(reptdata)+ geom_abline(intercept = 0, slope = 1, col = "dark gray", linetype = "dashed",  size = 0.8)+
#   theme_bw(base_size=18) +
#   geom_line(data=df_r,aes(Mean_m*1000, Mean_i_pred*1000),color="#E69F00", size = 1)+
#   theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+xlab("Mass mainland (g))")+ ylab("Mass island (g))")+
#   scale_x_continuous(breaks=seq(0,500,50), limits= c(0,1)) +
#   scale_y_continuous(breaks=seq(0,500,50), limits= c(0,1)) +
#   int_text + slo_text +
#   labs(tag = "c")
# Ralt_bt
# ggplotly(Ralt_bt)

#amphibians ####
phylocor<-list(Binomial= amph_phylo_cor)
# metaamph_alt <- rma.mv(logMean_i ~ logMean_m, V= var_i,  data=amphdata, random= RE,
#                        R = phylocor, method = "REML")
# summary(metaamph_alt) 

# Compute t-student H0: slope=1.Seen in:  https://stackoverflow.com/questions/33060601/test-if-the-slope-in-simple-linear-regression-equals-to-a-given-constant-in-r  
tstats <- (1-metaamph_alt$beta[2,1])/metaamph_alt$se[2]
# Calculates two tailed probability
pval<- 2 * pt(abs(tstats), df = df.residual(metaamph_alt), lower.tail = FALSE)
print(pval) # t-test:   0.9615736 , p: 0.3375756

# Use LRT, or Wald-type test from car package
linearHypothesis(metaamph_alt, "logMean_m = 1") # Chisq: 1.1505     0.2834

coef<-round(data.frame(b =metaamph_alt$b, lci = metaamph_alt$ci.lb, uci =  metaamph_alt$ci.ub), digits = 3)
int <- paste0("a = ", coef$b[1], " (", coef$lci[1], " - ", coef$uci[1], ")")
slo <- paste0("b = ", coef$b[2], " (", coef$lci[2], " - ", coef$uci[2], ")")

int_text<-annotate(geom="text", x= 2.8, y= -1.2, label= int, size = 4)
slo_text<-annotate(geom="text", x= 2.8, y= -1.5, label= slo, size = 4)

logMean_m <- seq(from = min(amphdata$logMean_m), to = max(amphdata$logMean_m), length.out = 1000)

#predict for vector of mass
df_a<-predict(metaamph_alt, newmods = cbind(logMean_m), addx=TRUE)
df_a<-data.frame(df_a)
df_a$logMean_m<-df_a$X.logMean_m
df_a$Mean_m<-exp(df_a$logMean_m)
df_a$Mean_i_pred<-exp(df_a$pred)

# import silhouette
# raster format
sil_A<- readPNG("Silhouettes/Will-Booker.Hyla-versicolor_CC0.1.0.png")

Aalt<-ggplot(reptdata)+ geom_abline(intercept = 0, slope = 1, col = "dark gray", linetype = "dashed",  size = 0.8)+ 
  annotation_custom(rasterGrob(sil_A,
                               x = unit(0.14, "npc"),
                               y = unit(0.85, "npc"),
                               width = unit(0.16,"npc"),
                               height = unit(0.18,"npc")),
                              -Inf, Inf, -Inf, Inf) +theme_bw(base_size=18) +
  geom_line(data=df_a,aes(logMean_m, pred),color="#009E73", size = 1.2)+
  theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+xlab("ln(mass mainland (g))")+ ylab("ln(mass island (g))")+ 
  scale_x_continuous(breaks=seq(-2,4,1),limits= c(-1.5,4)) +
  scale_y_continuous(breaks=seq(-2,4,1),limits= c(-1.5,4)) +
  int_text + slo_text +
  labs(tag = "d")
# Aalt

# back-transformed
# int_text<-annotate(geom="text", x= 20, y= 3, label= int, size = 4)
# slo_text<-annotate(geom="text", x= 20, y= 1, label= slo, size = 4)
# 
# #calculate size of points
# wi    <- 1/sqrt(reptdata$var_i)
# size  <- 2 + 20.0 * (wi - min(wi))/(max(wi) - min(wi))
# 
# Aalt<-ggplot(amphdata)+ geom_abline(intercept = 0, slope = 1, col = "dark gray", linetype = "dashed",  size = 0.8)+ 
#   theme_bw(base_size=18) +
#   geom_line(data=df_a,aes(Mean_m, Mean_i_pred),color="#009E73", size = 1)+
#   theme(element_blank(), axis.text=element_text(size=18, colour ="black"))+xlab("Mass mainland (g))")+ ylab("Mass island (g))")+ 
#   scale_x_continuous(breaks=seq(0,120,20), limits= c(0,120)) +
#   scale_y_continuous(breaks=seq(0,120,20), limits= c(0,120)) +
#   int_text + slo_text +
#   labs(tag = "d")
# Aalt
# ggplotly(Aalt)

multi_alt<-ggarrange(Malt, Balt, Ralt, Aalt, ncol = 2, nrow = 2, align = "hv")
# tiff('Results/Figures/ED_Fig2.tif', res=300, width=3000, height=3000, compression = "lzw")
# multi_alt
# dev.off()

pdf('Results/Figures/ED_Fig 2.pdf', width=10, height=9.5)
multi_alt
dev.off()

# saveRDS(metamam_alt, file = "Data/Final data/metamam_alt.Rdata")
# saveRDS(metabird_alt, file = "Data/Final data/metabird_alt.Rdata")
# saveRDS(metarept_alt, file = "Data/Final data/metarept_alt.Rdata")
# saveRDS(metaamph_alt, file = "Data/Final data/metaamph_alt.Rdata")

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/alternative_models_R_session.txt")
sessionInfo()
sink()

### End of script ####
