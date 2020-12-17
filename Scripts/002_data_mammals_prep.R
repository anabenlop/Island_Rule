##############################################################
# Authors: 
# Ana Benítez-López (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 17th of September 2019
# Modified August 2020, and November 2020

##############################################################
# Description of script and instructions
##############################################################

# This script is to build prepare the dataset for mammalis, including adding a common control,  
# and imputing the SD that are missing for the paper:

# Benitez-Lopez et al. The island rule explains consistent patterns of 
# body size evolution

##############################################################
# Packages needed
##############################################################
#load libraries
library(metagear)
library(metafor)
library(ggplot2)
library(plotly)
library(janitor)
library(dplyr)

# Clear memory
 rm(list=ls())

##############################################################
# Importing datasets
##############################################################

#load data
mamdata<-read.csv("Data/mamdata.csv", header = TRUE, stringsAsFactors = FALSE)
allometry<-read.csv("Data/allometric_relationships.csv", header = TRUE, stringsAsFactors = FALSE)
diet<-read.csv("Data/M_Traits_guild.csv", stringsAsFactors = FALSE) #Mammal traits Elton

# create observation level identifier
mamdata$ID<-paste0("ES",1:nrow(mamdata))

##############################################################
# Converting to mass equivalents using allometric relationships
##############################################################

p <- mamdata %>% 
  mutate(id = 1:nrow(mamdata)) %>% 
  inner_join(allometry, by = "Allometry") %>% 
  distinct(id, .keep_all = TRUE) %>%
  transmute(ID,
            Mean_m, Mean_i,sd_m, sd_i, 
            a,b, 
            Allometry, 
            Mass_m_allom2 = case_when(
              Measure.x == "Weight" ~ Mean_m,
              log == "log" ~ exp(a + b *log(Mean_m)),
              log == "log10" & BM_units == "g" ~ 10^(a +b * log10(Mean_m)),
              log == "log10" & BM_units == "kg" ~ 10^(a +b*log10(Mean_m)) * 1000,
              TRUE ~ 99999
            ),
            Mass_i_allom2 = case_when(
              Measure.x == "Weight" ~ Mean_i,
              log == "log" ~ exp(a + b *log(Mean_i)),
              log == "log10" & BM_units == "g" ~ 10^(a +b * log10(Mean_i)),
              log == "log10" & BM_units == "kg" ~ 10^(a +b*log10(Mean_i)) * 1000,
              TRUE ~ 99999

            )) %>% 
  ggplot() +
  aes(x = log10(Mass_m_allom2), y = log10(Mass_i_allom2)) +
  geom_point() + geom_smooth(method= "lm") + geom_abline(intercept= 0, slope = 1, colour ="red")
p
ggplotly(p)

#join mass allometry with original data
mamdata<-left_join(mamdata,p$data[,c("ID", "Mass_m_allom2",  "Mass_i_allom2")], by ="ID")
mamdata$Mass_m_allom2<-ifelse(mamdata$Measure == "Weight", mamdata$Mean_m, mamdata$Mass_m_allom2)
mamdata$Mass_i_allom2<-ifelse(mamdata$Measure == "Weight", mamdata$Mean_i, mamdata$Mass_i_allom2)

#log-transform mainland mass
mamdata$logmass<-log10(mamdata$Mass_m_allom2)

#calculate sd_allom based on multiplicative factor
mamdata$allom_factor_m<-mamdata$Mass_m_allom2/mamdata$Mean_m
mamdata$allom_factor_i<-mamdata$Mass_i_allom2/mamdata$Mean_i

mamdata$sd_m_allom<-mamdata$sd_m*mamdata$allom_factor_m
mamdata$sd_i_allom<-mamdata$sd_i*mamdata$allom_factor_i

#calculate effect size and sampling variance for allometric data
mamdata$RR_allom<-log(mamdata$Mass_i_allom2/mamdata$Mass_m_allom2)

mamdata$var_allom<-ifelse(is.na(mamdata$sd_m_allom) | is.na(mamdata$sd_i_allom),NA, 
                         (mamdata$sd_m_allom)^2/(mamdata$N_m*(mamdata$Mass_m_allom2)^2) +
                           (mamdata$sd_i_allom)^2/(mamdata$N_i*(mamdata$Mass_i_allom2)^2)) 

##################################
##Data imputation             ####
##################################
#impute SD
impute_missingness(mamdata) #20.8% sd_m 21.7% sd_i
mamdata$imputed <- ifelse(is.na(mamdata$sd_m_allom) | is.na(mamdata$sd_i_allom), "Yes", "No")
data_imp<-impute_SD(mamdata,columnSDnames= c("sd_m_allom", "sd_i_allom"),columnXnames=c("Mass_m_allom2", "Mass_i_allom2"), method="Bracken1992")

# summary(data_imp$sd_i_allom)
# summary(data_imp$sd_m_allom)

data_imp$var_imp<-ifelse(c(is.na(data_imp$sd_m_allom) | is.na(data_imp$sd_i_allom)),NA,
                        (data_imp$sd_m_allom)^2/(data_imp$N_m*(data_imp$Mass_m_allom2)^2) +
                        (data_imp$sd_i_allom)^2/(data_imp$N_i*(data_imp$Mass_i_allom2)^2))  
impute_missingness(data_imp)

#create common control
data_imp$Shared_control<-paste0(data_imp$Reference,"_",data_imp$Mainland,"_", data_imp$Binomial)
uni_shared<-data.frame(Shared_control = unique(data_imp$Shared_control), CommonControl = paste0("CC",1:length(unique(data_imp$Shared_control))))
uni_shared$Shared_control <-as.character(uni_shared$Shared_control)
mamdata_temp<-inner_join(data_imp,uni_shared, by = "Shared_control")

#Join diet database
mamdata_temp<-left_join(mamdata_temp, diet[,c("Species","guild", "BM")], by=c("Binomial" = "Species"))
nrow(mamdata_temp[is.na(mamdata_temp$guild),]) #3 rows species with diet not assigned
mamdata_temp[mamdata_temp$Binomial == "Oncifelis guigna", "guild"] <- "Carn" # Leopardus guigna in EltonTraits
mamdata_temp[mamdata_temp$Binomial == "Galeopterus variegatus", "guild"] <- "Herb" # Galeopterus variegates in EltonTraits 
mamdata_temp[mamdata_temp$Binomial == "Oncifelis guigna", "BM"] <- 5.157939 # Leopardus guigna in EltonTraits
mamdata_temp[mamdata_temp$Binomial == "Galeopterus variegatus", "BM"] <- 1.1122 # Galeopterus variegates in EltonTraits
mamdata_temp[mamdata_temp$guild=="Frug","guild"] <- "Herb" # assign frugivores to herbivores category

#keep only what we need
mamdata_temp$RR<-mamdata_temp$RR_allom #change names to simplify stuff
mamdata_temp$var<-mamdata_temp$var_imp #change names to simplify stuff
mamdata_temp$Mean_m<-mamdata_temp$Mass_m_allom2 #change names to simplify stuff
mamdata_temp$Mean_i<-mamdata_temp$Mass_i_allom2 #change names to simplify stuff
mamdata_temp$sd_m<-mamdata_temp$sd_m_allom #change names to simplify stuff
mamdata_temp$sd_i<-mamdata_temp$sd_i_allom #change names to simplify stuff

mamdata_def<-mamdata_temp[,c("Reference", "ID","CommonControl", "Mainland","Island", "Class", "Order","Family",  
                             "Binomial","Species_main","Species_island", "guild", "Sex", "Measure",
                             "Mean_m","Mean_i","sd_m","sd_i","N_m", "N_i", 
                             "RR","var", "Long_i", "Lat_i", "logmass", "Island_km2", 
                             "Dist_near_mainland", "NDVI", "SDNDVI", "tmean", "tseas", "prec","Archipielago",
                             "Phylogeny", "Data_source_type", "imputed")]

# save data
write.csv(mamdata_def,file= "Data/mamdata_def.csv", row.names = FALSE) 

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/data_mamprep_R_session.txt")
sessionInfo()
sink()

# End of script ####