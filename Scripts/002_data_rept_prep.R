##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estacion Biologica de Donana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 17th of September 2019
# Modified July 2020

##############################################################
# Description of script and instructions
##############################################################

# This script is to prepare the dataset for reptiles, including adding a common control,  
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

#load data
reptdata<-read.csv("Data/reptdata.csv", header = TRUE, stringsAsFactors = F)
allometry<-read.csv("Data/allometric_relationships.csv", header = TRUE, stringsAsFactors = F)
diet1<-read.csv("Data/R_Traits_1.csv", stringsAsFactors = F) #Reptile traits 
diet2<-read.csv("Data/R_Traits_2.csv", stringsAsFactors = F) #Reptile traits 

# create observation level identifier
reptdata$ID<-paste0("ES",1:nrow(reptdata))

p <- reptdata %>% 
  mutate(id = 1:nrow(reptdata)) %>% 
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
            )
  ) %>% 
  ggplot() +
  aes(x = log10(Mass_m_allom2), y = log10(Mass_i_allom2)) +
  geom_point() + geom_smooth(method= "lm") + geom_abline(intercept= 0, slope = 1, colour ="red")
p
ggplotly(p)

#join mass allometry with original data
reptdata<-left_join(reptdata,p$data[,c("ID", "Mass_m_allom2",  "Mass_i_allom2")], by ="ID")
reptdata$Mass_m_allom2<-ifelse(reptdata$Measure == "Weight", reptdata$Mean_m, reptdata$Mass_m_allom2)
reptdata$Mass_i_allom2<-ifelse(reptdata$Measure == "Weight", reptdata$Mean_i, reptdata$Mass_i_allom2)

#log-transform mainland mass
reptdata$logmass<-log10(reptdata$Mass_m_allom2)

#calculate sd_allom based on multiplicative factor
reptdata$allom_factor_m<-reptdata$Mass_m_allom2/reptdata$Mean_m
reptdata$allom_factor_i<-reptdata$Mass_i_allom2/reptdata$Mean_i

reptdata$sd_m_allom<-reptdata$sd_m*reptdata$allom_factor_m
reptdata$sd_i_allom<-reptdata$sd_i*reptdata$allom_factor_i

#calculate effect size and sampling variance for allometric data
reptdata$RR_allom<-log(reptdata$Mass_i_allom2/reptdata$Mass_m_allom2)
reptdata$var_allom<-ifelse(is.na(reptdata$sd_m_allom) | is.na(reptdata$sd_i_allom),NA, 
                          (reptdata$sd_m_allom)^2/(reptdata$N_m*(reptdata$Mass_m_allom2)^2) +
                          (reptdata$sd_i_allom)^2/(reptdata$N_i*(reptdata$Mass_i_allom2)^2)) 

##DATA IMPUTATION####
#impute SD
impute_missingness(reptdata) #11.5 sd_m 10.94 sd_i % missing
data_imp<-impute_SD(reptdata,columnSDnames= c("sd_m_allom", "sd_i_allom"),columnXnames=c("Mass_m_allom2", "Mass_i_allom2"), method="Bracken1992")

summary(data_imp$sd_i_allom)
summary(data_imp$sd_m_allom)

data_imp$var_imp<-ifelse(c(is.na(data_imp$sd_m_allom) | is.na(data_imp$sd_i_allom)),NA,
                         (data_imp$sd_m_allom)^2/(data_imp$N_m*(data_imp$Mass_m_allom2)^2) +
                          (data_imp$sd_i_allom)^2/(data_imp$N_i*(data_imp$Mass_i_allom2)^2))  

impute_missingness(data_imp)

#create common control
data_imp$Shared_control<-paste0(data_imp$Reference,"_",data_imp$Mainland,"_", data_imp$Binomial)
uni_shared<-data.frame(Shared_control = unique(data_imp$Shared_control), CommonControl = paste0("CC",1:length(unique(data_imp$Shared_control))))
uni_shared$Shared_control <-as.character(uni_shared$Shared_control)
reptdata_temp<-inner_join(data_imp,uni_shared, by = "Shared_control")

#Join diet database
reptdata_temp$Binomial<-as.character(reptdata_temp$Binomial)

diet1<-diet1[,c("species","diet")]
diet2<-diet2[,c("Binomial","diet")]
diet1$species<-as.character(diet1$species) #first source
diet2$Binomial<-as.character(diet2$Binomial) #second source 

diet<-full_join(diet2,diet1, by = c("Binomial" = "species"))

#remove species for which diet is not available from any of the datasets
diet_temp<-diet[!c(is.na(diet$diet.x) & is.na(diet$diet.y)),] #3722 species

reptdata_temp<-left_join(reptdata_temp, diet_temp, by="Binomial") #join diet data
reptdata_temp$Diet<-as.character(reptdata_temp$Diet) #need to convert to character so that the ifelse function works properly
reptdata_temp$diet.x<-as.character(reptdata_temp$diet.x)
reptdata_temp$diet.y<-as.character(reptdata_temp$diet.y)

reptdata_temp$Diet<-ifelse(reptdata_temp$Diet !="",reptdata_temp$Diet,ifelse(is.na(reptdata_temp$diet.x), reptdata_temp$diet.y, reptdata_temp$diet.x)) #assign a single diet

# nrow(reptdata_temp[is.na(reptdata_temp$Diet),]) #19 species without diet assigned
# View(reptdata_temp[is.na(reptdata_temp$Diet),]) #19 species without diet assigned

###Assign diet based on literature (http://reptile-database.reptarium.cz/)
reptdata_temp[reptdata_temp$Binomial=="Alligator mississippiensis", "Diet"] <- "Carnivorous"
reptdata_temp[reptdata_temp$Binomial=="Aspidoscelis costata", "Diet"] <- "Carnivorous" #Aspidoscelis costatus
reptdata_temp[reptdata_temp$Binomial=="Hierophis viridiflavus", "Diet"] <- "Carnivorous"
reptdata_temp[reptdata_temp$Binomial=="Drysdalia coronoides", "Diet"] <- "Carnivorous"
reptdata_temp[reptdata_temp$Binomial=="Elapognathus coronatus", "Diet"] <- "Carnivorous"
reptdata_temp[reptdata_temp$Binomial=="Podarcis taurica", "Diet"] <- "Carnivorous" #Podarcis tauricus
reptdata_temp[reptdata_temp$Binomial=="Eumeces okadae", "Diet"] <- "Carnivorous"
reptdata_temp[reptdata_temp$Binomial=="Microlophus albermalensis", "Diet"] <- "Carnivorous"
reptdata_temp[reptdata_temp$Binomial=="Trimeresurus stejnegeri", "Diet"] <- "Carnivorous"
reptdata_temp[reptdata_temp$Binomial=="Hemidactylus minutus", "Diet"] <- "Carnivorous"
reptdata_temp[reptdata_temp$Binomial=="Hemidactylus masirahensis", "Diet"] <- "Carnivorous"
reptdata_temp[reptdata_temp$Binomial=="Oocatochus rufodorsatus", "Diet"] <- "Carnivorous"
reptdata_temp[reptdata_temp$Binomial=="Ptyas korros", "Diet"] <- "Carnivorous"
reptdata_temp[reptdata_temp$Binomial=="Gloydius saxatilis", "Diet"] <- "Carnivorous"

reptdata_temp$guild<-reptdata_temp$Diet #change name to common term across databases
                                       

#keep only what we need
reptdata_temp$RR<-reptdata_temp$RR_allom #change names to simplify stuff
reptdata_temp$var<-reptdata_temp$var_imp #change names to simplify stuff
reptdata_temp$Mean_m<-reptdata_temp$Mass_m_allom2 #change names to simplify stuff
reptdata_temp$Mean_i<-reptdata_temp$Mass_i_allom2 #change names to simplify stuff
reptdata_temp$sd_m<-reptdata_temp$sd_m_allom #change names to simplify stuff
reptdata_temp$sd_i<-reptdata_temp$sd_i_allom #change names to simplify stuff

reptdata_def<-reptdata_temp[,c("Reference", "ID","CommonControl", "Mainland","Island", "Class", "Order","Family",  
                               "Binomial","Species_main","Species_island", "guild", "Sex", "Measure",
                               "Mean_m","Mean_i","sd_m","sd_i","N_m", "N_i", 
                               "RR","var", "Long_i", "Lat_i", "logmass", "Island_km2", 
                               "Dist_near_mainland", "NDVI", "SDNDVI", "tmean", "tseas", "prec", "Phylogeny", "Data_source_type")]


write.csv(reptdata_def,file= "Data/reptdata_def.csv", row.names = FALSE) #455 rows

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/data_prep_rept_R_session.txt")
sessionInfo()
sink()

