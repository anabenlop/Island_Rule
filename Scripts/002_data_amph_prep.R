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

# This script is to prepare the dataset for amphibians, including adding a common control,  
# and imputing the SD that are missing for the paper:

# Benitez-Lopez et al. The island rule explains consistent patterns of 
# body size evolution across terrestrial vertebrates. 


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
# rm(list=ls())

#load data
amphdata<-read.csv("~/New projects/Island rule/Data/amphdata.csv", header = TRUE, stringsAsFactors = FALSE)
allometry<-read.csv("~/New projects/Island rule/Data/allometric_relationships.csv", header = TRUE, stringsAsFactors = FALSE)
#no diet data

amphdata$ID<-paste0("ES",1:nrow(amphdata))

p <- amphdata %>% 
  mutate(id = 1:nrow(amphdata)) %>% 
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
amphdata<-left_join(amphdata,p$data[,c("ID", "Mass_m_allom2",  "Mass_i_allom2")], by ="ID")
amphdata$Mass_m_allom2<-ifelse(amphdata$Measure == "Weight", amphdata$Mean_m, amphdata$Mass_m_allom2)
amphdata$Mass_i_allom2<-ifelse(amphdata$Measure == "Weight", amphdata$Mean_i, amphdata$Mass_i_allom2)

#log-transform mainland mass
amphdata$logmass<-log10(amphdata$Mass_m_allom2)

#calculate sd_allom based on multiplicative factor
amphdata$allom_factor_m<-amphdata$Mass_m_allom2/amphdata$Mean_m
amphdata$allom_factor_i<-amphdata$Mass_i_allom2/amphdata$Mean_i

amphdata$sd_m_allom<-amphdata$sd_m*amphdata$allom_factor_m
amphdata$sd_i_allom<-amphdata$sd_i*amphdata$allom_factor_i

#calculate effect size and sampling variance for allometric data
amphdata$RR_allom<-log(amphdata$Mass_i_allom2/amphdata$Mass_m_allom2)
amphdata$var_allom<-ifelse(is.na(amphdata$sd_m_allom) | is.na(amphdata$sd_i_allom),NA, 
                           (amphdata$sd_m_allom)^2/(amphdata$N_m*(amphdata$Mass_m_allom2)^2) +
                             (amphdata$sd_i_allom)^2/(amphdata$N_i*(amphdata$Mass_i_allom2)^2))

##DATA IMPUTATION####
#impute SD
impute_missingness(amphdata) #7.26% missing sd_i
data_imp<-impute_SD(amphdata,columnSDnames= c("sd_m_allom", "sd_i_allom"),columnXnames=c("Mass_m_allom2", "Mass_i_allom2"), method="Bracken1992")

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
amphdata_temp<-inner_join(data_imp,uni_shared, by = "Shared_control")

#keep only what we need
amphdata_temp$RR<-amphdata_temp$RR_allom #change names to simplify stuff
amphdata_temp$var<-amphdata_temp$var_imp #change names to simplify stuff
amphdata_temp$Mean_orig_m<-amphdata_temp$Mean_m #change names to simplify stuff
amphdata_temp$Mean_orig_i<-amphdata_temp$Mean_i #change names to simplify stuff
amphdata_temp$sd_orig_m<-amphdata_temp$sd_m #change names to simplify stuff
amphdata_temp$sd_orig_i<-amphdata_temp$sd_i #change names to simplify stuff
amphdata_temp$Mean_m<-amphdata_temp$Mass_m_allom2 #change names to simplify stuff
amphdata_temp$Mean_i<-amphdata_temp$Mass_i_allom2 #change names to simplify stuff
amphdata_temp$sd_m<-amphdata_temp$sd_m_allom #change names to simplify stuff
amphdata_temp$sd_i<-amphdata_temp$sd_i_allom #change names to simplify stuff

amphdata_def<-amphdata_temp[,c("Reference", "ID","CommonControl", "Mainland","Island", "Class", "Order","Family",  
                               "Binomial","Species_main","Species_island",  "Sex", "Measure",
                               "Mean_m","Mean_i","sd_m","sd_i","N_m", "N_i", 
                               "RR","var", "Long_i", "Lat_i", "logmass", "Island_km2", 
                               "Dist_near_mainland", "NDVI", "SDNDVI", "tmean", "tseas", "prec", "Phylogeny", "Data_source_type")] 

write.csv(amphdata_def,file= "~/New projects/Island rule/Data/amphdata_def.csv", row.names = FALSE)

# saving session information with all packages versions for reproducibility purposes
sink("~/New projects/Island rule/Data/Final data/data_prep_amph_R_session.txt")
sessionInfo()
sink()
