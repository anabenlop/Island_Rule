##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estacion Biologica de Donana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 17th of September 2019
# Modified August 2020

##############################################################
# Description of script and instructions
##############################################################

# This script is to prepare the dataset for birds, including adding a common control,  
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
birddata<-read.csv("Data/birddata.csv", header = TRUE, stringsAsFactors = FALSE)
allometry<-read.csv("Data/allometric_relationships.csv", header = TRUE, stringsAsFactors = FALSE)
diet<-read.csv("Data/B_traits_guild.csv", stringsAsFactors = FALSE) #Bird traits Elton
mig_status<-read.csv("Data/SpeciesList3_1_migbehav_v2_0.csv", stringsAsFactors = FALSE) #Bird migratory status:  Eyres et al. 2017 https://onlinelibrary.wiley.com/doi/pdf/10.1111/jav.01308
mig_status$Binomial <- paste0(mig_status$IOC3_1_Genus," ", mig_status$IOC3_1_Species)

# create observation level identifier
birddata$ID<-paste0("ES",1:nrow(birddata))

p <- birddata %>% 
  mutate(id = 1:nrow(birddata)) %>% 
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
birddata<-left_join(birddata,p$data[,c("ID", "Mass_m_allom2",  "Mass_i_allom2")], by ="ID")
birddata$Mass_m_allom2<-ifelse(birddata$Measure == "Weight", birddata$Mean_m, birddata$Mass_m_allom2)
birddata$Mass_i_allom2<-ifelse(birddata$Measure == "Weight", birddata$Mean_i, birddata$Mass_i_allom2)

#log-transform mainland mass
birddata$logmass<-log10(birddata$Mass_m_allom2)

#calculate sd_allom based on multiplicative factor
birddata$allom_factor_m<-birddata$Mass_m_allom2/birddata$Mean_m
birddata$allom_factor_i<-birddata$Mass_i_allom2/birddata$Mean_i

birddata$sd_m_allom<-birddata$sd_m*birddata$allom_factor_m
birddata$sd_i_allom<-birddata$sd_i*birddata$allom_factor_i

#calculate effect size and sampling variance for allometric data
birddata$RR_allom<-log(birddata$Mass_i_allom2/birddata$Mass_m_allom2)
birddata$var_allom<-ifelse(is.na(birddata$sd_m_allom) | is.na(birddata$sd_i_allom),NA, 
                          (birddata$sd_m_allom)^2/(birddata$N_m*(birddata$Mass_m_allom2)^2) +
                          (birddata$sd_i_allom)^2/(birddata$N_i*(birddata$Mass_i_allom2)^2))


##DATA IMPUTATION####
#impute SD
impute_missingness(birddata) #0.93 % sd_m and 1.1% sd_i missing
data_imp<-impute_SD(birddata,columnSDnames= c("sd_m_allom", "sd_i_allom"),columnXnames=c("Mass_m_allom2", "Mass_i_allom2"), method="Bracken1992")

summary(data_imp$sd_i_allom)
summary(data_imp$sd_m_allom)

data_imp$var_imp<-ifelse(c(is.na(data_imp$sd_m_allom) | is.na(data_imp$sd_i_allom)),NA,
                  (data_imp$sd_m_allom)^2/(data_imp$N_m*(data_imp$Mass_m_allom2)^2) +
                  (data_imp$sd_i_allom)^2/(data_imp$N_i*(data_imp$Mass_i_allom2)^2))   

impute_missingness(data_imp)

#remove species for which both tarsus and wing has been measured, remove only the tarsus measures except for Apteryx
data_imp_temp<-data_imp[!(data_imp$Reference == "Pigot et al. 2020" & data_imp$Measure == "Tarsus Length"),]
data_imp_temp<-rbind(data_imp_temp, data_imp[data_imp$Binomial == "Apteryx australis",])
data_imp<-data_imp_temp

#create common control
data_imp$Shared_control<-paste0(data_imp$Reference,"_",data_imp$Mainland,"_", data_imp$Binomial)
uni_shared<-data.frame(Shared_control = unique(data_imp$Shared_control), CommonControl = paste0("CC",1:length(unique(data_imp$Shared_control))))
uni_shared$Shared_control <-as.character(uni_shared$Shared_control)
birddata_temp<-inner_join(data_imp,uni_shared, by = "Shared_control")

#Join diet database
diet <- diet[!duplicated(diet$Species),] #remove duplicates
birddata_temp<-left_join(birddata_temp, diet[,c("Species","Diet_cat")], by=c("Binomial" = "Species"))
birddata_temp[is.na(birddata_temp$Diet_cat),] #all species with diet assigned
birddata_temp$guild<-birddata_temp$Diet_cat #change name to common term across databases

#Join migratory status database
birddata_temp <- left_join(birddata_temp, mig_status[,c("Binomial","Migratory_status", "Migratory_status_3")], by=c("Binomial" = "Binomial"))
nrow(birddata_temp[is.na(birddata_temp$Migratory_status),]) #30 species without mig status assigned
nrow(birddata_temp[is.na(birddata_temp$Migratory_status_3),]) #30 species without mig status assigned

# fix species without mig status
birddata_temp[birddata_temp$Binomial == "Aethopyga latouchii", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") #based on christinae, which includes latouchii
birddata_temp[birddata_temp$Binomial == "Alcedo azurea", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # recorded as Ceyx azureus 
birddata_temp[birddata_temp$Binomial == "Alcippe cinereiceps", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # All members in genus Alcippe are resident
birddata_temp[birddata_temp$Binomial == "Psittacara leucophthalmus", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # recorded as Aratinga leucophtalma
birddata_temp[birddata_temp$Binomial == "Eupsittula astec", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Includes Aratinga nana/astec
birddata_temp[birddata_temp$Binomial == "Zanda funerea", c("Migratory_status", "Migratory_status_3")] <- c("resident","partial resident") # Recorded as Calyptorhynchus funereus
birddata_temp[birddata_temp$Binomial == "Synoicus ypsilophorus", c("Migratory_status", "Migratory_status_3")] <- c("resident","partial resident") # Recorded as Coturnix ypsilophora
birddata_temp[birddata_temp$Binomial == "Dryobates scalaris", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Picoides scalaris
birddata_temp[birddata_temp$Binomial == "Dyaphorophyia castanea", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Platysteira castanea
birddata_temp[birddata_temp$Binomial == "Dyaphorophyia chalybea", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Platysteira chalybea
birddata_temp[birddata_temp$Binomial == "Garrulax elliotii", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Trochalopteron elliotii
birddata_temp[birddata_temp$Binomial == "Haemorhous mexicanus", c("Migratory_status", "Migratory_status_3")] <- c("directional migratory", "partial directional migrant") # Recorded as Carpodacus mexicanus
birddata_temp[birddata_temp$Binomial == "Hemicircus sordidus", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as	Hemicircus concretus/sordidus
birddata_temp[birddata_temp$Binomial == "Nesoptilotis leucotis", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Lichenostomus leucotis
birddata_temp[birddata_temp$Binomial == "Macronous kelleyi", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Macronus kelleyi
birddata_temp[birddata_temp$Binomial == "Macronous ptilosus", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Macronus ptilosus
birddata_temp[birddata_temp$Binomial == "Meiglyptes grammithorax", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Meiglyptes tristis
birddata_temp[birddata_temp$Binomial == "Pardaliparus venustulus", c("Migratory_status", "Migratory_status_3")] <- c("dispersive migratory", "partial dispersive migrant") # Recorded as Periparus venustulus
birddata_temp[birddata_temp$Binomial == "Zoothera interpres", c("Migratory_status", "Migratory_status_3")] <- c("resident","partial resident") # Recorded as Geockila interpres
birddata_temp[birddata_temp$Binomial == "Turdus libonyanus", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Turdus libonyana
birddata_temp[birddata_temp$Binomial == "Trogon ambiguus", c("Migratory_status", "Migratory_status_3")] <- c("directional migratory", "partial directional migrant") # Included in Trogon elegans
birddata_temp[birddata_temp$Binomial == "Passer italiae", c("Migratory_status", "Migratory_status_3")] <- c("directional migratory", "partial directional migrant") # Included in Passer hispaniolensis
birddata_temp[birddata_temp$Binomial == "Sephanoides sephaniodes", c("Migratory_status", "Migratory_status_3")] <- c("directional migratory", "partial directional migrant") # Recorded as Sephanoides sephanoides
birddata_temp[birddata_temp$Binomial == "Spinus psaltria", c("Migratory_status", "Migratory_status_3")] <- c("directional migratory", "partial directional migrant") # Recorded as Carduelis psaltria
birddata_temp[birddata_temp$Binomial == "Trichastoma tickelli", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Recorded as Pellorneum tickelli
birddata_temp[birddata_temp$Binomial == "Trichoglossus meyeri", c("Migratory_status", "Migratory_status_3")] <- c("resident","full resident") # Included in Trichoglossus flavoviridis


# table(birddata_temp$Migratory_status)
# table(birddata_temp$Migratory_status_3)

# Remove marine species
birddata_temp <- birddata_temp[birddata_temp$Migratory_status_3 != "marine",]
# Remove full migratory species
birddata_temp <- birddata_temp[birddata_temp$Migratory_status_3 != "full directional migrant",]

#keep only what we need
birddata_temp$RR<-birddata_temp$RR_allom #change names to simplify stuff
birddata_temp$var<-birddata_temp$var_imp #change names to simplify stuff
birddata_temp$Mean_m<-birddata_temp$Mass_m_allom2 #change names to simplify stuff
birddata_temp$Mean_i<-birddata_temp$Mass_i_allom2 #change names to simplify stuff
birddata_temp$sd_m<-birddata_temp$sd_m_allom #change names to simplify stuff
birddata_temp$sd_i<-birddata_temp$sd_i_allom #change names to simplify stuff

birddata_def<-birddata_temp[,c("Reference", "ID","CommonControl", "Mainland","Island", "Class", "Order","Family",  
                               "Binomial","Species_main","Species_island", "guild","Mean_m","Mean_i","sd_m","sd_i","N_m", "N_i", 
                               "RR","var", "Long_i", "Lat_i", "logmass", "Island_km2", 
                               "Dist_near_mainland", "NDVI", "SDNDVI", "tmean", "tseas", "prec", "Migratory_status", "Migratory_status_3",
                               "Phylogeny", "Data_source_type")]

write.csv(birddata_def,file= "Data/birddata_def.csv", row.names = FALSE) #727 rows

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/data_birdprep_R_session.txt")
sessionInfo()
sink()

# End of script  ####
