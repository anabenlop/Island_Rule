##############################################################
# Authors: 
# Ana Benítez-López (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estación Biológica de Doñana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 18th of September 2019, adapted from Sánchez-Tojar et al. 2019 (https://osf.io/yjua8/)

##############################################################
# Description of script and instructions
##############################################################

# This script is to build the phylogeny for reptiles and estimate the 
# phylogenetic relatedness among the species included in:

# Benítez-López et al. Assessing the generality of the island rule in tetrapods.


##############################################################
# Packages needed
##############################################################

#load libraries
library(dplyr)
library(ape)
library(treebase)
library(rotl)
library(diagram)

# Clear memory
 rm(list=ls())

##############################################################
# Importing datasets
##############################################################

# final and clean database 
reptdata<-read.csv("Data/reptdata_def.csv", header = TRUE, stringsAsFactors = FALSE)

# generating list of species
species <- sort(unique(as.character(reptdata$Binomial)))

##############################################################
# Formatting species data
##############################################################

# obtaining dataframe listing the Open Tree identifiers potentially 
# matching our list of species.

taxa <- tnrs_match_names(names = species)

# according to the `approximate_match` column, there might be 
# 0 typo in the species list
nrow(taxa[taxa$approximate_match==TRUE,]) # 0 species
taxa[taxa$approximate_match==TRUE,]  ### All good!!

# exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
taxa[taxa$number_matches != 1,]
ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa, ott_id = ott_id_tocheck[i]))
}

# Seems ok, all of them are synonyms

# check synonyms and change name accordingly
taxa[taxa$is_synonym==TRUE,] # 8 species returned

species[species=="Aspidoscelis costata"] <- "Aspidoscelis costatus" #they were an unique species in the past
species[species=="Aspidoscelis lineattissima"] <- "Aspidoscelis lineattissimus" #they were an unique species in the past
species[species=="Aspidoscelis sexlineata"] <- "Aspidoscelis sexlineatus" #they were an unique species in the past
species[species=="Eumeces okadae"] <- "Plestiodon latiscutatus" #synonym
species[species=="Japalura polygonata"] <- "Diploderma polygonatum (species in Holozoa)" #synonym
# species[species=="Japalura polygonata"] <- "Diploderma polygonatum" #synonym
species[species=="Japalura swinhonis"] <- "Diploderma swinhonis" #synonym
species[species=="Mabuya macrorhyncha"] <- "Psychosaura macrorhyncha" #synonym
species[species=="Podarcis raffoneae"] <- "Podarcis raffonei" #synonym
species[species=="Podarcis taurica"] <- "Podarcis tauricus" #synonym

reptdata[reptdata$Binomial=="Aspidoscelis costata","Binomial"] <- "Aspidoscelis costatus"
reptdata[reptdata$Binomial=="Aspidoscelis lineattissima","Binomial"] <- "Aspidoscelis lineattissimus"
reptdata[reptdata$Binomial=="Aspidoscelis sexlineata","Binomial"] <- "Aspidoscelis sexlineatus"
reptdata[reptdata$Binomial=="Eumeces okadae","Binomial"] <- "Plestiodon latiscutatus"
reptdata[reptdata$Binomial=="Japalura polygonata","Binomial"] <-  "Diploderma polygonatum (species in Holozoa)"
# reptdata[reptdata$Binomial=="Japalura polygonata","Binomial"] <-  "Diploderma polygonatum"
reptdata[reptdata$Binomial=="Japalura swinhonis","Binomial"] <- "Diploderma swinhonis"
reptdata[reptdata$Binomial=="Mabuya macrorhyncha","Binomial"] <-  "Psychosaura macrorhyncha"
reptdata[reptdata$Binomial=="Podarcis raffoneae","Binomial"] <-  "Podarcis raffonei"
reptdata[reptdata$Binomial=="Podarcis taurica","Binomial"] <-  "Podarcis tauricus"

# rerun 
taxa.c<- tnrs_match_names(names = species)

taxa.c[taxa.c$approximate_match==TRUE,] # no species returned!!
taxa.c[taxa.c$is_synonym==TRUE,] # no species returned!!

##############################################################
# Retrieving phylogenetic relationships
##############################################################

# retrieving phylogenetic relationships among taxa in the form 
# of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa.c[["ott_id"]], label_format = "name")
plot(tree, cex=.5, label.offset =.1, no.margin = TRUE)

# Notice that the species names shown in the tree are not exactly 
# the same as the species names that we had in our list. This is 
# because those names had synonyms in the tree of life database, 
# and we are using those names for the plot.


##############################################################
# Dealing with polytomies
##############################################################

# we can check for the existence of polytomies by running the 
# following code. If polytomies exist, the output will be 
# `FALSE`, and vice versa.

is.binary.tree(tree) # there are some polytomies

# to take care of these polytomies, we are going to use a 
# randomization approach
set.seed(111) #making it replicable, at least for this version of R (i.e. v.3.5.1), I AM USING THE SAME! 
tree_random <- multi2di(tree,random=TRUE)
is.binary.tree(tree_random)


##############################################################
# Final checks
##############################################################

# exploring whether our tree covers all the species we wanted 
# it to include, and making sure that the species names in our 
# database match those in the tree. We use the following code.

tree_random$tip.label <- gsub("_"," ", tree_random$tip.label)
intersect(as.character(tree_random$tip.label), as.character(species))

species[!species %in% as.character(tree_random$tip.label)] #listed in our database but not in the tree, none!
tree_random$tip.label[!as.character(tree_random$tip.label) %in% species] # listed in the tree but not in our database, none!

# they are the same species: Diploderma polygonatum and Diploderma polygonatum (species in Holozoa) 
# we fix it here
# tree_random$tip.label[tree_random$tip.label =="Diploderma polygonatum (species in Holozoa)"] <-"Diploderma polygonatum"

tiff("Results/Reptiles/Phylogeny/rept_phylogenetic_tree_pruned.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_random, cex=.5, label.offset =.1, no.margin = TRUE)

dev.off()

# we can now save the tree
# save(tree_random, file = "Data/Final data/rept_tree_random.Rdata")


##############################################################
# Computing branch lengths
##############################################################

# we are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# before we need to make sure that tree labels and database
# use the same nomenclature
setdiff(reptdata$Binomial, as.character(tree_random$tip.label))
setdiff(as.character(tree_random$tip.label),reptdata$Binomial)
# all good!

# exclude species in the tree that are not in dataset, none in this case
drops <- tree_random$tip.label[!tree_random$tip.label %in% reptdata$Binomial]
rept.tree_random.fixed <- drop.tip(tree_random, drops)

# save the new tree
write.tree(rept.tree_random.fixed, file = "Data/Final data/rept.tree_random.fixed.tre")

# compute branch lengths of tree
phylo_branch <- compute.brlen(rept.tree_random.fixed, method = "Grafen", power = 1)

# check tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

##############################################################
# Phylogenetic matrix
##############################################################

# matrix to be included in the models
rept_phylo_cor <- vcv(phylo_branch, cor = T)

# remove rows not in correlation matrix
reptdata_ph <- reptdata[which(reptdata$Binomial %in% rownames(rept_phylo_cor)),] # we haven't lost any species, yei!

##create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID<-data.frame(Binomial = unique(reptdata_ph$Binomial), SPID = paste0("SP",1:length(unique(reptdata_ph$Binomial))))
SpID$Binomial<-as.character(SpID$Binomial)
reptdata_ph<-inner_join(reptdata_ph,SpID, by = "Binomial")

# finally, save matrix for future analyses
save(rept_phylo_cor, file = "Data/Final data/rept_phylo_cor.Rdata")


# exporting fixed dataset for analyses (NOT NEEDED NOW)
write.csv(reptdata_ph, 
           "Data/Final data/reptdata_ph.csv", row.names = FALSE)

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/rept_phylogeny_R_session.txt")
sessionInfo()
sink()

# End of script ####