##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estacion Biologica de Donana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 18th of September 2019
# Modified July 2020

##############################################################
# Description of script and instructions
##############################################################

# This script is to build the phylogeny for amphibians and estimate the 
# phylogenetic relatedness among the species included in:

# Benitez-Lopez et al. The island rule explains consistent patterns of 
# body size evolution across terrestrial vertebrates. 

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
amphdata<-read.csv("Data/amphdata_def.csv", header = TRUE, stringsAsFactors = FALSE)

# generating list of species
species <- sort(unique(as.character(amphdata$Binomial)))


##############################################################
# Formatting species data
##############################################################

# obtaining dataframe listing the Open Tree identifiers potentially 
# matching our list of species.

taxa <- tnrs_match_names(names = species)

# according to the `approximate_match` column, there might be 
# 0 typo in the species list
nrow(taxa[taxa$approximate_match==TRUE,])
taxa[taxa$approximate_match==TRUE,]  #### zero species to fix

# exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
taxa[taxa$number_matches != 1,]  #NONE!
#ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"] 

#for(i in 1:length(ott_id_tocheck)){
#  print(inspect(taxa, ott_id = ott_id_tocheck[i]))
#}

# No species return more than one match

# check synonyms and change name accordingly
taxa[taxa$is_synonym==TRUE,] # 6 species returned

species[species=="Bufo boulengeri"] <- "Bufotes boulengeri" #synonym
species[species=="Bufo viridis"] <- "Bufotes viridis" #synonym
species[species=="Hypogeophis brevis"] <- "Grandisonia brevis" #synonym
species[species=="Hypsiboas albomarginatus"] <- "Boana albomarginata" #synonym
species[species=="Leptodactylus marmoratus"] <- "Adenomera marmorata" #synonym
species[species=="Ololygon trapicheiroi"] <- "Scinax trapicheiroi" #synonym
species[species=="Rana chalconota"] <- "Hylarana chalconota" #synonym

amphdata[amphdata$Binomial=="Bufo boulengeri","Binomial"] <- "Bufotes boulengeri"
amphdata[amphdata$Binomial=="Bufo viridis","Binomial"] <- "Bufotes viridis"
amphdata[amphdata$Binomial=="Hypogeophis brevis","Binomial"] <-  "Grandisonia brevis"
amphdata[amphdata$Binomial=="Hypsiboas albomarginatus","Binomial"] <-  "Boana albomarginata"
amphdata[amphdata$Binomial=="Leptodactylus marmoratus","Binomial"] <- "Adenomera marmorata"
amphdata[amphdata$Binomial=="Ololygon trapicheiroi","Binomial"] <-  "Scinax trapicheiroi"
amphdata[amphdata$Binomial=="Rana chalconota","Binomial"] <- "Hylarana chalconota"

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

species[!species %in% as.character(tree_random$tip.label)] #listed in our database but not in the tree, 0 species!
tree_random$tip.label[!as.character(tree_random$tip.label) %in% species] # listed in the tree but not in our database. 0 species!

tiff("Results/Amphibians/Phylogeny/amph_phylogenetic_tree_pruned.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_random, cex=.5, label.offset =.1, no.margin = TRUE)

dev.off()

# we can now save the tree
save(tree_random, file = "Data/Final data/amph_tree_random.Rdata")


##############################################################
# Computing branch lengths
##############################################################

# we are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# before we need to make sure that tree labels and database
# use the same nomenclature
setdiff(amphdata$Binomial, as.character(tree_random$tip.label))
setdiff(as.character(tree_random$tip.label),amphdata$Binomial)
# all good!

# exclude species in the tree that are not in dataset, none in this case
drops <- tree_random$tip.label[!tree_random$tip.label %in% amphdata$Binomial]
amph.tree_random.fixed <- drop.tip(tree_random, drops)

# save the new tree
write.tree(amph.tree_random.fixed, file = "Data/Final data/amph.tree_random.fixed.tre")

# compute branch lengths of tree
phylo_branch <- compute.brlen(amph.tree_random.fixed, method = "Grafen", power = 1)

# check tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

##############################################################
# Phylogenetic matrix
##############################################################

# matrix to be included in the models
amph_phylo_cor <- vcv(phylo_branch, cor = T)

# remove rows not in correlation matrix
amphdata_ph <- amphdata[which(amphdata$Binomial %in% rownames(amph_phylo_cor)),] # we haven't lost any species, yei!

##create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID<-data.frame(Binomial = unique(amphdata_ph$Binomial), SPID = paste0("SP",1:length(unique(amphdata_ph$Binomial))))
SpID$Binomial<-as.character(SpID$Binomial)
amphdata_ph<-inner_join(amphdata_ph,SpID, by = "Binomial")


# finally, save matrix for future analyses
save(amph_phylo_cor, file = "Data/Final data/amph_phylo_cor.Rdata")


# exporting fixed dataset for analyses
write.csv(amphdata_ph, 
          "Data/Final data/amphdata_ph.csv", row.names = FALSE)

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/amph_phylogeny_R_session.txt")
sessionInfo()
sink()

# End of script ####
