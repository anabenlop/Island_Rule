##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estacion Biologica de Donana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 18th of September 2019, adapted from Sánchez-Tojar et al. 2019 (https://osf.io/yjua8/)

##############################################################
# Description of script and instructions
##############################################################

# This script is used to build the phylogeny for mammals and estimate the 
# phylogenetic relatedness among the species included in:

# Benitez-Lopez et al. The island rule explains consistent patterns of 
# body size evolution


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
mamdata<-read.csv("Data/mamdata_def.csv", header = TRUE, stringsAsFactors = FALSE) #1051

# generating list of species
species <- sort(unique(as.character(mamdata$Binomial))) #218 species
# species_is <- sort(unique(as.character(mamdata$Species_island))) #328 species

##############################################################
# Formatting species data
##############################################################

# obtaining dataframe listing the Open Tree identifiers potentially 
# matching our list of species.

taxa <- tnrs_match_names(species) 

# according to the `approximate_match` column, there might be 
# 0 typos in the species list 
# nrow(taxa[taxa$approximate_match==TRUE,])
taxa[taxa$approximate_match==TRUE,] # no species returned

# fixing those typos (example in case they were unresolved matches)
#species[species=="Crocidura attenuatta"] <- "Crocidura attenuata"

#mamdata$Binomial<-as.character(mamdata$Binomial)
#mamdata[mamdata$Binomial=="Crocidura attenuatta","Binomial"] <- "Crocidura attenuata"

# rerun
# taxa.c <- tnrs_match_names(names = species)
# 
# taxa.c[taxa.c$approximate_match==TRUE,] # no species returned

# exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
taxa.c <- taxa
taxa.c[taxa.c$number_matches != 1,]
ott_id_tocheck <- taxa.c[taxa.c$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa.c, ott_id = ott_id_tocheck[i]))
}

# Seems ok, they are synonyms

# check synonyms and change name accordingly
taxa.c[taxa.c$is_synonym==TRUE,] # 9 species returned

species[species=="Aonyx cinerea"] <- "Aonyx cinereus"
species[species=="Galago zanzibaricus"] <- "Galagoides zanzibaricus"
species[species=="Oncifelis guigna"] <- "Leopardus guigna"
species[species=="Herpestes javanicus"] <- "Urva javanica"
species[species=="Lasiurus cinereus"] <- "Aeorestes cinereus"
species[species=="Lissonycteris angolensis"] <- "Myonycteris angolensis"
species[species=="Spermophilus beecheyi"] <- "Otospermophilus beecheyi"
species[species=="Spermophilus parryii"] <- "Urocitellus parryii"
species[species=="Tarsius bancanus"] <- "Cephalopachus bancanus"
species[species=="Tarsius syrichta"] <- "Carlito syrichta"

mamdata[mamdata$Binomial=="Aonyx cinerea", "Binomial"] <- "Aonyx cinereus"
mamdata[mamdata$Binomial=="Galago zanzibaricus","Binomial"] <- "Galagoides zanzibaricus"
mamdata[mamdata$Binomial=="Oncifelis guigna","Binomial"] <- "Leopardus guigna"
mamdata[mamdata$Binomial=="Herpestes javanicus", "Binomial"] <- "Urva javanica"
mamdata[mamdata$Binomial=="Lasiurus cinereus","Binomial"] <- "Aeorestes cinereus"
mamdata[mamdata$Binomial=="Lissonycteris angolensis","Binomial"] <- "Myonycteris angolensis"
mamdata[mamdata$Binomial=="Spermophilus beecheyi","Binomial"] <- "Otospermophilus beecheyi"
mamdata[mamdata$Binomial=="Spermophilus parryii","Binomial"] <- "Urocitellus parryii"
mamdata[mamdata$Binomial=="Tarsius bancanus","Binomial"] <- "Cephalopachus bancanus"
mamdata[mamdata$Binomial=="Tarsius syrichta","Binomial"] <- "Carlito syrichta"

# rerun 2
taxa.c2 <- tnrs_match_names(names = species)

taxa.c2[taxa.c2$approximate_match==TRUE,] # no species returned
taxa.c2[taxa.c2$is_synonym==TRUE,] # no species returned


##############################################################
# Retrieving phylogenetic relationships
##############################################################

# retrieving phylogenetic relationships among taxa in the form 
# of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa.c2[["ott_id"]], label_format = "name")
plot(tree, cex=.5, label.offset =.1, no.margin = TRUE)

##############################################################
# Dealing with polytomies
##############################################################

# we run the following function to check for the existence of polytomies.
# If polytomies exist, the output will be `FALSE`, and vice versa.

is.binary.tree(tree) # there are some polytomies


# to take care of these polytomies, we are going to use a 
# randomization approach
set.seed(111) #making it replicable, at least for this version of R (i.e. v.3.5.3)
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

species[!species %in% as.character(tree_random$tip.label)] #listed in our database but not in the tree
tree_random$tip.label[!as.character(tree_random$tip.label) %in% species] # listed in the tree but not in our database

##I have a problem with a species labeled "mrcaott319315ott366155" in the tree
## Also, Eulemur fulvus and "Miniopterus schreibersii" are in my data but not in the tree

# try to see which species is that

test<-tnrs_match_names(names = c("Eulemur fulvus","Miniopterus schreibersii"))
tree_test <- tol_induced_subtree(ott_ids = c("394961",  "238425"), label_format = "name")

# Miniopterus schreibersii is mrcaott319315ott366155 and Eulumur fulvus is Eulemur, which we cannot fix
# we fix Miniopterus here
tree_random$tip.label[tree_random$tip.label =="mrcaott319315ott366155"] <-"Miniopterus schreibersii"

tiff("Results/Mammals/Phylogeny/mam_phylogenetic_tree_pruned.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_random, cex=.5, label.offset =.1, no.margin = TRUE)

dev.off()

# we can now save the tree
save(tree_random, file = "Data/Final data/mam_tree_random.Rdata")


##############################################################
# Computing branch lengths
##############################################################

# we are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# before we need to make sure that tree labels and database
# use the same nomenclature
setdiff(mamdata$Binomial, as.character(tree_random$tip.label)) # Eulemur fulvus still missing
setdiff(as.character(tree_random$tip.label),mamdata$Binomial)

# exclude species in the tree that are not in dataset
  drops <- tree_random$tip.label[!tree_random$tip.label %in% mamdata$Binomial]
  mam.tree_random.fixed <- drop.tip(tree_random, drops)

# save the new tree
write.tree(mam.tree_random.fixed, file = "Data/Final data/mam.tree_random.fixed.tre")

# compute branch lengths of tree
phylo_branch <- compute.brlen(mam.tree_random.fixed, method = "Grafen", power = 1)

# check if tree is ultrametric
is.ultrametric(phylo_branch) # TRUE


##############################################################
# Phylogenetic matrix
##############################################################

# matrix to be included in the models
mam_phylo_cor <- vcv(phylo_branch, cor = T)

# remove rows not in correlation matrix
mamdata_ph <- mamdata[which(mamdata$Binomial %in% rownames(mam_phylo_cor)),] # we lose 1 datapoint (Eulemur fulvus), 1046 rows

##create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID<-data.frame(Binomial = unique(mamdata_ph$Binomial), SPID = paste0("SP",1:length(unique(mamdata_ph$Binomial))))
SpID$Binomial<-as.character(SpID$Binomial)
mamdata_ph<-inner_join(mamdata_ph,SpID, by = "Binomial")

# finally, save matrix for future analyses
save(mam_phylo_cor, file = "Data/Final data/mam_phylo_cor.Rdata")


# exporting fixed dataset for analyses
write.csv(mamdata_ph, 
           "Data/Final data/mamdata_ph.csv", row.names = FALSE)

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/mam_phylogeny_R_session.txt")
sessionInfo()
sink()

# End of script ####