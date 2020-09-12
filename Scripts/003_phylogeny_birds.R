##############################################################
# Authors: 
# Ana Benitez-Lopez (@anabenlop)
# Scholar Profile: https://scholar.google.com/citations?user=HC_j51sAAAAJ&hl=es
# Department of Integrative Ecology, Estacion Biologica de Donana (EBD-CSIC, ESP) 
# Email: abenitez81@gmail.com

# Script first created on the 18th of September 2019

##############################################################
# Description of script and instructions
##############################################################

# This script is to build the phylogeny for birds and estimate the 
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
birddata<-read.csv("Data/birddata_def.csv", header = TRUE, stringsAsFactors = FALSE) #use birddatadef2 for tarsus
birddata <- birddata[birddata$Binomial != "Emberiza calandra", ] #Remove Emberiza calandra, not available in OpenTreeofLife

# generating list of species
species <- sort(unique(as.character(birddata$Binomial)))

##############################################################
# Formatting species data
##############################################################

# obtaining dataframe listing the Open Tree identifiers potentially 
# matching our list of species.

taxa <- tnrs_match_names(names = species)

# it does not match  Aethopyga latouchii, Emberiza calandra, Hemicircus sordidus, Trichastoma tickelli and  Meiglyptes grammithorax, 
# we try by changing the species name to previous name:  Phaiopicus grammithorax

birddata[birddata$Binomial=="Aethopyga latouchii","Binomial"] <- "Aethopyga christinae" #previously clumped
birddata[birddata$Binomial=="Meiglyptes grammithorax","Binomial"] <- "Meiglyptes tristis" #previously clumped
birddata[birddata$Binomial=="Trichastoma tickelli","Binomial"] <- "Pellorneum tickelli" #synonym
birddata[birddata$Binomial=="Hemicircus sordidus","Binomial"] <- "Hemicircus concretus" #conspecific
birddata[birddata$Binomial=="Zanda funerea","Binomial"] <- "Calyptorhynchus funereus" #conspecific


# try again
species.r <- sort(unique(as.character(birddata$Binomial)))
taxa.r <- tnrs_match_names(names = species.r)

# according to the `approximate_match` column, there might be 
# 0 typo in the species list
nrow(taxa.r[taxa.r$approximate_match==TRUE,])
taxa.r[taxa.r$approximate_match==TRUE,]

# exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
taxa.r[taxa.r$number_matches != 1,]

# taxa.r<-taxa.r[!is.na(taxa.r$unique_name),] #remove NA

ott_id_tocheck <- taxa.r[taxa.r$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa.r, ott_id = ott_id_tocheck[i]))
}

# Seems ok, they are synonyms

# check synonyms and change name accordingly
taxa.r[taxa.r$is_synonym==TRUE,] # 23 species returned

species.r[species.r=="Actinodura souliei"] <- "Sibia souliei"
species.r[species.r=="Alcedo azurea"] <- "Ceyx azureus"
species.r[species.r=="Ceyx pusillus"] <- "Alcedo pusilla"
species.r[species.r=="Copsychus malabaricus"] <- "Kittacincla malabarica"
species.r[species.r=="Cormobates leucophaea"] <- "Cormobates leucophaeus"
species.r[species.r=="Dryobates scalaris"] <- "Picoides scalaris"
species.r[species.r=="Garrulax elliotii"] <- "Trochalopteron elliotii"
species.r[species.r=="Liocichla ripponi"] <- "Liocichla phoenicea ripponi"
species.r[species.r=="Macronous kelleyi"] <- "Macronus kelleyi"
species.r[species.r=="Napothera brevicaudata"] <- "Turdinus brevicaudatus"
species.r[species.r=="Nigrita canicapillus"] <- "Nigrita canicapilla"
species.r[species.r=="Nisaetus cirrhatus"] <- "Spizaetus cirrhatus"
species.r[species.r=="Peneoenanthe pulverulenta"] <- "Peneothello pulverulenta"
species.r[species.r=="Phaenicophaeus tristis"] <- "Rhopodytes tristis"
species.r[species.r=="Pheugopedius felix"] <- "Thryothorus felix"
species.r[species.r=="Phylidonyris pyrrhopterus"] <- "Phylidonyris pyrrhoptera"
species.r[species.r=="Sephanoides sephaniodes"] <- "Sephanoides sephanoides"
species.r[species.r=="Pomatorhinus swinhoei"] <- "Erythrogenys swinhoei"
species.r[species.r=="Saxicola torquatus"] <- "Saxicola torquata"
species.r[species.r=="Sericornis citreogularis"] <- "Neosericornis citreogularis"
species.r[species.r=="Sylvia melanocephala"] <- "Curruca melanocephala"
# species.r[species.r=="Tribonyx ventralis"] <- "Gallinula ventralis"
species.r[species.r=="Turdus libonyanus"] <- "Turdus libonyana"
species.r[species.r=="Zoothera interpres"] <- "Geokichla interpres"

birddata[birddata$Binomial=="Actinodura souliei","Binomial"] <- "Sibia souliei"
birddata[birddata$Binomial=="Alcedo azurea","Binomial"] <- "Ceyx azureus"
birddata[birddata$Binomial=="Ceyx pusillus", "Binomial"] <- "Alcedo pusilla"
birddata[birddata$Binomial=="Copsychus malabaricus","Binomial"] <- "Kittacincla malabarica"
birddata[birddata$Binomial=="Cormobates leucophaea","Binomial"] <-  "Cormobates leucophaeus"
birddata[birddata$Binomial=="Dryobates scalaris","Binomial"] <-  "Picoides scalaris"
birddata[birddata$Binomial=="Garrulax elliotii","Binomial"] <-  "Trochalopteron elliotii"
birddata[birddata$Binomial=="Liocichla ripponi","Binomial"] <- "Liocichla phoenicea ripponi"
birddata[birddata$Binomial=="Macronous kelleyi","Binomial"] <- "Macronus kelleyi"
birddata[birddata$Binomial=="Napothera brevicaudata","Binomial"] <- "Turdinus brevicaudatus"
birddata[birddata$Binomial=="Nigrita canicapillus","Binomial"] <- "Nigrita canicapilla"
birddata[birddata$Binomial=="Nisaetus cirrhatus", "Binomial"] <- "Spizaetus cirrhatus"
birddata[birddata$Binomial=="Peneoenanthe pulverulenta","Binomial"] <- "Peneothello pulverulenta"
birddata[birddata$Binomial=="Phaenicophaeus tristis","Binomial"] <- "Rhopodytes tristis"
birddata[birddata$Binomial=="Pheugopedius felix", "Binomial"] <- "Thryothorus felix"
birddata[birddata$Binomial=="Phylidonyris pyrrhopterus", "Binomial"] <- "Phylidonyris pyrrhoptera"
birddata[birddata$Binomial=="Sephanoides sephaniodes", "Binomial"] <- "Sephanoides sephanoides"
birddata[birddata$Binomial=="Pomatorhinus swinhoei", "Binomial"] <- "Erythrogenys swinhoei"
birddata[birddata$Binomial=="Saxicola torquatus", "Binomial"] <- "Saxicola torquata"
birddata[birddata$Binomial=="Sericornis citreogularis", "Binomial"] <- "Neosericornis citreogularis"
birddata[birddata$Binomial=="Sylvia melanocephala","Binomial"] <- "Curruca melanocephala"
# birddata[birddata$Binomial=="Tribonyx ventralis","Binomial"] <- "Gallinula ventralis"
birddata[birddata$Binomial=="Turdus libonyanus","Binomial"] <- "Turdus libonyana"
birddata[birddata$Binomial=="Zoothera interpres", "Binomial"] <- "Geokichla interpres"

# rerun 2 for syn
taxa.r <- tnrs_match_names(names = species.r)

taxa.r[taxa.r$approximate_match==TRUE,] # no species returned!! 
taxa.r[taxa.r$is_synonym==TRUE,] # no species returned

##############################################################
# Retrieving phylogenetic relationships
##############################################################

# retrieving phylogenetic relationships among taxa in the form 
# of a trimmed sub-tree
# taxa.r<-taxa.r[!is.na(taxa.r$is_synonym),] #remove NA
tree <- tol_induced_subtree(ott_ids = taxa.r[["ott_id"]], label_format = "name")
plot(tree, cex=.3, label.offset =.1, no.margin = TRUE)

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
intersect(as.character(tree_random$tip.label), as.character(species.r))

species.r[!species.r %in% as.character(tree_random$tip.label)] #listed in our database but not in the tree:  "Amazona ochrocephala"  "Miliaria calandra" 
tree_random$tip.label[!as.character(tree_random$tip.label) %in% species.r] #  "mrcaott165615ott165621"   
 
# try to see which species is that

#test<-tnrs_match_names(names = "Amazona ochrocephala")

# they are the same species: "Amazona ochrocephala" and "mrcaott165615ott165621" 
# we fix it here
tree_random$tip.label[tree_random$tip.label =="mrcaott165615ott165621"] <-"Amazona ochrocephala"

tiff("Results/Birds/Phylogeny/bird_phylogenetic_tree_pruned2.tiff",
     height=80, width=10,
     units='cm', compression="lzw", res=1200)

plot(tree_random, cex=.5, label.offset =.1, no.margin = TRUE)

dev.off()

# we can now save the tree
save(tree_random, file = "Data/Final data/bird_tree_random.Rdata")


##############################################################
# Computing branch lengths
##############################################################

# we are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# before we need to make sure that tree labels and database
# use the same nomenclature
setdiff(birddata$Binomial, as.character(tree_random$tip.label)) 
setdiff(as.character(tree_random$tip.label),birddata$Binomial)
# all good!

# exclude species in the tree that are not in dataset
drops <- tree_random$tip.label[!tree_random$tip.label %in% birddata$Binomial]
bird.tree_random.fixed <- drop.tip(tree_random, drops)

# save the new tree
write.tree(bird.tree_random.fixed, file = "Data/Final data/bird.tree_random.fixed.tre")

# compute branch lengths of tree
phylo_branch <- compute.brlen(bird.tree_random.fixed, method = "Grafen", power = 1)

# check tree is ultrametric
is.ultrametric(phylo_branch) # TRUE


##############################################################
# Phylogenetic matrix
##############################################################

# matrix to be included in the models
bird_phylo_cor <- vcv(phylo_branch, cor = T)

# remove rows not in correlation matrix
birddata_ph <- birddata[which(birddata$Binomial %in% rownames(bird_phylo_cor)),] # 706 we lose Emberiza calandra 

##create Species ID to distinguish later between variation explained by non-phylogenetic and phylogenetic effects
SpID<-data.frame(Binomial = unique(birddata_ph$Binomial), SPID = paste0("SP",1:length(unique(birddata_ph$Binomial))))
SpID$Binomial<-as.character(SpID$Binomial)
birddata_ph<-inner_join(birddata_ph,SpID, by = "Binomial") #706 rows

# finally, save matrix for future analyses
save(bird_phylo_cor, file = "Data/Final data/bird_phylo_cor.Rdata")

# exporting fixed dataset for analyses
write.csv(birddata_ph, 
           "Data/Final data/birddata_ph.csv", row.names = FALSE)

# saving session information with all packages versions for reproducibility purposes
sink("Data/Final data/bird_phylogeny_R_session.txt")
sessionInfo()
sink()

# End of script ####
