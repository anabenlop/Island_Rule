##### Island_Rule ######
 Data and scrips to run the analyses for the paper: 
 
 Benítez-López et al. The island rule explains consistent patterns of body size evolution in terrestrial vertebrates.
 
 All data are available at: https://figshare.com/account/home#/projects/89102 and in this GitHub repository.
 
The Supplementary Dataset 1 in figshare contains the raw data in xls format. The dataset corresponds to the islandrule_clean.csv file in the GitHub repository: https://github.com/anabenlop/Island_Rule
Following the sequence of scripts in the repository from 001_data_prep.R to 005_metaregressions.R, these data is merged with spatial and diet data, allometric models are applied to convert
original size measurements to body mass equivalents, phylogenetic correlation matrixes are calculated and, finally, the models are run. Scripts from 006 until 015 examine heterogeneity and run 
sensitivity analyses to assess the influence of data imputation, publication bias, the use of archipelagos or small sample bias. We also run alternative models to assess the robustness of the island rule (008_Alternativeodels.R).
All scripts are commented.
 
The final datasets used in the analyses are stored in the folder Final data in the github repository. These files are: mamdata_ph.csv, birddata_ph.csv, reptdata_ph.csv, amphdata_ph.csv.
birddata_ph_tarsus.csv corresponds to the dataset using tarsus length as an alternative metric for double checking the consistency in our results.

Some intermediate files are used in the workflow:

# Allometric models
The table with all info to convert linear dimensions to body mass equivalents based on allometric models is included in the file allometric_relationships.csv. This file is loaded and used in the scripts
002_data_mammals_prep.R, 002_data_birds_prep.R, 002_data_rept_prep.R and 002_data_amph_prep.R 

# Migratory behaviour
SpeciesList3_1_migbehav_v2_0.csv corresponds to the dataset on bird migratrory behaviour by Eyres et al. (2017), used to classify bird species as migratory or not. Fully migratory species were removed from the
analysis using the script 002_data_birds_prep.R

# Diet
B_traits_guild.csv and M_traits_guild.csv are based on the EltonTraits database (Wilman, H. et al. EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals: Ecological Archives E095-178. Ecology 95, 2027-2027 (2014)) and used to classify species as carnivores or not. 
For reptiles we used diet information from Scharf, I. et al. Late bloomers and baby boomers: ecological drivers of longevity in squamates and the tuatara. Global Ecol. Biogeogr. 24, 396-405 (2015). and
Meiri, S. Traits of lizards of the world: Variation around a successful evolutionary design. Global Ecol. Biogeogr. 27, 1168-1172 (2018).
