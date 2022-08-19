###########################################################################################
# Script: Source trait data
# AFE
# August 2022
###########################################################################################

# Libreries:-------------------------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(readxl)
library(hablar)
library(rfishbase)
library(data.table)

rm(list=ls())
myd <- getwd()

# Read data:-------------------------------------------------------------------------------
FMorph <- read.csv2("C:/Users/Usuario/Documents/PHD/Datasets/FishMORPH/FISHMORPH_Database.csv", h=T)
load(paste0(myd, "/HNC.RData"))
load(paste0(myd, "/HNB.RData"))
load(paste0(myd, "/HE.RData"))
load(paste0(myd, "/ContN.RData"))
load(paste0(myd, "/ContE.RData"))

###########################################################################################
# 1) Retrieve measurements avilable in FishMorph: -----------------------------------------

# Species for which data are available:
vecsps <- as.vector(sort(unique(c(names(HNB[,-c(1:2)]),names(ContN[,-c(1:2)]), names(ContE[,-c(1:2)]))))) # 98
vecsps <- plyr::revalue(vecsps, c("Tampichthys dichroma"="Tampichthys dichromus"))
genus_sp <- c("Astyanax sp","Gila sp", "Pseudoxiphophorus sp", "Chirostoma sp", "Oreochromis sp", "Poecilia sp","Poeciliopsis sp")
vecsps <- setdiff(vecsps, genus_sp) # 91 (- 7) OK

FMorphSubs <- subset(FMorph, FMorph$Genus.species %in% vecsps)           # 69

# Try other updated & non updated names instead:
NP <- as.vector(setdiff(vecsps, FMorphSubs$Genus.species))               # 22
NP

# Retrieve the list of data for which some measurements are unavailable: 
sapply(FMorphSubs[,c(6:15)], function(x) sum(is.na(x)))
# MBl BEl VEp REs OGp RMl BLs PFv PFs CPt 
# 9   0   0   0   0   2   0   1   9   0 

FMorphFill <- data.frame(matrix(ncol=17, nrow=22))
names(FMorphFill)<- names(FMorphSubs)
FMorphFill[,5] <- NP

F_GenusSp <- data.frame(matrix(ncol=17, nrow=7))
names(F_GenusSp)<- names(FMorphSubs)
F_GenusSp[,5] <- genus_sp

69 + 22 + 7 # 98, OK. Sourced (with some needing updates) + (to fill) + "genus sp"

# Save fishtrait data:
FMorphSubs$Genus.species <- plyr::revalue(FMorphSubs$Genus.species, c("Tampichthys dichromus"="Tampichthys dichroma"))
fishtrait <- rbind(FMorphSubs, FMorphFill, F_GenusSp)
rownames(fishtrait) <- 1:nrow(fishtrait)
save(fishtrait, file=paste0(myd,"/fishtrait.RData"))              

# NOTES:
# No info for 29 records. 
# 7 are "Genus sp". I'll take the mean across Genus members present in the survey
# Measurement of pictures for the other 22 records

###########################################################################################
# 2) Fill out missing trait data & update values: -----------------------------------------


# 2A) MBl:
fishtrait$Genus.species <- plyr::revalue(fishtrait$Genus.species, c("Agonostomus monticola"="Dajaus monticola"))

FB_SpsT <- distinct(species(fishtrait$Genus.species))        
FB_EstimateT <- distinct(estimate(fishtrait$Genus.species))  

# To visualize the differences between measures and the vals reported in FMorph:
BS <- data.frame(matrix(ncol=8, nrow=98))                               
names(BS) <- c("Species", "MaxLength", "Type","EstMaxLengthTL", "EstMaxLengthSL", "FB_Val", "OE", "FMorph")

BS[,1] <- fishtrait$Genus.species
BS[,2] <- FB_SpsT$Length[match(FB_SpsT$Species, BS$Species)]
BS[,3] <- FB_SpsT$LTypeMaxM[match(FB_SpsT$Species, BS$Species)]
BS[,4] <- round(FB_EstimateT$MaxLengthTL[match(FB_EstimateT$Species, BS$Species)], 2)
BS[,5] <- round(FB_EstimateT$MaxLengthSL[match(FB_EstimateT$Species, BS$Species)], 2)
BS[,6] <- ifelse(BS$Type=="SL", BS$MaxLength, BS$EstMaxLengthSL)
BS[,7] <- ifelse(BS$Type=="SL", "Observed_FB", "Estimated_FB")
BS[,8] <- fishtrait$MBl[match(fishtrait$Genus.species, BS$Species)]

# Save:
write.csv2(BS, file=paste0(myd, "/BS.csv"))
save(BS, file=paste0(myd, "/BS.RData"))

# NOTES: 
# Values are a mix of SL & TL in FishBase:
# Keep those that are SL, and assign a observed/estimated classification.
# Go through all vals reported in Miller, 2005 and update those that are NAs in the
# database or that are larger than those in FMorph.
# Standard length is synonym for "longitud patron, LP" in Miller, 2005


# 2B) Update BS measurements:
BSU <- read.csv2(paste0(myd, "/BSUpdated.csv"), h=T, sep=";")

setdiff(sort(unique(BSU$Species)), sort(unique(fishtrait$Genus.species))) # "Chirostoma aculeatum" "Chirostoma charari"
BSU <- BSU[!BSU$Species %in% c("Chirostoma aculeatum","Chirostoma charari"),]
identical(sort(unique(BSU$Species)), sort(unique(fishtrait$Genus.species))) #TRUE

# Update MBl values in fishtrait:
fishtrait$MBl <- BSU$Final[match(as.character(fishtrait$Genus.species), as.character(BSU$Species))]


###########################################################################################
# 3) Fill out missing species data:

# TBC
###########################################################################################
# End of script ###########################################################################

