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
vecsps <- as.vector(sort(unique(c(names(HNB[,-c(1:2)]),names(ContN[,-c(1:2)]), names(ContE[,-c(1:2)]))))) # 100
vecsps <- plyr::revalue(vecsps, c("Tampichthys dichroma"="Tampichthys dichromus"))
genus_sp <- c("Astyanax sp","Gila sp", "Pseudoxiphophorus sp", "Chirostoma sp", "Oreochromis sp", "Poecilia sp","Poeciliopsis sp")
vecsps <- setdiff(vecsps, genus_sp) # 93 (- 7) OK

FMorphSubs <- subset(FMorph, FMorph$Genus.species %in% vecsps)           # 70

# Try other updated & non updated names instead:
NP <- as.vector(setdiff(vecsps, FMorphSubs$Genus.species))               # 23
NP

# Retrieve the list of data for which some measurements are unavailable: 
sapply(FMorphSubs[,c(6:15)], function(x) sum(is.na(x)))
# MBl BEl VEp REs OGp RMl BLs PFv PFs CPt 
# 10   0   0   0   0   2   0   1   9   0 

FMorphFill <- data.frame(matrix(ncol=17, nrow=23))
names(FMorphFill)<- names(FMorphSubs)
FMorphFill[,5] <- NP

F_GenusSp <- data.frame(matrix(ncol=17, nrow=7))
names(F_GenusSp)<- names(FMorphSubs)
F_GenusSp[,5] <- genus_sp

70 + 23 + 7 # 98, OK. Sourced (with some needing updates) + (to fill) + "genus sp"

# Save fishtrait data:
FMorphSubs$Genus.species <- plyr::revalue(FMorphSubs$Genus.species, c("Tampichthys dichromus"="Tampichthys dichroma"))
fishtrait <- rbind(FMorphSubs, FMorphFill, F_GenusSp)
rownames(fishtrait) <- 1:nrow(fishtrait)
save(fishtrait, file=paste0(myd,"/fishtrait.RData"))              

# NOTES:
# No info for 30 records. 
# 7 are "Genus sp". I'll take the mean across Genus members present in the survey
# Measurement of pictures for the other 23 records (+ 3 of Gila sp in the area)

###########################################################################################
# 2) Fill out missing trait data & update values: -----------------------------------------

# 2A) MBl:
fishtrait$Genus.species <- plyr::revalue(fishtrait$Genus.species, c("Agonostomus monticola"="Dajaus monticola"))

FB_SpsT <- distinct(species(fishtrait$Genus.species))        
FB_EstimateT <- distinct(estimate(fishtrait$Genus.species))  

# To visualize the differences between measures and the vals reported in FMorph:
BS <- data.frame(matrix(ncol=8, nrow=100))                               
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

identical(sort(unique(BSU$Species)), sort(unique(fishtrait$Genus.species))) # TRUE

# Update MBl values in fishtrait:
fishtrait$MBl <- BSU$Final[match(as.character(fishtrait$Genus.species), as.character(BSU$Species))]


###########################################################################################
# 3) Fill out missing species data:
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/Traits/Measurements")
file_names <- dir() # 26
mes <- do.call(rbind,lapply(file_names,read.csv))   # 312 (12*26)                                # All data
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2")   # Go back to main dir

# Retrieve mesurements:
file_names <- str_split_fixed(file_names, "_", 2)[,2]
file_names <- str_split_fixed(file_names, "\\.", 2)[,1]
file_names <- gsub("_", " ", file_names)

mes$S <- rep(file_names, each=12)
mes$M <- rep(c("Bl", "Bd", "Hd", "Ed", "CPd", "CFd", "Eh", "Mo", "Jl", "PFh", "PFi", "PFiII"), times=26)

sum(mes$Length==0)  # 6
sum(mes$Length <0)  # 0
mes$Length[mes$Length==0] <- NA # Measurements that could not be taken
sum(is.na(mes$Length))

mes <- mes[c("S", "M", "Length")]
names(mes) <- c("Species", "Measure", "Length")

mes_mat <- spread(mes, key="Measure", value="Length")

# Compute ratios:
fishtrait_fill <- subset(fishtrait, fishtrait$Genus.species %in% file_names) #23, OK
names(fishtrait_fill)

gila_mes_mat <- subset(mes_mat, mes_mat$Species %in% c("Gila conspersa", "Gila minaceae", "Gila pulchra"))
mes_mat <- mes_mat[! mes_mat$Species %in% c("Gila conspersa", "Gila minaceae", "Gila pulchra"),]

fishtrait_fill$BEl <- mes_mat$Bl/mes_mat$Bd
fishtrait_fill$VEp <- mes_mat$Eh/mes_mat$Bd
fishtrait_fill$REs <- mes_mat$Ed/mes_mat$Hd
fishtrait_fill$OGp <- mes_mat$Mo/mes_mat$Bd
fishtrait_fill$RMl <- mes_mat$Jl/mes_mat$Hd
fishtrait_fill$BLs <- mes_mat$Hd/mes_mat$Bd
fishtrait_fill$PFv <- mes_mat$PFiII/mes_mat$Bd #Or PFi
fishtrait_fill$PFs <- mes_mat$PFh/mes_mat$Bl
fishtrait_fill$CPt <- mes_mat$CFd/mes_mat$CPd

save(fishtrait_fill, file="fishtrait_fill.RData")

fishtrait_fill$Reference <- rep("Miller, R.R. (2009) Freshwater Fishes of México (Peces Dulceacuícolas de México). 1st edn.", nrow(fishtrait_fill))
fishtrait_fill$Type.of.illustration <- rep("Picture", nrow(fishtrait_fill))
fishtrait_fill$Type.of.illustration <- ifelse(fishtrait_fill$Genus.species %in% c("Allotoca zacapuensis",
                                                                                  "Ictalurus mexicanus",
                                                                                  "Pseudoxiphophorus jonesii",
                                                                                  "Sicydium multipunctatum",
                                                                                  "Yuriria chapalae"), "Drawing",fishtrait_fill$Type.of.illustration)
fishtrait <- fishtrait[! fishtrait$Genus.species %in% fishtrait_fill$Genus.species,]
fishtrait <- retype(fishtrait)
fishtrait_fill <- retype(fishtrait_fill)
fishtrait <- bind_rows(fishtrait, fishtrait_fill)

# Genus sp records:
fishtrait$Genus <- rep(NA, nrow(fishtrait))
fishtrait$Genus <- str_split_fixed(fishtrait$Genus.species, " ", 2)[,1]

Ore <- subset(fishtrait, fishtrait$Genus=="Oreochromis")  # Average between O. mossambicus, O. aureus & O. niloticus
Ore_m <- apply(Ore[,c(6:15)], 2, mean, na.rm=TRUE) 
Ore_sp <- c(NA, NA, NA, NA,"Oreochromis sp", Ore_m, NA, NA, NA) 

Heter <- subset(fishtrait, fishtrait$Genus=="Pseudoxiphophorus") # Average between P. jonesii & P. bimaculatus
Heter_m <- apply(Heter[,c(6:15)], 2, mean, na.rm=TRUE) 
Heter_sp <- c(NA, NA, NA, NA,"Pseudoxiphophorus sp", Heter_m, NA, NA, NA)

Ast <- subset(fishtrait, fishtrait$Genus=="Astyanax") # Average between A. aeneus & A. mexicanus
Ast_m <- apply(Ast[,c(6:15)], 2, mean, na.rm=TRUE)
Ast_sp <- c(NA, NA, NA, NA,"Astyanax sp", Ast_m, NA, NA, NA)

Poe <- subset(fishtrait, fishtrait$Genus=="Poecilia") # Average between P. butleri, P. mexicana, P. reticulata, P. sphenops & P. chica
Poe_m <- apply(Poe[,c(6:15)], 2, mean, na.rm=TRUE) 
Poe_sp <- c(NA, NA, NA, NA,"Poecilia sp", Poe_m, NA, NA, NA)

Poel <- subset(fishtrait, fishtrait$Genus=="Poeciliopsis") 
# Average between P. balsas, P. gracilis, P. infans, P. turrubarensis, P. viriosa, P. baenschi, P. turneri
Poel_m <- apply(Poel[,c(6:15)], 2, mean, na.rm=TRUE) 
Poel_sp <- c(NA, NA, NA, NA,"Poeciliopsis sp", Poel_m, NA, NA, NA)

Chir <- subset(fishtrait, fishtrait$Genus=="Chirostoma")
# Average between C. aculeatum, C. arge, C. chapalae, C. riojai, C. charari, C. consocium, C. humboldtianum, C. jordani, C. labarcae
Chir_m <- apply(Chir[,c(6:15)], 2, mean, na.rm=TRUE) 
Chir_sp <- c(NA, NA, NA, NA,"Chirostoma sp", Chir_m, NA, NA, NA)

# Average results of all the possible native species in the genus Gila:
# According to Miller (2005) the ones that are more likely to be are: Gila conspersa, Gila minaceae, Gila modesta & Gila pulchra.
Gila_modesta <- subset(FMorph, FMorph$Genus.species %in% "Gila modesta")
Gila_modesta$MBl <- 10.2  # MBl update (Miller, 2005)

Gila <- as.data.frame(matrix(ncol=18, nrow=3))
names(Gila) <- names(fishtrait)

Gila$Genus.species <- gila_mes_mat$Species
Gila$MBl <- c(16.5, 60, 15.4) # Miller (2005)
Gila$BEl <- gila_mes_mat$Bl/gila_mes_mat$Bd
Gila$VEp <- gila_mes_mat$Eh/gila_mes_mat$Bd
Gila$REs <- gila_mes_mat$Ed/gila_mes_mat$Hd
Gila$OGp <- gila_mes_mat$Mo/gila_mes_mat$Bd
Gila$RMl <- gila_mes_mat$Jl/gila_mes_mat$Hd
Gila$BLs <- gila_mes_mat$Hd/gila_mes_mat$Bd
Gila$PFv <- gila_mes_mat$PFiII/gila_mes_mat$Bd
Gila$PFs <- gila_mes_mat$PFh/gila_mes_mat$Bl
Gila$CPt <- gila_mes_mat$CFd/gila_mes_mat$CPd

Gila_modesta$Genus <- rep("Gila", nrow(Gila_modesta))
Gila <- rbind(Gila_modesta, Gila)
Gila <- retype(Gila)
Gil_m <- apply(Gila[,c(6:15)], 2, mean, na.rm=TRUE) 
Gil_sp <- c(NA, NA, NA, NA, "Gila sp", Gil_m, NA, NA, NA)


fishtrait <- fishtrait[! fishtrait$Genus.species %in% c("Astyanax sp", "Poecilia sp",
                                                        "Gila sp", "Poeciliopsis sp",
                                                        "Oreochromis sp", "Chirostoma sp",
                                                        "Pseudoxiphophorus sp"),]
fishtrait_complete <- rbind(fishtrait, Ore_sp, Ast_sp, Poe_sp, Poel_sp, Heter_sp, Chir_sp, Gil_sp)
fishtrait_complete <- retype(fishtrait_complete)

save(fishtrait_complete, file="fishtrait_complete.RData")

###########################################################################################
# End of script ###########################################################################

