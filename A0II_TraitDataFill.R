#===============================================================================
# Script: Trait Data
# AFE
# update August 2024
#===============================================================================




#===============================================================================
# Libraries:--------------------------------------------------------------------
#===============================================================================
library(dplyr)
library(tidyverse)
library(readxl)
library(hablar)
library(rfishbase)
library(data.table)


rm(list=ls())


path_lists <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists"
path_traits <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Traits"




#===============================================================================
# Read FishMORPH: --------------------------------------------------------------
#===============================================================================

FMorph <- read.csv2("C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_GeneralMethods/FishMORPH/FISHMORPH_Database.csv")
load(paste0(path_lists, "/l67.RData")) # subset of 67 sites
load(paste0(path_lists, "/l53.RData")) # subset of 53 sites




#===============================================================================
# Retrieve measurements available in FishMORPH: --------------------------------
#===============================================================================

## Species for which data are available:----------------------------------------
vecsps <- as.vector(sort(unique(names(l67)[-c(1:3)]))) # 102
genus_sp <- c("Astyanax sp","Gila sp", "Pseudoxiphophorus sp", "Chirostoma sp", "Oreochromis sp", "Poecilia sp","Poeciliopsis sp")
vecsps <- setdiff(vecsps, genus_sp)                    # 95 (- 7) 
vecsps <- plyr::revalue(vecsps, c("Tampichthys dichroma"="Tampichthys dichromus")) # as T. dichromus in FISHMORPH
vecsps <- plyr::revalue(vecsps, c("Dajaus monticola"="Agonostomus monticola"))     # as A. monticola in FISHMORPH
vecsps <- plyr::revalue(vecsps, c("Amphilophus istlanus"="Cichlasoma istlanum"))   # as C. istlanum in FISHMORPH
vecsps <- plyr::revalue(vecsps, c("Herichthys bartoni"="Nosferatu bartoni"))       # as N. bartoni in FISHMORPH
vecsps <- plyr::revalue(vecsps, c("Herichthys labridens"="Nosferatu labridens"))   # as N. labridens in FISHMORPH


## Retrieve: -------------------------------------------------------------------
FMorphSubs <- subset(FMorph, FMorph$Genus.species %in% vecsps) # 70

# Try other updated & non updated names instead:
NP <- as.vector(setdiff(vecsps, FMorphSubs$Genus.species))     # 24
NP
FMorphSubsMissing <- NP

save(FMorphSubs, file=paste0(path_traits, "/FMorphSubs.RData"))
save(FMorphSubsMissing, file=paste0(path_traits, "/FMorphSubsMissing.RData"))


## Check unavailable measurements: ---------------------------------------------
sapply(FMorphSubs[,c(6:15)], function(x) sum(is.na(x)))
# MBl BEl VEp REs OGp RMl BLs PFv PFs CPt 
# 10   0   0   0   0   2   0   1   9   0 




#===============================================================================
# Retrieve data (main FMorphSubsMissing): --------------------------------------
#===============================================================================
setwd("C:/Users/afe1/Dropbox/CentralMexicoFish/0_FinalAug2024") # dir with fish trait measurements

## Read Excel sheets from FIJI: ------------------------------------------------
filenames <- list.files(pattern=".csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)
names(ldf) <- str_split_fixed(str_split_fixed(filenames, "_",2)[,2], "\\.",2)[,1]
ldf <- lapply(ldf, function(x) {x<- within(x, rm(X, Area, Mean, Min, Max))})
ldf <- lapply(ldf, function(x) {x$M <- c("TL", "SL", "BD", "HD", "EDv", "EDh", "EP", "HL",
                                         "PED_D", "CAU_D", "SNL", "Mo", "Jl", "PFp", "PFl");x})
ldf2 <- do.call(rbind, ldf)
ldf2$Pic <- rep(names(ldf), each=15)

ldf2 <- ldf2[,!names(ldf2) == "Angle"]
rownames(ldf2) <- 1:nrow(ldf2)
ldf2 <- spread(ldf2, key="M", value="Length")
ldf2 <- hablar::retype(ldf2)


## Read corrected measurements: ------------------------------------------------
setwd("C:/Users/afe1/Dropbox/CentralMexicoFish/0_FinalAug2024_Corr")
filenamesC <- list.files(pattern=".csv", full.names=TRUE)
ldfC <- lapply(filenamesC, read.csv)
names(ldfC) <- str_split_fixed(str_split_fixed(filenamesC, "_",2)[,2], "\\.",2)[,1]
ldfC <- lapply(ldfC, function(x) {x<- within(x, rm(X, Area, Mean, Min, Max, Angle))})


## Add corrections: ------------------------------------------------------------
names(ldfC)
ldf2$Jl[ldf2$Pic=="Allodontichthys_tamazulae_MaleM"] <- ldfC$Allodontichthys_tamazulae_MaleM_JlPFl$Length[1]    # Jl A tamazulae
ldf2$PFl[ldf2$Pic=="Allodontichthys_tamazulae_MaleM"] <- ldfC$Allodontichthys_tamazulae_MaleM_JlPFl$Length[2]   # PFl A tamazulae

ldf2$PFl[ldf2$Pic=="Allotoca_zacapuensis_MaleM"] <- ldfC$Allotoca_zacapuensis_MaleM_PFl$Length[1]               # PFl A zacapuensis

ldf2$Jl[ldf2$Pic=="Chirostoma_consocia_AdultM"] <- ldfC$Chirostoma_consocia_AdultM_Jl$Length[1]                 # Jl C consocia

ldf2$Jl[ldf2$Pic=="Chirostoma_humboldtianum_AdultM"] <- ldfC$Chirostoma_humboldtianum_AdultM_Jl_PFp_PFl$Length[1]    # Jl C humboldtianum
ldf2$PFp[ldf2$Pic=="Chirostoma_humboldtianum_AdultM"] <- ldfC$Chirostoma_humboldtianum_AdultM_Jl_PFp_PFl$Length[2]   # PFp C humboldtianum
ldf2$PFl[ldf2$Pic=="Chirostoma_humboldtianum_AdultM"] <- ldfC$Chirostoma_humboldtianum_AdultM_Jl_PFp_PFl$Length[3]   # PFl C humboldtianum

ldf2$Jl[ldf2$Pic=="Chirostoma_jordani_AdultM"] <- ldfC$Chirostoma_jordani_AdultM_Jl$Length[1]                   # Add Jl of C jordani

ldf2$PFl[ldf2$Pic=="Chirostoma_labarcae_AdultM"] <- ldfC$Chirostoma_labarcae_AdultM_PFl$Length[1]               # PFl of C labarcae

ldf2$PFl[ldf2$Pic=="Chirostoma_mezquital_AdultM"] <- ldfC$Chirostoma_mezquital_AdultM_PFl$Length[1]             # PFl of C mezquital

ldf2$Jl[ldf2$Pic=="Heterandria_jonesii_MaleM"] <- ldfC$Heterandria_jonesii_MaleM_Jl$Length[1]                   # Jl H jonesii

ldf2$Jl[ldf2$Pic=="Poecilia_chica_MaleM"] <- ldfC$Poecilia_chica_MaleM_Jl$Length[1]                             # Jl P chica

ldf2$Jl[ldf2$Pic=="Poeciliopsis_turneri_MaleHolotipeM"] <- ldfC$Poeciliopsis_turneri_MaleHolotipeM_Jl$Length[1] # Jl P turneri

ldf2$Jl[ldf2$Pic=="Skiffia_lermae_MaleM"] <- ldfC$Skiffia_lermae_MaleM_Jl$Length[1]                             # Jl S lermae


## Add 0s & NAs: ---------------------------------------------------------------
ldf2$CAU_D[ldf2$Pic=="Atherinella_balsana_AdultHolotipeM"] <- NA                                                # broken caudal A. balsana
ldf2$CAU_D[ldf2$Pic=="Hybopsis_boucardi_MaleM"] <- NA                                                           # folded caudal H. boucardi
ldf2$CAU_D[ldf2$Pic=="Gila_pulchra_MaleM"] <- NA                                                                # folded caudal G. pulchra

setwd("C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2") # re-set wd




#===============================================================================
# Append data: -----------------------------------------------------------------
#===============================================================================
Sps27 <- data.frame(matrix(nrow=27, ncol=17))

## Meta:------------------------------------------------------------------------
names(Sps27) <- c("SpecCode", "Suporder", "Order", "Family",
                  "Genus.species", "MBl",
                  "BEl", "VEp", "REs", "OGp",
                  "RMl", "BLs","PFv","PFs","CPt",
                  "Type.of.illustration", "Reference")
Sps27$Genus.species <- paste0(str_split_fixed(ldf2$Pic, "_", 3)[,1], " ", str_split_fixed(ldf2$Pic, "_", 3)[,2]) 
# rest of meta is NA, I will make tables for supplementary

## Ratios:----------------------------------------------------------------------
# traits in FM
Sps27$BEl <- ldf2$SL/ldf2$BD        # BEl: hydrodynamism
Sps27$VEp <- ldf2$EP/ldf2$BD        # VEp: position of fish and/or its prey in the water column
Sps27$REs <- ldf2$EDv/ldf2$HD       # REs: visual acuity
Sps27$OGp <- ldf2$Mo/ldf2$BD        # OGp: feeding position in the water column
Sps27$RMl <- ldf2$Jl/ldf2$HD        # RMl: size of mouth and strength of jaw
Sps27$BLs <- ldf2$HD/ldf2$BD        # BLs: hydrodynamism and head size
Sps27$PFv <- ldf2$PFp/ldf2$BD       # PFv: pectoral fin use for swimming
Sps27$PFs <- ldf2$PFl/ldf2$SL       # PFs: pectoral fin use for swimming
Sps27$CPt <- ldf2$CAU_D/ldf2$PED_D  # CPt: caudal propulsion efficiency through reduction of drag


## Append: ---------------------------------------------------------------------
identical(names(FMorphSubs), names(Sps27))    # TRUE

dt <- as.data.frame(rbind(FMorphSubs, Sps27)) # 97
97-3 # 94, OK
dt <- dt %>% hablar::retype()
str(dt)

## Re-name to match l67 & l53: -------------------------------------------------
setdiff(sort(unique(dt$Genus.species)), sort(unique(names(l67))))
setdiff(sort(unique(names(l67))), sort(unique(dt$Genus.species)))

dt$Genus.species[dt$Genus.species=="Agonostomus monticola"] <- "Dajaus monticola"
dt$Genus.species[dt$Genus.species=="Chirostoma consocia"] <- "Chirostoma consocium"
dt$Genus.species[dt$Genus.species=="Cichlasoma istlanum"] <- "Amphilophus istlanus"
dt$Genus.species[dt$Genus.species=="Heterandria jonesii"] <- "Pseudoxiphophorus jonesii"
dt$Genus.species[dt$Genus.species=="Hybopsis boucardi"] <- "Notropis boucardi"
dt$Genus.species[dt$Genus.species=="Nosferatu bartoni"] <- "Herichthys bartoni"
dt$Genus.species[dt$Genus.species=="Nosferatu labridens"] <- "Herichthys labridens"
dt$Genus.species[dt$Genus.species=="Scartomyzon austrinus"] <- "Moxostoma austrinum"
dt$Genus.species[dt$Genus.species=="Tampichthys dichromus"] <- "Tampichthys dichroma"

setdiff(sort(unique(dt$Genus.species)), sort(unique(names(l67)))) # OK
setdiff(sort(unique(names(l67))), sort(unique(dt$Genus.species))) # OK


## Unidentified records: -------------------------------------------------------
Av7 <- data.frame(matrix(nrow=7, ncol=17))
names(Av7) <- c("SpecCode", "Suporder", "Order", "Family",
                  "Genus.species", "MBl",
                  "BEl", "VEp", "REs", "OGp",
                  "RMl", "BLs","PFv","PFs","CPt",
                  "Type.of.illustration", "Reference")

Av7$Genus.species <- c("Astyanax sp","Gila sp", "Pseudoxiphophorus sp", "Chirostoma sp", "Oreochromis sp", "Poecilia sp","Poeciliopsis sp")

Av7[Av7$Genus.species %like% "Astyanax", c(7:15)] <- dt[dt$Genus.species == "Astyanax mexicanus", c(7:15)] 
# Astyanax sp idem Astyanax mexicanus

Av7[Av7$Genus.species %like% "Gila", c(7:15)] <- as.vector(apply(dt[dt$Genus.species %like% "Gila", c(7:15)], 2, mean, na.rm = TRUE))  
# Gila sp av three sps

#Av7$BLs[Av7$Genus.species %like% "Gila"]
#mean(c(dt$BLs[dt$Genus.species =="Gila pulchra"], dt$BLs[dt$Genus.species =="Gila conspersa"], dt$BLs[dt$Genus.species =="Gila minaceae"])) # check, OK

Av7[Av7$Genus.species %like% "Oreochromis", c(7:15)] <- dt[dt$Genus.species == "Oreochromis aureus", c(7:15)]      
# Oreochromis sp idem Oreochromis aureus

Av7[Av7$Genus.species %like% "Chirostoma", c(7:15)] <- dt[dt$Genus.species == "Chirostoma riojai", c(7:15)]      
# Chirostoma sp idem Chirostoma riojai (in one site, and to C. humboldtiana in another) [FIXED IN A0II_Community Data & downstream]

Av7[Av7$Genus.species %like% "Pseudoxiphophorus", c(7:15)] <- dt[dt$Genus.species == "Pseudoxiphophorus bimaculatus", c(7:15)]      
# Pseudoxiphophorus sp idem	Pseudoxiphophorus bimaculatus

Av7[Av7$Genus.species %like% "Poeciliopsis", c(7:15)] <- dt[dt$Genus.species == "Poeciliopsis gracilis", c(7:15)]      
# Poeciliopsis sp idem Poeciliopsis gracilis

Av7[Av7$Genus.species %like% "Poecilia", c(7:15)] <- dt[dt$Genus.species == "Poecilia sphenops", c(7:15)]      
# Poecilia sp idem Poecilia sphenops (P sp found in el Rincon - not in 67 - should be P. butleri)


## Append unidentified records: ------------------------------------------------
identical(names(dt), names(Av7))    # TRUE

dt <- as.data.frame(rbind(dt, Av7)) # 104
dt <- dt[! dt$Genus.species %in% c("Gila conspersa", "Gila minaceae", "Gila pulchra"),] # 101
str(dt)
names(dt)
dt <- dt[,names(dt) %in% c("Genus.species", "MBl", "BEl", "VEp", "REs", "OGp", "RMl", "BLs",
                          "PFv", "PFs", "CPt")]




#===============================================================================
# Body size: -------------------------------------------------------------------
#===============================================================================
# NOTE: length of some species only available on one of the two sources 
# (FishBase or Miller)
BS <- data.frame(matrix(nrow=101, ncol=4))
names(BS) <- c("Genus.species", "MaxLengthSL_FB", "MaxLP_Miller2009", "Final")
BS$Genus.species <- dt$Genus.species


## Fishbase (SL, max observed or estimated):------------------------------------
mbls <- rfishbase::estimate(dt$Genus.species)
gilas <- rfishbase::estimate(c("Gila conspersa", "Gila pulchra", "Gila minaceae")) # no minaceae, & Miller av. estimate larger anyway
setdiff(sort(unique(dt$Genus.species)), sort(unique(mbls$Species))) # missing "sp" & C. mezquital
BS$MaxLengthSL_FB <- mbls$MaxLengthSL[match(BS$Genus.species, mbls$Species)]
sum(is.na(BS$MaxLengthSL_FB)) # 8
BS$MaxLengthSL_FB[is.na(BS$MaxLengthSL_FB)] <- 0


## Miller 2009 reported SLs:----------------------------------------------------
#BSMiller2009Template <- BS[names(BS) %in% c("Genus.species","MaxLP_Miller2009")]
#write.csv(BSMiller2009Template,
#          row.names = F,
#          file="C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Traits/BSMiller2009Template.csv")
BSMiller2009 <- read.csv("C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Traits/BSMiller2009TemplateComplete.csv")
# NOTE: despite being native, no LP reported for S. austrinus in Miller 2009
# NOTE: av for three sps of Gila already in Excel value, and vals for other unidentified records too.

BSMiller2009 <- BSMiller2009[!BSMiller2009$Genus.species == "Herichthys cyanoguttatus",] # forgot to remove in first round of revision Aug 2024
BS$MaxLP_Miller2009 <- BSMiller2009$MaxLP_Miller2009[match(BS$Genus.species, BSMiller2009$Genus.species)]
sum(is.na(BS$MaxLP_Miller2009)) # 10
BS$MaxLP_Miller2009[is.na(BS$MaxLP_Miller2009)] <- 0
BS$Final <- ifelse(BS$MaxLengthSL_FB >= BS$MaxLP_Miller2009, BS$MaxLengthSL_FB, BS$MaxLP_Miller2009)
sum(BS$Final==0) # 1 (O sp)

BS$Final[BS$Genus.species=="Oreochromis sp"] <- BS$Final[BS$Genus.species=="Oreochromis aureus"]      # final from FB
BS$Final[BS$Genus.species=="Poeciliopsis sp"] <- BS$Final[BS$Genus.species=="Poeciliopsis gracilis"]  # FB estimate is larger
BS$Final[BS$Genus.species=="Chirostoma sp"] <- BS$Final[BS$Genus.species=="Chirostoma riojai"]        # FB estimate is larger
# The three above are changed because FB had not read "sp", but the FB estimate for 
# the species used as proxy is larger than the estimate given by Miller 2009


## Append body size data: ------------------------------------------------------
identical(sort(unique(dt$Genus.species)), sort(unique(BS$Genus.species))) # TRUE
dt$MBl <- BS$Final[match(dt$Genus.species, BS$Genus.species)]
sapply(dt[,c(2:11)], function(x) sum(is.na(x))) # mostly PFs...
# MBl BEl VEp REs OGp RMl BLs PFv PFs CPt 
# 0   0   0   0   0   2   0   2  12   2
(2+2+12+2)/(101*10)*100 # 1.78%



#===============================================================================
# Some sps treated effectively the same...--------------------------------------
#===============================================================================
# NOTE: this could cause issues in trait space..., so rm:

unique(dt$Genus.species)
dt <- dt[!dt$Genus.species %in% c("Astyanax sp", "Pseudoxiphophorus sp", 
                                  "Chirostoma sp", "Oreochromis sp", 
                                  "Poecilia sp", "Poeciliopsis sp"),] 
101-6 # 95 OK
#NOTE: rm all but Gila sp because they have the exact same traits as other sps.


sapply(dt[,c(2:11)], function(x) sum(is.na(x))) # mostly PFs...
# MBl BEl VEp REs OGp RMl BLs PFv PFs CPt 
# 0   0   0   0   0   2   0   1  9   2
(2+1+9+2)/(95*10)*100 # 1.5%           

#===============================================================================
# Save data: -------------------------------------------------------------------
#===============================================================================
fishtrait_complete <- dt
save(fishtrait_complete, file=paste0(path_traits, "/fishtrait_complete.RData")) 



# End of script ################################################################
