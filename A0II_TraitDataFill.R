###########################################################################################
# Script: Trait data fill
# AFE
# Aug-Nov 2022
###########################################################################################

# Libraries:-------------------------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(readxl)
library(hablar)
library(rfishbase)
library(data.table)

rm(list=ls())
myd <- getwd()

# 1) Fill out missing species data:--------------------------------------------------------
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/Traits/Measurements")
file_names <- dir() # 26
mes <- do.call(rbind,lapply(file_names,read.csv))   # 312                                # All data
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2")   # Go back to main dir

file_names <- str_split_fixed(file_names, "_", 2)[,2]
file_names <- str_split_fixed(file_names, "\\.", 2)[,1]
file_names <- gsub("_", " ", file_names)

mes$S <- rep(file_names, each=12) # Retrieve measurements
mes$M <- rep(c("Bl", "Bd", "Hd", "Ed", "CPd", "CFd", "Eh", "Mo", "Jl", "PFh", "PFi", "PFiII"), times=26)

sum(mes$Length==0)  # 6
sum(mes$Length <0)  # 0
mes$Length[mes$Length==0] <- NA # Measurements that could not be taken
sum(is.na(mes$Length))

mes <- mes[c("S", "M", "Length")]
names(mes) <- c("Species", "Measure", "Length")

mes_mat <- spread(mes, key="Measure", value="Length")

# 2) Corrected measurements Nov (by trait):------------------------------------------------
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/Traits/Measurements_Revised")
file_namesC <- dir() 
mesC <- do.call(rbind,lapply(file_namesC,read.csv))
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2") 

file_namesC <- gsub("_", " ", str_split_fixed(file_namesC, "_", 2)[,2])
file_namesC <- str_split_fixed(file_namesC, "\\.", 2)[,1]

mesC$S <- paste0(str_split_fixed(file_namesC, " ", 3)[,1], " ", str_split_fixed(file_namesC, " ", 3)[,2])
mesC$M <- str_split_fixed(file_namesC, " ", 3)[,3]

sort(unique(mesC$M)) # "CFd" "CPd" "Ed"  "Eh"  "Hd"  "Jl"  "PFl"


# CFd:------------------------------------------------------------------------------------
unique(mesC$S[mesC$M=="CFd"]) # P baenschi
mes_mat$CFd[mes_mat$Species=="Poeciliopsis baenschi"] <- mesC$Length[mesC$S=="P baenschi" & mesC$M=="CFd"]

# CPd:------------------------------------------------------------------------------------
unique(mesC$S[mesC$M=="CPd"])
# "I mexicanus" "X resolanae"
mes_mat$CPd[mes_mat$Species=="Ictalurus mexicanus"] <- mesC$Length[mesC$S=="I mexicanus" & mesC$M=="CPd"]
mes_mat$CPd[mes_mat$Species=="Xenotaenia resolanae"] <- mesC$Length[mesC$S=="X resolanae" & mesC$M=="CPd"]

# Ed:--------------------------------------------------------------------------------------
unique(mesC$S[mesC$M=="Ed"])
#"A balsana"       "C consocium"     "C humboldtianum" "C jordani" 
#"C labarcae"      "P jonesii"       "X variatus"

mes_mat$Ed[mes_mat$Species=="Atherinella balsana"] <- mesC$Length[mesC$S=="A balsana" & mesC$M=="Ed"]
mes_mat$Ed[mes_mat$Species=="Chirostoma consocium"] <- mesC$Length[mesC$S=="C consocium" & mesC$M=="Ed"]
mes_mat$Ed[mes_mat$Species=="Chirostoma humboldtianum"] <- mesC$Length[mesC$S=="C humboldtianum" & mesC$M=="Ed"]
mes_mat$Ed[mes_mat$Species=="Chirostoma jordani"] <- mesC$Length[mesC$S=="C jordani" & mesC$M=="Ed"]
mes_mat$Ed[mes_mat$Species=="Chirostoma labarcae"] <- mesC$Length[mesC$S=="C labarcae" & mesC$M=="Ed"]
mes_mat$Ed[mes_mat$Species=="Pseudoxiphophorus jonesii"] <- mesC$Length[mesC$S=="P jonesii" & mesC$M=="Ed"]
mes_mat$Ed[mes_mat$Species=="Xiphophorus variatus"] <- mesC$Length[mesC$S=="X variatus" & mesC$M=="Ed"]

# Eh:--------------------------------------------------------------------------------------
unique(mesC$S[mesC$M=="Eh"])
# "A balsana"       "A tamazulae"     "C charari"       "C consocium"    
# "C humboldtianum" "C jordani"       "C labarcae"      "G conspersa"    
# "G pulchra"       "N boucardi"      "P baenschi"      "P chica"        
# "P jonesii"       "P turneri"       "S lermae"        "X resolanae"    
# "X variatus"     

mes_mat$Eh[mes_mat$Species=="Atherinella balsana"] <- mesC$Length[mesC$S=="A balsana" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Allodontichthys tamazulae"] <- mesC$Length[mesC$S=="A tamazulae" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Chirostoma charari"] <- mesC$Length[mesC$S=="C charari" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Chirostoma consocium"] <- mesC$Length[mesC$S=="C consocium" & mesC$M=="Eh"]

mes_mat$Eh[mes_mat$Species=="Chirostoma humboldtianum"] <- mesC$Length[mesC$S=="C humboldtianum" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Chirostoma jordani"] <- mesC$Length[mesC$S=="C jordani" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Chirostoma labarcae"] <- mesC$Length[mesC$S=="C labarcae" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Gila conspersa"] <- mesC$Length[mesC$S=="G conspersa" & mesC$M=="Eh"]

mes_mat$Eh[mes_mat$Species=="Gila pulchra"] <- mesC$Length[mesC$S=="G pulchra" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Notropis boucardi"] <- mesC$Length[mesC$S=="N boucardi" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Poeciliopsis baenschi"] <- mesC$Length[mesC$S=="P baenschi" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Poecilia chica"] <- mesC$Length[mesC$S=="P chica" & mesC$M=="Eh"]

mes_mat$Eh[mes_mat$Species=="Pseudoxiphophorus jonesii"] <- mesC$Length[mesC$S=="P jonesii" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Poeciliopsis turneri"] <- mesC$Length[mesC$S=="P turneri" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Skiffia lermae"] <- mesC$Length[mesC$S=="S lermae" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Xenotaenia resolanae"] <- mesC$Length[mesC$S=="X resolanae" & mesC$M=="Eh"]
mes_mat$Eh[mes_mat$Species=="Xiphophorus variatus"] <- mesC$Length[mesC$S=="X variatus" & mesC$M=="Eh"]

# Hd:--------------------------------------------------------------------------------------
unique(mesC$S[mesC$M=="Hd"])
# "A balsana"       "A tamazulae"     "C charari"       "C consocium"    
# "C humboldtianum" "C jordani"       "C labarcae"      "G conspersa"    
# "G pulchra"       "N boucardi"      "P baenschi"      "P chica"        
# "P jonesii"       "P turneri"       "S lermae"        "X resolanae"    
# "X variatus"

mes_mat$Hd[mes_mat$Species=="Atherinella balsana"] <- mesC$Length[mesC$S=="A balsana" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Allodontichthys tamazulae"] <- mesC$Length[mesC$S=="A tamazulae" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Chirostoma charari"] <- mesC$Length[mesC$S=="C charari" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Chirostoma consocium"] <- mesC$Length[mesC$S=="C consocium" & mesC$M=="Hd"]

mes_mat$Hd[mes_mat$Species=="Chirostoma humboldtianum"] <- mesC$Length[mesC$S=="C humboldtianum" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Chirostoma jordani"] <- mesC$Length[mesC$S=="C jordani" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Chirostoma labarcae"] <- mesC$Length[mesC$S=="C labarcae" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Gila conspersa"] <- mesC$Length[mesC$S=="G conspersa" & mesC$M=="Hd"]

mes_mat$Hd[mes_mat$Species=="Gila pulchra"] <- mesC$Length[mesC$S=="G pulchra" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Notropis boucardi"] <- mesC$Length[mesC$S=="N boucardi" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Poeciliopsis baenschi"] <- mesC$Length[mesC$S=="P baenschi" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Poecilia chica"] <- mesC$Length[mesC$S=="P chica" & mesC$M=="Hd"]

mes_mat$Hd[mes_mat$Species=="Pseudoxiphophorus jonesii"] <- mesC$Length[mesC$S=="P jonesii" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Poeciliopsis turneri"] <- mesC$Length[mesC$S=="P turneri" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Skiffia lermae"] <- mesC$Length[mesC$S=="S lermae" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Xenotaenia resolanae"] <- mesC$Length[mesC$S=="X resolanae" & mesC$M=="Hd"]
mes_mat$Hd[mes_mat$Species=="Xiphophorus variatus"] <- mesC$Length[mesC$S=="X variatus" & mesC$M=="Hd"]

# Jl:--------------------------------------------------------------------------------------
unique(mesC$S[mesC$M=="Jl"])
# "A balsana"        "A zacapuensis"    "Aztecula sallaei" "C audax"         
# "C labarcae"       "G conspersa"      "G minaceae"       "G pulchra"       
# "X resolanae"

mes_mat$Jl[mes_mat$Species=="Atherinella balsana"] <- mesC$Length[mesC$S=="A balsana" & mesC$M=="Jl"]
mes_mat$Jl[mes_mat$Species=="Allotoca zacapuensis"] <- mesC$Length[mesC$S=="A zacapuensis" & mesC$M=="Jl"]
mes_mat$Jl[mes_mat$Species=="Aztecula sallaei"] <- mesC$Length[mesC$S=="Aztecula sallaei" & mesC$M=="Jl"]
mes_mat$Jl[mes_mat$Species=="Characodon audax"] <- mesC$Length[mesC$S=="C audax" & mesC$M=="Jl"]

mes_mat$Jl[mes_mat$Species=="Chirostoma labarcae"] <- mesC$Length[mesC$S=="C labarcae" & mesC$M=="Jl"]
mes_mat$Jl[mes_mat$Species=="Gila conspersa"] <- mesC$Length[mesC$S=="G conspersa" & mesC$M=="Jl"]
mes_mat$Jl[mes_mat$Species=="Gila minaceae"] <- mesC$Length[mesC$S=="G minaceae" & mesC$M=="Jl"]
mes_mat$Jl[mes_mat$Species=="Gila pulchra"] <- mesC$Length[mesC$S=="G pulchra" & mesC$M=="Jl"]

mes_mat$Jl[mes_mat$Species=="Xenotaenia resolanae"] <- mesC$Length[mesC$S=="X resolanae" & mesC$M=="Jl"]

# PFl:-------------------------------------------------------------------------------------
unique(mesC$S[mesC$M=="PFl"])
# "A zacapuensis"    "Aztecula sallaei" "C audax"          "G conspersa"     
# "G minaceae"       "I mexicanus"      "P chica"          "P jonesii"       
# "P turneri"        "S lermae"         "S multipunctatum" "X variatus"      
# "Y chapalae"  

mes_mat$PFiII[mes_mat$Species=="Allotoca zacapuensis"] <- mesC$Length[mesC$S=="A zacapuensis" & mesC$M=="PFl"]
mes_mat$PFiII[mes_mat$Species=="Aztecula sallaei"] <- mesC$Length[mesC$S=="Aztecula sallaei" & mesC$M=="PFl"]
mes_mat$PFiII[mes_mat$Species=="Characodon audax"] <- mesC$Length[mesC$S=="C audax" & mesC$M=="PFl"]
mes_mat$PFiII[mes_mat$Species=="Gila conspersa"] <- mesC$Length[mesC$S=="G conspersa" & mesC$M=="PFl"]

mes_mat$PFiII[mes_mat$Species=="Gila minaceae"] <- mesC$Length[mesC$S=="G minaceae" & mesC$M=="PFl"]
mes_mat$PFiII[mes_mat$Species=="Ictalurus mexicanus"] <- mesC$Length[mesC$S=="I mexicanus" & mesC$M=="PFl"]
mes_mat$PFiII[mes_mat$Species=="Poecilia chica"] <- mesC$Length[mesC$S=="P chica" & mesC$M=="PFl"]
mes_mat$PFiII[mes_mat$Species=="Pseudoxiphophorus jonesii"] <- mesC$Length[mesC$S=="P jonesii" & mesC$M=="PFl"]

mes_mat$PFiII[mes_mat$Species=="Poeciliopsis turneri"] <- mesC$Length[mesC$S=="P turneri" & mesC$M=="PFl"]
mes_mat$PFiII[mes_mat$Species=="Skiffia lermae"] <- mesC$Length[mesC$S=="S lermae" & mesC$M=="PFl"]
mes_mat$PFiII[mes_mat$Species=="Sicydium multipunctatum"] <- mesC$Length[mesC$S=="S multipunctatum" & mesC$M=="PFl"]
mes_mat$PFiII[mes_mat$Species=="Xiphophorus variatus"] <- mesC$Length[mesC$S=="X variatus" & mesC$M=="PFl"]
mes_mat$PFiII[mes_mat$Species=="Yuriria chapalae"] <- mesC$Length[mesC$S=="Y chapalae" & mesC$M=="PFl"]


# Compute ratios:--------------------------------------------------------------------------

sps26fill <- data.frame(matrix(nrow=26, ncol=17))
names(sps26fill) <- c("SpecCode", "Suporder", "Order", "Family",
                      "Genus.species", "MBl",
                      "BEl", "VEp", "REs", "OGp",
                      "RMl", "BLs","PFv","PFs","CPt",
                      "Type.of.illustration", "Reference")

sps26fill$Genus.species <- mes_mat$Species
sps26fill$Reference <- rep("Miller, R.R. (2009) Freshwater Fishes of Mexico (Peces Dulceacuicolas de Mexico). 1st edn.", nrow(sps26fill))
sps26fill$Type.of.illustration <- rep("Picture", nrow(sps26fill))
sps26fill$Type.of.illustration <- ifelse(sps26fill$Genus.species %in% c("Allotoca zacapuensis",
                                                                        "Ictalurus mexicanus",
                                                                        "Pseudoxiphophorus jonesii",
                                                                        "Sicydium multipunctatum",
                                                                        "Yuriria chapalae",
                                                                        "Gila minaceae"), "Drawing",sps26fill$Type.of.illustration)

sps26fill$BEl <- mes_mat$Bl/mes_mat$Bd   #Body elongation
sps26fill$VEp <- mes_mat$Eh/mes_mat$Bd   #Verical eye position
sps26fill$REs <- mes_mat$Ed/mes_mat$Hd   #Relative eye size
sps26fill$OGp <- mes_mat$Mo/mes_mat$Bd   #Oral gape position
sps26fill$RMl <- mes_mat$Jl/mes_mat$Hd   #Relative maxillary length
sps26fill$BLs <- mes_mat$Hd/mes_mat$Bd   #Body lateral shape
sps26fill$PFv <- mes_mat$PFh/mes_mat$Bd  #Pectoral fin vertical position
sps26fill$PFs <- mes_mat$PFiII/mes_mat$Bl #Pectoral fin size
sps26fill$CPt <- mes_mat$CFd/mes_mat$CPd  #Caudal peduncle throttling

sum(is.na(sps26fill[,c(7:15)])) # 1
save(sps26fill, file="sps26fill.RData")

