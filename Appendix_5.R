###########################################################################################
# Script: Sensitivity Analysis/ Observer test II
# AFE
# Nov 2022
###########################################################################################

# Libraries:------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(missForest)
library(factoextra)
rm(list=ls())

# Read Data:------------------------------------------------------------------------------
FMorph <- read.csv2("C:/Users/Usuario/Documents/PHD/Datasets/FishMORPH/FISHMORPH_Database.csv", h=T)
FMorph$Genus <- str_split_fixed(FMorph$Genus, " ", 2)[,1]
F_Morph_subset <- subset(FMorph, FMorph$Genus %in% c("Algansea",
                                                             "Allodontichthys",
                                                             "Allotoca",
                                                             "Atherinella",
                                                             "Aztecula",
                                                             "Characodon",
                                                             "Chirostoma",
                                                             "Gila",
                                                             "Ictalurus",
                                                             "Moxostoma",
                                                             "Notropis",
                                                             "Poecilia",
                                                             "Poeciliopsis",
                                                             "Pseudoxiphophorus",
                                                             "Sicydium",
                                                             "Skiffia",
                                                             "Xenotaenia",
                                                             "Xiphophorus",
                                                             "Yuriria"))
# NOTE: All found but Xenotaenia.

##########################################################################################
# Principal component analysis:-----------------------------------------------------------
str(F_Morph_subset)
#F_Morph_subset <- subset(F_Morph_subset, !F_Morph_subset$Genus %in% "Gila")
picture_traits <- F_Morph_subset[,c(7:15)]
picture_traits <- hablar::retype(picture_traits)
rownames(picture_traits) <- F_Morph_subset$Genus.species
str(picture_traits)
picture_traitsII <- missForest(as.matrix(picture_traits))

PCAFM <- prcomp(picture_traitsII$ximp, scale=T)

(bI_PCAFM <- fviz_pca_ind(PCAFM,
                          title = "",
                          label = "none", # hide individual labels
                          habillage = F_Morph_subset$Genus, # color by groups
                          addEllipses = TRUE, # TRUE for concentration ellipses
                          legend.title="Genus",
                          ellipse.type = "convex"
)) 

(bII_PCAFM <- fviz_pca_ind(PCAFM,
                          title = "",
                          axes = c(3,4),
                          label = "none", # hide individual labels
                          habillage = F_Morph_subset$Genus, # color by groups
                          addEllipses = TRUE, # TRUE for concentration ellipses
                          legend.title="Genus",
                          ellipse.type = "convex"
)) 

##########################################################################################
# Principal component analysisII:-----------------------------------------------------------
load("sps26fill.RData")

F_Morph_subset <- hablar::retype(F_Morph_subset)
F_Morph_subset$Observer <- "FM"
str(F_Morph_subset)

sps26fill$Genus <- str_split_fixed(sps26fill$Genus.species, " ", 2)[,1]
sps26fill$Observer <- "AFE"
str(sps26fill)

FMorph_AFE <- bind_rows(F_Morph_subset, sps26fill) # combine both datasets
picture_traits_combined <- FMorph_AFE[,c(7:15)]
rownames(picture_traits_combined) <- FMorph_AFE$Genus.species
picture_traits_combined_II <- missForest(as.matrix(picture_traits_combined))

PCAFMAFE <- prcomp(picture_traits_combined_II$ximp, scale=T) # PCA combined

(screep <- fviz_eig(PCAFMAFE, addlabels=TRUE, hjust = -0.3))

(bI_PCAFMAFE <- fviz_pca_ind(PCAFMAFE,
                                title = "",
                                label="none",
                                habillage = FMorph_AFE$Genus, # color by groups
                                addEllipses = TRUE, # TRUE for concentration ellipses
                                legend.title="Genus",
                                ellipse.type = "convex",
                                alpha.ind = 0.5))
observer_groups <- FMorph_AFE$Observer
observer_groups[observer_groups=="FM"] <- ""
(bI_PCAFMAFE <- bI_PCAFMAFE + geom_text(aes(label = observer_groups),
                                        alpha = 0.5, size = 2, nudge_y = 0.1, show.legend = FALSE))

(bI_PCAFMAFE_3_4 <- fviz_pca_ind(PCAFMAFE,
                             axes = c(3,4),
                             title = "",
                             label="none",
                             habillage = FMorph_AFE$Genus, # color by groups
                             addEllipses = TRUE, # TRUE for concentration ellipses
                             legend.title="Genus",
                             ellipse.type = "convex",
                             alpha.ind = 0.5))
observer_groups <- FMorph_AFE$Observer
observer_groups[observer_groups=="FM"] <- ""
(bI_PCAFMAFE_3_4 <- bI_PCAFMAFE_3_4 + geom_text(aes(label = observer_groups),
                                      alpha = 0.5, size = 2, nudge_y = 0.1, show.legend = FALSE))

##########################################################################################
# Final plot:-----------------------------------------------------------------------------

(finalp <- fviz_pca_ind(PCAFMAFE,
                             title = "",
                             label="none",
                             habillage = FMorph_AFE$Genus, # color by groups
                             addEllipses = TRUE, # TRUE for concentration ellipses
                             legend.title="Genus",
                             ellipse.type = "convex",
                             alpha.ind = 0.5,
                             repel=TRUE, ellipse.alpha = 0.02))
(finalp <- finalp + theme(legend.position = "none"))

observer_groups <- FMorph_AFE$Observer
observer_groups[observer_groups=="FM"] <- ""
observer_groups[observer_groups=="AFE"] <- "M"
(finalp <- finalp + geom_text(aes(label = observer_groups),
                                        alpha = 0.5, size = 3, nudge_y = 0.1, show.legend = FALSE))

plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots" # Dir to save main plots
ggsave(finalp, filename = paste0(plot_dir, "/SM/Ms/TestObserver_App5.jpg"), width=10, height=8)

