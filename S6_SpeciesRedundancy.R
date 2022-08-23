###########################################################################################
# Script: Compute species distinctiveness and uniqueness
# AFE
# August 2022
###########################################################################################

# Libraries:------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(mFD)
library(funrar)

rm(list=ls())
myd <- getwd()
plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots" # Dir to save main plots

# Community data:-------------------------------------------------------------------------
load(paste0(myd, "/HNC.RData"))     # Conservative Historical Native Assemblage
load(paste0(myd, "/HNB.RData"))     # Broad Historical Native Assemblage
load(paste0(myd, "/ContN.RData"))   # Contemporary Native Assemblage
load(paste0(myd, "/ContE.RData"))   # Contemporary Native Assemblage
load(paste0(myd, "/ContAll.RData")) # Contemporary Exotic assemblage
load(paste0(myd, "/All.RData"))     # All Data

load(paste0(myd, "/dist_mat1.RData"))  # Trait distance matrix
load(paste0(myd, "/tt.RData"))

# Format community data: -----------------------------------------------------------------
# Site-species matrix
l <- list("HNB"=HNB, "ContAll"=ContAll)
l <- lapply(l, function(x) {rownames(x) <- x$SiteNameE;x
                            within(x, rm(SiteNameE, DrainageBasinE))
                            as.matrix(x)}) #OK

rownames(All) <- paste0(All$SiteNameE,"_", All$Period)
rownames(All) <- gsub(" ", "_", rownames(All))
All <- within(All, rm(SiteNameE, DrainageBasinE, Period))
All <- as.matrix(All)


colnames(All)[colnames(All) == "Agonostomus monticola"] <- "Dajaus monticola"
All <- All[,order(colnames(All))]
dist_mat1 <- dist_mat1[order(rownames(dist_mat1)), order(colnames(dist_mat1))]


tt <- tt[order(rownames(tt)),]
tt <- as.matrix(tt)
identical(colnames(All), colnames(dist_mat1))  # TRUE
identical(colnames(All), rownames(tt))         # TRUE


# Uniqueness (regional):------------------------------------------------------------------
Uii <- uniqueness(All, dist_matrix = dist_mat1)

#Ui_dim <- uniqueness_dimensions(t(All), 
#                                tt,
#                                metric="euclidean")


# Native/Introduced/Remaining
nat <- sort(unique(names(ContN[,-c(1:2)]))) # 48
nat[nat=="Agonostomus monticola"] <-"Dajaus monticola"
ext <- sort(unique(setdiff(names(HNB), names(ContN)))) # 30
int <- sort(unique(setdiff(names(ContE), names(HNB)))) # 22 (translocated species not taken into account)

Uii$Status <- rep(NA, nrow(Uii))
Uii$Status <- ifelse(Uii$species %in% nat, "Native Remaining", Uii$Status)
Uii$Status <- ifelse(Uii$species %in% ext, "Extirpated", Uii$Status)
Uii$Status <- ifelse(Uii$species %in% int, "Introduced", Uii$Status)
sum(is.na(Uii$Status))
str(Uii)
Uii$Status <- as.factor(Uii$Status)

Uii$StatusII <- rep("Native", nrow(Uii))
Uii$StatusII <- ifelse(Uii$Status %in% "Introduced", "Introduced", Uii$StatusII)

Uii$StatusIII <- rep("Others", nrow(Uii))
Uii$StatusIII <- ifelse(Uii$Status %in% "Extirpated", "Extirpated", Uii$StatusIII)

# Plot

(p1 <- ggplot(Uii, aes(x=Ui, y=reorder(species, Ui))) + 
  geom_point(stat='identity', aes(col=Status), size=2, alpha=0.5)+
  labs(x="Uniqueness", y="Species")+
  scale_color_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
  theme(axis.text = element_text(size = 20))+
  theme_bw())
(p2 <- ggplot(Uii, aes(x=Ui, y=reorder(species, Ui))) + 
    geom_point(stat='identity', aes(col=StatusII), size=2, alpha=0.5)+
    labs(x="Uniqueness", y="Species")+
    scale_color_manual(values=c("#F0E442", "Darkgray"))+
    theme(axis.text = element_text(size = 20))+
    theme_bw())
(p3 <- ggplot(Uii, aes(x=Ui, y=reorder(species, Ui))) + 
    geom_point(stat='identity', aes(col=StatusIII), size=2, alpha=0.5)+
    labs(x="Uniqueness", y="Species")+
    scale_color_manual(values=c("#0072B2", "Darkgray"))+
    theme(axis.text = element_text(size = 20))+
    theme_bw())


