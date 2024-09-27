###########################################################################################
# Script: Sensitivity Analysis / Individual Trait effects on trait diversity
# AFE
# Nov 2022
###########################################################################################

# Libraries: ------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyverse)

rm(list=ls())
myd <- getwd()
plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots" # Dir to save main plots

# Community data:-------------------------------------------------------------------------
load(paste0(myd, "/tt.RData"))    # Trait data matrix
load("lsS5.RData")                # List of sites where to compute FD (already vector, but names in the same order)

# Functions:------------------------------------------------------------------------------
# FD and TD (Hill numbers) use the same units and TD is always greater than or equal to FD
# For q=0, Redundancy is S-FD0
hill_dir <- "C:/Users/Usuario/Documents/PHD/Functions/FunD-master/"

fskt <- read.table(paste0(hill_dir,"Plant_Abundance.txt"))
fskt = list(fs=fskt$Fushan, kt=fskt$Kenting)
dij_fskt = read.table(paste0(hill_dir,"Plant_Distance_Matrix.txt"))
dij_fskt = as.matrix(dij_fskt) # All needed for the code to run

source("C:/Users/Usuario/Documents/PHD/Functions/FunD_Rcode.R") # Warnings don't affect FD_MLE()

# Generate 10 trait matrices:------------------------------------------------------------
rownames(tt)[rownames(tt)=="Dajaus monticola"] <- "Agonostomus monticola"
tt <- tt[order(rownames(tt)),]
columns <- names(tt)
list_mats <- lapply(seq_along(columns), function(x) tt[-which(names(tt) %in% columns[x])])
str(tt[[5]])
list_mats[[11]] <- tt
names(list_mats) <- c(names(tt), "All")

list_dists <- lapply(list_mats, function(x) {as.matrix(daisy(x, metric = "euclidean"))})
list_distsn <- lapply(list_dists, function(x) {(x - min(x)) / (max(x) - min(x))})

# Compute FD:-----------------------------------------------------------------------------
class(lsS5[[1]])
list_sens <- list()
for (i in 1:length(list_distsn)){
  dist_mat <- list_distsn[[i]]
  FD <- lapply(lsS5, function(x) {FD_MLE(x, dist_mat, mean(dist_mat[dist_mat>0]), 0)}) # trait diversity
  list_sens[[i]] <- as.data.frame(do.call(rbind, FD))
}

# Correlations:---------------------------------------------------------------------------
dt <- as.data.frame(do.call(cbind, list_sens))
names(dt) <- names(list_distsn)

min(cor(dt, method = "pearson"))  # 0.9867384
min(cor(dt, method = "spearman")) # 0.9822769
