###########################################################################################
# Script: Trait distance matrices
# AFE
# August 2022
###########################################################################################

# Libraries:-------------------------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(hablar)
library(ggplot2)
library(factoextra)
library(cluster)
library(corrplot)
library(mFD)

rm(list=ls())
myd <- getwd()

# Read data:-------------------------------------------------------------------------------
load(paste0(myd, "/HNB.RData"))   # Broad Historical Native Assemblage
load(paste0(myd, "/ContN.RData")) # Contemporary Native Assemblage
load(paste0(myd, "/ContE.RData")) # Contemporary exotic assemblage

load("fishtrait_complete_imp.RData")
tt <- fishtrait_complete_imp
class(tt)

# Correlations & exploration: -------------------------------------------------------------
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
} # custom function (http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram#data-for-correlation-analysis)

C <- cor(tt)
p.mat <- cor.mtest(tt)
corrplot(C, method="number", type = "lower", order = "hclust")

file_path <- paste0(myd, "/Plots/Correlation matrix.png")
png(height=700, width=700, file=file_path)

corrplot(C, method="number", type = "lower", order = "hclust", 
         p.mat = p.mat, sig.level = 0.01, insig = "blank")

dev.off()

# Most correlations are not significant.
# The highest is 0.60 between BEl & PFv.


# Create trait distance matrix:------------------------------------------------------------
round(diag(cov(tt)), 2)                 # Need to scale variances.
tt <- as.data.frame(tt)
tt[] <- data.frame(apply(tt, 2, scale)) # Defaults are TRUE center & scale
round(diag(cov(tt)), 2) 

dist_mat1 <- as.matrix(daisy(tt, metric = "euclidean"))
dist_matS <- (dist_mat1 - min(dist_mat1)) / (max(dist_mat1) - min(dist_mat1)) # Re-scale between 0 & 1

save(dist_mat1, file="dist_mat1.RData")

qfsp <- quality.fspaces(as.dist(dist_mat1), deviation_weighting = "absolute", maxdim_pcoa = 6)
#qfsp <- quality.fspaces(as.dist(dist_mat1), deviation_weighting = "squarred", maxdim_pcoa = 10)
round(qfsp$quality_fspaces, 3) 
quality.fspaces.plot(qfsp, "mad", fspaces_plot = rownames(qfsp$quality_fspaces))
coords_qfsp <- qfsp$details_fspaces$sp_pc_coord #(same output of prcomp)

# NOTE: Quality increases consistently with increasing n? of traits.


# Mapping exotics / natives extirpated & remaining:----------------------------------------
nat <- sort(unique(names(ContN[,-c(1:2)]))) # 48
nat[nat=="Agonostomus monticola"] <-"Dajaus monticola"
ext <- sort(unique(setdiff(names(HNB), names(ContN)))) # 30
int <- sort(unique(setdiff(names(ContE), names(HNB)))) # 22 (translocated species not taken into account)

tt$Status <- rep(NA, nrow(tt))
tt$Status <- ifelse(rownames(tt) %in% nat, "N", tt$Status)
tt$Status <- ifelse(rownames(tt) %in% ext, "E", tt$Status)
tt$Status <- ifelse(rownames(tt) %in% int, "I", tt$Status)
sum(is.na(tt$Status))


# PCA:-------------------------------------------------------------------------------------
PCA <- prcomp(tt[,c(1:10)])
coords <- PCA$x
save(coords, file="coords.RData")

# Facto extra (exploration of trait space)
fviz_eig(PCA)
#26.6 + 21.5 = 3 48.1 (First two components)

groups <- as.factor(tt$Status)

# Individuals:
(iI <- fviz_pca_ind(PCA,
                   axes = c(1,2),
                   col.ind = "cos2", # Color by the quality of representation
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                   repel = TRUE
))

(iII <- fviz_pca_ind(PCA,
                    axes = c(3,4),
                    col.ind = "cos2", # Color by the quality of representation
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                    repel = TRUE     # Avoid text overlapping
))


# Variables:
(vI <- fviz_pca_var(PCA,
                   col.var = "contrib", # Color by contributions to the PC
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE     # Avoid text overlapping
))
(vII <- fviz_pca_var(PCA,
                    axes = c(3,4),
                    col.var = "contrib", # Color by contributions to the PC
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE     # Avoid text overlapping
))

# Biplots:
(bI <- fviz_pca_biplot(PCA, repel = TRUE,
                      col.var = "contrib", # Variables color
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      col.ind = "#696969",  # Individuals color
                      label = "var",
                      addEllipses = TRUE,
                      ellipse.alpha = 0.1,
                      ellipse.type = "convex"
))
(bII <- fviz_pca_biplot(PCA, 
                       axes = c(3,4),
                       repel = TRUE,
                       col.var = "contrib", # Variables color
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                       col.ind = "#696969",  # Individuals color
                       label = "var",
                       addEllipses = TRUE,
                       ellipse.alpha = 0.1,
                       ellipse.type = "convex"
))

# Eigenvalues:
eig.val <- get_eigenvalue(PCA)
eig.val # Reach > 80% var explained with first five components

# Results for Variables
PCA <- get_pca_var(PCA)
PCA$coord          # Coordinates
PCA$contrib        # Contributions to the PCs
PCA$cos2           # Quality of representation 

# NOTES:
fviz_contrib(PCA, choice = "var", axes = 1, top = 10) # Body elongation (hydrodynamism) and pectoral fin use for swimming
# Dim 1:BEl, PFv,  PFs
# Hydrodynamism and use of pectoral fins
fviz_contrib(PCA, choice = "var", axes = 2, top = 10) # Visual acuity, maximum body length but also size of mouth and strengh of jaw and hydrodynamism and head size
# Dim 2:REs, MBl,  RMl, BLs
# ~Trophic level (predator-prey)
fviz_contrib(PCA, choice = "var", axes = 3, top = 10)
# Dim 3:VEp, OGp,  PFs
# Feeding strategy in the water column
fviz_contrib(PCA, choice = "var", axes = 4, top = 10)
# Dim 4:CPt(mainly), MBl
# caudal propulsion efficiency through reduction of drag, swimming mode (also somewhat affected by body size)

# Results for individuals
PCA <- get_pca_ind(PCA)
PCA$coord          # Coordinates
PCA$contrib        # Contributions to the PCs
PCA$cos2           # Quality of representation 

# PLot Individuals (Status):
(iIII <- fviz_pca_ind(PCA,
                      label = "none", # hide individual labels
                      habillage = tt$Status, # color by groups
                      palette = c("#0072B2", "#F0E442", "Darkgray"),
                      addEllipses = FALSE # TRUE for concentration ellipses
))

# NOTES:
# Dim 1: seems to be a distinction in hydrodynamism between introduced and extirpated fish, with remaining natives showing more middle values
# Dim 2: here natives appear to have slightly more different values than extirpated and introduced fish.
# PERMANOVAs??

(iIV <- fviz_pca_ind(PCA,
                     axes = c(3,4),
                     label = "none", # hide individual labels
                     habillage = tt$Status, # color by groups
                     palette = c("#0072B2", "#F0E442", "Darkgray"),
                     addEllipses = FALSE # TRUE for concentration ellipses
))

# NOTES:
# Dim 3: seems to be a distinction in feeding strategy in the water column between introduced and extirpated fish, with remaining natives showing more middle values
# Dim 4: clear distinction between groups not observed along this dimension
# PERMANOVAs??

# Save plots: -----------------------------------------------------------------------------
# TBC

###########################################################################################
# End of script ###########################################################################

