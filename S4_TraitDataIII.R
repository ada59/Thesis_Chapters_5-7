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
library(gridExtra)

rm(list=ls())
myd <- getwd()
plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots" # Dir to save main plots

# Read data:-------------------------------------------------------------------------------
load(paste0(myd, "/HNB.RData"))   # Broad Historical Native Assemblage
load(paste0(myd, "/ContN.RData")) # Contemporary Native Assemblage
load(paste0(myd, "/ContE.RData")) # Contemporary exotic assemblage
load(paste0(myd, "/NDT67.RData")) # Double-check translocated species

load("fishtrait_complete_imp.RData")
tt <- fishtrait_complete_imp
class(tt)

# Correlations & exploration: -------------------------------------------------------------
#cor.mtest <- function(mat, ...) {
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
# Other : https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
C <- cor(tt)
p.mat <- cor.mtest(tt)

file_path <- paste0(myd, "/Plots/Correlation matrix.png")
png(height=700, width=700, file=file_path)

corrplot(C, method="color", type = "lower", order = "hclust", 
         p.mat = p.mat$p, sig.level = 0.05, insig = "blank" ,diag=FALSE)$corrPos -> p1
text(p1$x, p1$y, round(p1$corr, 2))

dev.off()

# The highest positives are 0.51 between BEl & PFv and OGp & REs.
# The highest negative is -0-48 between MBl and REs.
# All traits kept.


# Create trait distance matrix:------------------------------------------------------------
round(diag(cov(tt)), 2)                 # Need to scale variances.
tt <- as.data.frame(tt)
tt[] <- data.frame(apply(tt, 2, scale)) # Defaults are TRUE center & scale
round(diag(cov(tt)), 2) 

save(tt, file="tt.RData")

dist_mat1 <- as.matrix(daisy(tt, metric = "euclidean"))
dist_mat1 <- (dist_mat1 - min(dist_mat1)) / (max(dist_mat1) - min(dist_mat1)) # Re-scale between 0 & 1
range(dist_mat1) #OK

save(dist_mat1, file="dist_mat1.RData")

# Assess quality of trait space
qfsp <- quality.fspaces(as.dist(dist_mat1), deviation_weighting = "absolute", maxdim_pcoa = 10)
#qfsp <- quality.fspaces(as.dist(dist_mat1), deviation_weighting = "squarred", maxdim_pcoa = 10)
round(qfsp$quality_fspaces, 3) 
quality.fspaces.plot(qfsp, "mad", fspaces_plot = rownames(qfsp$quality_fspaces))
coords_qfsp <- qfsp$details_fspaces$sp_pc_coord #(same output of prcomp)

# NOTE: Quality increases consistently with increasing number of traits.

# Classifying species in exotics/natives_extirpated/natives_remaining (REGIONAL):----------
nat <- sort(unique(names(ContN[,-c(1:2)]))) # 48 
nat[nat=="Agonostomus monticola"] <-"Dajaus monticola"
ext <- sort(unique(setdiff(names(HNB), names(ContN)))) # 30
int <- sort(unique(setdiff(names(ContE), names(HNB)))) # 22
48+22+30 # 100, OK

write.csv2(ext, file="ext.csv", row.names = FALSE)

# Any species translocated?
intersect(names(ContE), names(HNB))

#View(NDT67) [these are in the exotics column]
# "Chapalichthys encaustus"
# "Goodea atripinnis"
# "Poecilia butleri"
# "Poecilia sphenops"

intersect(names(ContE), names(HNB)) %in% nat
intersect(names(ContE), names(HNB)) %in% ext # Chapalichthys encaustus translocate din one site but extipated from its native range?

tt$Status <- rep(NA, nrow(tt))
tt$Status <- ifelse(rownames(tt) %in% nat, "Native", tt$Status)
tt$Status <- ifelse(rownames(tt) %in% ext, "Extirpated", tt$Status)
tt$Status <- ifelse(rownames(tt) %in% int, "Introduced", tt$Status)
tt$Status <- ifelse(rownames(tt) == "Chapalichthys encaustus", "Native", tt$Status) 
# Chapalichthys encaustus counted as native remaining instead of extirpated since it seems to be translocated in one site
sum(is.na(tt$Status))
sum(tt$Status=="Native") #49
sum(tt$Status=="Extirpated") #29

# PCA:-------------------------------------------------------------------------------------
PCA <- prcomp(tt[,c(1:10)])
coords <- PCA$x
groups <- as.factor(tt$Status)
save(coords, file="coords.RData")

#dist_mat2 <- as.matrix(cluster::daisy(coords, metric = "euclidean"))
#dist_mat2 <- (dist_mat2 - min(dist_mat2)) / (max(dist_mat2) - min(dist_mat2))
#dist1 and dist2 are the same, OK

# Facto extra (exploration of trait space):------------------------------------------------
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
                   title = "Variables-PCA 1-2",
                   col.var = "contrib", # Color by contributions to the PC
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE     # Avoid text overlapping
))
(vII <- fviz_pca_var(PCA,
                    title = "Variables-PCA 3-4",
                    axes = c(3,4),
                    col.var = "contrib", # Color by contributions to the PC
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE     # Avoid text overlapping
))

# Biplots:
(bI <- fviz_pca_biplot(PCA, repel = TRUE,
                      title = "PCA-Biplot 1-2",
                      col.var = "contrib", # Variables color
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      col.ind = "#696969",  # Individuals color
                      label = "var",
                      addEllipses = TRUE,
                      ellipse.alpha = 0.1,
                      ellipse.type = "convex"
))
(bII <- fviz_pca_biplot(PCA,
                       title = "PCA-Biplot 3-4",
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
eig.val 

# NOTES:
# Dim 1 (25.865711) and Dim 2 (20.990517) = 46.87 (~46.9)
# Reach > 80% var explained with first five components

# Results for Variables
PCAv <- get_pca_var(PCA)
PCAv$coord          # Coordinates
PCAv$contrib        # Contributions to the PCs
PCAv$cos2           # Quality of representation 

# NOTES:
fviz_contrib(PCA, choice = "var", axes = 1, top = 10) 
# Dim 1: RES, PFv, BEl,  OGp (Visual acuity, use of fin for swimming, hydrodinamism and oral gape position with some effect of body size)
# Predatory/ non predatory division?

fviz_contrib(PCA, choice = "var", axes = 2, top = 10) 
# Dim 2: RMl, VEp, BLs (Strength of jaw, vertical eye position and shape of body)
# Could also be another axis of predatory/ non predatory behaviour

fviz_contrib(PCA, choice = "var", axes = 3, top = 10)
# Dim 3:PFs, VEp, OGp

# Feeding strategy in the water column
fviz_contrib(PCA, choice = "var", axes = 4, top = 10)
# Dim 4:CPt(mainly), BLs and MBl

# Results for individuals
PCAi <- get_pca_ind(PCA)
PCAi$coord          # Coordinates
PCAi$contrib        # Contributions to the PCs
PCAi$cos2           # Quality of representation 

# PLot Individuals (Status):
(iIII <- fviz_pca_ind(PCA,
                      title = "Individuals-PCA 1-2 (by Status in 2005)",
                      label = "none", # hide individual labels
                      habillage = tt$Status, # color by groups
                      palette = c("#0072B2", "#F0E442", "Darkgray"),
                      addEllipses = TRUE, # TRUE for concentration ellipses
                      ellipse.type = "convex"
))


(iIV <- fviz_pca_ind(PCA,
                     title = "Individuals-PCA 3-4 (by Status in 2005)",
                     axes = c(3,4),
                     label = "none", # hide individual labels
                     habillage = tt$Status, # color by groups
                     palette = c("#0072B2", "#F0E442", "Darkgray"),
                     addEllipses = TRUE, # TRUE for concentration ellipses
                     ellipse.type = "convex"
))

# PERMANOVAs??

# Save plots: -----------------------------------------------------------------------------
axes12 <- grid.arrange(vI, iIII, ncol=2)
axes34 <- grid.arrange(vII, iIV, ncol=2)
ggsave(axes12, file= paste0(plot_dir, "/biplot1-2.jpg"), width = 15, height = 7) 
ggsave(axes34, file= paste0(plot_dir, "/biplot3-4.jpg"), width = 15, height = 7) 

###########################################################################################
# End of script ###########################################################################


