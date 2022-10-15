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
library(fishtree)
library(vegan)
library(pairwiseAdonis)
library(ggpubr)
library(grid)
library(ggforce)

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
#  mat <- as.matrix(mat)
#  n <- ncol(mat)
#  p.mat<- matrix(NA, n, n)
#  diag(p.mat) <- 0
#  for (i in 1:(n - 1)) {
#    for (j in (i + 1):n) {
#      tmp <- cor.test(mat[, i], mat[, j], ...)
#      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
#    }
#  }
#  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
#  p.mat
#} 

C <- cor(tt)
p.mat <- cor.mtest(tt)

file_path <- paste0(plot_dir, "/SM/Ms/Correlationmatrix.png")
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

# NOTE: Quality increases consistently with increasing number of summary traits.

# Classifying species in exotics/natives_extirpated/natives_remaining (REGIONAL):----------
nat <- sort(unique(names(ContN[,-c(1:2)]))) # 48 
nat[nat=="Agonostomus monticola"] <-"Dajaus monticola"
ext <- sort(unique(setdiff(names(HNB), names(ContN)))) # 30
int <- sort(unique(setdiff(names(ContE), names(HNB)))) # 22
48+22+30 # 100, OK

# Any species translocated?
intersect(names(ContE), names(HNB))
"Poecilia mexicana" %in% names(HNB)
# In NDT67: only as exotic in manantial puerta del río

#View(NDT67) [these are in the exotics column]
# "Chapalichthys encaustus"
# "Goodea atripinnis"
# "Poecilia butleri"
# "Poecilia sphenops"
# In paper, Poecilia mexicana & Chirostoma chapalae also mentioned.

intersect(names(ContE), names(HNB)) %in% nat
intersect(names(ContE), names(HNB)) %in% ext 
# Chapalichthys encaustus translocated in one site but extipated from its native range

tt$Status <- rep(NA, nrow(tt))
tt$Status <- ifelse(rownames(tt) %in% nat, "Native Remaining", tt$Status)
tt$Status <- ifelse(rownames(tt) %in% ext, "Extirpated", tt$Status)
tt$Status <- ifelse(rownames(tt) %in% int, "Introduced", tt$Status)
tt$Status <- ifelse(rownames(tt) == "Chapalichthys encaustus", "Native Remaining", tt$Status) 
# Chapalichthys encaustus counted as native remaining instead of extirpated since it seems to be translocated in one site
sum(is.na(tt$Status))
sum(tt$Status=="Native Remaining") #49
sum(tt$Status=="Extirpated") #29

# Classifying introduced species according to source:--------------------------------------
# int (= exotic component)
# Sources from Gesundheit & Macias Garcia, 2018.

# 1: Aquaculture
# 2: Sportfishing
# 3: Aquarium trade
# 4: Contaminants in aquaculture or other source

tt$Source <- rep(NA, nrow(tt))
tt$Source <- ifelse(rownames(tt) %in% c("Oreochromis sp", "Oreochromis mossambicus",
                                        "Oreochromis aureus", "Oreochromis niloticus"), "1", tt$Source)
tt$Source <- ifelse(rownames(tt) %in% c("Cyprinus carpio", "Lepomis macrochirus",
                                        "Micropterus salmoides"), "1/2", tt$Source)
tt$Source <- ifelse(rownames(tt) %in% c("Pomoxis nigromaculatus"), "2", tt$Source)
tt$Source <- ifelse(rownames(tt) %in% c("Amatitlania nigrofasciata",
                                        "Astatotilapia burtoni"), "3", tt$Source)
tt$Source <- ifelse(rownames(tt) %in% c("Pseudoxiphophorus bimaculatus",
                                        "Pseudoxiphophorus jonesii",
                                        "Gambusia yucatana",
                                        "Poeciliopsis gracilis",
                                        "Pseudoxiphophorus sp",
                                        "Poeciliopsis sp",
                                        "Poecilia mexicana"), "4", tt$Source)
tt$Source <- ifelse(rownames(tt) %in% c("Xiphophorus hellerii",
                                        "Poecilia reticulata",
                                        "Xiphophorus maculatus",
                                        "Xiphophorus variatus",
                                        "Poecilia sp"), "3/4", tt$Source)
# Poeciliopsis sp" /"Poecilia sp" /"Pseudoxiphophorus sp"  added following sps classification
sum(is.na(tt$Source))
tt$Source <- ifelse(is.na(tt$Source), tt$Status, tt$Source)

tt$SourceII <- rep(NA, nrow(tt))
tt$SourceII <- ifelse(tt$Source %in% c("1", "2", "1/2"), "IA", tt$SourceII)
tt$SourceII <- ifelse(tt$Source %in% c("3", "4", "3/4"), "IB", tt$SourceII)
sum(is.na(tt$SourceII))
tt$SourceII <- ifelse(is.na(tt$SourceII), tt$Status, tt$SourceII)
tt$SourceII[tt$SourceII=="Extirpated"] <- "E"
tt$SourceII[tt$SourceII=="Native Remaining"] <- "NR"
tt2 <- tt
save(tt2, file="tt2.RData")
###########################################################################################
# PCA:-------------------------------------------------------------------------------------
PCA <- prcomp(tt[,c(1:10)])
coords <- PCA$x
groups <- as.factor(tt$Status)
save(coords, file="coords.RData")

#dist_mat2 <- as.matrix(cluster::daisy(coords, metric = "euclidean"))
#dist_mat2 <- (dist_mat2 - min(dist_mat2)) / (max(dist_mat2) - min(dist_mat2))
#dist1 and dist2 are the same, OK

###########################################################################################
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
                    ylim=c(-1,1),
                    xlim=c(-1,1),
                    title = "Variables-PCA 1-2",
                    col.var = "contrib", # Color by contributions to the PC
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE     # Avoid text overlapping
))
(vI <- vI + geom_circle(aes(x0=0, y0=0, r = 1), col="grey70"))
(vII <- fviz_pca_var(PCA,
                     ylim=c(-1,1),
                     xlim=c(-1,1),
                     title = "Variables-PCA 3-4",
                     axes = c(3,4),
                     col.var = "contrib", # Color by contributions to the PC
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                     col.circle = "grey70",
                     repel = TRUE     # Avoid text overlapping
))
(vII <- vII + geom_circle(aes(x0=0, y0=0, r = 1), col="grey70"))
(panels_v <- ggarrange(vI,
                       vII,
                       common.legend = TRUE, 
                       legend="bottom",
                       align="hv",
                       font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))

ggsave(panels_v, filename = paste0(plot_dir, "/SM/Ms/panels_v.jpg"), width=9, height=5)

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

(panels_b <- ggarrange(bI,
                       bII,
                       common.legend = TRUE, 
                       legend="bottom",
                       align="hv",
                       font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))

ggsave(panels_b, filename = paste0(plot_dir, "/SM/Ms/panels_b.jpg"), width=9, height=5)

# Eigenvalues:
eig.val <- get_eigenvalue(PCA)
eig.val 

# NOTES:
# Dim 1 (25.865711) and Dim 2 (20.990517) = 46.87 (~46.9)
# Reach > 70% variance explained with first four components

# Results for Variables
PCAv <- get_pca_var(PCA)
PCAv$coord          # Coordinates
PCAv$contrib        # Contributions to the PCs
PCAv$cos2           # Quality of representation 

# NOTES:
c1 <- fviz_contrib(PCA, choice = "var", axes = 1, top = 10) 
# Dim 1: RES, PFv, BEl,  OGp (Visual acuity, use of fin for swimming, hydrodinamism and oral gape position)

c2 <- fviz_contrib(PCA, choice = "var", axes = 2, top = 10) 
# Dim 2: RMl, VEp, BLs (Strength of jaw, vertical eye position and shape of body)

c3 <- fviz_contrib(PCA, choice = "var", axes = 3, top = 10)
# Dim 3:PFs, VEp, OGp

# Feeding strategy in the water column
c4 <- fviz_contrib(PCA, choice = "var", axes = 4, top = 10)
# Dim 4:CPt(mainly), BLs and MBl

contrib <- grid.arrange(c1, c2, c3, c4)
ggsave(contrib, filename=paste0(plot_dir, "/SM/Ms/ContributionsTraits.jpg"), width=10, height = 8)

# Results for individuals
PCAi <- get_pca_ind(PCA)
PCAi$coord          # Coordinates
PCAi$contrib        # Contributions to the PCs
PCAi$cos2           # Quality of representation 

###########################################################################################
# Plots to save : -------------------------------------------------------------------------
(Main <- fviz_pca_ind(PCA,
                      title = "",
                      label = "none", # hide individual labels
                      habillage = tt$SourceII, # color by groups
                      palette = c("#0072B2","#D55E00","#E69F00","darkgray"),
                      #palette = c("#0072B2", "#F0E442", "Darkgray"),
                      addEllipses = TRUE, # TRUE for concentration ellipses
                      legend.title="Status in 2005",
                      ellipse.type = "convex"
))

(Main_biplotB <- fviz_pca_biplot(PCA, repel = TRUE,
                                title = "PCA-Biplot 1-2",
                                col.var = "grey31", # Variables color
                                habillage = tt$SourceII, # color by groups
                                palette = c("#0072B2","#D55E00","#E69F00","darkgray"),
                                label = "var",
                                addEllipses = TRUE,
                                ellipse.alpha = 0.2,
                                ellipse.type = "convex",
                                legend.title = "Status"
))

(Main_biplotB34 <- fviz_pca_biplot(PCA, repel = TRUE,
                                 title = "PCA-Biplot 3-4",
                                 axes = c(3,4),
                                 col.var = "grey31", # Variables color
                                 habillage = tt$SourceII, # color by groups
                                 palette = c("#0072B2","#D55E00","#E69F00","darkgray"),
                                 label = "var",
                                 addEllipses = TRUE,
                                 ellipse.alpha = 0.2,
                                 ellipse.type = "convex",
                                 legend.title = "Status"
))

(panels_bSource <- ggarrange(Main_biplotB,
                             Main_biplotB34,
                             common.legend = TRUE, 
                             legend="bottom",
                             align="hv",
                             font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))


# Save plots: -----------------------------------------------------------------------------
ggsave(Main, filename = paste0(plot_dir, "/SM/Ms/TraitSpacePanel.jpg"), width=8, height=6)
ggsave(panels_bSource, filename = paste0(plot_dir, "/SM/Ms/panels_bSource.jpg"), width=10, height=6)
save(Main, file="Main.RData")
###########################################################################################
# End of script ###########################################################################


