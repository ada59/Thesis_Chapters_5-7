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
library(taxize)
library(vegan)
library(pairwiseAdonis)


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
# custom function (http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram#data-for-correlation-analysis)
# Other : https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

C <- cor(tt)
p.mat <- cor.mtest(tt)

file_path <- paste0(plot_dir, "/SM/Correlation matrix.png")
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
# In NDT67:
# only as exotic in manantial puerta del río

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
tt$Status <- ifelse(rownames(tt) %in% nat, "Native", tt$Status)
tt$Status <- ifelse(rownames(tt) %in% ext, "Extirpated", tt$Status)
tt$Status <- ifelse(rownames(tt) %in% int, "Introduced", tt$Status)
tt$Status <- ifelse(rownames(tt) == "Chapalichthys encaustus", "Native", tt$Status) 
# Chapalichthys encaustus counted as native remaining instead of extirpated since it seems to be translocated in one site
sum(is.na(tt$Status))
sum(tt$Status=="Native") #49
sum(tt$Status=="Extirpated") #29

# Classifying introduced species according to source:--------------------------------------
# int (= exotic component)
# Sources from Gesundheit & Macias Garcia, 2018.

# 1 : Aquaculture
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
tt$SourceII <- ifelse(tt$Source %in% c("1", "2", "1/2"), "Aquaculture&Sportfishing", tt$SourceII)
tt$SourceII <- ifelse(tt$Source %in% c("3", "4", "3/4"), "Aquarium&Contaminant", tt$SourceII)
sum(is.na(tt$SourceII))
tt$SourceII <- ifelse(is.na(tt$SourceII), tt$Status, tt$SourceII)

# Classifying species according to their phylogeny:----------------------------------------
example <- c("Poecilia sphenops", "Cyprinus carpio")
exampled_df <- tax_name(example, get = c('family', "genus"), db = 'ncbi')

tt$Family <- rep(NA, nrow(tt))
tt$Family <- tax_name(rownames(tt), get = "family", db = "ncbi")$family

sum(is.na(tt$Family)) # 14
rownames(tt)[is.na(tt$Family)]

tt$Family <- ifelse(rownames(tt) %in% c("Chirostoma aculeatum",
                                        "Chirostoma charari",
                                        "Chirostoma sp"), "Atherinopsidae", tt$Family) # FishBase
tt$Family <- ifelse(rownames(tt) %in% c("Sicydium multipunctatum"), "Gobiidae", tt$Family) # FishBase
tt$Family <- ifelse(rownames(tt) %in% c("Gila sp", "Tampichthys dichroma", "Algansea popoche"), "Leuciscidae", tt$Family) # FishBase
tt$Family <- ifelse(rownames(tt) %in% c("Gambusia senilis",
                                        "Poecilia sp",
                                        "Poeciliopsis sp",
                                        "Pseudoxiphophorus sp"), "Poeciliidae", tt$Family) # FishBase
tt$Family <- ifelse(rownames(tt) %in% c("Oreochromis sp", "Nosferatu labridens"), "Cichlidae", tt$Family)
tt$Family <- ifelse(rownames(tt) %in% c("Astyanax sp"), "Characidae", tt$Family)
sum(is.na(tt$Family))

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

ggsave(vI, filename = paste0(plot_dir, "/SM/FactoExtra/variables_12.jpg"), width=8, height=8)
ggsave(vII, filename = paste0(plot_dir, "/SM/FactoExtra/variables_34.jpg"), width=8, height=8)

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

ggsave(bI, filename = paste0(plot_dir, "/SM/FactoExtra/biplot_12.jpg"), width=8, height=8)
ggsave(bII, filename = paste0(plot_dir, "/SM/FactoExtra/biplot_34.jpg"), width=8, height=8)


# Eigenvalues:
eig.val <- get_eigenvalue(PCA)
eig.val 

# NOTES:
# Dim 1 (25.865711) and Dim 2 (20.990517) = 46.87 (~46.9)
# Reach > 70% variance explained with first four components
# Reach > 80% var explained with first five components

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
ggsave(contrib, filename=paste0(plot_dir, "/SM/FactoExtra/ContributionsTraits.jpg"), width=10, height = 8)

# Results for individuals
PCAi <- get_pca_ind(PCA)
PCAi$coord          # Coordinates
PCAi$contrib        # Contributions to the PCs
PCAi$cos2           # Quality of representation 

###########################################################################################
# Plots to save : -------------------------------------------------------------------------
(iA <- fviz_pca_ind(PCA,
                      title = "Individuals-PCA 1-2 (by Status in 2005)",
                      label = "none", # hide individual labels
                      habillage = tt$Status, # color by groups
                      palette = c("#0072B2", "#F0E442", "Darkgray"),
                      addEllipses = TRUE, # TRUE for concentration ellipses
                      ellipse.type = "convex"
))
(iB <- fviz_pca_ind(PCA,
                      title = "Individuals-PCA 1-2 (by Source in 2005)",
                      label = "none", # hide individual labels
                      habillage = tt$SourceII, # color by groups
                      #palette = c("#0072B2", "#F0E442", "Darkgray"),
                      addEllipses = TRUE, # TRUE for concentration ellipses
                      ellipse.type = "convex"
))
(iC <- fviz_pca_ind(PCA,
                    title = "Individuals-PCA 1-2 (by Family in 2005)",
                    label = "none", # hide individual labels
                    habillage = tt$Family, # color by groups
                    #palette = c("#0072B2", "#F0E442", "Darkgray"),
                    addEllipses = TRUE, # TRUE for concentration ellipses
                    ellipse.type = "convex"
))
(iD <- fviz_pca_ind(PCA,
                    title = "Individuals-PCA 1-2 (by SourceII in 2005)",
                    label = "none", # hide individual labels
                    habillage = tt$Source, # color by groups
                    #palette = c("#0072B2", "#F0E442", "Darkgray"),
                    addEllipses = TRUE, # TRUE for concentration ellipses
                    ellipse.type = "convex"
))

(i2A <- fviz_pca_ind(PCA,
                     title = "Individuals-PCA 3-4 (by Status in 2005)",
                     axes = c(3,4),
                     label = "none", # hide individual labels
                     habillage = tt$Status, # color by groups
                     palette = c("#0072B2", "#F0E442", "Darkgray"),
                     addEllipses = TRUE, # TRUE for concentration ellipses
                     ellipse.type = "convex"
))
(i2B <- fviz_pca_ind(PCA,
                     title = "Individuals-PCA 3-4 (by Source in 2005)",
                     axes = c(3,4),
                     label = "none", # hide individual labels
                     habillage = tt$SourceII, # color by groups
                     #palette = c("#0072B2", "#F0E442", "Darkgray"),
                     addEllipses = TRUE, # TRUE for concentration ellipses
                     ellipse.type = "convex"
))
(i2C <- fviz_pca_ind(PCA,
                     title = "Individuals-PCA 3-4 (by Family in 2005)",
                     axes = c(3,4),
                     label = "none", # hide individual labels
                     habillage = tt$Family, # color by groups
                     #palette = c("#0072B2", "#F0E442", "Darkgray"),
                     addEllipses = TRUE, # TRUE for concentration ellipses
                     ellipse.type = "convex"
))
(i2D <- fviz_pca_ind(PCA,
                     title = "Individuals-PCA 3-4 (by SourceII in 2005)",
                     axes = c(3,4),
                     label = "none", # hide individual labels
                     habillage = tt$Source, # color by groups
                     #palette = c("#0072B2", "#F0E442", "Darkgray"),
                     addEllipses = TRUE, # TRUE for concentration ellipses
                     ellipse.type = "convex"
))

# Save plots: -----------------------------------------------------------------------------
groups12 <- grid.arrange(iA, iB, ncol=2)
ggsave(groups12, filename= paste0(plot_dir, "/SM/FactoExtra/habillageStatus&Source12.jpg"), width = 16, height = 8) 

groups34 <- grid.arrange(i2A, i2B, ncol=2)
ggsave(groups34, filename= paste0(plot_dir, "/SM/FactoExtra/habillageStatus&Source34.jpg"), width = 16, height = 8) 

ggsave(iC, filename= paste0(plot_dir, "/SM/FactoExtra/habillageFamily12.jpg"), width = 8, height = 8) 

ggsave(i2C, filename= paste0(plot_dir, "/SM/FactoExtra/habillageFamily34.jpg"), width = 8, height = 8) 

###########################################################################################
# PERMANOVA:-------------------------------------------------------------------------------
# Used to compare groups of objects and test the null hypothesis 
# that the centroids and dispersion of the groups as defined 
# by measure space are equivalent for all groups.

identical(colnames(dist_mat1), rownames(tt)) #TRUE
grouping <- as.factor(tt$Status)

adonis(dist_mat1~grouping) # adonis 2?
adonis(coords~grouping, method = "euclidean")

# R-squared of 0.089*100= 8.9%
pairwise.adonis(coords, grouping, sim.function='vegdist',sim.method='euclidian')


grouping2 <- as.factor(tt$SourceII)

adonis(dist_mat1~grouping2) # adonis 2?
adonis(coords~grouping2, method = "euclidean")

# R-squared of 0.17439*100= 17.43%
pairwise.adonis(coords, grouping2, sim.function='vegdist',sim.method='euclidian')
# OK, so overall we see that the group Aquaculture&Sportfishing is the most
# different in composition. (But sample size)

# ANOVA:------------------------------------------------------------------------------------
coords2 <- data.frame(coords, tt$Status, tt$SourceII)
coords2$tt.Status <- as.factor(coords2$tt.Status)
coords2$tt.SourceII <- as.factor(coords2$tt.SourceII)

mod1_dim1 <- lm(PC1~tt$Status, data=coords2)
mod2_dim1 <- aov(PC1~tt$Status, data=coords2)
summary(mod1_dim1)
summary(mod2_dim1)
anova(mod1_dim1)
anova(mod2_dim1)
TukeyHSD(mod2_dim1)

mod1_dim2 <- lm(PC2~tt$Status, data=coords2)
mod2_dim2 <- aov(PC2~tt$Status, data=coords2)
summary(mod1_dim2)
summary(mod2_dim2)
anova(mod1_dim2)
anova(mod2_dim2)
TukeyHSD(mod2_dim2)

mod1_dim3 <- lm(PC3~tt$Status, data=coords2)
mod2_dim3 <- aov(PC3~tt$Status, data=coords2)
summary(mod1_dim3) #  0.1564 
summary(mod2_dim3)
anova(mod1_dim3)
anova(mod2_dim3)
TukeyHSD(mod2_dim3)

dim4 <- lm(PC4~tt$Status, data=coords2)
plot(dim1)
summary(dim4)

###########################################################################################
# End of script ###########################################################################


