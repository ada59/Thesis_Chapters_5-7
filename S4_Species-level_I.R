#===============================================================================
# Script: Trait distance matrices
# AFE
# August 2022
#===============================================================================




#===============================================================================
# Libraries:--------------------------------------------------------------------
#===============================================================================
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
library(ggpubr)
library(grid)
library(ggforce)
library(taxize)

rm(list=ls())
myd <- getwd()
path_lists <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists"   # Raw data lists
path_traits <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Traits" # Dir to trait objects


path_plots5 <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots/Chapter5"   # Dir to save plots Chapter 5
path_plots6 <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots/Chapter6"   # Dir to save plots Chapter 6




#===============================================================================
# Read data:--------------------------------------------------------------------
#===============================================================================
load(paste0(path_lists,"/l67.RData"))     # main list with community data matrices (67 sites)
load(paste0(path_lists, "/l53.RData"))    # main list with community data matrices (53 sites)
load(paste0(path_lists, "/NDT67.RData"))  # List for checks on raw data
load(paste0(path_lists, "/NDT83.RData"))  # List for checks on raw data


load(paste0(path_traits,"/fishtrait_complete_imp.RData")) # trait dataset with no missing values
tt <- fishtrait_complete_imp
class(tt)

head(tt)
hist(tt[,1])  # checked all. MBl only non-normally distributed
range(tt[,1]) #  4.14 100.00




#===============================================================================
# Correlations : ---------------------------------------------------------------
#===============================================================================
C1 <- cor(tt, method = "spearman")

file_path <- paste0(path_plots5, "/S4_PearsonCorrelationmatrix_NoMBlLog.png")
png(width = 700, height = 700, file=file_path)
corrplot(C1, method="number", type = "lower", order = "hclust", number.cex = 1.5, tl.cex = 1.5)
dev.off()



tt[,1] <- log(tt[,1])
colnames(tt)[colnames(tt)=="MBl"] <- "logMBl" # needed for PCA
C2 <- cor(tt) # Pearson

file_path <- paste0(path_plots6, "/S4_PearsonCorrelationmatrix.png")
png(width = 700, height = 700, file=file_path)
corrplot(C2, method="number", type = "lower", order = "hclust", number.cex = 1.5, tl.cex = 1.5)
dev.off()




#===============================================================================
# Distance Matrix: -------------------------------------------------------------
#===============================================================================
round(diag(cov(tt)), 2) # need to scale variances.
tt <- as.data.frame(tt)
tt[] <- data.frame(apply(tt, 2, scale)) # defaults are TRUE center & scale
round(diag(cov(tt)), 2) # equal to 1
save(tt, file=paste0(path_traits, "/tt.RData"))

dist_mat1 <- as.matrix(daisy(tt, metric = "euclidean"))
dist_mat1 <- (dist_mat1 - min(dist_mat1)) / (max(dist_mat1) - min(dist_mat1)) # re-scale between 0 & 1
range(dist_mat1) # OK
save(dist_mat1, file=paste0(path_traits, "/dist_mat1.RData"))




#===============================================================================
# Assign REGIONAL status groups: -----------------------------------------------
#===============================================================================
HNC <- subset(l67, l67$Period == "HNC")
HNC <- HNC[, colSums(HNC != 0) > 0]        # 77

HNB <- subset(l67, l67$Period == "HNB")
HNB <- HNB[, colSums(HNB != 0) > 0]        # 80

ContN <- subset(l67, l67$Period == "ContN")
ContN <- ContN[, colSums(ContN != 0) > 0]  # 67

ContE <- subset(l67, l67$Period == "ContE")
ContE <- ContE[, colSums(ContE != 0) > 0]  # 37

nat <- sort(unique(names(ContN[,-c(1:3)]))) # 47 natives remaining in the area

ext_conservative <- sort(unique(setdiff(names(HNC), names(ContN)))) # 29 extirpated natives (conservative)
ext_conservative

ext <- sort(unique(setdiff(names(HNB), names(ContN)))) # 30 extirpated (broad)
ext

setdiff(ext, ext_conservative) # A zacapuensis (OK, checked case in NDT67)

int <- sort(unique(setdiff(names(ContE), names(HNB)))) # 18
47+18+30 # 95, OK


## Translocated in some sites: -------------------------------------------------
intersect(names(ContE), names(ContN))
intersect(names(ContE), nat) # idem ContN but without Site, Basin & Period
# "Goodea atripinnis"
# "Poecilia butleri"
# "Poecilia sphenops"
# The three above are native in some sites, and introduced in others (some translocations)

intersect(names(ContE), names(HNB)) %in% nat
intersect(names(ContE), names(HNB)) %in% ext 

intersect(names(ContE), names(HNB))
# Previous 3 +
# Chapalichthys encaustus (extirpated in native locs, but introduced in Ameca River at La Vega)


## Particular cases: -----------------------------------------------------------
#"Poecilia mexicana" %in% names(HNB) # F
# Poecilia mexicana in NDT83 as native in the Panuco River
# In NDT67: only as exotic in manantial puerta del río
# C. chapalae: in NDT67, only as native in Grande de Santiago River
# C. chapalae: in NDT83, as translocated invLake Cuitzeo



# ORIGIN INTRODUCED: -----------------------------------------------------------
# C. carpio (Europe & Asia)

# O aureus, O. mossambicus, O. niloticus, O sp (Africa)
# A burtoni (Africa)

# Lepomis macrochirus (North America)
# Micropterus salmoides (North America)
# Pomoxis nigromaculatus (North America)
# Herichthys cyanoguttatus (North Mexico) [removed, see email thread "Respuestas a dudas" 30/11/2024]

# A. nigrofasciata (Central America)
# Gambusia yucatana (Central America, tropical climate provinces in south of Mexico)
# Poeciliopsis gracilis: (South of Mexico)
# Poeciliopsis sp (proxy P gracilis): (South of Mexico)

# Poecilia reticulata (South America)

# Pseudoxiphophorus bimaculatus: (Pacific coast Mexico)
# Pseudoxiphophorus jonesii: (Pacific coast Mexico)
# Pseudoxiphophorus sp (proxy P bimaculatus): (Pacific coast Mexico)
# Xiphophorus hellerii: (Pacific coast Mexico)
# Xiphophorus maculatus: (Pacific coast Mexico)
# Xiphophorus variatus: (Pacific coast Mexico)

# Poecilia sp (proxy P. sphenops): (overlap with area sampled)
# P. mexicana (overlap with area sampled)
# C. encaustus (overlap with area sampled)
# G. atripinnis, P. butleri, & P. cf sphenops (overlap with area sampled)


tt$Status <- rep(NA, nrow(tt))
tt$Status <- ifelse(rownames(tt) %in% nat, "Native Remaining", tt$Status)
tt$Status <- ifelse(rownames(tt) %in% ext, "Extirpated", tt$Status)
tt$Status <- ifelse(rownames(tt) %in% int, "Introduced", tt$Status)
tt$Status <- ifelse(rownames(tt) == "Chapalichthys encaustus", "Native Remaining", tt$Status) 
# Chapalichthys encaustus counted as native remaining
# P butleri, G atripinnis & P cf sphenops corretly classified as Native Remaning 


sum(is.na(tt$Status))
sum(tt$Status=="Native Remaining") # 49
sum(tt$Status=="Extirpated")       # 28 (Historical Broad)
sum(tt$Status=="Introduced")       # 18

sort(unique(rownames(tt)[tt$Status=="Native Remaining"]))

## Any that would change category if considering group of 83?:------------------

# There are a few sps classified as Extirpated in 67, that were found in other
# locations in 83. These are:

# Chirostoma chapalae in NDT83 as translocated in P. Cointzio
# Herychtys labridens present as native in NDT83 in Laguna de la Media Luna
# P viriosa found as native in NDT83 in rio en 6 de Enero
# Notropis boucardi found as native in NDT83 in rio Nexapa


# There is also a species that in the 67 susbet would be classified as 
# introduced that is native of a location in the 83 subset.

# Poecilia mexicana in NDT83 as native in the Panuco River


# NOTE #########################################################################
# Consider the above in a sensitivity analysis

tt$StatusSen <- tt$Status
tt$StatusSen[rownames(tt) %in% c("Chirostoma chapalae",
                                 "Herichthys labridens", "Poeciliopsis viriosa", 
                                 "Notropis boucardi", "Poecilia mexicana")] <- "Native Remaining"
 
sum(is.na(tt$StatusSen))
sum(tt$StatusSen=="Native Remaining") # 54 (49 + 5)
sum(tt$StatusSen=="Extirpated")       # 24 (Historical Broad)
sum(tt$StatusSen=="Introduced")       # 17, OK
################################################################################


## Classifying introduced species according to source:---------------------------
# int (= exotic component)
# Sources from Gesundheit & Macias Garcia, 2018.
# Removed "sp" except for Gila after addendum!

# 1: Aquaculture
# 2: Sportfishing
# 3: Aquarium trade
# 4: Contaminants in aquaculture or other source

tt$Source <- rep(NA, nrow(tt))
tt$Source <- ifelse(rownames(tt) %in% c("Oreochromis mossambicus",
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
                                        "Poecilia mexicana"), "4", tt$Source)
tt$Source <- ifelse(rownames(tt) %in% c("Xiphophorus hellerii",
                                        "Poecilia reticulata",
                                        "Xiphophorus maculatus",
                                        "Xiphophorus variatus"), "3/4", tt$Source)

sum(is.na(tt$Source))
tt$Source <- ifelse(is.na(tt$Source), tt$Status, tt$Source)

tt$SourceII <- rep(NA, nrow(tt))
tt$SourceII <- ifelse(tt$Source %in% c("1", "2", "1/2"), "IA", tt$SourceII)
tt$SourceII <- ifelse(tt$Source %in% c("3", "4", "3/4"), "IB", tt$SourceII)
sum(is.na(tt$SourceII))
tt$SourceII <- ifelse(is.na(tt$SourceII), tt$Status, tt$SourceII)
tt$SourceII[tt$SourceII=="Extirpated"] <- "E"
tt$SourceII[tt$SourceII=="Native Remaining"] <- "NR"
table(tt$SourceII)


## SourceIISen: ----------------------------------------------------------------
tt$SourceIISen <- tt$SourceII
tt$SourceIISen[rownames(tt) %in% c("Chirostoma chapalae",
                                   "Herichthys labridens", "Poeciliopsis viriosa", 
                                   "Notropis boucardi", "Poecilia mexicana")] <- "NR"


## SourceBasic: ----------------------------------------------------------------
tt$SourceIIBasic <- tt$SourceII
tt$SourceIIBasic[tt$SourceIIBasic %in% c("IA", "IB")] <- "I"
table(tt$SourceIIBasic)


## SourceBasicSen: -------------------------------------------------------------
tt$SourceIIBasicSen <- tt$SourceIIBasic
tt$SourceIIBasicSen[rownames(tt) %in% c("Chirostoma chapalae",
                                        "Herichthys labridens", "Poeciliopsis viriosa", 
                                        "Notropis boucardi", "Poecilia mexicana")] <- "NR"
table(tt$SourceIIBasicSen)

tt2 <- tt
save(tt2, file=paste0(path_traits, "/tt2.RData"))




#===============================================================================
# Quality of trait space:-------------------------------------------------------
#===============================================================================
# https://cmlmagneville.github.io/mFD/articles/mFD_general_workflow.html
# https://cmlmagneville.github.io/mFD/articles/Continuous_traits_framework.html #compute-the-functional-space
#qfsp <- quality.fspaces(as.dist(dist_mat1), deviation_weighting = c("absolute", "squared"), maxdim_pcoa = 10)
#round(qfsp$quality_fspaces, 3) # quality increases consistently with increasing number of summary traits.

#qualityspace <- quality.fspaces.plot(
#  fspaces_quality            = qfsp,
#  quality_metric             = "mad",
#  fspaces_plot               = row.names(qfsp$quality_fspaces),
#  name_file                  = NULL,
#  range_dist                 = NULL,
#  range_dev                  = NULL,
#  range_qdev                 = NULL,
#  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
#  gradient_deviation_quality = c(low = "yellow", high = "red"),
#  x_lab                      = "Trait-based distance")

#ggsave(qualityspace, filename = paste0(path_plots5, "/S4_qualityspace.jpg"), width=15, height=15)

#coords_qfsp <- qfsp$details_fspaces$sp_pc_coord #(same output of prcomp)

# Above: OLD, mFD has been updated. For all-continuous traits, now use: 
# function tr.cont.fspace()

fspace <- tr.cont.fspace(
  sp_tr        = tt[, c(1:10)], 
  pca          = TRUE, 
  nb_dim       = 10, 
  scaling      = "scale_center",
  compute_corr = "pearson") # ok, identical to prcomp (qfsp performs PCoA, so not identical)
#View(fspace$sp_faxes_coord)
round(fspace$quality_metrics, 3)




#===============================================================================
# PCA:--------------------------------------------------------------------------
#===============================================================================
rownames(tt) <- paste0(substr(rownames(tt), 1, 1), ".", str_split_fixed(rownames(tt), " ", 2)[,2])


PCA <- prcomp(tt[,c(1:10)])
coords <- PCA$x
save(coords, file=paste0(path_traits,"/coords.RData"))


(fviz_eig(PCA)) # Elbow between 2 & 3




#===============================================================================
# Facto extra (exploration of trait space Chapter 5):---------------------------
#===============================================================================

## Individuals:-----------------------------------------------------------------
(iI <- fviz_pca_ind(PCA,
                    axes = c(1,2),
                    labelsize=3,
                    col.ind = "cos2", # Color by the quality of representation
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                    repel = TRUE))
# geom_text(aes(label = custom_labels), size = 3) [to control label size with ggplot instead & label custom within fviz_pca_ind()]
(iI <- iI + theme(
  axis.title = element_text(size = 14),  
  axis.text = element_text(size = 12),   
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12), 
  plot.title = element_text(size = 16),  
))
(iII <- fviz_pca_ind(PCA,
                     axes = c(3,4),
                     labelsize=3,
                     col.ind = "cos2", # Color by the quality of representation
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                     repel = TRUE     # Avoid text overlapping
))
(iII <- iII + theme(
  axis.title = element_text(size = 14),  
  axis.text = element_text(size = 12),   
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12), 
  plot.title = element_text(size = 16),  
))


ggsave(iI, filename = paste0(path_plots5, "/S4_2D_individuals.jpg"), width=9, height=9)
ggsave(iII, filename = paste0(path_plots5, "/S4_2D_34_individuals.jpg"), width=9, height=9)


## Variables:-------------------------------------------------------------------
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
                     col.var = "contrib", 
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                     col.circle = "grey70",
                     repel = TRUE     
))
(vII <- vII + geom_circle(aes(x0=0, y0=0, r = 1), col="grey70"))
(panels_v <- ggarrange(vI,
                       vII,
                       common.legend = TRUE, 
                       legend="bottom",
                       align="hv",
                       font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))

ggsave(panels_v, filename = paste0(path_plots5, "/S4_panels_variables.jpg"), width=9, height=5)
ggsave(vI, filename = paste0(path_plots5, "/S4_variables_1_2.jpg"), width=6, height=5)

## Biplots:---------------------------------------------------------------------
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
ggsave(panels_b, filename = paste0(path_plots5, "/S4_panels_biplots.jpg"), width=9, height=5)


## Eigenvalues:-----------------------------------------------------------------
eig.val <- get_eigenvalue(PCA)
eig.val # Dim 1 (25.8) and Dim 2 (23) = 48.8%


## Contributions:---------------------------------------------------------------
round(PCA$rotation*100)


(c1 <- fviz_contrib(PCA, choice = "var", axes = 1, top = 10)) 
(c2 <- fviz_contrib(PCA, choice = "var", axes = 2, top = 10)) 
(c3 <- fviz_contrib(PCA, choice = "var", axes = 3, top = 10))
(c4 <- fviz_contrib(PCA, choice = "var", axes = 4, top = 10))


contrib <- grid.arrange(c1, c2, c3, c4)
ggsave(contrib, filename=paste0(path_plots5, "/S4_panels_contributionsDim1-4.jpg"), width=10, height = 8)
ggsave(c1, filename = paste0(path_plots5, "/S4_contributions_1.jpg"), width=5, height=3)
ggsave(c2, filename = paste0(path_plots5, "/S4_contributions_2.jpg"), width=5, height=3)


## Variables & inviduals results: ----------------------------------------------
PCAv <- get_pca_var(PCA) # variable results
PCAv$coord               # coordinates
PCAv$contrib             # contributions to the PCs
PCAv$cos2                # quality of representation 

PCAi <- get_pca_ind(PCA) # individuals results
PCAi$coord               # Coordinates
PCAi$contrib             # Contributions to the PCs
PCAi$cos2                # Quality of representation 



## Family-grouping plot --------------------------------------------------------
taxa <- read.csv(file=paste0(path_lists, "/SM_Chapter5/taxa.csv"))
taxa$abbrev <- paste0(substr(taxa$Authority, 1, 1), ".", str_split_fixed(taxa$Authority, " ", 3)[,2])
taxa$abbrev[taxa$abbrev=="T.dichromus"] <- "T.dichroma"
taxa$abbrev[taxa$abbrev=="G.Baird"] <- "G.sp"
tt$Family <- taxa$Family[match(rownames(tt), taxa$abbrev)]

sum(is.na(tt$Family)) # 0

(iF <- fviz_pca_ind(PCA,
                    title = "Individuals-PCA 1-2 (by Family)",
                    label = "none", # hide individual labels
                    habillage = tt$Family, # color by groups
                    addEllipses = TRUE, # TRUE for concentration ellipses
                   ellipse.type = "convex"
))
(iF34 <- fviz_pca_ind(PCA,
                    title = "Individuals-PCA 1-2 (by Family)",
                    axes = c(3,4),
                    label = "none", 
                    habillage = tt$Family, 
                    addEllipses = TRUE, 
                    ellipse.type = "convex"
))

ggsave(iF, filename = paste0(path_plots5, "/S4_family_12.jpg"), width=6, height=5)
ggsave(iF34, filename = paste0(path_plots5, "/S4_family_34.jpg"), width=6, height=5)



#===============================================================================
# Chapter 6 Plots: -------------------------------------------------------------
#===============================================================================
(main_c6 <- fviz_pca_ind(PCA,
                        title = "",
                        label = "none", # hide individual labels
                        habillage = tt$SourceII, # color by groups
                        palette = c("#0072B2","#D55E00","#E69F00","darkgray"),
                        #palette = c("#0072B2", "#F0E442", "Darkgray"),
                        addEllipses = TRUE, # TRUE for concentration ellipses
                        legend.title="Status in 2005",
                        ellipse.type = "convex"
) + scale_shape_manual(values=c(19,19,19,19)))

(main_biplotc6 <- fviz_pca_biplot(PCA, repel = TRUE,
                                title = "PCA-Biplot 1-2",
                                col.var = "grey31", # Variables color
                                habillage = tt$SourceII, # color by groups
                                palette = c("#0072B2","#D55E00","#E69F00","darkgray"),
                                label = "var",
                                addEllipses = TRUE,
                                ellipse.alpha = 0.2,
                                ellipse.type = "convex",
                                legend.title = "Status"
)+ scale_shape_manual(values=c(19,19,19,19)))

(main_c6_34 <- fviz_pca_ind(PCA,
                         title = "",
                         axes = c(3,4),
                         label = "none", # hide individual labels
                         habillage = tt$SourceII, # color by groups
                         palette = c("#0072B2","#D55E00","#E69F00","darkgray"),
                         #palette = c("#0072B2", "#F0E442", "Darkgray"),
                         addEllipses = TRUE, # TRUE for concentration ellipses
                         legend.title="Status in 2005",
                         ellipse.type = "convex"
) + scale_shape_manual(values=c(19,19,19,19)))

(main_biplotc6_34 <- fviz_pca_biplot(PCA, repel = TRUE,
                                  title = "PCA-Biplot 1-2",
                                  axes = c(3,4),
                                  col.var = "grey31", # Variables color
                                  habillage = tt$SourceII, # color by groups
                                  palette = c("#0072B2","#D55E00","#E69F00","darkgray"),
                                  label = "var",
                                  addEllipses = TRUE,
                                  ellipse.alpha = 0.2,
                                  ellipse.type = "convex",
                                  legend.title = "Status"
)+ scale_shape_manual(values=c(19,19,19,19)))


ggsave(main_c6, filename = paste0(path_plots6, "/S4_TraitSpacePanel_12.jpg"), width=8, height=6)
ggsave(main_biplotc6, filename = paste0(path_plots6, "/S4_TraitSpacePanelVariables_12.jpg"), width=8, height=6)
ggsave(main_c6_34, filename = paste0(path_plots6, "/S4_TraitSpacePanel_34.jpg"), width=8, height=6)
ggsave(main_biplotc6_34, filename = paste0(path_plots6, "/S4_TraitSpacePanelVariables_34.jpg"), width=8, height=6)

save(main_c6, file="trait_space_12_status.RData")
save(main_biplotc6, file="trait_space_12_status_biplot.RData")




#===============================================================================
# Chapter 6 Sensitivity
#===============================================================================
path_plots_sen <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots/Chapter6/Sensitivity"

(main_c6_sen1 <- fviz_pca_ind(PCA,
                         title = "",
                         label = "none", # hide individual labels
                         habillage = tt$SourceIISen, # color by groups
                         palette = c("#0072B2","#D55E00","#E69F00","darkgray"),
                         #palette = c("#0072B2", "#F0E442", "Darkgray"),
                         addEllipses = TRUE, # TRUE for concentration ellipses
                         legend.title="Status in 2005",
                         ellipse.type = "convex"
) + scale_shape_manual(values=c(19,19,19,19)))

# H. labridens is one vertex 
# P mexicana, [TBC]

(main_c6_sen2 <- fviz_pca_ind(PCA,
                              title = "",
                              label = "none", # hide individual labels
                              habillage = tt$SourceIIBasic, # color by groups
                              palette = c("#0072B2","#D55E00","#E69F00","darkgray"),
                              #palette = c("#0072B2", "#F0E442", "Darkgray"),
                              addEllipses = TRUE, # TRUE for concentration ellipses
                              legend.title="Status in 2005",
                              ellipse.type = "convex"
) + scale_shape_manual(values=c(19,19,19,19)))
(main_c6_sen2 <- fviz_pca_ind(PCA,
                              title = "",
                              label = "none", # hide individual labels
                              habillage = tt$SourceIIBasicSen, # color by groups
                              palette = c("#0072B2","#D55E00","#E69F00","darkgray"),
                              #palette = c("#0072B2", "#F0E442", "Darkgray"),
                              addEllipses = TRUE, # TRUE for concentration ellipses
                              legend.title="Status in 2005",
                              ellipse.type = "convex"
) + scale_shape_manual(values=c(19,19,19,19)))



#===============================================================================
# Statistical differences Dim1 & Dim2:------------------------------------------
#===============================================================================

# In theory, KW H should be valid despite the largely uneven sample sizes.
dt_tests <- data.frame("Dim1"=PCA$x[,1], "Dim2"=PCA$x[,2], "SourceII"=tt$SourceII)
str(dt_tests)
dt_tests$SourceII <- as.factor(dt_tests$SourceII)


## AOV: ------------------------------------------------------------------------
Dim1fit <- aov(Dim1 ~ SourceII, data=dt_tests)
summary(Dim1fit)
shapiro.test(residuals(Dim1fit))

Dim2fit <- aov(Dim2 ~ SourceII, data=dt_tests)
summary(Dim2fit)
shapiro.test(residuals(Dim2fit)) 


## Kruskal: --------------------------------------------------------------------
Dim1fit_KW <- kruskal.test(Dim1 ~ SourceII, data=dt_tests)
Dim1fit_KW
pairwise.wilcox.test(dt_tests$Dim1,dt_tests$SourceII,
                     p.adjust.method = "BH") 
# IA with the rest of groups

Dim2fit_KW <- kruskal.test(Dim2 ~ SourceII, data=dt_tests)
Dim2fit_KW
pairwise.wilcox.test(dt_tests$Dim2,dt_tests$SourceII,
                     p.adjust.method = "BH") 
# Both IA and E are different from the other groups 


## OR: --------------------------------------------------------------------------
pairwise.wilcox.test(dt_tests$Dim1,dt_tests$SourceII,
                     p.adjust.method = "none") 
pairwise.wilcox.test(dt_tests$Dim2,dt_tests$SourceII,
                     p.adjust.method = "none") 

p2 <- c(2.3e-05, 0.5501, 0.3823, 0.0019, 7.2e-06, 0.9160,
        0.01228, 0.01378, 0.00252, 6.3e-05, 0.00045, 0.85537)
p.adjust(p2, method = "BH") # no difference!





#===============================================================================
# Statistical differences Dim1 & Dim2 (Sensitivity):----------------------------
#===============================================================================
dt_tests <- data.frame("Dim1"=PCA$x[,1], "Dim2"=PCA$x[,2], 
                           "SourceII"=tt$SourceIIBasic) # for SourceIISen (idem Source II)

## Kruskal: --------------------------------------------------------------------
Dim1fit_KW <- kruskal.test(Dim1 ~ SourceII, data=dt_tests)
Dim1fit_KW
pairwise.wilcox.test(dt_tests$Dim1,dt_tests$SourceII,
                     p.adjust.method = "BH") 
# Introduced & extirpated are different

Dim2fit_KW <- kruskal.test(Dim2 ~ SourceII, data=dt_tests)
Dim2fit_KW
pairwise.wilcox.test(dt_tests$Dim2,dt_tests$SourceII,
                     p.adjust.method = "BH") 
# NR are different from I & E 





# End of script ################################################################


