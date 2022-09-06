###########################################################################################
# Script: Sensitivity & Additional Analyses:
# Additional measures of trait diversity.
# Are the results of the community-based analysis strongly affected by the measure of trait diversity used?
# AFE
# September 2022
###########################################################################################


# Libraries:------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(ggrepel)
library(mFD)
library(data.table)


rm(list=ls())
myd <- getwd()
plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots" # Dir to save main plots

# Community data:-------------------------------------------------------------------------
load(paste0(myd, "/All.RData"))        # All surveys
load(paste0(myd, "/dist_mat1.RData"))  # Species x species trait distance matrix
load(paste0(myd, "/tt.RData"))         # Species x trait matrix
load(paste0(myd, "/coords.RData"))     # Coordinates in trait space

names(All)[names(All)=="Agonostomus monticola"] <- "Dajaus monticola"
rownames(All) <- paste0(All$SiteNameE, "_", All$Period)
All <- within(All, rm(SiteNameE, Period, DrainageBasinE))
All <- All[,order(names(All))]
coords <- coords[order(rownames(coords)),]

identical(rownames(coords), names(All)) #TRUE

coords <- as.matrix(coords)
All <- as.matrix(All)

All <- All[rowSums(All)>0,] # remove coms with no species

###########################################################################################
# Observed: -------------------------------------------------------------------------------

alpha_fd_indices10 <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = coords,
  asb_sp_w         = All,
  ind_vect         = c("fdis","fori","fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE) # All axes

alpha_fd_indices4 <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = coords[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = All,
  ind_vect         = c("fdis","fori","fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE) # All axes

# Correlations:
cor(alpha_fd_indices10$functional_diversity_indices[,c(2:7)],
    alpha_fd_indices4$functional_diversity_indices[,c(2:7)]) # Dispersion is highly correlated (0.98)

div <- as.data.frame(alpha_fd_indices10$functional_diversity_indices[,c(1:8)])

###########################################################################################
# Shape of the curve (FDis) ---------------------------------------------------------------

hnc <- subset(div, rownames(div) %like% "HNC")
hnb <- subset(div, rownames(div) %like% "HNB")
contN <- subset(div, rownames(div) %like% "ContN")
contAll <- subset(div, rownames(div) %like% "ContAll")

# Historical conservative:
hnc_1 <- lm(fdis ~ sp_richn, data=hnc)
hnc_2 <- lm(fdis ~ log(sp_richn+1), data=hnc)

AIC(hnc_1, hnc_2) # Curve

# Historical broad:
hnb_1 <- lm(fdis ~ sp_richn, data=hnb)
hnb_2 <- lm(fdis ~ log(sp_richn+1), data=hnb)

AIC(hnb_1, hnb_2) # Curve

# Current native:
contN_1 <- lm(fdis ~ sp_richn, data=contN)
contN_2 <- lm(fdis ~ log(sp_richn+1), data=contN)

AIC(contN_1, contN_2) # Curve

# Current All:
contAll_1 <- lm(fdis ~ sp_richn, data=contAll)
contAll_2 <- lm(fdis ~ log(sp_richn+1), data=contAll)

AIC(contAll_1, contAll_2) # Curve

# NOTE:
# A curve is a better fit in all cases.

div$Period <- rep(NA, nrow(div))
div$Period <- ifelse(rownames(div) %like% "HNC", "A) Historical conservative", div$Period)
div$Period <- ifelse(rownames(div) %like% "HNB", "B) Historical broad", div$Period)
div$Period <- ifelse(rownames(div) %like% "ContN", "C) Current Native", div$Period)
div$Period <- ifelse(rownames(div) %like% "ContAll", "D) Current Native + Exotic", div$Period)

(oFDis <- ggplot(div, aes(x=sp_richn, y=fdis, color=Period))+
    geom_point(aes(color=Period))+
    geom_smooth(method=lm, formula = y ~ log(x+1))+
    theme_classic()+
    facet_wrap(~Period))

(oFDisII <- ggplot(div, aes(x=Period, y=fdis, fill=Period))+
    geom_violin(aes(fill=Period))+
    geom_boxplot(aes(fill=Period), width=0.1, position=position_dodge(1))+
    theme_classic())


###########################################################################################
# F Div and identity by period ------------------------------------------------------------
div <- gather(div, key="Metric", value="Value", -Period)
divI <- subset(div, div$Metric  %in% c("fdis", "fori", "fspe"))
divII <- subset(div, div$Metric %in% c("fide_PC1", "fide_PC2", "fide_PC3", "fide_PC4"))

(fmetricsI <- ggplot(divI, aes(y = Value, x = Period, color=Period)) +
   geom_point(alpha=0.05, position=position_dodge(1))+
   geom_violin(alpha=0.4) +
   geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
   theme_bw()+
   facet_wrap(~Metric)) # The current communities are the least original.

(fmetricsI <- ggplot(divII, aes(y = Value, x = Period, color=Period)) +
    geom_point(alpha=0.05, position=position_dodge(1))+
    geom_violin(alpha=0.4) +
    geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
    theme_bw()+
    facet_wrap(~Metric)) # The current communities are the least original.

###########################################################################################
# End of script ###########################################################################

