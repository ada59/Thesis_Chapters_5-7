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
  scaling          = TRUE, # TRUE????
  check_input      = TRUE,
  details_returned = TRUE) # All axes

#alpha_fd_indices4 <- mFD::alpha.fd.multidim(
#  sp_faxes_coord   = coords[ , c("PC1", "PC2", "PC3", "PC4")],
#  asb_sp_w         = All,
#  ind_vect         = c("fdis","fori","fspe", "fide"),
#  scaling          = TRUE,
#  check_input      = TRUE,
#  details_returned = TRUE) # All axes

###########################################################################################
# Correlations:----------------------------------------------------------------------------
divDis <- data.frame(alpha_fd_indices10$functional_diversity_indices[,c(1:2)])
divDis$SiteName <-  str_split_fixed(rownames(divDis),"_",2)[,1]
divDis$Period <- str_split_fixed(row.names(divDis),"_",2)[,2]
load("div3.RData")
new <- merge(divDis, div3, by=c("SiteName", "Period"))
cor(new$F0, new$fdis, method = "spearman")

###########################################################################################
# Shape of the curve (FDis) ---------------------------------------------------------------
divDis <- divDis[divDis$sp_richn > 1,]

hnc <- subset(divDis, rownames(divDis) %like% "HNC")
hnb <- subset(divDis, rownames(divDis) %like% "HNB")
contN <- subset(divDis, rownames(divDis) %like% "ContN")
contAll <- subset(divDis, rownames(divDis) %like% "ContAll")

# Historical conservative:
str(hnc)

hnc_1 <- lm(fdis ~ sp_richn, data=hnc)
hnc_2 <- lm(log(fdis) ~ log(sp_richn), data=hnc)
summary(hnc_1)
summary(hnc_2)

AIC(hnc_1, hnc_2) 

# Historical broad:
hnb_1 <- lm(fdis ~ sp_richn, data=hnb)
hnb_2 <- lm(log(fdis) ~ log(sp_richn), data=hnb)
summary(hnb_1)
summary(hnb_2)

AIC(hnb_1, hnb_2) 

# Current native:
contN_1 <- lm(fdis ~ sp_richn, data=contN)
contN_2 <- lm(log(fdis) ~ log(sp_richn), data=contN)
summary(contN_1)
summary(contN_2)

AIC(contN_1, contN_2) 

# Current All:
contAll_1 <- lm(fdis ~ sp_richn, data=contAll)
contAll_2 <- lm(log(fdis) ~ log(sp_richn), data=contAll)

AIC(contAll_1, contAll_2)


divDis$Period <- rep(NA, nrow(divDis))
divDis$Period <- ifelse(rownames(divDis) %like% "HNC", "A) Historical conservative", divDis$Period)
divDis$Period <- ifelse(rownames(divDis) %like% "HNB", "B) Historical broad", divDis$Period)
divDis$Period <- ifelse(rownames(divDis) %like% "ContN", "C) Current Native", divDis$Period)
divDis$Period <- ifelse(rownames(divDis) %like% "ContAll", "D) Current Native + Exotic", divDis$Period)

(oFDis <- ggplot(divDis, aes(x=sp_richn, y=fdis, color=Period))+
    geom_point(aes(color=Period))+
    geom_smooth()+
    theme_classic()+
    facet_wrap(~Period))

###########################################################################################
# Null models -----------------------------------------------------------------------------
load("rand.RData")
rand <- rand[, order(names(rand))]
coords <- coords[order(rownames(coords)),]
identical(colnames(dist_mat1), rownames(coords)) # FALSE

rand_l <- split(rand, f=rownames(rand))

F0null <- list()
for (i in 1:length(rand_l)){
  assemblage <- as.matrix(t(rand_l[[i]]))
  T0_n <- colSums(assemblage)
  F0_n <- FD_MLE(assemblage, dist_mat1, mean(dist_mat1[dist_mat1>0]), q=0)
  R0_n <- (T0_n - F0_n)
  F0null[[i]] <- cbind(T0_n, F0_n, R0_n)
}

###########################################################################################
# End of script ###########################################################################

