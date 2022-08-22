###########################################################################################
# Script: Compute redundancies and test for temporal change
# AFE
# August 2022
###########################################################################################

# Libraries:------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(raster)
library(gridExtra)
library(cluster)
library(factoextra)
library(ggrepel)
library(adiv)
library(mFD)

rm(list=ls())
myd <- getwd()
plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots" # Dir to save main plots

# Community data:-------------------------------------------------------------------------
load(paste0(myd, "/HNC.RData"))   # Conservative Historical Native Assemblage
load(paste0(myd, "/HNB.RData"))   # Broad Historical Native Assemblage
load(paste0(myd, "/ContN.RData")) # Contemporary Native Assemblage
load(paste0(myd, "/ContE.RData")) # Contemporary Exotic assemblage

load(paste0(myd, "/dist_mat1.RData"))  # Trait distance matrix

# Functions:------------------------------------------------------------------------------
# FD and TD (Hill numbers) use the same units and TD is always greater than or equal to FD
# For q=0, Redundancy is S-FD0
hill_dir <- "C:/Users/Usuario/Documents/PHD/Functions/FunD-master/"

fskt <- read.table(paste0(hill_dir,"Plant_Abundance.txt"))
fskt = list(fs=fskt$Fushan, kt=fskt$Kenting)
dij_fskt = read.table(paste0(hill_dir,"Plant_Distance_Matrix.txt"))
dij_fskt = as.matrix(dij_fskt) # All needed for the code to run

source("C:/Users/Usuario/Documents/PHD/Functions/FunD_Rcode.R")

# Data formatting:------------------------------------------------------------------------
ContAll <- full_join(ContN, ContE, by=c("SiteNameE", "DrainageBasinE"))
names(ContAll)

ContAll$`Goodea atripinnis`<- ContAll$`Goodea atripinnis.x`+ ContAll$`Goodea atripinnis.y`
ContAll$`Poecilia butleri`<- ContAll$`Poecilia butleri.x`+ ContAll$`Poecilia butleri.y`
ContAll$`Poecilia sphenops`<- ContAll$`Poecilia sphenops.x`+ ContAll$`Poecilia sphenops.y`
ContAll <- ContAll[, ! names(ContAll) %in% c("Goodea atripinnis.x", "Goodea atripinnis.y",
                                            "Poecilia butleri.x", "Poecilia butleri.y",
                                            "Poecilia sphenops.x", "Poecilia sphenops.y")]

l <- bind_rows(HNC, HNB, ContN, ContAll)
l[is.na(l)] <- 0
sum(is.na(l)) # 0

l$Period <- rep(c("HNC", "HNB", "ContN", "ContAll"), each=67)
rownames(l) <- paste0(l$SiteNameE, "_", l$Period)
l <- within(l, rm(SiteNameE, DrainageBasinE, Period))
l <- l[, order(names(l))]

str(l)

ls <- split(l, rownames(l)) # 67*4

setdiff(colnames(dist_mat1), names(ls[[1]]))

colnames(dist_mat1)[colnames(dist_mat1)=="Dajaus monticola"] <-"Agonostomus monticola"
rownames(dist_mat1)[rownames(dist_mat1)=="Dajaus monticola"] <-"Agonostomus monticola"

dist_mat1 <- dist_mat1[, order(colnames(dist_mat1))]
dist_mat1 <- dist_mat1[order(rownames(dist_mat1)),]

identical(colnames(dist_mat1), names(ls[[1]])) #TRUE

ls <- lapply(ls, function(x) {as.vector(as.matrix(x))})

##########################################################################################
# Observed TD and FD: --------------------------------------------------------------------
divT0 <- lapply(ls, function(x) {FD_MLE(x, dist_mat1, min(dist_mat1[dist_mat1>0]), 0)})
divF0 <- lapply(ls, function(x) {FD_MLE(x, dist_mat1, mean(dist_mat1[dist_mat1>0]), 0)})

divT0d <- as.data.frame(do.call(rbind, divT0))
divT0d$SiteName <- str_split_fixed(rownames(divT0d), "_", 2)[,1]
divT0d$Period <- str_split_fixed(rownames(divT0d), "_", 2)[,2]

divF0d <- as.data.frame(do.call(rbind, divF0))
divF0d$SiteName <- str_split_fixed(rownames(divF0d), "_", 2)[,1]
divF0d$Period <- str_split_fixed(rownames(divF0d), "_", 2)[,2]

div <- left_join(divT0d, divF0d, by=c("SiteName", "Period"))
names(div) <- c("T0", "SiteName", "Period", "F0")

# Observed trends: -----------------------------------------------------------------------
# Subsets:
hnc <- subset(div, div$Period=="HNC")
hnb <- subset(div, div$Period=="HNB")
contN <- subset(div, div$Period=="ContN")
contAll <- subset(div, div$Period=="ContAll")

str(hnc)

sum(hnc$T0==0)
sum(hnb$T0==0)
sum(contN$T0==0)   #21
sum(contAll$T0==0) #14 (67-53, OK)

contN <- contN[!contN$T0==0,]       # 46
contAll <- contAll[!contAll$T0==0,] # 53

# Historical conservative:
hnc_1 <- lm(F0 ~ T0, data=hnc)
hnc_2 <- lm(F0 ~ log(T0+1), data=hnc)
hnc_3 <- lm(F0 ~ log(T0), data=hnc)

AIC(hnc_1, hnc_2, hnc_3) # Curve

# Historical broad:
hnb_1 <- lm(F0 ~ T0, data=hnb)
hnb_2 <- lm(F0 ~ log(T0+1), data=hnb)
hnb_3 <- lm(F0 ~ log(T0), data=hnb)

AIC(hnb_1, hnb_2, hnb_3)

# Current native:
contN_1 <- lm(F0 ~ T0, data=contN)
contN_2 <- lm(F0 ~ log(T0+1), data=contN)
contN_3 <- lm(F0 ~ log(T0), data=contN)

AIC(contN_1, contN_2, contN_3)

# Current All:
contAll_1 <- lm(F0 ~ T0, data=contAll)
contAll_2 <- lm(F0 ~ log(T0+1), data=contAll)
contAll_3 <- lm(F0 ~ log(T0), data=contAll)

AIC(contAll_1, contAll_2, contAll_3)

# A curve is a better fit in all cases.

(p1 <- ggplot(div, aes(x=T0, y=F0, color=Period))+
  geom_point(aes(color=Period))+
  geom_smooth()+
  theme_bw()+
  facet_wrap(~Period))

(p2 <- ggplot(div, aes(x=T0, y=F0, color=Period))+
    geom_point(aes(color=Period))+
    geom_smooth(method=lm, formula = y ~ log(x+1))+
    theme_bw()+
    facet_wrap(~Period))

# Remove locality with richness = 17 in the historical period.

hnc66 <- hnc[!hnc$T0 ==17,]
hnb66 <- hnb[!hnb$T0 ==17,]

hnc66_1 <- lm(F0 ~ T0, data=hnc66)
hnc66_2 <- lm(F0 ~ log(T0+1), data=hnc66)

AIC(hnc66_1, hnc66_2)

hnb66_1 <- lm(F0 ~ T0, data=hnb66)
hnb66_2 <- lm(F0 ~ log(T0+1), data=hnb66)

AIC(hnb66_1, hnb66_2)

# The locality with more species doesn't drive the trend observed 
# (asymptote is a better fit)

##########################################################################################
# Null TD and FD: ------------------------------------------------------------------------



###########################################################################################
# End of script ###########################################################################
