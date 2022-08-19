###########################################################################################
# Script: Imputatitions
# AFE
# August 2022
###########################################################################################

# Libraries:-------------------------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(hablar)
library(missForest)

rm(list=ls())
myd <- getwd()

# Read data:-------------------------------------------------------------------------------
load("fishtrait_complete.RData")
str(fishtrait_complete)

# Imputations: ----------------------------------------------------------------------------
rownames(fishtrait_complete) <- fishtrait_complete$Genus.species
fishtrait_complete <- fishtrait_complete[c("MBl", "BEl", "VEp", "REs", "OGp", "RMl","BLs","PFv","PFs","CPt")]


sum(is.na(fishtrait_complete))  # 15
(15/(100*10))*100     # 1.5 %

ntimes <- 100
t_imp <- list()
t_err <- list()
for (i in 1:ntimes){
  mat <- as.matrix(fishtrait_complete)
  miss <- missForest(mat)
  mat_imp <- miss$ximp
  mat_err <- miss$OOBerror
  t_imp[[i]] <- mat_imp
  t_err[[i]] <- mat_err
}
t_impm <- do.call(cbind, t_imp)
t_impm <- array(t_impm, dim=c(dim(t_imp[[1]]), length(t_imp)))
t_impm <- as.matrix(apply(t_impm, c(1, 2), mean))              # Mean across the 100 objects
rownames(t_impm) <- rownames(fishtrait_complete)
colnames(t_impm) <- colnames(fishtrait_complete)

fishtrait_complete_imp <- t_impm

save(t_imp, file="t_imp.RData")
save(t_err, file="t_err.RData")
save(fishtrait_complete_imp, file="fishtrait_complete_imp.RData")

###########################################################################################
# End of script ###########################################################################

