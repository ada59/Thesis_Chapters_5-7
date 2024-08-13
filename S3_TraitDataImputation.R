#===============================================================================
# Script: Imputation of missing values
# AFE
# August 2022
#===============================================================================



#===============================================================================
# Libraries:--------------------------------------------------------------------
#===============================================================================
library(dplyr)
library(tidyverse)
library(hablar)
library(missForest)

rm(list=ls())
myd <- getwd()




#===============================================================================
# Read data:--------------------------------------------------------------------
#===============================================================================
path_traits <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Traits"
load(paste0(path_traits, "/fishtrait_complete.RData"))
str(fishtrait_complete)




#===============================================================================
# Imputations: -----------------------------------------------------------------
#===============================================================================
fishtrait_completeII <- fishtrait_complete[c("MBl", "BEl", "VEp", "REs", "OGp", "RMl","BLs","PFv","PFs","CPt")]
rownames(fishtrait_completeII) <- fishtrait_complete$Genus.species

sum(is.na(fishtrait_completeII))  # 14
(14/(95*10))*100     # 1.5 %

ntimes <- 100
t_imp <- list()
t_err <- list()
for (i in 1:ntimes){
  mat <- as.matrix(fishtrait_completeII)
  miss <- missForest(mat)
  mat_imp <- miss$ximp
  mat_err <- miss$OOBerror
  t_imp[[i]] <- mat_imp
  t_err[[i]] <- mat_err
}
t_impm <- do.call(cbind, t_imp)
t_impm <- array(t_impm, dim=c(dim(t_imp[[1]]), length(t_imp)))
t_impm <- as.matrix(apply(t_impm, c(1, 2), mean))              # Mean across the 100 objects
rownames(t_impm) <- rownames(fishtrait_completeII)
colnames(t_impm) <- colnames(fishtrait_completeII)

fishtrait_complete_imp <- t_impm

save(t_imp, file=paste0(path_traits, "/t_imp.RData"))
save(t_err, file=paste0(path_traits, "/t_err.RData"))
save(fishtrait_complete_imp, file=paste0(path_traits, "/fishtrait_complete_imp.RData"))




# End of script ################################################################

