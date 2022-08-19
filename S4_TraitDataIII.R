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

rm(list=ls())
myd <- getwd()

# Read data:-------------------------------------------------------------------------------
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
corrplot(C, method="circle", type = "lower", order = "hclust")

file_path <- paste0(myd, "/Plots/Correlation matrix.png")
png(height=700, width=700, file=file_path)

corrplot(C, method="circle", type = "lower", order = "hclust", 
         p.mat = p.mat, sig.level = 0.01, insig = "blank")

dev.off()


