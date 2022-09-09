###########################################################################################
# Script: Model assumptions
# AFE
# August 2022
###########################################################################################

library(dplyr)
library(tidyverse)

rm(list=ls())
myd <- getwd()
sm_plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots/SM" # Dir to save main plots


load("model_l.RData")


# Histograms:

par(mfrow = c(1, 2))
hist(model_l$hnc_linear$residuals, main= "Linear")
hist(model_l$hnc_loglinear$residuals, main="Log-Linear")

par(mfrow = c(1, 2))
hist(model_l$hnb_linear$residuals, main= "Linear")
hist(model_l$hnb_loglinear$residuals, main="Log-Linear")

par(mfrow = c(1, 2))
hist(model_l$contN_linear$residuals, main= "Linear")
hist(model_l$contN_loglinear$residuals, main="Log-Linear")

par(mfrow = c(1, 2))
hist(model_l$contAll_linear$residuals, main= "Linear")
hist(model_l$contAll_loglinear$residuals, main="Log-Linear")
dev.off()

# Plots:

par(mfrow = c(2, 2))
plot(model_l$hnc_linear)

par(mfrow = c(2, 2))
plot(model_l$hnc_loglinear)

par(mfrow = c(2, 2))
plot(model_l$hnb_linear)

par(mfrow = c(2, 2))
plot(model_l$hnb_loglinear)

par(mfrow = c(2, 2))
plot(model_l$contN_linear)

par(mfrow = c(2, 2))
plot(model_l$contN_loglinear)

par(mfrow = c(2, 2))
plot(model_l$contAll_linear)

par(mfrow = c(2, 2))
plot(model_l$contAll_loglinear)
dev.off()