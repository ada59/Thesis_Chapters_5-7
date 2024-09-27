#===============================================================================
# Script: Model the relationship between SR and T0
# AFE
# August/Sept 2022 & August 2024
#===============================================================================




#===============================================================================
# Libraries:--------------------------------------------------------------------
#===============================================================================
library(Rmisc)
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
library(AICcmodavg)
library(ggpubr)
library(PairedData)
library(grid)
library(lmtest)
library(coefplot)
library(stats)
library(investr)

rm(list=ls())

path_lists <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists"   # Raw data lists
path_traits <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Traits" # Dir to trait objects
path_plots6 <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots/Chapter6"   # Dir to save plots Chapter 6




#===============================================================================
# read community data:----------------------------------------------------------
#===============================================================================
load(paste0(path_lists,"/l67.RData"))     # main list with community data matrices (67 sites)
load(paste0(path_lists, "/l53.RData"))    # main list with community data matrices (53 sites)
load(paste0(path_lists, "/NDT67.RData"))  # List for checks on raw data
load(paste0(path_lists, "/NDT83.RData"))  # List for checks on raw data
load(paste0(path_traits, "/dist_mat1.RData"))

l67$SiteNameE_Period <- paste0(l67$SiteNameE, "_", l67$Period)
l67 <- within(l67, rm(SiteNameE, DrainageBasinE, Period))

l67 <- l67[,order(colnames(l67))]
dist_mat1 <- dist_mat1[order(rownames(dist_mat1)), order(colnames(dist_mat1))]

lassemb <- split(l67, f=l67$SiteNameE_Period)
lassemb <- lapply(lassemb, function(x) {within(x, rm(SiteNameE_Period))}) # always 95, OK
identical(colnames(dist_mat1), names(lassemb[[1]])) # TRUE, crucial to obtain correct results below
lassemb <- lapply(lassemb, function(x) {as.vector(as.matrix(x))}) # list of assemblages (wrapping of as vector needed for FD_MLE)
      



#===============================================================================
# read hill functions:----------------------------------------------------------
#===============================================================================
# FD and TD (Hill numbers) use the same units and TD is always greater than or equal to FD
# For q=0, Redundancy is S-FD0
hill_dir <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/Functions/FunD-master"

fskt <- read.table(paste0(hill_dir,"/Plant_Abundance.txt"))
fskt = list(fs=fskt$Fushan, kt=fskt$Kenting)
dij_fskt = read.table(paste0(hill_dir,"/Plant_Distance_Matrix.txt"))
dij_fskt = as.matrix(dij_fskt) # All needed for the code to run

source("C:/Users/afe1/OneDrive - University of St Andrews/PHD/Functions/FunD_Rcode.R") 
# NOTE: Warnings don't affect FD_MLE(), Error is also OK and does not affect FD_MLE()




#===============================================================================
# Observed TD and FD: ----------------------------------------------------------
#===============================================================================
divT0 <- lapply(lassemb, function(x) {FD_MLE(x, dist_mat1, min(dist_mat1[dist_mat1>0]), 0)})  # taxonomic diversity
divF0 <- lapply(lassemb, function(x) {FD_MLE(x, dist_mat1, mean(dist_mat1[dist_mat1>0]), 0)}) # trait diversity

divT0d <- as.data.frame(do.call(rbind, divT0))
divT0d$SiteName <- str_split_fixed(rownames(divT0d), "_", 2)[,1]
divT0d$Period <- str_split_fixed(rownames(divT0d), "_", 2)[,2]

divF0d <- as.data.frame(do.call(rbind, divF0))
divF0d$SiteName <- str_split_fixed(rownames(divF0d), "_", 2)[,1]
divF0d$Period <- str_split_fixed(rownames(divF0d), "_", 2)[,2]

div <- left_join(divT0d, divF0d, by=c("SiteName", "Period"))
names(div) <- c("T0", "SiteName", "Period", "F0")
div <- div %>% relocate(T0, .after = Period)
str(div)

div$SiteName <- as.factor(div$SiteName)
div$Period <- as.factor(div$Period)

div$R0 <- div$T0-div$F0 # trait redundancy


save(lassemb, file="lassemb.RData")                       # used in other CM Chapter
save(divT0, file=paste0(path_traits, "/divT0.RData"))     # used in other CM Chapter
save(divF0, file=paste0(path_traits, "/divF0.RData"))     # used in other CM Chapter




#===============================================================================
# Data subsets: ----------------------------------------------------------------
#===============================================================================

## Considering all localities with data for each period: -----------------------
hist(div$F0)
div2 <- div[!div$T0==0,]   # 270                    [KEEPING SR = 1]
hist(div2$F0)
div3 <- div2[!div2$T0==1,] # 210 (MAIN ANALYSIS)    [NOT KEEPING SR = 1]
hist(div3$F0)
#save(div2, file = "div2.RData")  # these include contE, which won't be used in analyses
save(div3, file = "div3.RData")  # these include contE, which won't be used in analyses


div3 <- div3   # change to div2 for sensitivity


par(mfrow=c(1,2))
hist(div2$T0)
hist(div2$F0)
dev.off()
par(mfrow=c(1,2))
hist(div3$T0)
hist(div3$F0) # More normally distributed
dev.off()


#hnc2 <- subset(div2, div2$Period=="HNC")           # 67
#hnb2 <- subset(div2, div2$Period=="HNB")           # 67
#contN2 <- subset(div2, div2$Period=="ContN")       # 46
#contAll2 <- subset(div2, div2$Period=="ContAll")   # 53

hnc3 <- subset(div3, div3$Period=="HNC")           # 48
hnb3 <- subset(div3, div3$Period=="HNB")           # 60
contN3 <- subset(div3, div3$Period=="ContN")       # 33
contAll3 <- subset(div3, div3$Period=="ContAll")   # 49


## Using only sites in contN3:--------------------------------------------------
#hnc3 <- subset(hnc3, hnc3$SiteName %in% contN3$SiteName)
#hnb3 <- subset(hnb3, hnb3$SiteName %in% contN3$SiteName)

# Remove locality with richness = 17 in the historical period:
#hnc3 <- hnc3[!hnc3$T0 ==17,]
#hnb3 <- hnb3[!hnb3$T0 ==17,]




#===============================================================================
# Relationship between T0 & Period:---------------------------------------------
#===============================================================================
str(div2)
str(div3)
div3_1T0 <- lm(T0 ~ Period, data=div3)
summary(div3_1T0)

plot(T0 ~ Period, data=div3)
dev.off()




#===============================================================================
# Model sources: ---------------------------------------------------------------
#===============================================================================
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0026735 [Guillemot 2011]
# https://tuos-bio-data-skills.github.io/intro-stats-book/non-linear-regression-in-R.html
# https://statisticsbyjim.com/regression/r-squared-invalid-nonlinear-regression/ [ON WHY NOT USE R2]
# https://www.datacamp.com/tutorial/introduction-to-non-linear-model-and-insights-using-r [ok to fit linear with NLS]
# https://stackoverflow.com/questions/61341287/how-to-calculate-confidence-intervals-for-nonlinear-least-squares-in-r





#===============================================================================
# Linear regression:------------------------------------------------------------
#===============================================================================

l3 <- list("HNC"=hnc3, "HNB"=hnb3, "ContN"=contN3, "ContAll"=contAll3)
save(l3, file = "l3.RData")  


## lm: -------------------------------------------------------------------------
fit_linear_lm_3 <- lapply(l3, function(x) {lm(F0 ~ T0, data=x)})       # linear OLS
summlinear3 <- lapply(fit_linear_lm_3, function(x) {summary(x)})       # summaries



## Plotting & assumptions (lm): ------------------------------------------------

#### Historical Conservative ####
summlinear3[[1]]
plot(F0 ~ T0, data=hnc3, main="Historical Conservative")
curve(predict(fit_linear_lm_3[[1]], newdata = data.frame(T0 = x)), col = "green", lty=1, add = TRUE)
par(mfrow=c(2,2))
plot(fit_linear_lm_3[[1]])
dev.off()
plot(fit_linear_lm_3[[1]]$fitted,fit_linear_lm_3[[1]]$residuals)
plot(cooks.distance(fit_linear_lm_3[[1]])) # Index 16 (assemblage with 17 sps)


#### Historical Broad ####
summlinear3[[2]]
plot(F0 ~ T0, data=hnb3, main="Historical Broad")
curve(predict(fit_linear_lm_3[[2]], newdata = data.frame(T0 = x)), col = "green", lty=1, add = TRUE)
par(mfrow=c(2,2))
plot(fit_linear_lm_3[[2]])
dev.off()
plot(fit_linear_lm_3[[1]]$fitted,fit_linear_lm_3[[1]]$residuals)
plot(cooks.distance(fit_linear_lm_3[[2]])) # Index 16 (assemblage with 17 sps)
hnb3[c(17,33),] # Sites with higher richness.


#### Contemporary Native ####
summlinear3[[3]]
plot(F0 ~ T0, data=hnb3, main="Contemporary Native")
curve(predict(fit_linear_lm_3[[3]], newdata = data.frame(T0 = x)), col = "green", lty=1, add = TRUE)
par(mfrow=c(2,2))
plot(fit_linear_lm_3[[3]])
dev.off()
#plot(fit_linear_lm_3[[1]]$fitted,fit_linear_lm_3[[1]]$residuals)
plot(cooks.distance(fit_linear_lm_3[[3]])) # Index 16 (assemblage with 17 sps)
contN3[c(19,32),] # Sites with relatively higher richness.


#### Contemporary Native + Introduced ####
summlinear3[[4]]
plot(F0 ~ T0, data=hnb3, main="Contemporary Native + Introduced")
curve(predict(fit_linear_lm_3[[4]], newdata = data.frame(T0 = x)), col = "green", lty=1, add = TRUE)
par(mfrow=c(2,2))
plot(fit_linear_lm_3[[4]])
dev.off()
#plot(fit_linear_lm_3[[1]]$fitted,fit_linear_lm_3[[1]]$residuals)
plot(cooks.distance(fit_linear_lm_3[[4]])) # Index 16 (assemblage with 17 sps)
contAll3[c(18,26,47,49),] 

# NOTES: #######################################################################
# Outliers tend to be at large values of species richness?
################################################################################


## nls: ------------------------------------------------------------------------

fit_linear_nls_hnc_3 <- nls(F0 ~ b0 + b*T0, data = l3[[1]], 
                            start = list(b0 = 0, b = 1)) # historical conservative
summary(fit_linear_nls_hnc_3)
summlinear3[[1]] # idem estimates

fit_linear_nls_hnb_3 <- nls(F0 ~ b0 + b*T0, data = l3[[2]], 
                            start = list(b0 = 0, b = 1)) # historical broad
summary(fit_linear_nls_hnb_3)
summlinear3[[2]] # idem estimates

fit_linear_nls_contN_3 <- nls(F0 ~ b0 + b*T0, data = l3[[3]], 
                            start = list(b0 = 0, b = 1)) # contemporary native
summary(fit_linear_nls_contN_3)
summlinear3[[3]] # idem estimates

fit_linear_nls_contAll_3 <- nls(F0 ~ b0 + b*T0, data = l3[[4]], 
                              start = list(b0 = 0, b = 1)) # contemporary native
summary(fit_linear_nls_contAll_3)
summlinear3[[4]] # idem estimates




#===============================================================================
## Power regression (log-log & nls):--------------------------------------------
#===============================================================================

## power regression with lm: ---------------------------------------------------
lapply(l3, function(x) {range(x$F0)}) # No need to do log + 1

fit_power_lm_3 <- lapply(l3, function(x) {lm(log(F0) ~ log(T0), data=x)}) #log(F0 + 1) is keeping 1s
summloglog3 <- lapply(fit_power_lm_3, function(x) {summary(x)})


# Power regression with nls: (different error family, additive)-----------------

fit_power_nls_hnc_3 <- nls(F0 ~ a * T0^b, data = l3[[1]], 
                 start = list(a = exp(coef(fit_power_lm_3[[1]])[1]), b = coef(fit_power_lm_3[[1]])[2]))
fit_power_nls_hnb_3 <- nls(F0 ~ a * T0^b, data = l3[[2]], 
                 start = list(a = exp(coef(fit_power_lm_3[[2]])[1]), b = coef(fit_power_lm_3[[2]])[2]))
fit_power_nls_contN_3 <- nls(F0 ~ a * T0^b, data = l3[[3]], 
                 start = list(a = exp(coef(fit_power_lm_3[[3]])[1]), b = coef(fit_power_lm_3[[3]])[2]))
fit_power_nls_contAll_3 <- nls(F0 ~ a * T0^b, data = l3[[4]], 
                 start = list(a = exp(coef(fit_power_lm_3[[4]])[1]), b = coef(fit_power_lm_3[[4]])[2]))



## Plotting & assumptions: -----------------------------------------------------

#### Historical Conservative ####
plot(F0 ~ T0, data=hnc3, main="Historical Conservative")
curve(exp(predict(fit_power_lm_3[[1]], newdata = data.frame(T0 = x))), col = "black", lty=2, add = TRUE)
curve(predict(fit_power_nls_hnc_3, newdata = data.frame(T0 = x)), col = "red", lty=3, add = TRUE)
legend("bottomright", legend = c("Log-Log", "Nls"),
       col = c("black", "red"),
       lty = c(2,3), cex=0.5)

par(mfrow=c(2,2))
plot(fit_power_nls_hnc_3) 
dev.off()

file_path <- paste0(path_plots6, "/S5_AssumptionsHNC.png") # Assumptions HNC
png(width = 600, height = 300, file=file_path)
par(mfrow=c(1,2))
qqnorm(resid(fit_power_nls_hnc_3)/sd(resid(fit_power_nls_hnc_3)))
abline(0,1) 
plot(predict(fit_power_nls_hnc_3), resid(fit_power_nls_hnc_3), xlab="Fitted values", ylab="Residuals", main="Residuals vs Fitted")
abline(h=0)
dev.off()


#### Historical Broad ####
plot(F0 ~ T0, data=hnb3, main="Historical Broad")
curve(exp(predict(fit_power_lm_3[[2]], newdata = data.frame(T0 = x))), col = "black", lty=2, add = TRUE)
curve(predict(fit_power_nls_hnb_3, newdata = data.frame(T0 = x)), col = "red", lty=3, add = TRUE)
legend("bottomright", legend = c("Log-Log", "Nls"),
       col = c("black", "red"),
       lty = c(2,3), cex=0.5)

par(mfrow=c(2,2))
plot(fit_power_lm_3[[2]]) # Cooks distance large max 0.14
dev.off()
plot(predict(fit_power_nls_hnb_3), resid(fit_power_nls_hnb_3), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit_power_nls_hnb_3)/sd(resid(fit_power_nls_hnb_3)))
abline(0,1) #OK (the power model seems good)



#### Contemporary Native ####
plot(F0 ~ T0, data=contN3)
curve(exp(predict(fit_power_lm_3[[3]], newdata = data.frame(T0 = x))), col = "black", lty=2, add = TRUE)
curve(predict(fit_power_nls_contN_3, newdata = data.frame(T0 = x)), col = "red", lty=3, add = TRUE)
legend("bottomright", legend = c("Log-Log", "Nls"),
       col = c("black", "red"),
       lty = c(2, 3), cex=0.5)

par(mfrow=c(2,2))
plot(fit_power_lm_3[[3]])
dev.off()
plot(predict(fit_power_nls_contN_3), resid(fit_power_nls_contN_3), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit_power_nls_contN_3)/sd(resid(fit_power_nls_contN_3)))
abline(0,1) #OK (less good than in the historical period)



#### Contemporary Native + Introduced ####
plot(F0 ~ T0, data=contAll3)
curve(exp(predict(fit_power_lm_3[[4]], newdata = data.frame(T0 = x))), col = "black", lty=2, add = TRUE)
curve(predict(fit_power_nls_contAll_3, newdata = data.frame(T0 = x)), col = "red", lty=3, add = TRUE)
legend("bottomright", legend = c("Log-Log", "Nls"),
       col = c("black", "red"),
       lty = c(2, 3), cex=0.5)

par(mfrow=c(2,2))
plot(fit_power_lm_3[[4]])
dev.off()
plot(predict(fit_power_nls_contAll_3), resid(fit_power_nls_contAll_3), xlab="Fitted values", ylab="Residuals")
abline(h=0) # At low & high values we observe all dots under 0
qqnorm(resid(fit_power_nls_contAll_3)/sd(resid(fit_power_nls_contAll_3)))
abline(0,1) #OK (less good than in the historical period) 




#===============================================================================
# Asymptote (nls):--------------------------------------------------------------
#===============================================================================
dev.off()

#### Historical Conservative ####
#  When x = 0, y = 0 When x is huge, y ∼ k
plot(F0 ~ T0, data=hnc3) # xlim=c(1, 30) with predict x2 below to confirm 5 is the asymptote

x=c(2:17)
x2=c(2:30) # to confirm asymptote by increasing x axis
my.k <- 6  # 6. parameter for model, value at which model reaches asymptote
my.b <- 1  # 1. parameter for model, linear increase at the beginning
# NOTE: when using nls() if estimates are in a reasonable range, 
# the output is the same, nls finds parameters that minimise the sum 
# of squared residuals.

# my.k*(1-exp((-1)*my.b*x)) can be used to help determine k & b

fit4_hnc3 <- nls(F0~k*(1-exp((-1)*b*T0)), data=hnc3, start=list(k=my.k,b=my.b))
lines(x, predict(fit4_hnc3,newdata=data.frame(T0=x)),lwd=2,col='blue') # plot the model
summary(fit4_hnc3) # Parameters seem to be in the correct range.


#### Historical Broad ####
plot(F0 ~ T0, data=hnb3)
x=c(2:17)
my.k <- 4.5  # parameter for model
my.b <- 1    # parameter for model

fit4_hnb3 <- nls(F0~k*(1-exp((-1)*b*T0)), data=hnb3, start=list(k=my.k,b=my.b))
lines(x, predict(fit4_hnb3,newdata=data.frame(T0=x)),lwd=2,col='blue') # plot the model
summary(fit4_hnb3) # Parameters seem to be in the correct range.


#### Contemporary Native ####
plot(F0 ~ T0, data=contN3)
x=c(2:9)
my.k <- 3  # parameter for model
my.b <- 1    # parameter for model

fit4_contN3 <- nls(F0~k*(1-exp((-1)*b*T0)), data=contN3, start=list(k=my.k,b=my.b))
lines(x, predict(fit4_contN3,newdata=data.frame(T0=x)),lwd=2,col='blue') # plot the model
summary(fit4_contN3) # Parameters seem to be in the correct range.

#### Contemporary Native + Introduced ####
plot(F0 ~ T0, data=contAll3)
x=c(2:9)
my.k <- 3  # parameter for model
my.b <- 1    # parameter for model

fit4_contAll3 <- nls(F0~k*(1-exp((-1)*b*T0)), data=contAll3, start=list(k=my.k,b=my.b))
lines(x, predict(fit4_contAll3,newdata=data.frame(T0=x)),lwd=2,col='blue') # plot the model
summary(fit4_contAll3) # Parameters seem to be in the correct range.



## Plotting & assumptions: -----------------------------------------------------
dev.off()
plot(predict(fit4_hnc3), resid(fit4_hnc3), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit4_hnc3)/sd(resid(fit4_hnc3)))
abline(0,1) # pretty good


file_path <- paste0(path_plots6, "/S5_AssumptionsHNB.png")
png(width = 600, height = 300, file=file_path)
par(mfrow=c(1,2))
qqnorm(resid(fit4_hnb3)/sd(resid(fit4_hnb3)))
abline(0,1) 
plot(predict(fit4_hnb3), resid(fit4_hnb3), xlab="Fitted values", ylab="Residuals", main="Residuals vs Fitted")
abline(h=0) # passable
dev.off()


file_path <- paste0(path_plots6, "/S5_AssumptionsContN.png")
png(width = 600, height = 300, file=file_path)
par(mfrow=c(1,2))
qqnorm(resid(fit4_contN3)/sd(resid(fit4_contN3)))
abline(0,1) # not great
plot(predict(fit4_contN3), resid(fit4_contN3), xlab="Fitted values", ylab="Residuals", main="Residuals vs Fitted")
abline(h=0) # not great...
dev.off()


file_path <- paste0(path_plots6, "/S5_AssumptionsContAll.png")
png(width = 600, height = 300, file=file_path)
par(mfrow=c(1,2))
qqnorm(resid(fit4_contAll3)/sd(resid(fit4_contAll3)))
abline(0,1) #OK (not great) /  over-predicting?
plot(predict(fit4_contAll3), resid(fit4_contAll3), xlab="Fitted values", ylab="Residuals", main="Residuals vs Fitted")
abline(h=0) # At low & high values we observe all dots under 0
dev.off()



#===============================================================================
# Logistic (nls): --------------------------------------------------------------
#===============================================================================
# SSlogis == y~a/(1 + exp(-b * (x-c)) 
# [Logistic but slightly different from GUillemot 2011]
#(a/(1+b*exp(-c*T0))) -(a/(1+b)) # formula sigmoidal in Guillemot 2011


#### Historical Conservative ####
plot(F0 ~ T0, data=hnc3) # a priori no indication of sigmoidal curve
#my.a <- 5.5
#my.b <- 1
#my.c <- 1
#fit5_hnc <- nls(F0~(a/(1+b*exp(-c*T0))) -(a/(1+b)), data=hnc3, start=list(a=my.a,b=my.b, c=my.c))
# tried multiple vals above, it seems to not perform well (over-parametrization)

fit5_hnc <- nls(F0 ~ SSlogis(T0, Asym, xmid, scal), data = hnc3) # SSlogis
summary(fit5_hnc)

#### Historical Broad ####
plot(F0 ~ T0, data=hnb3) 
#fit5_hnb <- nls(F0~(a/(1+b*exp(-c*T0))) -(a/(1+b)), data=hnb3, start=list(a=my.a,b=my.b, c=my.c))
fit5_hnb <- nls(F0 ~ SSlogis(T0, Asym, xmid, scal), data = hnb3)

#### Contemporary Native ####
plot(F0 ~ T0, data=contN3) 
#fit5_contN <- nls(F0~(a/(1+b*exp(-c*T0))) -(a/(1+b)), data=contN3, start=list(a=my.a,b=my.b, c=my.c))
fit5_contN <- nls(F0 ~ SSlogis(T0, Asym, xmid, scal), data = contN3)

#### Contemporary Native + Introduced ####
plot(F0 ~ T0, data=contAll3) 
#fit5_contAll <- nls(F0~(a/(1+b*exp(-c*T0))) -(a/(1+b)), data=contAll3, start=list(a=my.a,b=my.b, c=my.c))
fit5_contAll <- nls(F0 ~ SSlogis(T0, Asym, xmid, scal), data = contAll3)


## Plots & Assumptions: --------------------------------------------------------
plot(F0 ~ T0, data=hnc3)
x=c(2:17)
lines(x, predict(fit5_hnc,newdata=data.frame(T0=x)),lwd=2,col='blue') # plot the model
plot(predict(fit5_hnc), resid(fit5_hnc), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit5_hnc)/sd(resid(fit5_hnc)))
abline(0,1)

plot(F0 ~ T0, data=hnb3)
x=c(2:17)
lines(x, predict(fit5_hnb,newdata=data.frame(T0=x)),lwd=2,col='blue') # plot the model
plot(predict(fit5_hnb), resid(fit5_hnb), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit5_hnb)/sd(resid(fit5_hnb)))
abline(0,1)

plot(F0 ~ T0, data=contN3)
x=c(2:9)
lines(x, predict(fit5_contN,newdata=data.frame(T0=x)),lwd=2,col='blue') # plot the model
plot(predict(fit5_contN), resid(fit5_contN), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit5_contN)/sd(resid(fit5_contN)))
abline(0,1) # not great at all

plot(F0 ~ T0, data=contAll3)
x=c(2:9)
lines(x, predict(fit5_contAll,newdata=data.frame(T0=x)),lwd=2,col='blue') # plot the model
plot(predict(fit5_contAll), resid(fit5_contAll), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit5_contAll)/sd(resid(fit5_contAll)))
abline(0,1) # not great at all





#===============================================================================
# Comparisons: -----------------------------------------------------------------
#===============================================================================

fun_comp <- function(x){
  akaike_comp <- as.data.frame(aictab(x))
  akaike_comp <- within(akaike_comp, rm("ModelLik", "AICcWt", "Cum.Wt"))
  bayesian_comp <- as.data.frame(bictab(x))  # this will pensalise complex models more strongly than AIC!
  bayesian_comp <- within(bayesian_comp, rm("K", "ModelLik", "BICWt", "Cum.Wt", "LL"))
  comparison <- left_join(akaike_comp, bayesian_comp, by="Modnames")
  comparison <- comparison %>% relocate(c(BIC, Delta_BIC), .before=LL)
  return(comparison)
}
#NOTE: remember you can compare lm & linearised power because F0 is log-transformed 
# (i.e., response variables are different)


#### Historical conservative ####
list_nls_hnc <- list("Linear"=fit_linear_nls_hnc_3,  # Linear
                     "Power"=fit_power_nls_hnc_3,    # Power
                     "Asympt"=fit4_hnc3,             # Asymptote
                     "Sigmoidal"=fit5_hnc)           # Sigmoidal
hnc_comp <- fun_comp(list_nls_hnc)
hnc_comp$Period <- "Historical Conservative"
##### NOTE: POWER. // ASYMPTOTE if considering 1s


#### Historical Broad ####
list_nls_hnb <- list("Linear"=fit_linear_nls_hnb_3, 
                     "Power"=fit_power_nls_hnb_3,
                     "Asympt"=fit4_hnb3, 
                     "Sigmoidal"=fit5_hnb)
hnb_comp <- fun_comp(list_nls_hnb)
hnb_comp$Period <- "Historical Broad"
##### NOTE: ASYMPTOTE. // ASYMPTOTE


#### Contemporary Native ####
list_nls_contN <- list("Linear"=fit_linear_nls_contN_3, 
                     "Power"=fit_power_nls_contN_3,
                     "Asympt"=fit4_contN3, 
                     "Sigmoidal"=fit5_contN)
contN_comp <- fun_comp(list_nls_contN)
contN_comp$Period <- "Contemporary Native"
##### NOTE: ASYMPTOTE. // ASYMPTOTE


#### Contemporary Native + Introduced ####
list_nls_contAll <- list("Linear"=fit_linear_nls_contAll_3,
                         "Power"=fit_power_nls_contAll_3,
                         "Asympt"=fit4_contAll3, 
                         "Sigmoidal"=fit5_contAll)
contAll_comp <- fun_comp(list_nls_contAll)
contAll_comp$Period <- "Contemporary Native + Introduced"
##### NOTE: ASYMPTOTE. // ASYMPTOTE


TableComp <- as.data.frame(rbind(hnc_comp, hnb_comp, contN_comp, contAll_comp))
TableComp <- TableComp %>% relocate(Period, .before = Modnames)

write.csv(TableComp, file="TableModelComparisonsMain.csv", row.names=F)




#===============================================================================
# Data Ggplots: ----------------------------------------------------------------
#===============================================================================
range(hnc3$T0)      # 2 17
range(hnb3$T0)
range(contN3$T0)    # 2 9
range(contAll3$T0)

historical_x <- c(2:17)
contemporary_x <- c(2:9)


# power:
pred_hnc <- as.data.frame(cbind("T0"=historical_x, predFit(fit_power_nls_hnc_3, data.frame("T0"=historical_x), interval="confidence")))
pred_hnc$Period <- "Historical Conservative"


# asymptotes:
pred_hnb <- as.data.frame(cbind("T0"=historical_x, predFit(fit4_hnb3, data.frame("T0"=historical_x), interval="confidence")))
pred_hnb$Period <- "Historical Broad"
pred_contN <- as.data.frame(cbind("T0"=contemporary_x, predFit(fit4_contN3, data.frame("T0"=contemporary_x), interval="confidence")))
pred_contN$Period <- "Contemporary Native"
pred_contAll <- as.data.frame(cbind("T0"=contemporary_x, predFit(fit4_contAll3, data.frame("T0"=contemporary_x), interval="confidence")))
pred_contAll$Period <- "Contemporary Native + Introduced"

pred_redundancy <- as.data.frame(rbind(pred_hnc, pred_hnb, pred_contN, pred_contAll))
pred_redundancy$Period <- factor(pred_redundancy$Period, levels = c("Historical Conservative",
                                                                    "Historical Broad",
                                                                    "Contemporary Native", 
                                                                    "Contemporary Native + Introduced"))        

save(pred_redundancy, file="pred_redundancy.RData")




#===============================================================================
# Plots: -----------------------------------------------------------------------
#===============================================================================
div3 <- div3[!div3$Period=="ContE",]   # rm Cont E from rawdata
div3$Period <- recode_factor(div3$Period, 
                            HNC = "Historical Conservative", 
                            HNB = "Historical Broad",
                            ContN = "Contemporary Native",
                            ContAll = "Contemporary Native + Introduced")
div3$Period <- factor(div3$Period, levels = c("Historical Conservative",
                                              "Historical Broad",
                                              "Contemporary Native", 
                                              "Contemporary Native + Introduced"))

div3_subset_main <- droplevels(subset(div3, div3$Period %in% c("Historical Broad", "Contemporary Native",
                                                 "Contemporary Native + Introduced")))
pred_subset_main <- subset(pred_redundancy, pred_redundancy$Period %in% c("Historical Broad", "Contemporary Native",
                                                                          "Contemporary Native + Introduced"))
colorBlind4   <- c("#28E2E5", "#61D04F", "#F5C710", "#CD0BBC")
colorBlind3   <- c("#61D04F", "#F5C710", "#CD0BBC")


## Free fit loess: -------------------------------------------------------------
(p1 <- ggplot(div3, aes(x=T0, y=F0, color=Period))+
  geom_point(aes(color=Period))+
  geom_smooth()+
  theme_classic()+
  scale_color_manual(values=colorBlind4)+
  facet_wrap(~Period)) # free fit
ggsave(p1, file= paste0(path_plots6, "/S5_F0vsT0FreeFit4.jpg"), width = 8, height = 6) 



## Plots (4 subsets): ----------------------------------------------------------
(p2 <- ggplot() + 
   geom_point(data=div3, aes(x = T0, y = F0, col=Period), 
              size = 2, alpha = 0.7, shape = 19) +
   geom_line(data = pred_redundancy, 
             aes(x = T0, y = fit, col=Period), 
             linetype = 1, size = 1) +
   geom_ribbon(data = pred_redundancy, 
               aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.1) +
   scale_color_manual(values=colorBlind4, 
                      breaks=c("Historical Conservative", "Historical Broad", "Contemporary Native", "Contemporary Native + Introduced"), 
                      labels=c("Historical conservative (n=48)", "Historical broad (n=60)", 
                               "Contemporary Native (n=33)", "Contemporary Native + Introduced (n=49)"))+ 
   theme_classic()+
   facet_wrap(~Period)
) 
ggsave(p2, file= paste0(path_plots6, "/S5_F0vsT0BestModel.jpg"), width = 8, height = 6) # for supplementary




#===============================================================================
# Plots All: -------------------------------------------------------------------
#===============================================================================

#plot(F0 ~ T0, data=hnc3, main="Historical Narrow")
#curve(predict(fit_linear_lm_3[[1]], newdata = data.frame(T0 = x)), col = "black", lty=1, add = TRUE)
#curve(predict(fit_power_nls_hnc_3, newdata = data.frame(T0 = x)), col = "red", lty=2, add = TRUE)
#x=c(2:17)
#lines(x, predict(fit4_hnc3,newdata=data.frame(T0=x)),lwd=1, lty=3, col='blue') # plot the model
#lines(x, predict(fit5_hnc,newdata=data.frame(T0=x)),lwd=1, lty=4, col='green') # plot the model
#legend("bottomright", legend = c("Linear", "Power", "Asymptote", "Logistic"),
#       col = c("black", "red", "blue", "green"),
#       lty = c(1,2,3,4), cex=0.5)

#plot(F0 ~ T0, data=hnb3, main="Historical Broad")
#curve(predict(fit_linear_lm_3[[2]], newdata = data.frame(T0 = x)), col = "black", lty=1, add = TRUE)
#curve(predict(fit_power_nls_hnb_3, newdata = data.frame(T0 = x)), col = "red", lty=2, add = TRUE)
#x=c(2:17)
#lines(x, predict(fit4_hnb3,newdata=data.frame(T0=x)),lwd=1, lty=3, col='blue') # plot the model
#lines(x, predict(fit5_hnb,newdata=data.frame(T0=x)),lwd=1, lty=4, col='green') # plot the model
#legend("bottomright", legend = c("Linear", "Power", "Asymptote", "Logistic"),
#       col = c("black", "red", "blue", "green"),
#       lty = c(1,2,3,4), cex=0.5)

#plot(F0 ~ T0, data=contN3, main="Contemporary Native")
#curve(predict(fit_linear_lm_3[[3]], newdata = data.frame(T0 = x)), col = "black", lty=1, add = TRUE)
#curve(predict(fit_power_nls_contN_3, newdata = data.frame(T0 = x)), col = "red", lty=2, add = TRUE)
#x=c(2:9)
#lines(x, predict(fit4_contN3,newdata=data.frame(T0=x)),lwd=1, lty=3, col='blue') # plot the model
#lines(x, predict(fit5_contN,newdata=data.frame(T0=x)),lwd=1, lty=4, col='green') # plot the model
#legend("bottomright", legend = c("Linear", "Power", "Asymptote", "Logistic"),
#       col = c("black", "red", "blue", "green"),
#       lty = c(1,2,3,4), cex=0.5)

#plot(F0 ~ T0, data=contAll3, main="Contemporary Native + Introduced")
#curve(predict(fit_linear_lm_3[[4]], newdata = data.frame(T0 = x)), col = "black", lty=1, add = TRUE)
#curve(predict(fit_power_nls_contAll_3, newdata = data.frame(T0 = x)), col = "red", lty=2, add = TRUE)
#x=c(2:9)
#lines(x, predict(fit4_contAll3,newdata=data.frame(T0=x)),lwd=1, lty=3, col='blue') # plot the model
#lines(x, predict(fit5_contAll,newdata=data.frame(T0=x)),lwd=1, lty=4, col='green') # plot the model
#legend("bottomright", legend = c("Linear", "Power", "Asymptote", "Logistic"),
#       col = c("black", "red", "blue", "green"),
#       lty = c(1,2,3,4), cex=0.5)


# End of script ################################################################

