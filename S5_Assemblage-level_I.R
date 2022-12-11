###########################################################################################
# Script: Model the relationship between SR and T0
# AFE
# August/Sept 2022
###########################################################################################

# Libraries:------------------------------------------------------------------------------
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

rm(list=ls())
myd <- getwd()
plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots" # Dir to save main plots

# Community data:-------------------------------------------------------------------------
load(paste0(myd, "/HNC.RData"))        # Conservative Historical Native Assemblage
load(paste0(myd, "/HNB.RData"))        # Broad Historical Native Assemblage
load(paste0(myd, "/ContN.RData"))      # Contemporary Native Assemblage
load(paste0(myd, "/ContE.RData"))      # Contemporary Exotic assemblage
load(paste0(myd, "/dist_mat1.RData"))  # Trait distance matrix
load(paste0(myd, "/ls.RData"))         # List of assemblages

# Functions:------------------------------------------------------------------------------
# FD and TD (Hill numbers) use the same units and TD is always greater than or equal to FD
# For q=0, Redundancy is S-FD0
hill_dir <- "C:/Users/Usuario/Documents/PHD/Functions/FunD-master/"

fskt <- read.table(paste0(hill_dir,"Plant_Abundance.txt"))
fskt = list(fs=fskt$Fushan, kt=fskt$Kenting)
dij_fskt = read.table(paste0(hill_dir,"Plant_Distance_Matrix.txt"))
dij_fskt = as.matrix(dij_fskt) # All needed for the code to run

source("C:/Users/Usuario/Documents/PHD/Functions/FunD_Rcode.R") # Warnings don't affect FD_MLE()

# Final formatting: -----------------------------------------------------------------------
setdiff(colnames(dist_mat1), names(ls[[1]]))

colnames(dist_mat1)[colnames(dist_mat1)=="Dajaus monticola"] <-"Agonostomus monticola"
rownames(dist_mat1)[rownames(dist_mat1)=="Dajaus monticola"] <-"Agonostomus monticola"

dist_mat1 <- dist_mat1[,order(colnames(dist_mat1))]
dist_mat1 <- dist_mat1[order(rownames(dist_mat1)),]

identical(colnames(dist_mat1), names(ls[[1]])) #TRUE

ls <- lapply(ls, function(x) {as.vector(as.matrix(x))})
lsS5 <- ls
save(lsS5, file="lsS5.RData")

##########################################################################################
# Observed TD and FD: --------------------------------------------------------------------
divT0 <- lapply(lsS5, function(x) {FD_MLE(x, dist_mat1, min(dist_mat1[dist_mat1>0]), 0)})  # taxonomic diversity
divF0 <- lapply(lsS5, function(x) {FD_MLE(x, dist_mat1, mean(dist_mat1[dist_mat1>0]), 0)}) # trait diversity

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

##########################################################################################
# Observed trends: -----------------------------------------------------------------------
div2 <- div[!div$T0==0,]   # 233
div3 <- div2[!div2$T0==1,] # 190 (main analysis)

save(div2, file = "div2.RData")
save(div3, file = "div3.RData")

par(mfrow=c(1,2))
hist(div2$T0)
hist(div2$F0)

par(mfrow=c(1,2))
hist(div3$T0)
hist(div3$F0) # More normally distributed

str(div3)
hnc2 <- subset(div2, div2$Period=="HNC")           # 67
hnb2 <- subset(div2, div2$Period=="HNB")           # 67
contN2 <- subset(div2, div2$Period=="ContN")       # 46
contAll2 <- subset(div2, div2$Period=="ContAll")   # 53

hnc3 <- subset(div3, div3$Period=="HNC")           # 48
hnb3 <- subset(div3, div3$Period=="HNB")           # 60
contN3 <- subset(div3, div3$Period=="ContN")       # 33
contAll3 <- subset(div3, div3$Period=="ContAll")   # 49

# Tests: ---------------------------------------------------------------------------------
# Use only sites in contN or cont All:
#hnc3 <- subset(hnc3, hnc3$SiteName %in% contN3$SiteName)
#hnb3 <- subset(hnb3, hnb3$SiteName %in% contN3$SiteName)

#hnc3 <- subset(hnc3, hnc3$SiteName %in% contAll3$SiteName)
#hnb3 <- subset(hnb3, hnb3$SiteName %in% contAll3$SiteName)

# Remove locality with richness = 17 in the historical period:
#hnc3 <- hnc3[!hnc3$T0 ==17,]
#hnb3 <- hnb3[!hnb3$T0 ==17,]

##########################################################################################
# Relationship between T0 & Period:-------------------------------------------------------
str(div2)
str(div3)
div3_1T0 <- lm(T0 ~ Period, data=div3)
summary(div3_1T0)
dev.off()

plot(T0 ~ Period, data=div3)

# Relationship between T0 & F0:-----------------------------------------------------------

l2 <- list(hnc2, hnb2, contN2, contAll2)
l3 <- list(hnc3, hnb3, contN3, contAll3)

save(l2, file = "l2.RData")
save(l3, file = "l3.RData")

##########################################################################################
# 1) Linear regression:-------------------------------------------------------------------

#fit1_2 <- lapply(l2, function(x) {lm(F0 ~ T0, data=x)})
#summlinear2 <- lapply(fit1_2, function(x) {summary(x)})
#print(summlinear2[[1]])

fit1_3 <- lapply(l3, function(x) {lm(F0 ~ T0, data=x)})
summlinear3 <- lapply(fit1_3, function(x) {summary(x)})
AIC(fit1_3[[4]])
BIC(fit1_3[[4]])

# Plotting & assumptions: ----------------------------------------------------------------
summlinear3[[1]]
plot(F0 ~ T0, data=hnc3, main="Historical Conservative")
curve(predict(fit1_3[[1]], newdata = data.frame(T0 = x)), col = "green", lty=1, add = TRUE)
par(mfrow=c(2,2))
plot(fit1_3[[1]])
dev.off()
plot(fit1_3[[1]]$fitted,fit1_3[[1]]$residuals)
plot(cooks.distance(fit1_3[[1]])) # Index 16 (assemblage with 17 sps)

summlinear3[[2]]
plot(F0 ~ T0, data=hnb3, main="Historical Broad")
curve(predict(fit1_3[[2]], newdata = data.frame(T0 = x)), col = "green", lty=1, add = TRUE)
par(mfrow=c(2,2))
plot(fit1_3[[2]])
dev.off()
plot(fit1_3[[1]]$fitted,fit1_3[[1]]$residuals)
plot(cooks.distance(fit1_3[[2]])) # Index 16 (assemblage with 17 sps)
hnb3[c(17,33),] # Sites with higher richness.

summlinear3[[3]]
plot(F0 ~ T0, data=hnb3, main="Contemporary natives")
curve(predict(fit1_3[[3]], newdata = data.frame(T0 = x)), col = "green", lty=1, add = TRUE)
par(mfrow=c(2,2))
plot(fit1_3[[3]])
dev.off()
#plot(fit1_3[[1]]$fitted,fit1_3[[1]]$residuals)
plot(cooks.distance(fit1_3[[3]])) # Index 16 (assemblage with 17 sps)
contN3[c(19,32),] # Sites with relatively higher richness.

summlinear3[[4]]
plot(F0 ~ T0, data=hnb3, main="Contemporary natives + exotics")
curve(predict(fit1_3[[4]], newdata = data.frame(T0 = x)), col = "green", lty=1, add = TRUE)
par(mfrow=c(2,2))
plot(fit1_3[[4]])
dev.off()
#plot(fit1_3[[1]]$fitted,fit1_3[[1]]$residuals)
plot(cooks.distance(fit1_3[[4]])) # Index 16 (assemblage with 17 sps)
contAll3[c(18,26,47,49),] 

# Notes: Linear models don't seem to be a great fit
# especially for larger values of SR.
# Outliers tend to be at large SR's.


##########################################################################################
# 2) Power regression (log-log & nls):----------------------------------------------------
lapply(l3, function(x) {range(x$F0)}) # No need to do log + 1

#fit2_2 <- lapply(l2, function(x) {lm(log(F0+1) ~ log(T0+1), data=x)})
#summloglog2 <- lapply(fit2_2, function(x) {summary(x)})

fit2_3 <- lapply(l3, function(x) {lm(log(F0) ~ log(T0), data=x)})
summloglog3 <- lapply(fit2_3, function(x) {summary(x)})
AIC(fit2_3[[4]])
BIC(fit2_3[[4]])

# Power regression nls: (different error family, additive):-------------------------------

#fit3_2hnc <- nls((F0+1) ~ a * (T0+1)^b, data = l2[[1]], 
#                 start = list(a = exp(coef(fit2_2[[1]])[1]), b = coef(fit2_2[[1]])[2]))

#fit3_2hnb <- nls((F0+1) ~ a * (T0+1)^b, data = l2[[2]], 
#                start = list(a = exp(coef(fit2_2[[2]])[1]), b = coef(fit2_2[[2]])[2]))

#fit3_2contN <- nls((F0+1) ~ a * (T0+1)^b, data = l2[[3]], 
#                 start = list(a = exp(coef(fit2_2[[3]])[1]), b = coef(fit2_2[[3]])[2]))

#fit3_2contAll <- nls((F0+1) ~ a * (T0+1)^b, data = l2[[4]], 
#                 start = list(a = exp(coef(fit2_2[[4]])[1]), b = coef(fit2_2[[4]])[2]))

fit3_3hnc <- nls(F0 ~ a * T0^b, data = l3[[1]], 
                 start = list(a = exp(coef(fit2_3[[1]])[1]), b = coef(fit2_3[[1]])[2]))
fit3_3hnb <- nls(F0 ~ a * T0^b, data = l3[[2]], 
                 start = list(a = exp(coef(fit2_3[[2]])[1]), b = coef(fit2_3[[2]])[2]))
fit3_3contN <- nls(F0 ~ a * T0^b, data = l3[[3]], 
                 start = list(a = exp(coef(fit2_3[[3]])[1]), b = coef(fit2_3[[3]])[2]))
fit3_3contAll <- nls(F0 ~ a * T0^b, data = l3[[4]], 
                 start = list(a = exp(coef(fit2_3[[4]])[1]), b = coef(fit2_3[[4]])[2]))

# Plotting & assumptions: ----------------------------------------------------------------

plot(F0 ~ T0, data=hnc3, main="Historical Conservative")
curve(exp(predict(fit2_3[[1]], newdata = data.frame(T0 = x))), col = "black", lty=2, add = TRUE)
curve(predict(fit3_3hnc, newdata = data.frame(T0 = x)), col = "red", lty=3, add = TRUE)
legend("bottomright", legend = c("Log-Log", "Nls"),
       col = c("black", "red"),
       lty = c(2,3), cex=0.5)

par(mfrow=c(2,2))
plot(fit2_3[[1]]) # Cooks distance large max 0.2
dev.off()
plot(predict(fit3_3hnc), resid(fit3_3hnc), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit3_3hnc)/sd(resid(fit3_3hnc)))
abline(0,1) #OK (the power model seems good)


plot(F0 ~ T0, data=hnb3, main="Historical Broad")
curve(exp(predict(fit2_3[[2]], newdata = data.frame(T0 = x))), col = "black", lty=2, add = TRUE)
curve(predict(fit3_3hnb, newdata = data.frame(T0 = x)), col = "red", lty=3, add = TRUE)
legend("bottomright", legend = c("Log-Log", "Nls"),
       col = c("black", "red"),
       lty = c(2,3), cex=0.5)

par(mfrow=c(2,2))
plot(fit2_3[[2]]) # Cooks distance large max 0.14
dev.off()
plot(predict(fit3_3hnb), resid(fit3_3hnb), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit3_3hnb)/sd(resid(fit3_3hnb)))
abline(0,1) #OK (the power model seems good)

plot(F0 ~ T0, data=contN3)
curve(exp(predict(fit2_3[[3]], newdata = data.frame(T0 = x))), col = "black", lty=2, add = TRUE)
curve(predict(fit3_3contN, newdata = data.frame(T0 = x)), col = "red", lty=3, add = TRUE)
legend("bottomright", legend = c("Log-Log", "Nls"),
       col = c("black", "red"),
       lty = c(2, 3), cex=0.5)

par(mfrow=c(2,2))
plot(fit2_3[[3]]) # Cooks distance large max 0.20
dev.off()
plot(predict(fit3_3contN), resid(fit3_3contN), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit3_3contN)/sd(resid(fit3_3contN)))
abline(0,1) #OK (less good than in the historical period)


plot(F0 ~ T0, data=contAll3)
curve(exp(predict(fit2_3[[4]], newdata = data.frame(T0 = x))), col = "black", lty=2, add = TRUE)
curve(predict(fit3_3contAll, newdata = data.frame(T0 = x)), col = "red", lty=3, add = TRUE)
legend("bottomright", legend = c("Log-Log", "Nls"),
       col = c("black", "red"),
       lty = c(2, 3), cex=0.5)

par(mfrow=c(2,2))
plot(fit2_3[[4]]) # Cooks distance large max 0.11
dev.off()
plot(predict(fit3_3contAll), resid(fit3_3contAll), xlab="Fitted values", ylab="Residuals")
abline(h=0) # At low & high values we observe all dots under 0
qqnorm(resid(fit3_3contAll)/sd(resid(fit3_3contAll)))
abline(0,1) #OK (less good than in the historical period) /  over-predicting?



##########################################################################################
# 3) Asymptote (nls):---------------------------------------------------------------------
#  When x = 0, y = 0 When x is huge, y ∼ k
plot(F0 ~ T0, data=hnc3)
x=c(2:17)
my.k <- 6  # parameter for model
my.b <- 1  # parameter for model

fit4_hnc3 <- nls(F0~k*(1-exp((-1)*b*T0)), data=hnc3, start=list(k=my.k,b=my.b))
lines(x, predict(fit4_hnc3,newdata=data.frame(T0=x)),lwd=2,col='blue') # plot the model
summary(fit4_hnc3) # Parameters seem to be in the correct range.


plot(F0 ~ T0, data=hnb3)
x=c(2:17)
my.k <- 4.5  # parameter for model
my.b <- 1    # parameter for model

fit4_hnb3 <- nls(F0~k*(1-exp((-1)*b*T0)), data=hnb3, start=list(k=my.k,b=my.b))
lines(x, predict(fit4_hnb3,newdata=data.frame(T0=x)),lwd=2,col='blue') # plot the model
summary(fit4_hnb3) # Parameters seem to be in the correct range.


plot(F0 ~ T0, data=contN3)
x=c(2:9)
my.k <- 3  # parameter for model
my.b <- 1    # parameter for model

fit4_contN3 <- nls(F0~k*(1-exp((-1)*b*T0)), data=contN3, start=list(k=my.k,b=my.b))
lines(x, predict(fit4_contN3,newdata=data.frame(T0=x)),lwd=2,col='blue') # plot the model
summary(fit4_contN3) # Parameters seem to be in the correct range.


plot(F0 ~ T0, data=contAll3)
x=c(2:9)
my.k <- 3  # parameter for model
my.b <- 1    # parameter for model

fit4_contAll3 <- nls(F0~k*(1-exp((-1)*b*T0)), data=contAll3, start=list(k=my.k,b=my.b))
lines(x, predict(fit4_contAll3,newdata=data.frame(T0=x)),lwd=2,col='blue') # plot the model
summary(fit4_contAll3) # Parameters seem to be in the correct range.

# Assumptions: ---------------------------------------------------------------------------
dev.off()
plot(predict(fit4_hnc3), resid(fit4_hnc3), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit4_hnc3)/sd(resid(fit4_hnc3)))
abline(0,1) # pretty good

plot(predict(fit4_hnb3), resid(fit4_hnb3), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit4_hnb3)/sd(resid(fit4_hnb3)))
abline(0,1) # not that good

plot(predict(fit4_contN3), resid(fit4_contN3), xlab="Fitted values", ylab="Residuals")
abline(h=0)
qqnorm(resid(fit4_contN3)/sd(resid(fit4_contN3)))
abline(0,1) # not great

plot(predict(fit4_contAll3), resid(fit4_contAll3), xlab="Fitted values", ylab="Residuals")
abline(h=0) # At low & high values we observe all dots under 0
qqnorm(resid(fit4_contAll3)/sd(resid(fit4_contAll3)))
abline(0,1) #OK (not great) /  over-predicting?

##########################################################################################
# 4) Logistic model: ---------------------------------------------------------------------
# SSlogis == y~a/(1 + exp(-b * (x-c))
fit5_hnc <- nls(F0 ~ SSlogis(T0, Asym, xmid, scal), data = hnc3)
fit5_hnb <- nls(F0 ~ SSlogis(T0, Asym, xmid, scal), data = hnb3)
fit5_contN <- nls(F0 ~ SSlogis(T0, Asym, xmid, scal), data = contN3)
fit5_contAll <- nls(F0 ~ SSlogis(T0, Asym, xmid, scal), data = contAll3)

# Plots & Assumptions: -------------------------------------------------------------------
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

##########################################################################################
# Comparisons: ---------------------------------------------------------------------------

# Historical conservative:----------------------------------------------------------------
sum(resid(fit1_3[[1]])^2) # 16.7281
sum(resid(fit3_3hnc)^2)   # 14.92375
sum(resid(fit4_hnc3)^2)   # 14.90888
sum(resid(fit5_hnc)^2)    # 14.85297

AIC(fit1_3[[1]]) # 91.62078
AIC(fit3_3hnc)   # 86.14224*
AIC(fit4_hnc3)   # 86.0944*
AIC(fit5_hnc)    # 87.91405

BIC(fit1_3[[1]]) # 97.23438
BIC(fit3_3hnc)   # 91.75585*
BIC(fit4_hnc3)   # 91.708*
BIC(fit5_hnc)    # 95.39885

# Both power & asymptote.
# Tests:
# 1) When considering only the contN sites --> asymptote or logistic.
# 2) When considering only the contAll sites --> power or asymptotic.
# 3) When removing the obs of 17: asymptotic

# Historical broad:-----------------------------------------------------------------------
sum(resid(fit1_3[[2]])^2) # 28.91686
sum(resid(fit3_3hnb)^2)   # 25.38108
sum(resid(fit4_hnb3)^2)   # 24.24942*
sum(resid(fit5_hnb)^2)    # 24.45509

AIC(fit1_3[[2]])          # 132.4774
AIC(fit3_3hnb)            # 124.6522
AIC(fit4_hnb3)            # 121.9155*
AIC(fit5_hnb)             # 124.4222

BIC(fit1_3[[2]])          # 138.7605
BIC(fit3_3hnb)            # 130.9352
BIC(fit4_hnb3)            # 128.1985*
BIC(fit5_hnb)             # 132.7996
# Asymptote.

# Tests:
# 1) When considering only the contN sites --> asymptote
# 2) When considering only the contAll sites --> asymptote.
# 3) When removing the obs of 17: asymptotic

# Contemporary natives:-------------------------------------------------------------------
sum(resid(fit1_3[[3]])^2)   # 10.13201
sum(resid(fit3_3contN)^2)   # 9.553247
sum(resid(fit4_contN3)^2)   # 9.244473
sum(resid(fit5_contN)^2)    # 9.276654

AIC(fit1_3[[3]]) # 60.68329
AIC(fit3_3contN) # 58.74227*
AIC(fit4_contN3) # 57.65805*
AIC(fit5_contN)  # 59.77273

BIC(fit1_3[[3]]) # 65.17281
BIC(fit3_3contN) # 63.23179*
BIC(fit4_contN3) # 62.14757*
BIC(fit5_contN)  # 65.75876
# Both power & asymptote.


# Contemporary natives + exotics:-----------------------------------------------------------
sum(resid(fit1_3[[4]])^2)     # 16.60955
sum(resid(fit3_3contAll)^2)   # 14.71436
sum(resid(fit4_contAll3)^2)   # 12.71858
sum(resid(fit5_contAll)^2)    # 12.07803

AIC(fit1_3[[4]])              # 92.0457
AIC(fit3_3contAll)            # 86.10915
AIC(fit4_contAll3)            # 78.96692*
AIC(fit5_contAll)             # 78.43479*
AIC(fit2_3[[4]])

BIC(fit1_3[[4]])              # 97.72117
BIC(fit3_3contAll)            # 91.78461
BIC(fit4_contAll3)            # 84.64238*
BIC(fit5_contAll)             # 86.00207*
# Asymptote & logistic


##########################################################################################
# Generate data frames for final plotting with ggplot: -----------------------------------

x_hnc <- data.frame("T0" = seq(from = min(hnc3$T0), to = max(hnc3$T0), length= 1000))
x_hnb <- data.frame("T0" = seq(from = min(hnb3$T0), to = max(hnb3$T0), length= 1000))
x_contN <- data.frame("T0" = seq(from = min(contN3$T0), to = max(contN3$T0), length= 1000))
x_contAll <- data.frame("T0" = seq(from = min(contAll3$T0), to = max(contAll3$T0), length= 1000))


hnc_f <- data.frame(x_hnc, exp(predict(fit2_3[[1]], x_hnc, interval = "confidence")))
hnc_f$Period <- rep("Historical Conservative", nrow(hnc_f))

hnb_f <- data.frame(x_hnb, exp(predict(fit2_3[[2]], x_hnb, interval = "confidence")))
hnb_f$Period <- rep("Historical Broad", nrow(hnb_f))

contN_f <- data.frame(x_contN, exp(predict(fit2_3[[3]], x_contN, interval = "confidence")))
contN_f$Period <- rep("Contemporary Natives", nrow(contN_f))

contAll_f <- data.frame(x_contAll, exp(predict(fit2_3[[4]], x_contAll, interval = "confidence")))
contAll_f$Period <- rep("Contemporary Natives + Exotics", nrow(contAll_f))

save(hnc_f, file = "hnc_f.RData")
save(hnb_f, file = "hnb_f.RData")
save(contN_f, file = "contN_f.RData")
save(contAll_f, file = "contAll_f.RData")


##########################################################################################
# Plot trends (final ggplot):-------------------------------------------------------------
div3$Period <- recode_factor(div3$Period, 
                            HNC = "Historical Conservative", 
                            HNB = "Historical Broad",
                            ContN = "Contemporary Natives",
                            ContAll = "Contemporary Natives + Exotics")

div3_subset_main <- droplevels(subset(div3, div3$Period %in% c("Historical Broad", "Contemporary Natives",
                                                 "Contemporary Natives + Exotics")))

colorBlind4   <- c("#F0E442", "#E69F00", "#009E73","#0072B2")
colorBlind3   <- c("#E69F00", "#009E73","#0072B2")

(p1 <- ggplot(div, aes(x=T0, y=F0, color=Period))+
  geom_point(aes(color=Period))+
  geom_smooth()+
  theme_classic()+
  scale_color_manual(values=colorBlind4)+
  facet_wrap(~Period)) # free fit
ggsave(p1, file= paste0(plot_dir, "/SM/Thesis/S5_F0vsT0FreeFit4.jpg"), width = 8, height = 6) 

(p2 <- ggplot() + 
  geom_point(data=div3, aes(x = T0, y = F0, col=Period), 
              size = 2, alpha = 0.7, shape = 19) +
  geom_line(data = hnc_f , aes(x = T0, y = fit), linetype = 1, size = 1, col="#F0E442") + 
  #geom_ribbon(data = hnc_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3)+ 
  
  geom_line(data = hnb_f , aes(x = T0, y = fit), linetype = 1, size = 1, col="#E69F00") + 
  #geom_ribbon(data = hnb_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3)+ 
  
  geom_line(data = contN_f , aes(x = T0, y = fit), linetype = 1, size = 1, col="#009E73") +
  #geom_ribbon(data = contN_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3)+ 
  
  geom_line(data = contAll_f , aes(x = T0, y = fit), linetype = 1, size = 1, col="#0072B2") +
  #geom_ribbon(data = contAll_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3)+ 
  labs(x = "Species Richness", y = "Trait diversity")+
  scale_color_manual(values=colorBlind4, 
                     breaks=c("Historical conservative", "Historical broad", "Contemporary Natives", "Contemporary Natives + Exotics"), 
                     labels=c("Historical conservative (n=48)", "Historical broad (n=60)", 
                              "Contemporary Natives (n=33)", "Contemporary Natives + Exotics (n=49)"))+
  theme_bw())

#ggsave(p2, file= paste0(plot_dir, "/SM/Thesis/S5_F0vsT0Power4.jpg"), width = 8, height = 6) 

(p3 <- ggplot() + 
    geom_point(data=div3[div3$Period=="Historical Conservative",], aes(x = T0, y = F0), 
               size = 2, shape = 19, col="#F0E442") +
    geom_line(data = hnc_f , aes(x = T0, y = fit), linetype = 1, size = 1, col="#F0E442") + 
    geom_ribbon(data = hnc_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3)+ 
    labs(x = "Species Richness", y = "Trait diversity")+
    theme_bw())

ggsave(p3, file= paste0(plot_dir, "/SM/Ms/S5_F0vsT0HNCPower.jpg"), width = 8, height = 6) 

df_fit <- rbind(hnb_f, contN_f, contAll_f)

(p4 <- ggplot() + 
    geom_point(data=div3_subset_main, aes(x = T0, y = F0, col=Period), 
              size = 2, alpha = 0.7, shape = 19) +
    geom_line(data = df_fit , aes(x = T0, y = fit, col=Period), linetype = 1, size = 1) + 
    geom_ribbon(data = df_fit, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3)+ 
    labs(x = "Species Richness", y = "Trait diversity")+
    scale_color_manual(values=colorBlind3, 
                       breaks=c("Historical Broad", "Contemporary Natives", "Contemporary Natives + Exotics"), 
                       labels=c("Historical Broad (n=60)","Contemporary Natives (n=33)", "Contemporary Natives + Exotics (n=49)"))+
   facet_wrap(~Period, ncol=2)+
   theme_bw())
ggsave(p4, file= paste0(plot_dir, "/SM/Ms/S5_F0vsT0Power3.jpg"), width = 8, height = 6) 

##########################################################################################
# Plots SM: ------------------------------------------------------------------------------
#par(mfrow=c(2,2))

plot(F0 ~ T0, data=hnc3, main="Historical Narrow")
curve(predict(fit1_3[[1]], newdata = data.frame(T0 = x)), col = "black", lty=1, add = TRUE)
curve(predict(fit3_3hnc, newdata = data.frame(T0 = x)), col = "red", lty=2, add = TRUE)
x=c(2:17)
lines(x, predict(fit4_hnc3,newdata=data.frame(T0=x)),lwd=1, lty=3, col='blue') # plot the model
lines(x, predict(fit5_hnc,newdata=data.frame(T0=x)),lwd=1, lty=4, col='green') # plot the model
legend("bottomright", legend = c("Linear", "Power", "Asymptote", "Logistic"),
       col = c("black", "red", "blue", "green"),
       lty = c(1,2,3,4), cex=0.5)

plot(F0 ~ T0, data=hnb3, main="Historical Broad")
curve(predict(fit1_3[[2]], newdata = data.frame(T0 = x)), col = "black", lty=1, add = TRUE)
curve(predict(fit3_3hnb, newdata = data.frame(T0 = x)), col = "red", lty=2, add = TRUE)
x=c(2:17)
lines(x, predict(fit4_hnb3,newdata=data.frame(T0=x)),lwd=1, lty=3, col='blue') # plot the model
lines(x, predict(fit5_hnb,newdata=data.frame(T0=x)),lwd=1, lty=4, col='green') # plot the model
legend("bottomright", legend = c("Linear", "Power", "Asymptote", "Logistic"),
       col = c("black", "red", "blue", "green"),
       lty = c(1,2,3,4), cex=0.5)

plot(F0 ~ T0, data=contN3, main="Contemporary Natives")
curve(predict(fit1_3[[3]], newdata = data.frame(T0 = x)), col = "black", lty=1, add = TRUE)
curve(predict(fit3_3contN, newdata = data.frame(T0 = x)), col = "red", lty=2, add = TRUE)
x=c(2:9)
lines(x, predict(fit4_contN3,newdata=data.frame(T0=x)),lwd=1, lty=3, col='blue') # plot the model
lines(x, predict(fit5_contN,newdata=data.frame(T0=x)),lwd=1, lty=4, col='green') # plot the model
legend("bottomright", legend = c("Linear", "Power", "Asymptote", "Logistic"),
       col = c("black", "red", "blue", "green"),
       lty = c(1,2,3,4), cex=0.5)

plot(F0 ~ T0, data=contAll3, main="Contemporary Natives + Exotics")
curve(predict(fit1_3[[4]], newdata = data.frame(T0 = x)), col = "black", lty=1, add = TRUE)
curve(predict(fit3_3contAll, newdata = data.frame(T0 = x)), col = "red", lty=2, add = TRUE)
x=c(2:9)
lines(x, predict(fit4_contAll3,newdata=data.frame(T0=x)),lwd=1, lty=3, col='blue') # plot the model
lines(x, predict(fit5_contAll,newdata=data.frame(T0=x)),lwd=1, lty=4, col='green') # plot the model
legend("bottomright", legend = c("Linear", "Power", "Asymptote", "Logistic"),
       col = c("black", "red", "blue", "green"),
       lty = c(1,2,3,4), cex=0.5)



##########################################################################################
# End of script ##########################################################################

