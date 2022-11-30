###########################################################################################
# Script: Observed vs null trait diversity & redundancy and temporal change
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

rm(list=ls())
myd <- getwd()
plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots" # Dir to save main plots

# Community data:-------------------------------------------------------------------------
load(paste0(myd, "/HNC.RData"))        # Conservative Historical Native Assemblage
load(paste0(myd, "/HNB.RData"))        # Broad Historical Native Assemblage
load(paste0(myd, "/ContN.RData"))      # Contemporary Native Assemblage
load(paste0(myd, "/ContE.RData"))      # Contemporary Exotic assemblage
load(paste0(myd, "/dist_mat1.RData"))  # Trait distance matrix

# Functions:------------------------------------------------------------------------------
# FD and TD (Hill numbers) use the same units and TD is always greater than or equal to FD
# For q=0, Redundancy is S-FD0
hill_dir <- "C:/Users/Usuario/Documents/PHD/Functions/FunD-master/"

fskt <- read.table(paste0(hill_dir,"Plant_Abundance.txt"))
fskt = list(fs=fskt$Fushan, kt=fskt$Kenting)
dij_fskt = read.table(paste0(hill_dir,"Plant_Distance_Matrix.txt"))
dij_fskt = as.matrix(dij_fskt) # All needed for the code to run

source("C:/Users/Usuario/Documents/PHD/Functions/FunD_Rcode.R") # Warnings don't affect FD_MLE()

# Data formatting:------------------------------------------------------------------------
ContAll <- full_join(ContN, ContE, by=c("SiteNameE", "DrainageBasinE"))
names(ContAll)

ContAll$`Goodea atripinnis`<- ContAll$`Goodea atripinnis.x`+ ContAll$`Goodea atripinnis.y`
ContAll$`Poecilia butleri`<- ContAll$`Poecilia butleri.x`+ ContAll$`Poecilia butleri.y`
ContAll$`Poecilia sphenops`<- ContAll$`Poecilia sphenops.x`+ ContAll$`Poecilia sphenops.y`
ContAll <- ContAll[, ! names(ContAll) %in% c("Goodea atripinnis.x", "Goodea atripinnis.y",
                                            "Poecilia butleri.x", "Poecilia butleri.y",
                                            "Poecilia sphenops.x", "Poecilia sphenops.y")]

sum(is.na(ContAll[,-c(1:2)]))
ContAll[is.na(ContAll)] <- 0
sum(is.na(ContAll[,-c(1:2)]))
save(ContAll, file="ContAll.RData")

l <- bind_rows(HNC, HNB, ContN, ContAll)
l[is.na(l)] <- 0
sum(is.na(l)) # 0

l$Period <- rep(c("HNC", "HNB", "ContN", "ContAll"), each=67)
All  <- l
save(All, file="All.RData")
rownames(l) <- paste0(l$SiteNameE, "_", l$Period)
l <- within(l, rm(SiteNameE, DrainageBasinE, Period))
l <- l[, order(names(l))]

str(l)

ls <- split(l, rownames(l)) # 67*4

setdiff(colnames(dist_mat1), names(ls[[1]]))

colnames(dist_mat1)[colnames(dist_mat1)=="Dajaus monticola"] <-"Agonostomus monticola"
rownames(dist_mat1)[rownames(dist_mat1)=="Dajaus monticola"] <-"Agonostomus monticola"

dist_mat1 <- dist_mat1[,order(colnames(dist_mat1))]
dist_mat1 <- dist_mat1[order(rownames(dist_mat1)),]

identical(colnames(dist_mat1), names(ls[[1]])) #TRUE

ls <- lapply(ls, function(x) {as.vector(as.matrix(x))})
save(ls, file="ls.RData")

##########################################################################################
# Observed TD and FD: --------------------------------------------------------------------
divT0 <- lapply(ls, function(x) {FD_MLE(x, dist_mat1, min(dist_mat1[dist_mat1>0]), 0)})  # taxonomic diversity
divF0 <- lapply(ls, function(x) {FD_MLE(x, dist_mat1, mean(dist_mat1[dist_mat1>0]), 0)}) # trait diversity

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
# Subsets:
div2 <- div[!div$T0==0,]   # 233
div3 <- div2[!div2$T0==1,] # 190

save(div2, file = "div2.RData")
save(div3, file = "div3.RData")

par(mfrow=c(1,2))
hist(div2$T0)
hist(div2$F0)

par(mfrow=c(1,2))
hist(div3$T0)
hist(div3$F0)

str(div3)
hnc2 <- subset(div2, div2$Period=="HNC")           # 67
hnb2 <- subset(div2, div2$Period=="HNB")           # 67
contN2 <- subset(div2, div2$Period=="ContN")       # 46
contAll2 <- subset(div2, div2$Period=="ContAll")   # 53

hnc3 <- subset(div3, div3$Period=="HNC")           # 48
hnb3 <- subset(div3, div3$Period=="HNB")           # 60
contN3 <- subset(div3, div3$Period=="ContN")       # 33
contAll3 <- subset(div3, div3$Period=="ContAll")   # 49

# Relationship between T0 & Period:-------------------------------------------------------
str(div2)
str(div3)
div3_1T0 <- lm(T0 ~ Period, data=div3)
summary(div3_1T0)
dev.off()

plot(T0 ~ Period, data=div3)

div3B <- droplevels(subset(div3, div3$Period %in% c("HNB", "ContAll")))
div3_2T0 <- lm(T0 ~ Period, data=div3B)
summary(div3_2T0)
plot(T0 ~ Period, data=div3B)

##########################################################################################
# Relationship between T0 & F0:-----------------------------------------------------------

l2 <- list(hnc2, hnb2, contN2, contAll2)
l3 <- list(hnc3, hnb3, contN3, contAll3)

# Linear regression:----------------------------------------------------------------------

fit1_2 <- lapply(l2, function(x) {lm(F0 ~ T0, data=x)})
summlinear2 <- lapply(fit1_2, function(x) {summary(x)})
print(summlinear2[[1]])

fit1_3 <- lapply(l3, function(x) {lm(F0 ~ T0, data=x)})
summlinear3 <- lapply(fit1_3, function(x) {summary(x)})

fit1_3[[1]]$df.residual
print(summlinear3[[1]])

# Power regression (log-log):-------------------------------------------------------------
lapply(l3, function(x) {range(x$F0)}) # No need to do log + 1

fit2_2 <- lapply(l2, function(x) {lm(log(F0+1) ~ log(T0+1), data=x)})
summloglog2 <- lapply(fit2_2, function(x) {summary(x)})

fit2_3 <- lapply(l3, function(x) {lm(log(F0) ~ log(T0), data=x)})
summloglog3 <- lapply(fit2_3, function(x) {summary(x)})

fit2_3[[1]]$df.residual
print(summloglog3[[1]])

# Power regression nls: (different error family, additive):-------------------------------

fit3_2hnc <- nls((F0+1) ~ a * (T0+1)^b, data = l2[[1]], 
                 start = list(a = exp(coef(fit2_2[[1]])[1]), b = coef(fit2_2[[1]])[2]))

fit3_2hnb <- nls((F0+1) ~ a * (T0+1)^b, data = l2[[2]], 
                 start = list(a = exp(coef(fit2_2[[2]])[1]), b = coef(fit2_2[[2]])[2]))

fit3_2contN <- nls((F0+1) ~ a * (T0+1)^b, data = l2[[3]], 
                 start = list(a = exp(coef(fit2_2[[3]])[1]), b = coef(fit2_2[[3]])[2]))

fit3_2contAll <- nls((F0+1) ~ a * (T0+1)^b, data = l2[[4]], 
                 start = list(a = exp(coef(fit2_2[[4]])[1]), b = coef(fit2_2[[4]])[2]))


fit3_3hnc <- nls(F0 ~ a * T0^b, data = l3[[1]], 
                 start = list(a = exp(coef(fit2_3[[1]])[1]), b = coef(fit2_3[[1]])[2]))
fit3_3hnb <- nls(F0 ~ a * T0^b, data = l3[[2]], 
                 start = list(a = exp(coef(fit2_3[[2]])[1]), b = coef(fit2_3[[2]])[2]))
fit3_3contN <- nls(F0 ~ a * T0^b, data = l3[[3]], 
                 start = list(a = exp(coef(fit2_3[[3]])[1]), b = coef(fit2_3[[3]])[2]))
fit3_3contAll <- nls(F0 ~ a * T0^b, data = l3[[4]], 
                 start = list(a = exp(coef(fit2_3[[4]])[1]), b = coef(fit2_3[[4]])[2]))

##########################################################################################
# Plotting:-------------------------------------------------------------------------------
# Historical conservative data:-----------------------------------------------------------
plot(F0 ~ T0, data=hnc3, main="Historical Conservative")
curve(predict(fit1_3[[1]], newdata = data.frame(T0 = x)), col = "green", lty=1, add = TRUE)
curve(exp(predict(fit2_3[[1]], newdata = data.frame(T0 = x))), col = "black", lty=2, add = TRUE)
curve(predict(fit3_3hnc, newdata = data.frame(T0 = x)), col = "red", lty=3, add = TRUE)
legend("bottomright", legend = c("Linear", "Log-Log", "Nls"),
       col = c("green", "black", "red"),
       lty = c(1, 2, 3), cex=0.5)

# Historical broad data:-----------------------------------------------------------------
plot(F0 ~ T0, data=hnb3)
curve(predict(fit1_3[[2]], newdata = data.frame(T0 = x)), col = "green", lty=1, add = TRUE)
curve(exp(predict(fit2_3[[2]], newdata = data.frame(T0 = x))), col = "black", lty=2, add = TRUE)
curve(predict(fit3_3hnb, newdata = data.frame(T0 = x)), col = "red", lty=3, add = TRUE)
legend("bottomright", legend = c("Linear", "Log-Log", "Nls"),
       col = c("green", "black", "red"),
       lty = c(1, 2, 3), cex=0.5)

# Contemporary native data:--------------------------------------------------------------
plot(F0 ~ T0, data=contN3)
curve(predict(fit1_3[[3]], newdata = data.frame(T0 = x)), col = "green", lty=1, add = TRUE)
curve(exp(predict(fit2_3[[3]], newdata = data.frame(T0 = x))), col = "black", lty=2, add = TRUE)
curve(predict(fit3_3contN, newdata = data.frame(T0 = x)), col = "red", lty=3, add = TRUE)
legend("bottomright", legend = c("Linear", "Log-Log", "Nls"),
       col = c("green", "black", "red"),
       lty = c(1, 2, 3), cex=0.5)

# Contemporary native + exotic data:------------------------------------------------------
plot(F0 ~ T0, data=contAll3)
curve(predict(fit1_3[[4]], newdata = data.frame(T0 = x)), col = "green", lty=1, add = TRUE)
curve(exp(predict(fit2_3[[4]], newdata = data.frame(T0 = x))), col = "black", lty=2, add = TRUE)
curve(predict(fit3_3contAll, newdata = data.frame(T0 = x)), col = "red", lty=3, add = TRUE)
legend("bottomright", legend = c("Linear", "Log-Log", "Nls"),
       col = c("green", "black", "red"),
       lty = c(1, 2, 3), cex=0.5)

##########################################################################################
# Assumptions:----------------------------------------------------------------------------
# HNC:
par(mfrow=c(2,2))
plot(fit1_3[[1]])

par(mfrow=c(2,2))
plot(fit2_3[[1]])

# HNB:
par(mfrow=c(2,2))
plot(fit1_3[[2]])

par(mfrow=c(2,2))
plot(fit2_3[[2]])

# Cont N:
par(mfrow=c(2,2))
plot(fit1_3[[3]])

par(mfrow=c(2,2))
plot(fit2_3[[3]])

# Cont All:
par(mfrow=c(2,2))
plot(fit1_3[[4]])

par(mfrow=c(2,2))
plot(fit2_3[[4]])


##########################################################################################
#library(lmtest)
#lrtest(fit1_3[[1]], fit2_3[[1]], fit3_3hnc)
#lrtest(fit1_3[[2]], fit2_3[[2]], fit3_3hnb)
#lrtest(fit1_3[[3]], fit2_3[[3]], fit3_3contN)
#lrtest(fit1_3[[4]], fit2_3[[4]], fit3_3contAll)

# AIC ? / BIC ?

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

##########################################################################################
# Tests: ---------------------------------------------------------------------------------
# Test the fit in the historical period with the same "contemporary" site subsets.
hnc_cN <- subset(hnc3, hnc3$SiteName %in% contN3$SiteName)
hnb_cN <- subset(hnb3, hnb3$SiteName %in% contN3$SiteName)

hnc_cAll <- subset(hnc3, hnc3$SiteName %in% contAll3$SiteName)
hnb_cAll <- subset(hnb3, hnb3$SiteName %in% contAll3$SiteName)

m1 <- lm(F0 ~ T0, data= hnc_cN)
m2 <- lm(log(F0) ~ log(T0), data= hnc_cN)
summary(m1)
summary(m2)

m3 <- lm(F0 ~ T0, data= hnc_cAll)
m4 <- lm(log(F0) ~ log(T0), data= hnc_cAll)
summary(m3)
summary(m4)

m5 <- lm(F0 ~ T0, data= hnb_cN)
m6 <- lm(log(F0) ~ log(T0), data= hnb_cN)
summary(m5)
summary(m6)

m7 <- lm(F0 ~ T0, data= hnb_cN)
m8 <- lm(log(F0) ~ log(T0), data= hnb_cN)
summary(m7)
summary(m8)

# Remove locality with richness = 17 in the historical period:
hnc66 <- hnc3[!hnc3$T0 ==17,]
hnb66 <- hnb3[!hnb3$T0 ==17,]

m9 <- lm(F0 ~ T0, data=hnc66)
m10 <- lm(log(F0) ~ log(T0), data=hnc66)

summary(m9)
summary(m10)


m11 <- lm(F0 ~ T0, data=hnb66)
m12 <- lm(log(F0) ~ log(T0), data=hnb66)

summary(m11)
summary(m12)

# Function extract model coefficients:-------------------------------------------------------
#listaic <- data.frame(rbind(hnb_aic[,c(1,3:4)],contN_aic[,c(1,3:4)], contAll_aic[,c(1,3:4)]))
#listModels <- list(hnb_2, hnb_1, contN_2, contN_1, contAll_2, contAll_1)
#extract_m <- function(mod=NULL) {
#  x <- mod
# df <- x$df.residual
#  fstat <- summary(x)$fstatistic[1]
# adjr2 <- summary(x)$adj.r.squared
# pval <- summary(x)$coefficients[2,4]
# df <- as.data.frame(cbind(df, fstat, adjr2, pval))
#  return(df)
#}                         
#l_ext_m <- lapply(listModels, extract_m)
#l_ext_m <- do.call(rbind, l_ext_m)

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

#ggsave(p3, file= paste0(plot_dir, "/SM/Ms/S5_F0vsT0HNCPower.jpg"), width = 8, height = 6) 

df_fit <- rbind(hnb_f, contN_f, contAll_f)

#(p4 <- ggplot() + 
#    geom_point(data=div3_subset_main, aes(x = T0, y = F0, col=Period), 
#              size = 2, alpha = 0.7, shape = 19) +
#    geom_line(data = df_fit , aes(x = T0, y = fit, col=Period), linetype = 1, size = 1) + 
#    geom_ribbon(data = df_fit, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3)+ 
#    labs(x = "Species Richness", y = "Trait diversity")+
#    scale_color_manual(values=colorBlind3, 
#                       breaks=c("Historical Broad", "Contemporary Natives", "Contemporary Natives + Exotics"), 
#                       labels=c("Historical Broad (n=60)","Contemporary Natives (n=33)", "Contemporary Natives + Exotics (n=49)"))+
#   facet_wrap(~Period, ncol=2)+
#   theme_bw())
#ggsave(p4, file= paste0(plot_dir, "/SM/Ms/S5_F0vsT0Power3.jpg"), width = 8, height = 6) 


##########################################################################################
# Null models: ---------------------------------------------------------------------------
species_list <- unique(c(names(HNB[,-c(1:2)]), names(ContAll[,-c(1:2)]))) 

rand <- list()
rand17 <- list() # 17 sps max observed in any site and period

for (j in 1:999){
  for(i in 1:17){
    df <- as.data.frame(matrix(ncol=100, nrow=1))
    names(df) <- species_list
    samp <- sample(species_list, size = i, replace = FALSE)
    df[as.character(samp)] <- 1
    df[is.na(df)] <- 0
    rand17[[i]] <- df
  } 
  rand[[j]] <- do.call(rbind, rand17)
} 
rand <- do.call(rbind, rand)
sum(rowSums(rand)==17) #999, OK
save(rand, file="rand.RData")

# Compute FD & Redundancy: ---------------------------------------------------------------------
rand <- rand[, order(names(rand))]
dist_mat1 <- dist_mat1[order(rownames(dist_mat1)), order(colnames(dist_mat1))]
identical(colnames(dist_mat1), names(rand)) # TRUE

rand_l <- split(rand, f=rownames(rand))

F0null <- list()
for (i in 1:length(rand_l)){
  assemblage <- as.matrix(t(rand_l[[i]]))
  T0_n <- colSums(assemblage)
  F0_n <- FD_MLE(assemblage, dist_mat1, mean(dist_mat1[dist_mat1>0]), q=0)
  R0_n <- (T0_n - F0_n)
  F0null[[i]] <- cbind(T0_n, F0_n, R0_n)
}
F0null <- as.data.frame(do.call(rbind,  F0null))
save(F0null, file="F0null.RData")

##########################################################################################
# Plots:----------------------------------------------------------------------------------
(nullF0 <- ggplot(data=F0null, aes(x=as.factor(T0_n), y=F0_n))+
    geom_boxplot(alpha=0.3, color="gray69")+
    stat_boxplot(geom ='errorbar', color="gray69") + 
    xlab("SR")+
    ylab("T0")+
    theme_bw()) # null SR vs T0 (17)
(nullF0II <- ggplot(data=F0null[F0null$T0_n<10,], aes(x=as.factor(T0_n), y=F0_n))+
    geom_boxplot(alpha=0.3, color="gray69")+
    stat_boxplot(geom ='errorbar', color="gray69") + 
    xlab("SR")+
    ylab("T0")+
    theme_bw()) # null SR vs T0 (9)

(p1nullF0_hnc <- nullF0 + # Historical Conservative
    geom_point(data = div3[div3$Period=="Historical Conservative",], 
               aes(x=T0, y=F0), size=2, alpha=1, col="#F0E442")+
    geom_line(data = hnc_f , aes(x = T0, y = fit), linetype = 1, size = 1, col="#F0E442")+
    geom_ribbon(data = hnc_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3))

(p1nullF0_hnb <- nullF0 + # Historical Broad
    geom_point(data = div3[div3$Period=="Historical Broad",], 
               aes(x=T0, y=F0), size=2, alpha=1, col="#E69F00")+
    geom_line(data = hnb_f , aes(x = T0, y = fit), linetype = 1, size = 1, col="#E69F00")+
    geom_ribbon(data = hnb_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3))

(p1nullF0_contN <- nullF0II + # Contemporary Natives
    geom_point(data = div3[div3$Period=="Contemporary Natives",], 
               aes(x=T0, y=F0), size=2, alpha=1, col="#009E73")+
    geom_line(data = contN_f, aes(x = T0, y = fit), linetype = 1, size = 1, col="#009E73")+
    geom_ribbon(data = contN_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3))

(p1nullF0_contAll <- nullF0II + # Contemporary Natives + Exotics
    geom_point(data = div3[div3$Period=="Contemporary Natives + Exotics",], 
               aes(x=T0, y=F0), size=2, alpha=1, col="#0072B2")+
    geom_line(data = contAll_f , aes(x = T0, y = fit), linetype = 1, size = 1, col="#0072B2")+
    geom_ribbon(data = contAll_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3))


ggsave(p1nullF0_hnc, file= paste0(plot_dir, "/SM/Ms/Null-ObservedHNC.jpg"), width = 8, height = 6) 
ggsave(p1nullF0_hnb, file= paste0(plot_dir, "/Null-ObservedHNB.jpg"), width = 8, height = 6) 
ggsave(p1nullF0_contN, file= paste0(plot_dir, "/Null-ObservedContN.jpg"), width = 8, height = 6) 
ggsave(p1nullF0_contAll, file= paste0(plot_dir, "/Null-ObservedContAll.jpg"), width = 8, height = 6) 

###########################################################################################
# SES:-------------------------------------------------------------------------------------
null_l <- split(F0null, f=F0null$T0_n)
null_l_mean <- lapply(null_l, function(x) {c(mean(x$F0_n),
                                             mean(x$R0_n))})
null_l_sd <- lapply(null_l, function(x) {c(sd(x$F0_n),
                                           sd(x$R0_n))})
ll <- lapply(l3, function(x) {split(x, f=x$T0)})
vecs <- lapply(ll, function(x) {as.vector(names(x))})

means <- lapply(vecs, function(x) {null_l_mean[x]})
sds <- lapply(vecs, function(x) {null_l_sd[x]})

ses_results_hnc <- mapply(function(x,y,z){((x$F0 - y[1])/z[1])}, ll[[1]], means[[1]], sds[[1]], SIMPLIFY = FALSE)
ses_results_hnb <- mapply(function(x,y,z){((x$F0 - y[1])/z[1])}, ll[[2]], means[[2]], sds[[2]], SIMPLIFY = FALSE)
ses_results_cn <- mapply(function(x,y,z){((x$F0 - y[1])/z[1])}, ll[[3]], means[[3]], sds[[3]], SIMPLIFY = FALSE)
ses_results_call <- mapply(function(x,y,z){((x$F0 - y[1])/z[1])}, ll[[4]], means[[4]], sds[[4]], SIMPLIFY = FALSE)

ses_results_hncR0 <- mapply(function(x,y,z){((x$R0 - y[2])/z[2])}, ll[[1]], means[[1]], sds[[1]], SIMPLIFY = FALSE)
ses_results_hnbR0 <- mapply(function(x,y,z){((x$R0 - y[2])/z[2])}, ll[[2]], means[[2]], sds[[2]], SIMPLIFY = FALSE)
ses_results_cnR0 <- mapply(function(x,y,z){((x$R0 - y[2])/z[2])}, ll[[3]], means[[3]], sds[[3]], SIMPLIFY = FALSE)
ses_results_callR0 <- mapply(function(x,y,z){((x$R0 - y[2])/z[2])}, ll[[4]], means[[4]], sds[[4]], SIMPLIFY = FALSE)

as_v <- function(x){as.vector(unlist(x))}
ses_resultsF0 <- c(as_v(ses_results_hnc),as_v(ses_results_hnb),
                   as_v(ses_results_cn),as_v(ses_results_call))
ses_resultsR0 <- c(as_v(ses_results_hncR0),as_v(ses_results_hnbR0),
                   as_v(ses_results_cnR0),as_v(ses_results_callR0))

period <- rep(c("Historical Conservative", "Historical Broad", "Contemporary Natives", "Contemporary Natives + Exotics"), times=c(48,60,33,49))
ses_dt <- data.frame("Period"=period, "SESF0"=ses_resultsF0, "SESR0"=ses_resultsR0)
ses_dt <- ses_dt[!is.na(ses_dt$SESF0),]

ses_dt %>% group_by(Period) %>% summarise(count = n(),
                                          countsAbv = sum(SESF0 > 0, na.rm = TRUE),
                                          countsUnd = sum(SESF0 < 0, na.rm = TRUE),
                                          countsAbvPer = (sum(SESF0 > 0, na.rm = TRUE)/count)*100,
                                          countsUndPer = (sum(SESF0 < 0, na.rm = TRUE)/count)*100)

summF0 <- summarySE(ses_dt, measurevar=c("SESF0"), groupvars=c("Period"), na.rm = TRUE)
summR0 <- summarySE(ses_dt, measurevar=c("SESR0"), groupvars=c("Period"), na.rm = TRUE)
summF0$Period <- recode_factor(summF0$Period,
                               "Historical Conservative"="HNC",
                               "Historical Broad"="HNB",
                               "Contemporary Natives"="CN",
                               "Contemporary Natives + Exotics"="CN+E")
ses_dt$Period <- recode_factor(ses_dt$Period,
                               "Historical Conservative"="HNC",
                               "Historical Broad"="HNB",
                               "Contemporary Natives"="CN",
                               "Contemporary Natives + Exotics"="CN+E")

(summF0hnc <- ggplot(ses_dt[ses_dt$Period=="HNC",], 
                    aes(x=Period, y=SESF0, colour=Period)) +
  geom_point(alpha=0.4, col="#F0E442") +
  geom_violin(alpha=0.1, col="#F0E442", fill="#F0E442") +
  #geom_boxplot(width=0.1, alpha=0.5) +
  geom_point(data=summF0[summF0$Period=="HNC",], 
             aes(x=Period, y=SESF0), size=2.5, col="#F0E442") +
  geom_errorbar(data=summF0[summF0$Period=="HNC",],
                aes(ymin=SESF0-ci, ymax=SESF0+ci),
                width=0.2, col="#F0E442")+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y="SES", x="")+
  geom_hline(yintercept = 0, linetype="dashed"))
  
(summF03p <- ggplot(ses_dt[!ses_dt$Period=="HNC",],
                    aes(x=Period, y=SESF0, colour=Period)) + 
    geom_point(alpha=0.5) +
    geom_violin(alpha=0.1, aes(fill=Period)) +
    geom_point(data=summF0[!summF0$Period=="HNC",], 
               aes(x=Period, y=SESF0, colour=Period), size=2.5)+
    geom_errorbar(data=summF0[!summF0$Period=="HNC",], aes(ymin=SESF0-ci, ymax=SESF0+ci),
                  width=0.2)+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_manual(values=colorBlind3, labels=c("HNB (n=60)", "CN (n=33)", "CN+E(n=49)"))+
    scale_fill_manual(values=colorBlind3, labels=c("HNB (n=60)", "CN (n=33)", "CN+E(n=49)"))+
    labs(y="SES", x="Period")+
    geom_hline(yintercept = 0, linetype="dashed")) # Same plot for redundancy.

ggsave(summF0hnc, file= paste0(plot_dir, "/SM/Ms/SES_F0hnc.jpg"), width = 6, height = 4) 
ggsave(summF03p, file= paste0(plot_dir, "/SES_F03.jpg"), width = 6, height = 4)

###########################################################################################
# One-tailed t test:-----------------------------------------------------------------------
ses_dt_n <- ses_dt %>%
  group_by(Period) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
ses_dtF0 <- ses_dt %>%
  group_by(Period) %>%
  summarise(res = list(t.test(SESF0, mu=0, alternative="less")))
ses_dtF0$res  # same p-vals when assessing redundancy

###########################################################################################
# Quantile scores:-------------------------------------------------------------------------
divrank <- div3[,-c(5)]
divrank <- split(divrank, f=divrank$T0)
vrank <- names(divrank)
nrows <- lapply(divrank, function(x) {nrow(x)})

null_lrank <- lapply(null_l, function(x) {x$F0_n})
null_lrank <- null_lrank[vrank]

null_lrank2 <- mapply(function(x,y) {as.data.frame(t(replicate(x, y)))}, nrows, null_lrank, SIMPLIFY = FALSE)
mat.rank.all <- mapply(function(x,y) {data.frame(cbind(x,y))}, divrank, null_lrank2, SIMPLIFY = FALSE)
mat.rank.all <- do.call(rbind, mat.rank.all)

l <- list()
for (i in 1:nrow(mat.rank.all)){
  ranks <- rank(as.vector(as.numeric(mat.rank.all[i,c(4:1003)])))[1]
  pvals <- ranks/1000
  l[[i]] <- data.frame(cbind(ranks, pvals))
}
l <- do.call(rbind, l)
l <- data.frame(mat.rank.all[,c(1:4)], l) 
l$Period <- recode_factor(l$Period, 
                          "Historical Conservative"="HNC",
                          "Historical Broad"="HNB",
                          "Contemporary Natives"="CN",
                          "Contemporary Natives + Exotics"="CN+E")
lSigN <- l[l$pvals<0.025,]
lSigP <- l[l$pvals>0.975,]

(pvalsp <- ggplot(l, aes(x=Period, y=pvals, colour=Period)) + 
    geom_violin(aes(fill=Period), alpha=0.1)+
    geom_boxplot(width=0.1, alpha=0.5)+
    theme_bw()+
    labs(y="P-values")+
    scale_color_manual(values=colorBlind4)+
    scale_fill_manual(values=colorBlind4)+
    theme(legend.position = "none")+
    geom_hline(yintercept = 0.025, linetype="dashed", color="red", alpha=0.5)+
    geom_hline(yintercept = 0.975, linetype="dashed", color="red", alpha=0.5)+
    geom_hline(yintercept = 0.75, linetype="dashed", color="lightgray")+
    geom_hline(yintercept = 0.25, linetype="dashed", color="lightgray")+
    geom_hline(yintercept = 0.5, linetype="dashed", color="lightgray"))

ggsave(pvalsp, file= paste0(plot_dir, "/SM/Thesis/Pvals_F0.jpg"), width = 6, height = 4) 
ggsave(pvalsp, file= paste0(plot_dir, "/SM/Ms/Pvals_F0.jpg"), width = 6, height = 4) 

##########################################################################################
# Figures (Ms) ---------------------------------------------------------------------------

# Main Figure 3: -------------------------------------------------------------------------
(Figure3 <- ggarrange(p1nullF0_hnb,
                      p1nullF0_contN,
                      p1nullF0_contAll,
                      summF03p,
                      labels=c("A)", "B)", "C)", "D)"),
                      common.legend = T,
                      legend="bottom",
                      ncol=2,
                      nrow=2,
                      font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))

ggsave(Figure3, file= paste0(plot_dir, "/S5_Figure3.jpg"), width =8, height = 8) 

# Supplementary material, results for HNC (yellow):---------------------------------------
(HNCres <- ggarrange(p1nullF0_hnc,
                     summF0hnc,
                     align="hv",
                     widths = c(1, 0.5),
                     labels=c("A)", "B)"),
                     nrow=1,
                     font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "bottom")))
ggsave(HNCres, file= paste0(plot_dir, "/SM/Ms/S5_HNC.jpg"), width =8, height = 5) 


##########################################################################################
# End of script ##########################################################################
