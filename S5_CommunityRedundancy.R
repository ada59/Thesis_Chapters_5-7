###########################################################################################
# Script: Observed trait redundancy and temporal change
# AFE
# August/Sept 2022
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
library(Rmisc)
library(AICcmodavg)
library(ggpubr)
library(PairedData)

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

sum(is.na(ContAll[,-c(1:2)]))
ContAll[is.na(ContAll)] <- 0
sum(is.na(ContAll[,-c(1:2)]))
# save(ContAll, file="ContAll.RData")

l <- bind_rows(HNC, HNB, ContN, ContAll)
l[is.na(l)] <- 0
sum(is.na(l)) # 0

l$Period <- rep(c("HNC", "HNB", "ContN", "ContAll"), each=67)
All  <- l
#save(All, file="All.RData")
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
div$SiteName <- as.factor(div$SiteName)
div$Period <- as.factor(div$Period)
names(div) <- c("T0", "SiteName", "Period", "F0")

##########################################################################################
# Observed trends: -----------------------------------------------------------------------
# Subsets:
div2 <- div[!div$T0==0,]
hnc <- subset(div2, div2$Period=="HNC")
hnb <- subset(div2, div2$Period=="HNB")
contN <- subset(div2, div2$Period=="ContN")
contAll <- subset(div2, div2$Period=="ContAll")

# Historical conservative:----------------------------------------------------------------
hnc_1 <- lm(F0 ~ T0, data=hnc)
hnc_2 <- lm(F0 ~ log(T0+1), data=hnc)

AIC(hnc_1, hnc_2) # Curve

# Historical broad:-----------------------------------------------------------------------
hnb_1 <- lm(F0 ~ T0, data=hnb)
hnb_2 <- lm(F0 ~ log(T0+1), data=hnb)

AIC(hnb_1, hnb_2)

# Current native:-------------------------------------------------------------------------
contN_1 <- lm(F0 ~ T0, data=contN)
contN_2 <- lm(F0 ~ log(T0+1), data=contN)

AIC(contN_1, contN_2)

# Current All:----------------------------------------------------------------------------
contAll_1 <- lm(F0 ~ T0, data=contAll)
contAll_2 <- lm(F0 ~ log(T0+1), data=contAll)

AIC(contAll_1, contAll_2)
# NOTE: A curve is a better fit in all cases.

##########################################################################################
# Tests: ---------------------------------------------------------------------------------
# Test the fit in the historical period with the same "contemporary" site subsets.
hnc2 <- subset(hnc, hnc$SiteName %in% contN$SiteName)
hnb2 <- subset(hnb, hnb$SiteName %in% contN$SiteName)

hnc3 <- subset(hnc, hnc$SiteName %in% contAll$SiteName)
hnb3 <- subset(hnb, hnb$SiteName %in% contAll$SiteName)

# Compared to the ContN subset:
hnc2_1 <- lm(F0 ~ T0, data=hnc2)
hnc2_2 <- lm(F0 ~ log(T0+1), data=hnc2)
AIC(hnc2_1, hnc2_2) # Curve is still a better fit

hnb2_1 <- lm(F0 ~ T0, data=hnb2)
hnb2_2 <- lm(F0 ~ log(T0+1), data=hnb2)
AIC(hnb2_1, hnb2_2) # Curve is still a better fit

# Compared to the ContAll subset:
hnc3_1 <- lm(F0 ~ T0, data=hnc3)
hnc3_2 <- lm(F0 ~ log(T0+1), data=hnc3)
AIC(hnc3_1, hnc3_2) # Curve is still a better fit

hnb3_1 <- lm(F0 ~ T0, data=hnb3)
hnb3_2 <- lm(F0 ~ log(T0+1), data=hnb3)
AIC(hnb3_1, hnb3_2) # Curve is still a better fit

# Remove locality with richness = 17 in the historical period:
hnc66 <- hnc[!hnc$T0 ==17,]
hnb66 <- hnb[!hnb$T0 ==17,]

hnc66_1 <- lm(F0 ~ T0, data=hnc66)
hnc66_2 <- lm(F0 ~ log(T0+1), data=hnc66)

AIC(hnc66_1, hnc66_2)

hnb66_1 <- lm(F0 ~ T0, data=hnb66)
hnb66_2 <- lm(F0 ~ log(T0+1), data=hnb66)

AIC(hnb66_1, hnb66_2)

# The locality with more species doesn't drive the trend observed 
# (asymptotes continue to be a better fit)

##########################################################################################
# Plot trends:----------------------------------------------------------------------------
div$Period <- recode_factor(div$Period, HNC  = "A) Historical conservative", 
                            HNB = "B) Historical broad",
                            ContN = "C) Contemporary Natives",
                            ContAll = "D) Contemporary Natives + Exotics")


(p1 <- ggplot(div, aes(x=T0, y=F0, color=Period))+
  geom_point(aes(color=Period))+
  geom_smooth()+
  theme_classic()+
  facet_wrap(~Period))

(p2 <- ggplot(div, aes(x=T0, y=F0, color=Period))+
    geom_point(aes(color=Period))+
    geom_smooth(method=lm, formula = y ~ log(x+1))+
    theme_classic()+
    facet_wrap(~Period))

#ggsave(p1, file= paste0(plot_dir, "/RedundancyFreeFit.jpg"), width = 12, height = 10) 
#ggsave(p2, file= paste0(plot_dir, "/RedundancyLogLinear.jpg"), width = 8, height = 6) 


##########################################################################################
# Effect of period: ----------------------------------------------------------------------
div2_1 <- lm(F0 ~ T0, data=div2)
summary(div2_1)
div2_2 <- lm(F0 ~ log(T0+1), data=div2)
summary(div2_2)

par(mfrow = c(1, 2))
hist(div2_1$residuals, main= "Linear")
hist(div2_2$residuals, main="Log-Linear")

par(mfrow = c(2, 2))
plot(div2_1)
par(mfrow = c(2, 2))
plot(div2_2) #Heterocedasticity

AIC(div2_1, div2_2) # A curve is a better fit

div2_3 <- lm(F0 ~ T0 + Period, data=div2)
summary(div2_3)
div2_4 <- lm(F0 ~ log(T0+1) + Period, data=div2)
summary(div2_4)

par(mfrow = c(1, 2))
hist(div2_3$residuals, main= "Linear")
hist(div2_4$residuals, main="Log-Linear")

par(mfrow = c(2, 2))
plot(div2_3)
par(mfrow = c(2, 2))
plot(div2_4) #Heterocedasticity

div2_5 <- lm(F0 ~ T0*Period, data=div2)
summary(div2_5)
div2_6 <- lm(F0 ~ log(T0+1)*Period, data=div2)
summary(div2_6)

par(mfrow = c(1, 2))
hist(div2_5$residuals, main= "Linear")
hist(div2_6$residuals, main="Log-Linear")

par(mfrow = c(2, 2))
plot(div2_5)
par(mfrow = c(2, 2))
plot(div2_6) #Heterocedasticity

aictab(list(div2_1, div2_2, div2_3, div2_4, div2_5, div2_6))
# Log-linear without accounting for period.

div2$logT0 <- log(div2$T0 + 1)
div2_7 <- lm(F0~logT0*Period, data=div2)
summary(div2_7)
summary(div2_6) # Same, ok

##########################################################################################
# Plot results:---------------------------------------------------------------------------
# Model div2_6
hnc$logT0 <- log(hnc$T0 + 1)
hnb$logT0 <- log(hnb$T0 + 1)
contN$logT0 <- log(contN$T0 + 1)
contAll$logT0 <- log(contAll$T0 + 1)

pred <- data.frame("logT0" = c((seq(min(hnc$logT0), max(hnc$logT0), length.out = 67)),
                            (seq(min(hnb$logT0), max(hnb$logT0), length.out = 67)),
                            (seq(min(contN$logT0), max(contN$logT0), length.out = 46)),
                            (seq(min(contAll$logT0), max(contAll$logT0), length.out = 53))),
                   'Period' = c(rep("HNC", each = 67), rep("HNB", each = 67), 
                                rep("ContN", each=46), rep("ContAll", each=53)))
fit <- data.frame(pred, predict(div2_7, pred, se.fit = T))
head(fit)

fhnc <- subset(fit, Period=="HNC")
fhnb <- subset(fit, Period=="HNB")
fcontN <- subset(fit, Period=="ContN")
fcontAll <- subset(fit, Period=="ContAll")

dev.off()
plot(F0 ~ logT0, data=div2, pch=as.numeric(Period), ylab="Trait diversity", xlab = "Species richness", cex=0.5)

# Plot predicted lines for hnc
lines(fhnc$fit ~ fhnc$logT0, col="#F8766D")
lines((fhnc$fit + fhnc$se.fit) ~ fhnc$logT0, col="#F8766D", lty=3)
lines((fhnc$fit - fhnc$se.fit) ~ fhnc$logT0, col="#F8766D", lty=3)


# Plot predicted lines for hnb
lines(fhnb$fit ~ fhnb$logT0, col="#00BFC4")
lines((fhnb$fit + fhnb$se.fit) ~ fhnb$logT0, col="#00BFC4", lty=3)
lines((fhnb$fit - fhnb$se.fit) ~ fhnb$logT0, col="#00BFC4", lty=3)

# Plot predicted lines for contN
lines(fcontN$fit ~ fcontN$logT0, col="#7CAE00")
lines((fcontN$fit + fcontN$se.fit) ~ fcontN$logT0, col="#7CAE00", lty=3)
lines((fcontN$fit - fcontN$se.fit) ~ fcontN$logT0, col="#7CAE00", lty=3)

# Plot predicted lines for contAll
lines(fcontAll$fit ~ fcontAll$logT0, col="#C77CFF")
lines((fcontAll$fit + fcontAll$se.fit) ~ fcontAll$logT0, col="#C77CFF", lty=3)
lines((fcontAll$fit - fcontAll$se.fit) ~ fcontAll$logT0, col="#C77CFF", lty=3)

# add legend
legend("topleft", legend=c("HNC", "HNB", "ContN", "ContAll"),
       col=c("#F8766D", "#00BFC4","#7CAE00", "#C77CFF"),lty=1:2, cex=0.5)

# add title
mtext(side=3, "Trait diversity vs Species Richness", line=2, cex=1)


############################################
# Fix the issue with the fit lines in plot!!!!!!!!!!!!!!!!!!!!!!!!
############################################

##########################################################################################
# Paired test (Trait diversity): ---------------------------------------------------------
# http://www.sthda.com/english/wiki/paired-samples-t-test-in-r

# Plot:
str(div)

div %>% 
group_by(Period) %>%
summarise(mean = mean(F0, na.rm = TRUE),sd = sd(F0, na.rm = TRUE))

div2 %>%
group_by(Period) %>%
  summarise(mean = mean(F0, na.rm = TRUE),sd = sd(F0, na.rm = TRUE))

ggboxplot(div, x = "Period", y = "F0", 
          color = "Period",
          ylab = "F0", xlab = "Period")

div2plot <- div2
div2plot$Period  <- recode_factor(div2plot$Period, HNC  = "A) Historical conservative", 
                                               HNB = "B) Historical broad",
                                               ContN = "C) Contemporary Natives",
                                               ContAll = "D) Contemporary Natives + Exotics")
ggboxplot(div2plot, x = "Period", y = "F0", 
          color = "Period",
          ylab = "F0", xlab = "Period")


hconservative <- div$F0[div$Period == "A) Historical conservative"]
hbroad <- div$F0[div$Period == "B) Historical broad"]
cnative <- div$F0[div$Period == "C) Contemporary Natives"]
call <- div$F0[div$Period == "D) Contemporary Natives + Exotics"]

pd_1 <- paired(hconservative, cnative) 
plot(pd_1, type = "profile") + theme_bw()

pd_2 <- paired(hconservative, call) 
plot(pd_2, type = "profile") + theme_bw()

pd_3 <- paired(hbroad, cnative) 
plot(pd_3, type = "profile") + theme_bw()

pd_4 <- paired(hbroad, call) 
plot(pd_4, type = "profile") + theme_bw()

# compute the difference
sites53 <- droplevels(div2$SiteName[div2$Period=="ContAll"])
d <- with(div, 
          F0[Period == "B) Historical broad" & SiteName %in% sites53] - F0[Period == "D) Contemporary Natives + Exotics" & SiteName %in% sites53])

# Shapiro-Wilk normality test for the differences
shapiro.test(d) 
# => p-value = 0.001474 we can't assume normality (67 obs)
# => p-value = 0.0002118 we can't assume normality (53 obs)

# Non-parametric test (67 & 53):

res1 <- wilcox.test(hconservative, cnative, paired = TRUE)
res1

res2 <- wilcox.test(hconservative, call, paired = TRUE)
res2 

res3 <- wilcox.test(hbroad, cnative, paired = TRUE)
res3

res4 <- wilcox.test(hbroad, call, paired = TRUE)
res4 


#################################
# Test also with subset of 53!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#################################


##########################################################################################
# Redundancy:-----------------------------------------------------------------------------
div$R0 <- div$T0 - div$F0 # at the assemblage level
div$R0I <- (div$T0 - div$F0)/div$T0 # at the assemblage level (standardized)

div_red <- subset(div, div$T0>0)

div_red$R0II <- div2_2$residuals
div_red$R0III <- div2_4$residuals
div_red$R0IV <- div2_6$residuals

View(div_red)
##########################################################################################
# Explore the relationship with T0:-------------------------------------------------------

(redp <- ggplot(div_red, aes(x=T0, y=R0, color=Period)) + 
  geom_smooth()+
  geom_point()+
  theme_bw()+
  labs(y="Redundancy")) # A linear relationship seems to be a good fit

trait <- lm(F0 ~ log(T0+1), data = div)
summary(trait) # Adjusted R-squared: 0.7694 

traitII <- lm(F0 ~ log(T0+1) + as.factor(Period), data = div)
summary(traitII) # Adjusted R-squared: 0.8565 (No effect of period)

aictab(cand.set = list(trait, traitII)) # Better fit without period

red <- lm(R0 ~ T0, data=div_red)
summary(red)   # Adjusted R-squared: 0.8815

redII <- lm(R0 ~ T0 + as.factor(Period), data=div_red)
summary(redII)   # Adjusted R-squared: 0.8821 

trait_residuals <- data.frame("Period"=div_red$Period, "R0Res"=trait$residuals)
red_residuals <- data.frame("Period"=div_red$Period, "R0ResII"=red$residuals)

(residuals_traitp <- ggplot(trait_residuals, aes(x=Period, y=R0Res, color=Period)) + 
    geom_violin()+
    theme_bw()+
    labs(y="Redundancy")) 

(residuals_redp <- ggplot(red_residuals, aes(x=Period, y=R0ResII, color=Period)) + 
    geom_boxplot()+
    theme_bw()+
    labs(y="Redundancy"))

