###########################################################################################
# Script: Observed trait redundancy and temporal change
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
#------------------------------------- TRAIT DIVERSITY -----------------------------------
##########################################################################################

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
# Relationship between T0 & Period:
str(div2)
div2_1T0 <- lm(T0 ~ Period, data=div2)
summary(div2_1T0)

div2B <- subset(div2, div2$Period %in% c("HNB", "ContAll"))
div2_2T0 <- lm(T0 ~ Period, data=div2B)
summary(div2_2T0)

# Paired tests : -------------------------------------------------------------------------
# http://www.sthda.com/english/wiki/paired-samples-t-test-in-r

# Data visualization:
str(div)
names(div)

div %>% group_by(Period) %>% dplyr::summarise(count =n(), 
                                              meanT0 = mean(T0), sdT0= sd(T0),
                                              meanF0 = mean(F0), sdF0 = sd(F0))

div2 %>% group_by(Period) %>% dplyr::summarise(count=n(),
                                               meanT0 = mean(T0), sdT0= sd(T0),
                                               meanF0 = mean(F0), sdF0 = sd(F0))

# Taxonomic diversity:
T0hnc <- div$T0[div$Period == "A) Historical conservative"]
T0hnb <- div$T0[div$Period == "B) Historical broad"]
T0cnat <- div$T0[div$Period == "C) Contemporary Natives"]
T0call <- div$T0[div$Period == "D) Contemporary Natives + Exotics"]

# Trait diversity:
F0hnc <- div$F0[div$Period == "A) Historical conservative"]
F0hnb <- div$F0[div$Period == "B) Historical broad"]
F0cnat <- div$F0[div$Period == "C) Contemporary Natives"]
F0call <- div$F0[div$Period == "D) Contemporary Natives + Exotics"]

# Plot trend over time
# pd <- paired(before, after) 
# plot(pd, type = "profile") + theme_bw()

# Compute the difference:
sites46 <- as.character(droplevels(div2$SiteName[div2$Period=="ContN"]))
sites53 <- as.character(droplevels(div2$SiteName[div2$Period=="ContAll"]))

T0hnb53 <- div$T0[div$Period == "B) Historical broad" & div$SiteName %in% sites53]
T0call53 <- div$T0[div$Period == "D) Contemporary Natives + Exotics" & div$SiteName %in% sites53]
F0hnb53 <- div$F0[div$Period == "B) Historical broad" & div$SiteName %in% sites53]
F0call53 <- div$F0[div$Period == "D) Contemporary Natives + Exotics" & div$SiteName %in% sites53]

dT0 <- with(div, 
            T0[Period == "B) Historical broad"] - T0[Period == "D) Contemporary Natives + Exotics"])
shapiro.test(dT0) # Normality test for the differences
# => p-value = 5.623e-06 we can't assume normality (67 obs)

dF0 <- with(div, 
          F0[Period == "B) Historical broad"] - F0[Period == "D) Contemporary Natives + Exotics"])
shapiro.test(dF0) 
# => p-value = 0.001474 we can't assume normality (67 obs)
# => p-value = 0.0002118 we can't assume normality (53 obs)

# Non-parametric tests:
res <- wilcox.test(F0hnb53, F0call53, paired = TRUE, alternative = "two.sided") # change vecs for results
res

# Broad vs All (67)
# Historical broad T0 and Contemporary T0 p-value = 3.76e-06
# Historical broad F0 and Contemporary F0 p-value = 0.0001099

# Broad vs All (53)
# Historical broad T0 and Contemporary T0 p-value = 0.002095
# Historical broad F0 and Contemporary F0 p-value = 0.05706

# Plot these trends:

divWTp <- subset(div, div$Period %in% c("B) Historical broad", "C) Contemporary Natives", "D) Contemporary Natives + Exotics"))
divWTp <- gather(divWTp, key="Metric", value="Value", -c(2:3))

divWTp$Period <- as.character(divWTp$Period)

divWTp$Metric[divWTp$Metric=="T0"] <- "A) Taxonomic"
divWTp$Metric[divWTp$Metric=="F0"] <- "B) Trait"

divWTp$Period[divWTp$Period=="B) Historical broad"] <- "A) Historical"
divWTp$Period[divWTp$Period=="C) Contemporary Natives"] <- "B) Contemporary Natives"
divWTp$Period[divWTp$Period=="D) Contemporary Natives + Exotics"] <- "C) Contemporary Natives + Exotics"

divWTp <- droplevels(divWTp)

divWTp <- subset(divWTp, divWTp$Metric %in% c("A) Taxonomic", "B) Trait"))

(T0F0p <- ggplot(divWTp, aes(Value, factor(Period))) +
  geom_violin(aes(fill=factor(Period)), alpha=0.4) +
  geom_boxplot(aes(fill=factor(Period)), alpha=0.5, width=0.1, position=position_dodge(1))+
  xlab("Diversity") +
  ylab("Period") +
  labs(fill="Period")+
  scale_fill_manual(values=c("#7CAE00","#00BFC4", "#C77CFF"))+
  coord_flip()+
  theme_bw()+
  facet_wrap(~Metric))

T0F0p <- T0F0p + theme(axis.text.x=element_text(angle=25, hjust=1))

ggsave(T0F0p, file= paste0(plot_dir, "/T0F0p.jpg"), width = 10, height = 6)


##########################################################################################
#------------------------------------- TRAIT REDUNDANCY ----------------------------------
##########################################################################################

##########################################################################################
# Mesures of Redundancy:------------------------------------------------------------------
div$R0 <- div$T0 - div$F0 # at the assemblage level
div$R0I <- (div$T0 - div$F0)/div$T0 # at the assemblage level (standardized)

div_red <- subset(div, div$T0>0)

div_red$R0II <- div2_2$residuals
div_red$R0III <- div2_4$residuals
div_red$R0IV <- div2_6$residuals

##########################################################################################
# Explore the relationship between R0/F0 & T0:--------------------------------------------

div_red_all <- gather(div_red, key="RedundancyM", value="Value", -c(1:4))
div_red_subs <- subset(div_red_all, div_red_all$RedundancyM %in% c("R0", "R0I", "R0II"))

(redT0 <- ggplot(div_red_subs, aes(x=T0, y=Value, color=Period)) + 
  geom_smooth(se=FALSE)+
  geom_point()+
  theme_bw()+
  labs(y="Redundancy")+
  facet_wrap(~RedundancyM)) 

#ggsave(redT0, file= paste0(plot_dir, "/T0vsRedundancy.jpg"), width = 12, height = 5) 

red <- lm(R0I ~ T0, data=div_red)
summary(red)   
redII <- lm(R0I ~ log(T0 + 1), data=div_red)
summary(redII)   # Adjusted R-squared:  0.7265

AIC(red, redII) # Curve (R0 is fitted by lm)

redIII <- lm(R0II ~ T0, data=div_red)
summary(redIII)   
redIV <- lm(R0II ~ log(T0 + 1), data=div_red)
summary(redIV)   

# NOTES:
# R0II*
# Residual is actual y value − predicted y value
# Res < 0, observed lower than predicted
# Res > 0, observed higher than predicted

##########################################################################################
# Relationship between R0 & Period:
red1 <- lm(R0 ~ Period, data=div_red)
red2 <- lm(R0I ~ Period, data=div_red)
red3 <- lm(R0II ~ Period, data=div_red)
#summary() # OK

##########################################################################################
# Paired tests : -------------------------------------------------------------------------
# Trait redundancy:
R0hnb <- div_red$R0[div_red$Period == "B) Historical broad"]
R0call <- div_red$R0[div_red$Period == "D) Contemporary Natives + Exotics"]

R0Ihnb <- div_red$R0I[div_red$Period == "B) Historical broad"]
R0Icall <- div_red$R0I[div_red$Period == "D) Contemporary Natives + Exotics"]

R0IIhnb <- div_red$R0II[div_red$Period == "B) Historical broad"]
R0IIcall <- div_red$R0II[div_red$Period == "D) Contemporary Natives + Exotics"]

# Plot trend over time
# pd <- paired(before, after) 
# plot(pd, type = "profile") + theme_bw()

# 53 observations:
R0hnb53 <- div_red$R0[div_red$Period == "B) Historical broad" & div_red$SiteName %in% sites53]
R0call53 <- div_red$R0[div_red$Period == "D) Contemporary Natives + Exotics" & div_red$SiteName %in% sites53]

R0Ihnb53 <- div_red$R0I[div_red$Period == "B) Historical broad" & div_red$SiteName %in% sites53]
R0Icall53 <- div_red$R0I[div_red$Period == "D) Contemporary Natives + Exotics" & div_red$SiteName %in% sites53]

R0Ihnb46 <- div_red$R0I[div_red$Period == "B) Historical broad" & div_red$SiteName %in% sites46]
R0Icall46 <- div_red$R0I[div_red$Period == "C) Contemporary Natives" & div_red$SiteName %in% sites46]

R0IIhnb53 <- div_red$R0II[div_red$Period == "B) Historical broad" & div_red$SiteName %in% sites53]
R0IIcall53 <- div_red$R0II[div_red$Period == "D) Contemporary Natives + Exotics" & div_red$SiteName %in% sites53]

#dT0 <- with(div, 
#            T0[Period == "B) Historical broad"] - T0[Period == "D) Contemporary Natives + Exotics"])
#shapiro.test(dT0) # Normality test for the differences

res <- wilcox.test(R0Ihnb46, R0Icall46, paired = TRUE, alternative = "two.sided")
res

# Broad vs All (67) R0 : ¿
# Broad vs All (67) R0I : ?
# Broad vs All (67) R0II : 0.5519

# Plot these trends:

div_redWTp <- subset(div_red, div_red$Period %in% c("B) Historical broad", "C) Contemporary Natives", "D) Contemporary Natives + Exotics"))
div_redWTp <- div_redWTp[,-c(1,4:6,8:9)]

div_redWTp <- gather(div_redWTp, key="Metric", value="Value", -c(1:2))

div_redWTp$Period <- as.character(div_redWTp$Period)

div_redWTp$Period[div_redWTp$Period=="B) Historical broad"] <- "A) Historical"
div_redWTp$Period[div_redWTp$Period=="C) Contemporary Natives"] <- "B) Contemporary Natives"
div_redWTp$Period[div_redWTp$Period=="D) Contemporary Natives + Exotics"] <- "C) Contemporary Natives + Exotics"

div_redWTp <- droplevels(div_redWTp)
dim(div_redWTp)

166/3
67+46+53


(R0IIp <- ggplot(div_redWTp, aes(Value, factor(Period))) +  # Or R0Ip
    geom_violin(aes(fill=factor(Period)), alpha=0.4) +
    geom_boxplot(aes(fill=factor(Period)), alpha=0.5, width=0.1, position=position_dodge(1))+
    xlab("Redundancy") +
    ylab("Period") +
    labs(fill="Period")+
    scale_fill_manual(values=c("#7CAE00","#00BFC4", "#C77CFF"))+
    coord_flip()+
    theme_bw())

ggsave(R0IIp, file= paste0(plot_dir, "/R0IIp.jpg"), width = 9, height = 5)


##########################################################################################
# End of script ##########################################################################
