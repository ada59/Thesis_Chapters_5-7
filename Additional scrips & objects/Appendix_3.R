###########################################################################################
# Script: Model F0 vs T0
# AFE
# Oct 2022
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

##########################################################################################
rm(list=ls())
myd <- getwd()

load("div2.RData")

##########################################################################################
# Effect of period: ----------------------------------------------------------------------
div2_1 <- lm(F0 ~ T0, data=div2)
summary(div2_1)
div2_2 <- lm(F0 ~ log(T0+1), data=div2)
summary(div2_2)

par(mfrow = c(1, 2))
hist(div2_1$residuals, main= "Linear")
hist(div2_2$residuals, main="Logarithmic")

par(mfrow = c(2, 2))
plot(div2_1)
par(mfrow = c(2, 2))
plot(div2_2) 

AIC(div2_1, div2_2) # A curve is a better fit

div2_3 <- lm(F0 ~ T0 + Period, data=div2)
summary(div2_3)
div2_4 <- lm(F0 ~ log(T0+1) + Period, data=div2)
summary(div2_4)

par(mfrow = c(1, 2))
hist(div2_3$residuals, main= "Linear")
hist(div2_4$residuals, main="Logarithmic")

par(mfrow = c(2, 2))
plot(div2_3)
par(mfrow = c(2, 2))
plot(div2_4) 

div2_5 <- lm(F0 ~ T0*Period, data=div2)
summary(div2_5)
div2_6 <- lm(F0 ~ log(T0+1)*Period, data=div2)
summary(div2_6)

par(mfrow = c(1, 2))
hist(div2_5$residuals, main= "Linear")
hist(div2_6$residuals, main="Logarithmic")

par(mfrow = c(2, 2))
plot(div2_5)
par(mfrow = c(2, 2))
plot(div2_6)

aictab(list(div2_1, div2_2, div2_3, div2_4, div2_5, div2_6))
# Logarithmic without accounting for period.

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

