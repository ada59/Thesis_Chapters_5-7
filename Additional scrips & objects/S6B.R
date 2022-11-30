###########################################################################################
# Script: S6 B
# AFE
# August 2022
###########################################################################################

# Libraries:------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(mFD)
library(funrar)
library(data.table)
library(car)
library(grid)
library(ggpubr)


rm(list=ls())
myd <- getwd()
plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots" # Dir to save main plots


load("Dii.RData")
load("tt2.RData")

Dii$Period <- rep(NA, nrow(Dii))
Dii$Period <- ifelse(rownames(Dii) %like% "HNB", "Historical", Dii$Period)
Dii$Period <- ifelse(rownames(Dii) %like% "ContAll", "Contemporary", Dii$Period)

Dii$Site <- rep(NA, nrow(Dii))
Dii$Site <- str_split_fixed(rownames(Dii), "_", 2)[,1]

Dii <- gather(Dii, key="Species", value="Di", -c("Period", "Site"))
Dii <- Dii[!is.na(Dii$Di),]

Dii$SourceII <- tt2$SourceII[match(Dii$Species, rownames(tt2))]

# Subsets:

DiiI <- subset(Dii, ! Dii$SourceII == "NR")
DiiII <- subset(Dii, Dii$SourceII == "NR")

DiiIIinsp <- DiiII %>%
  group_by(Site, Species) %>%
  summarise(n=n())

DiiII_1 <- subset(DiiIIinsp, DiiIIinsp$n ==1)
DiiII_2 <- subset(DiiIIinsp, DiiIIinsp$n ==2)

LocE <- merge(x=DiiII_1,y=DiiII,by=c("Site", "Species"), all.x=FALSE, all.y=FALSE)
LocNR <- merge(x=DiiII_2,y=DiiII,by=c("Site", "Species"), all.x=FALSE, all.y=FALSE)

sum(LocE$Period=="Contemporary") #11 (Translocated)
sum(LocE$Period=="Historical")   #142

LocE$SourceII <- ifelse(LocE$Period=="Historical", "E", LocE$SourceII)
LocE <- within(LocE, rm(n))
LocNR <- within(LocNR, rm(n))

Dii_v2 <- rbind(DiiI, LocNR, LocE) # 507
# The E category contains regional and local ext.
# Some species in both the E and NR categories.

save(Dii_v2, file="Dii_v2.RData")


violin_p <- function(data, x, y, pal, legend_title="legend", labx="labx", laby="laby") {
  plot <- ggplot(data, aes(x = x, y = y, fill=x)) +
    geom_point(alpha=0.3)+
    geom_violin(alpha=0.4) +
    geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
    xlab(labx) +
    ylab(laby) +
    labs(fill=legend_title)+
    scale_fill_manual(values=pal)+
    theme_bw()
  print(plot)
}

pal4 <- c("#0072B2","#D55E00","#E69F00","darkgray")
dviolinv2 <- violin_p(data=Dii_v2, x=Dii_v2$SourceII, y=Dii_v2$Di,
                     pal=pal4,
                     legend_title = "Status",
                     labx="Status", laby="Local distinctiveness")

ggsave(dviolinv2, filename= paste0(plot_dir, "/SM/Ms/S6B_LocalDisPanel_v2.jpg"), width = 7, height = 5) 
#save(dviolinv2, file="S6B_dviolin_v2.RData")
