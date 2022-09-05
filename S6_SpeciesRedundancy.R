###########################################################################################
# Script: Compute species distinctiveness and uniqueness
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

rm(list=ls())
myd <- getwd()
plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots" # Dir to save main plots

# Community data:-------------------------------------------------------------------------
load(paste0(myd, "/HNC.RData"))     # Conservative Historical Native Assemblage
load(paste0(myd, "/HNB.RData"))     # Broad Historical Native Assemblage
load(paste0(myd, "/ContN.RData"))   # Contemporary Native Assemblage
load(paste0(myd, "/ContE.RData"))   # Contemporary Native Assemblage
load(paste0(myd, "/ContAll.RData")) # Contemporary Exotic assemblage
load(paste0(myd, "/All.RData"))     # All Data

load(paste0(myd, "/dist_mat1.RData"))  # Trait distance matrix
load(paste0(myd, "/tt.RData"))

# Format community data: -----------------------------------------------------------------
# Site-species matrix
l <- list("HNB"=HNB, "ContAll"=ContAll)
l <- lapply(l, function(x) {rownames(x) <- x$SiteNameE;x
                            within(x, rm(SiteNameE, DrainageBasinE))
                            as.matrix(x)}) #OK

rownames(All) <- paste0(All$SiteNameE,"_", All$Period)
rownames(All) <- gsub(" ", "_", rownames(All))
All <- within(All, rm(SiteNameE, DrainageBasinE, Period))
All <- as.matrix(All)


colnames(All)[colnames(All) == "Agonostomus monticola"] <- "Dajaus monticola"
All <- All[,order(colnames(All))]
dist_mat1 <- dist_mat1[order(rownames(dist_mat1)), order(colnames(dist_mat1))]


tt <- tt[order(rownames(tt)),]
identical(colnames(All), colnames(dist_mat1))  # TRUE
identical(colnames(All), rownames(tt))         # TRUE

###########################################################################################
# Uniqueness (regional):------------------------------------------------------------------
Uii <- uniqueness(All, dist_matrix = dist_mat1)

#Ui_dim <- uniqueness_dimensions(All, 
#                                tt,
#                               metric="euclidean")


# Native/Introduced/Remaining
nat <- sort(unique(names(ContN[,-c(1:2)]))) # 48
nat[nat=="Agonostomus monticola"] <-"Dajaus monticola"
ext <- sort(unique(setdiff(names(HNB), names(ContN)))) # 30
int <- sort(unique(setdiff(names(ContE), names(HNB)))) # 22 (translocated species not taken into account)

nat <- c(nat, "Chapalichthys encaustus")
ext <- setdiff(ext, "Chapalichthys encaustus")

# All traits:
Uii$Status <- rep(NA, nrow(Uii))
Uii$Status <- ifelse(Uii$species %in% nat, "Native Remaining", Uii$Status)
Uii$Status <- ifelse(Uii$species %in% ext, "Extirpated", Uii$Status)
Uii$Status <- ifelse(Uii$species %in% int, "Introduced", Uii$Status)
sum(is.na(Uii$Status))
str(Uii)
Uii$Status <- as.factor(Uii$Status)
Uii$species <- paste0(substr(Uii$species,1,1), "_",str_split_fixed(Uii$species, " ", 2)[,2])


# Plots:
#(p1 <- ggplot(Uii, aes(x=Ui, y=reorder(species, Ui))) + 
# geom_point(stat='identity', aes(col=Status), size=4, alpha=0.7)+
#  labs(x="Uniqueness", y="Species")+
#  scale_color_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
#  theme(axis.text = element_text(size = 30))+
# theme_bw())

Uii$Ui <- (Uii$Ui-min(Uii$Ui))/(max(Uii$Ui)-min(Uii$Ui))

(p1bp <- ggplot(Uii,
                aes(x = Status, y = Ui, fill=Status)) +
                geom_violin(alpha=0.3) +
                geom_boxplot(alpha=0.5, width=0.1)+
                xlab("Status") +
                ylab("Uniqueness") +
                scale_fill_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
                labs(title = "Regional-level trait rarity (=Uniqueness)")+
                theme_bw())

ggsave(p1bp, file= paste0(plot_dir, "/UniquenessAllSpeciesBP.jpg"), width = 9, height = 6)

###########################################################################################
# Distinctiveness (local):-----------------------------------------------------------------
Dii <- distinctiveness(All, dist_matrix = dist_mat1)
Dii <- as.data.frame(Dii)

Dii$Period <- rep(NA, nrow(Dii))
Dii$Period <- ifelse(rownames(Dii) %like% "HNC", "A) Historical Conservative", Dii$Period)
Dii$Period <- ifelse(rownames(Dii) %like% "HNB", "B) Historical Broad", Dii$Period)
Dii$Period <- ifelse(rownames(Dii) %like% "ContN", "C) Current Natives", Dii$Period)
Dii$Period <- ifelse(rownames(Dii) %like% "ContAll", "D) Current Natives + Exotics", Dii$Period)

Dii <- gather(Dii, key="Species", value="Di", -Period)
Dii <- Dii[!is.na(Dii$Di),]

Dii$Status <- rep(NA, nrow(Dii))
Dii$Status <- ifelse(Dii$Species %in% nat, "Native Remaining", Dii$Status)
Dii$Status <- ifelse(Dii$Species %in% ext, "Extirpated", Dii$Status)
Dii$Status <- ifelse(Dii$Species %in% int, "Introduced", Dii$Status)

Dii$Species <- as.factor(Dii$Species)
Dii$Period <- as.factor(Dii$Period)
Dii$Status <- as.factor(Dii$Status)
Dii$Species <- paste0(substr(Dii$Species,1,1), "_",str_split_fixed(Dii$Species, " ", 2)[,2])

DiiII <- subset(Dii, Dii$Period %in% c("B) Historical Broad", "D) Current Natives + Exotics"))

(pDii1 <- ggplot(DiiII, aes(y = Di, x = Species, color=Period, fill=Period)) +
    geom_boxplot(aes(fill=Period))+
    xlab("Species") +
    ylab("Distinctiveness") +
    labs(title = "Local-level species' trait rarity (=Distinctiveness)")+
    theme_classic()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    scale_fill_manual(values=c("#7CAE00", "#C77CFF"))+
    scale_color_manual(values=c("#7CAE00", "#C77CFF"))+
    facet_wrap(~Status))

(pDii2 <- ggplot(Dii, aes(y = Di, x = Period, fill=Status, color=Status)) +
    geom_violin(alpha=0.3) +
    geom_boxplot(width=0.1, position=position_dodge(1))+
    xlab("Period") +
    ylab("Distinctiveness") +
    labs(title = "Local-level trait rarity (=Distinctiveness)")+
    scale_fill_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
    scale_color_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
    theme_bw())

ggsave(pDii1, file= paste0(plot_dir, "/DistinctivenessAllSpecies.jpg"), width = 12, height =7) 
ggsave(pDii2, file= paste0(plot_dir, "/DistinctivenessAllSpecies2.jpg"), width = 10, height = 7) 

###########################################################################################
# End of script ###########################################################################
