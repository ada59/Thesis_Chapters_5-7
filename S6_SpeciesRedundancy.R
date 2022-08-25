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

Ui_dim <- uniqueness_dimensions(All, 
                                tt,
                                metric="euclidean")


# Native/Introduced/Remaining
nat <- sort(unique(names(ContN[,-c(1:2)]))) # 48
nat[nat=="Agonostomus monticola"] <-"Dajaus monticola"
ext <- sort(unique(setdiff(names(HNB), names(ContN)))) # 30
int <- sort(unique(setdiff(names(ContE), names(HNB)))) # 22 (translocated species not taken into account)


# All traits:
Uii$Status <- rep(NA, nrow(Uii))
Uii$Status <- ifelse(Uii$species %in% nat, "Native Remaining", Uii$Status)
Uii$Status <- ifelse(Uii$species %in% ext, "Extirpated", Uii$Status)
Uii$Status <- ifelse(Uii$species %in% int, "Introduced", Uii$Status)
sum(is.na(Uii$Status))
str(Uii)
Uii$Status <- as.factor(Uii$Status)
Uii$species <- paste0(substr(Uii$species,1,1), "_",str_split_fixed(Uii$species, " ", 2)[,2])


# Trait dimensions:
Ui_dim$Status <- rep(NA, nrow(Ui_dim))
Ui_dim$Status <- ifelse(Ui_dim$species %in% nat, "Native Remaining", Ui_dim$Status)
Ui_dim$Status <- ifelse(Ui_dim$species %in% ext, "Extirpated", Ui_dim$Status)
Ui_dim$Status <- ifelse(Ui_dim$species %in% int, "Introduced", Ui_dim$Status)
sum(is.na(Ui_dim$Status))
Ui_dim$Status <- as.factor(Ui_dim$Status)
Ui_dim$species <- paste0(substr(Ui_dim$species,1,1), "_",str_split_fixed(Ui_dim$species, " ", 2)[,2])
Ui_dim <- gather(Ui_dim, key="Trait", value="Ui", -c(1,13))

#Uii$StatusII <- rep("Native", nrow(Uii))
#Uii$StatusII <- ifelse(Uii$Status %in% "Introduced", "Introduced", Uii$StatusII)

#Uii$StatusIII <- rep("Others", nrow(Uii))
#Uii$StatusIII <- ifelse(Uii$Status %in% "Extirpated", "Extirpated", Uii$StatusIII)

# Plots:
(p1 <- ggplot(Uii, aes(x=Ui, y=reorder(species, Ui))) + 
  geom_point(stat='identity', aes(col=Status), size=4, alpha=0.7)+
  labs(x="Uniqueness", y="Species")+
  scale_color_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
  theme(axis.text = element_text(size = 30))+
  theme_bw())

(p1bp <- ggplot(Uii,
                aes(x = Status, y = Ui, fill=Status)) +
                geom_boxplot() +
                xlab("Status") +
                ylab("Uniqueness") +
                scale_fill_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
                labs(title = "Regional-level trait rarity (=Uniqueness)")+
                theme_bw())

#(p2 <- ggplot(Ui_dim, aes(x=Ui, y=reorder(species, Ui))) + 
#    geom_point(stat='identity', aes(col=Status), size=4, alpha=0.7)+
#    labs(x="Uniqueness", y="Species")+
#    scale_color_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
#    theme(axis.text = element_text(size = 30))+
#    theme_bw()+
#    facet_wrap(~Trait))

Ui_dim <- subset(Ui_dim, ! Ui_dim$Trait=="Ui_all")
(p2bp <- ggplot(Ui_dim,
               aes(x = Status, y = Ui, fill=Status)) +
               geom_boxplot() +
               xlab("Status") +
               ylab("Uniqueness") +
               ylim(0,0.1)+
               scale_fill_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
               labs(title = "Regional-level trait rarity (=Uniqueness/Trait)")+
               theme_bw()+
               facet_wrap(~Trait))
# Removed 50 rows containing non-finite values (stat_boxplot).


ggsave(p1, file= paste0(plot_dir, "/UniquenessAllSpecies.jpg"), width = 15, height = 15)
ggsave(p1bp, file= paste0(plot_dir, "/UniquenessAllSpeciesBP.jpg"), width = 10, height = 8)
ggsave(p2bp, file= paste0(plot_dir, "/UniquenessTraitAllSpeciesBP.jpg"), width = 15, height = 8) # Some extreme values removed (just preliminary overview)


###########################################################################################
# Distinctiveness (local):-----------------------------------------------------------------
Dii <- distinctiveness(All, dist_matrix = dist_mat1)

hnc <- subset(Dii, rownames(Dii) %like% "HNC")
hnb <- subset(Dii, rownames(Dii) %like% "HNB")
contN <- subset(Dii, rownames(Dii) %like% "ContN")
contAll <- subset(Dii, rownames(Dii) %like% "ContAll")

hnb_d <- apply(hnb, 2, mean, na.rm=TRUE)
sum(is.nan(hnb_d)) # 22 exotics, OK.
contall_d <- apply(contAll, 2, mean, na.rm=TRUE)
sum(is.nan(contall_d)) # 30 extirpated, OK

hnb_d_dat <- as.data.frame(hnb_d)
hnb_d_dat$species <- rownames(hnb_d_dat)
hnb_d_dat$Period <- rep("A) Historical", nrow(hnb_d_dat))
names(hnb_d_dat) <- c("Di", "Species", "Period")
rownames(hnb_d_dat) <- 1:nrow(hnb_d_dat)

contall_d <- as.data.frame(contall_d)
contall_d$species <- rownames(contall_d)
contall_d$Period <- rep("B) Contemporary", nrow(contall_d))
names(contall_d) <- c("Di", "Species", "Period")
rownames(contall_d) <- 1:nrow(contall_d)

Dii <- as.data.frame(rbind(hnb_d_dat, contall_d))

Dii$Status <- rep(NA, nrow(Dii))
Dii$Status <- ifelse(Dii$Species %in% nat, "Native Remaining", Dii$Status)
Dii$Status <- ifelse(Dii$Species %in% ext, "Extirpated", Dii$Status)
Dii$Status <- ifelse(Dii$Species %in% int, "Introduced", Dii$Status)

Dii$Status <- as.factor(Dii$Status)
Dii$Species <- paste0(substr(Dii$Species,1,1), "_",str_split_fixed(Dii$Species, " ", 2)[,2])

sum(is.nan(Dii$Di))
Dii$Di[is.nan(Dii$Di)] <- 0

(p3bp <- ggplot(Dii,
    aes(x = Status, y = Di, fill=Status)) +
    geom_boxplot() +
    xlab("Status") +
    ylab("Distinctiveness") +
    scale_fill_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
    labs(title = "Local-level trait rarity (=Distinctiveness)")+
    theme_bw()+
    facet_wrap(~Period))

ggsave(p3bp, file= paste0(plot_dir, "/DistinctivenessAllSpeciesBP.jpg"), width = 12, height = 10) # Look C encaustus

###########################################################################################
# End of script ###########################################################################
