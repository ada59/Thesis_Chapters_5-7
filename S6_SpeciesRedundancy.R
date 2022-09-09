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
load(paste0(myd, "/ContE.RData"))   # Contemporary Exotic Assemblage
load(paste0(myd, "/ContAll.RData")) # Contemporary Native + Exotic assemblage
load(paste0(myd, "/All.RData"))     # All Data

load(paste0(myd, "/dist_mat1.RData"))  # Species x species distance matrix
load(paste0(myd, "/tt.RData"))         # Species x trait matrix

# Format community data: -----------------------------------------------------------------
HNB$Period <- rep("HNB", nrow(HNB))
ContAll$Period <- rep("ContAll", nrow(ContAll))

dt <- bind_rows(HNB, ContAll) # OK
dt[is.na(dt)] <- 0
rownames(dt) <- paste0(dt$SiteNameE,"_", dt$Period)
dt <- within(dt, rm(SiteNameE, DrainageBasinE, Period))
dt <- dt[rowSums(dt)>0,] # 120 (-14, OK)
dt <- as.matrix(dt)
colnames(dt)[colnames(dt) == "Agonostomus monticola"] <- "Dajaus monticola"
dt <- dt[,order(colnames(dt))]
dist_mat1 <- dist_mat1[order(rownames(dist_mat1)), order(colnames(dist_mat1))]


tt <- tt[order(rownames(tt)),]
identical(colnames(dt), colnames(dist_mat1))  # TRUE
identical(colnames(dt), rownames(tt))         # TRUE

###########################################################################################
# Trait rarity measures:-------------------------------------------------------------------

Uii <- uniqueness(dt, dist_matrix = dist_mat1) # uniqueness
Uii$Ui <- (Uii$Ui-min(Uii$Ui))/(max(Uii$Ui)-min(Uii$Ui)) # normalise
Ui_dim <- uniqueness_dimensions(dt, 
                                tt,
                                metric="euclidean") 


Dii <- distinctiveness(dt, dist_matrix = dist_mat1) #distinctivenss
Dii <- as.data.frame(Dii)
Di_dim <- distinctiveness_dimensions(dt, 
                                     tt,
                                     metric="euclidean") 

Dii$Period <- rep(NA, nrow(Dii))
Dii$Period <- ifelse(rownames(Dii) %like% "HNB", "Historical", Dii$Period)
Dii$Period <- ifelse(rownames(Dii) %like% "ContAll", "Current", Dii$Period)
sum(is.na(Dii$Period))
Dii <- gather(Dii, key="Species", value="Di", -Period)
Dii <- Dii[!is.na(Dii$Di),]
length(unique(Dii$Species)) #99
setdiff(unique(colnames(dt)), unique(Dii$Species))

# "Gambusia yucatana" (it means, this species only appears in coms of richness = 1)
# Distinctiveness of this taxon cannot be calculated.


###########################################################################################
# Period & status -------------------------------------------------------------------------
nat <- sort(unique(names(ContN[,-c(1:2)]))) # 48
nat[nat=="Agonostomus monticola"] <-"Dajaus monticola"
ext <- sort(unique(setdiff(names(HNB), names(ContN)))) # 30
int <- sort(unique(setdiff(names(ContE), names(HNB)))) # 22 (translocated species not taken into account)

nat <- c(nat, "Chapalichthys encaustus")
ext <- setdiff(ext, "Chapalichthys encaustus")

Uii$Status <- rep(NA, nrow(Uii))
Uii$Status <- ifelse(Uii$species %in% nat, "Native Remaining", Uii$Status)
Uii$Status <- ifelse(Uii$species %in% ext, "Extirpated", Uii$Status)
Uii$Status <- ifelse(Uii$species %in% int, "Introduced", Uii$Status)
sum(is.na(Uii$Status))
Uii$Status <- as.factor(Uii$Status)

Dii$Status <- rep(NA, nrow(Dii))
Dii$Status <- ifelse(Dii$Species %in% nat, "Native Remaining", Dii$Status)
Dii$Status <- ifelse(Dii$Species %in% ext, "Extirpated", Dii$Status)
Dii$Status <- ifelse(Dii$Species %in% int, "Introduced", Dii$Status)

#Dii$Period <- as.factor(Dii$Period)
#Dii$Status <- as.factor(Dii$Status)
#Dii$Species <- paste0(substr(Dii$Species,1,1), "_",str_split_fixed(Dii$Species, " ", 2)[,2])

###########################################################################################
# Plots -----------------------------------------------------------------------------------

# Violin plots:

(uviolin <- ggplot(Uii,aes(x = Status, y = Ui, fill=Status)) +
   geom_point(alpha=0.3)+
   geom_violin(alpha=0.4) +
   geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
   xlab("Status") +
   ylab("Uniqueness") +
   scale_fill_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
   theme_bw())
Dii$Period <- as.factor(Dii$Period)
Dii$Period <- recode_factor(Dii$Period, "Historical"="A) Historical", "Current"="B)Current")
(dviolin <- ggplot(Dii, aes(y = Di, x = Period, fill=Status, color=Status)) +
    geom_point(alpha=0.05, position=position_dodge(1))+
    geom_violin(alpha=0.4) +
    geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
    xlab("Period") +
    ylab("Distinctiveness") +
    scale_fill_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
    scale_color_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
    theme_bw())
ggsave(dviolin, filename=paste0(plot_dir, "/dviolin.jpg"), width=10, height=8)
Dii$StatusII <- ifelse((Dii$Period=="Historical" & Dii$Status =="Native Remaining"),
                     "A) Native Remaining Historical", 
                     Dii$Status)
Dii$StatusII <- recode_factor(as.factor(Dii$StatusI), 
                              "A) Native Remaining Historical"="A) Native Remaining Historical",
                              "Native Remaining" = "B) Native Remaining Current",
                              "Extirpated"= "C) Extirpated",
                              "Introduced" = "D) Introduced")

(dviolinII <- ggplot(Dii,aes(x = StatusII, y = Di, fill=StatusII)) +
    geom_point(alpha=0.3)+
    geom_violin(alpha=0.4) +
    geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
    xlab("StatusII") +
    ylab("Distinctiveness") +
    scale_fill_manual(values=c("Lightgray", "Darkgray", "#0072B2", "#F0E442"))+
    theme_bw())


# Point plots:

#(upoint <- ggplot(Uii, aes(x=Ui, y=reorder(species, Ui))) + 
# geom_point(stat='identity', aes(col=Status), size=2, alpha=0.7)+
# labs(x="Uniqueness", y="Species")+
# theme(axis.text.y = element_blank())+
# scale_color_manual(values=c("#0072B2", "#F0E442", "Darkgray"))+
# theme_bw())

sum(is.na(Dii$Di)) #0
length(sort(unique(Dii$Species))) # 99 


Dii_m <- Dii %>%
  group_by(Species, StatusII, Period) %>%
  summarise(Di_m=mean(Di))
names(Dii_m)[names(Dii_m)=="Species"] <- "species"
length(unique(Dii_m$species))

merged <- full_join(Dii_m, Uii, by=c("species"))

merged <- merged[!is.na(merged$Di_m),]
merged$Di_m <- (merged$Di_m - min(merged$Di_m))/(max(merged$Di_m)-min(merged$Di_m)) # normalise
range(merged$Di_m)


(UivsDi <- ggplot(merged, aes(x=Di_m , y=Ui, color=Status))+
  geom_point()+
  stat_smooth()+
  facet_wrap(~Period))

# Boxplot distinctiveness (color Period):

(PeriodDii <- ggplot(Dii, aes(x=Species , y=Di, color=Period))+
    geom_boxplot()+
    theme(axis.text.x = element_blank())+
    theme_classic()+
    facet_wrap(~Status))

###########################################################################################
# Save plots: -----------------------------------------------------------------------------



###########################################################################################
# End of script ###########################################################################
