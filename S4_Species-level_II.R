#===============================================================================
# Script: Compute species rarity measures
# AFE
# August 2022
#===============================================================================




#===============================================================================
# Libraries:--------------------------------------------------------------------
#===============================================================================
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
library(corrplot)

rm(list=ls())
myd <- getwd()
path_lists <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists"   # Raw data lists
path_traits <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Traits" # Dir to trait objects


path_plots5 <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots/Chapter5"   # Dir to save plots Chapter 5
path_plots6 <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots/Chapter6"   # Dir to save plots Chapter 6




#===============================================================================
# Read data:--------------------------------------------------------------------
#===============================================================================
load(paste0(path_lists, "/l67.RData"))          # data
load(paste0(path_traits, "/dist_mat1.RData"))   # Species x species distance matrix
load(paste0(path_traits, "/tt2.RData"))         # Species x trait matrix
load(paste0(path_lists, "/NDT67.RData"))  # List for checks on raw data
load(paste0(path_lists, "/NDT83.RData"))  # List for checks on raw data


## Format community data: ------------------------------------------------------
HNC <- subset(l67, l67$Period == "HNC")
HNC <- HNC[, colSums(HNC != 0) > 0]

HNB <- subset(l67, l67$Period == "HNB")
HNB <- HNB[, colSums(HNB != 0) > 0]

ContN <- subset(l67, l67$Period == "ContN")
ContN <- ContN[, colSums(ContN != 0) > 0]

ContE <- subset(l67, l67$Period == "ContE")
ContE <- ContE[, colSums(ContE != 0) > 0]

ContAll <- subset(l67, l67$Period == "ContAll")
ContAll <- ContAll[, colSums(ContAll != 0) > 0]




#===============================================================================
# Data for analyses:------------------------------------------------------------
#===============================================================================
dt <- bind_rows(HNB, ContAll) # OK (or change HNB to HNC)
dt[is.na(dt)] <- 0
rownames(dt) <- paste0(dt$SiteNameE,"_", dt$Period)
dt <- within(dt, rm(SiteNameE, DrainageBasinE, Period))
dt <- dt[rowSums(dt)>0,] # 120 (-14, OK)
dt <- as.matrix(dt)

dt <- dt[,order(colnames(dt))]
dist_mat1 <- dist_mat1[order(rownames(dist_mat1)), order(colnames(dist_mat1))]

tt2 <- tt2[order(rownames(tt2)),]
identical(colnames(dt), colnames(dist_mat1))  # TRUE
identical(colnames(dt), rownames(tt2))        # TRUE



#===============================================================================
# Trait rarity measures:--------------------------------------------------------
#===============================================================================
class(dist_mat1)
class(dt)


## Uniqueness: -----------------------------------------------------------------
Uii <- uniqueness(dt, dist_matrix = dist_mat1) # uniqueness (regional)
Uii$Ui_norm <- (Uii$Ui-min(Uii$Ui))/(max(Uii$Ui)-min(Uii$Ui)) # normalise
Ui_dim <- uniqueness_dimensions(dt, 
                                tt2[,c(1:10)],
                                metric="euclidean") 


## Local distinctivennss: ------------------------------------------------------
Dii <- distinctiveness(dt, dist_matrix = dist_mat1) # distinctiveness (local)
Dii <- as.data.frame(Dii)
warnings()
Di_dim <- distinctiveness_dimensions(dt, 
                                     tt2[,c(1:10)],
                                     metric="euclidean") 
warnings() # warnings when coms have 1 sps only, OK


## Global distinctiveness: -----------------------------------------------------
Diig <- distinctiveness_global(dist_mat1) # global distinctiveness (regional)
Diig$global_di_norm <- (Diig$global_di-min(Diig$global_di))/(max(Diig$global_di)-min(Diig$global_di)) # normalise
dt1 <- as.data.frame(t(dt[1,]))
dt1[1,] <- 1
Dig_dimg <- distinctiveness_dimensions(as.matrix(dt1), 
                                       tt2[,c(1:10)],
                                       metric="euclidean") 




#===============================================================================
# Period & Status --------------------------------------------------------------
#===============================================================================

## Uniqueness: -----------------------------------------------------------------
Uii$SourceII <- tt2$SourceII[match(Uii$species, rownames(tt2))]
Uii$SourceIISen <- tt2$SourceIISen[match(Uii$species, rownames(tt2))]


## Global distinctiveness:------------------------------------------------------
Diig$SourceII <- tt2$SourceII[match(Diig$species, rownames(tt2))]
Diig$SourceIISen <- tt2$SourceIISen[match(Diig$species, rownames(tt2))]


## Local distinctiveness:-------------------------------------------------------
Dii$Period <- rep(NA, nrow(Dii))
Dii$Period <- ifelse(rownames(Dii) %like% "HNB", "Historical", Dii$Period)
Dii$Period <- ifelse(rownames(Dii) %like% "ContAll", "Contemporary Natives + Exotics", Dii$Period)
sum(is.na(Dii$Period))

Dii$Site <- rep(NA, nrow(Dii))
Dii$Site <- str_split_fixed(rownames(Dii), "_", 2)[,1]
Dii <- gather(Dii, key="Species", value="Di", -c("Period", "Site"))
Dii[is.nan(Dii$Di),]  # Take into account when classifying LE, LR & IB
Dii <- Dii[!is.na(Dii$Di),]
Dii$SourceII <- tt2$SourceII[match(Dii$Species, rownames(tt2))]

DiiI <- subset(Dii, ! Dii$SourceII == "NR") # groups E, IA & IB remain the same
DiiII <- subset(Dii, Dii$SourceII == "NR")  # in regional NR, there might be some translocations, some local extirpations etc


## re-define NR at the local level: --------------------------------------------
DiiIIinsp <- DiiII %>%
  group_by(Site, Species) %>%
  summarise(n=n_distinct(Period)) # if 2, then NR (means sps in both periods)

DiiII_1 <- subset(DiiIIinsp, DiiIIinsp$n ==1) # locally extirpated (or a few cases of translocations)
DiiII_2 <- subset(DiiIIinsp, DiiIIinsp$n ==2) # locally remaining

LocE <- merge(x=DiiII_1,y=DiiII,by=c("Site", "Species"), all.x=FALSE, all.y=FALSE)
LocNR <- merge(x=DiiII_2,y=DiiII,by=c("Site", "Species"), all.x=FALSE, all.y=FALSE)

sum(LocE$Period=="Contemporary Natives + Exotics") # 11
sum(LocE$Period=="Historical")                     # 134 

LocE[LocE$Period=="Contemporary Natives + Exotics",]
# C. encaustus, tralocated in rio Gde. de Santiago en Ocotlan (IB)  [TRANSLOCATED]
# P. butleri (in Ameca River at Los Murillos, Ayutla River & Barraquitas River) (IB) [TRANSLOCATED]
# G atripinis at Lake Chalco (IB) [TRANSLOCATED]
# Xenoophorus captivus at Stream at Agua de Enmedio (only 1 native, so that's why it's missing) [NATIVE REMAINING]
# P sphenops Tamazula River is (IB) [TRANSLOCATED]
# P sphenops Turbio River close to El Nopal is "Poecilia sp" (IB) [TRANSLOCATED]
# G atripinis at Xochimilco Canals (IB) [TRANSLOCATED]
# G multiradiatus at Zempoala Lagoon is native, only 1 native, so that's why it's missing) [NATIVE REMAINING]
# P sphenops at Zinzimeo Canals is (IB) [TRANSLOCATED]

# Where is P gracilis in Xochimilco Canals, for instance? Double check

LocE$SourceII <- ifelse(LocE$Period=="Historical", "E", LocE$SourceII)
LocE$SourceII <- ifelse(LocE$Period=="Contemporary Natives + Exotics", "IB", LocE$SourceII)
LocE$SourceII <- ifelse(LocE$Species=="Girardinichthys multiradiatus" & LocE$Site=="Zempoala Lagoon", "NR", LocE$SourceII)
LocE$SourceII <- ifelse(LocE$Species=="Xenoophorus captivus" & LocE$Site=="Stream at Agua de Enmedio", "NR", LocE$SourceII)
LocE$SourceII <- ifelse(LocE$Species=="Chirostoma jordani" & LocE$Site=="Zumpango Lake", "NR", LocE$SourceII) # only sps remaining

LocE <- within(LocE, rm(n))
LocNR <- within(LocNR, rm(n))

Dii_v2 <- rbind(DiiI, LocNR, LocE) # 515 (diff from before)

#c2 <- Dii_v2 %>%
#  group_by(Species) %>%
#  distinct(SourceII) # for visual checks
#Dii_v2[Dii_v2$Species=="Characodon lateralis",] # for visual checks

length(unique(Dii_v2$Species))                        # 94    
setdiff(unique(colnames(dt)), unique(Dii_v2$Species)) # Gambusia yucatana

Dii_v2$SourceII[Dii_v2$SourceII=="E"] <-"LE"          # locally extirpated
Dii_v2$SourceII[Dii_v2$SourceII=="NR"] <-"L_NR"       # locally remaining
table(Dii_v2$SourceII)
# IA   IB L_NR   LE 
# 36   36  255  188 

save(Dii_v2, file="Dii_v2.RData")

# NOTES ########################################################################
# Some species have missing occurrences (when both historical and contemporary 
# community species richness = 1 or < 1) 
# For these, local distinctivenness cannot be calculated:
# C. lateralis as extirpated in Tunal River & Spring at Cerro Gordo
# "Gambusia yucatana" (it means, this species only appears in coms of richness 
# = 1). Local distinctiveness of this taxon cannot be computed anywhere.
# G multiradiatus at Rio Lerma cerca de Solis
# G viviparus at Texcoco 1 & Texcoco 2
# O. nilkoticus at La Colorada Lagoon
################################################################################




#===============================================================================
# Averaging Local Distinctiveness: ---------------------------------------------
#===============================================================================

## Checks: ---------------------------------------------------------------------
count_sps_dis <- Dii_v2 %>%
  group_by(Species) %>%
  summarise(n=n())
range(count_sps_dis$n) # 1 41 values per species

sum(count_sps_dis$n==1) # 20 species with 1 ocurrence only
sum(count_sps_dis$n==2) # 23 species with 2 ocurrences

count_sps_dis <- Dii_v2 %>%
  group_by(Species, SourceII) %>%
  summarise(n=n())
range(count_sps_dis$n) # 1 24 values per species & source


## Averages: -------------------------------------------------------------------
Dis_av <- Dii_v2 %>%
  group_by(Species, SourceII) %>% 
  summarise(av=mean(Di)) # this gives an average per species & per souce 
# (e.g. in some cases, the same species will have different av. distinctivenness 
# depending on whether it has the status of native remaining, IB, native extirpated])

Dis_av$av_norm <- (Dis_av$av-min(Dis_av$av))/(max(Dis_av$av)-min(Dis_av$av)) # normalise
Gy <- data.frame("Species"="Gambusia yucatana","av"=NA, "av_norm"=NA)
Dis_avx <- rbind(Dis_av, Gy)

identical(sort(unique(Dis_avx$Species)), sort(unique(Uii$species)))  # TRUE
identical(sort(unique(Dis_avx$Species)), sort(unique(Diig$species))) # TRUE

Dis_avSps <- Dii_v2 %>%
  group_by(Species) %>% 
  summarise(av=mean(Di)) # this gives an average per species & per souce 
Dis_avSps$av_norm <- (Dis_avSps$av-min(Dis_avSps$av))/(max(Dis_avSps$av)-min(Dis_avSps$av)) # normalise
Dis_avxSps <- rbind(Dis_avSps, Gy)




#===============================================================================
# Correlations: ----------------------------------------------------------------
#===============================================================================
shapiro.test(Dis_avxSps$av_norm)
shapiro.test(Uii$Ui_norm)
shapiro.test(Diig$global_di_norm)

cor_sing <- cor(data.frame("L.Dis"=Dis_avxSps$av, 
               "Uniq"=Uii$Ui, 
               "G.Dis"=Diig$global_di), method="spearman", use = "complete.obs")
cor(data.frame(Dis_avxSps$av_norm, Uii$Ui_norm, Diig$global_di_norm), method="spearman", use = "complete.obs") # identical, correct

file_path <- paste0(path_plots6, "/S4_SpearmanCorrelationRarity.png")
png(width = 400, height = 400, file=file_path)
corrplot(cor_sing, method="number", type = "lower", 
         order = "hclust", number.cex = 1.5, tl.cex = 1.5, diag = FALSE)
dev.off()


medianUii <- Uii %>%
  group_by(SourceII) %>%
  summarise(med = median(Ui_norm), sd= sd(Ui_norm))
medianUii

medianDii <- Dii_v2 %>%
  group_by(SourceII) %>%
  summarise(med = median(Di), sd= sd(Di))
medianDii

mediangbDii <- Diig %>%
  group_by(SourceII) %>%
  summarise(med = median(global_di_norm), sd= sd(global_di_norm))
mediangbDii



#===============================================================================
# Statistical differences Trait Singularity:------------------------------------
#===============================================================================

## Uniqueness:------------------------------------------------------------------
Uii$SourceII <- as.factor(Uii$SourceII)
sum(is.na(Uii))
str(Uii)

Uii_fit1 <- aov(Ui~SourceII, data=Uii)
summary(Uii_fit1)
par(mfrow = c(2, 2))
plot(Uii_fit1)
shapiro.test(residuals(Uii_fit1))

Uii_fit2 <- kruskal.test(Ui_norm~SourceII, data=Uii) # idem whther I use Ui_norm or not
Uii_fit2 # 0.03
pairwise.wilcox.test(Uii$Ui, Uii$SourceII,
                     p.adjust.method = "BH", exact=FALSE) # obs have ties (i.e., cases with two obs with same values, so cannot use exact=T)
# RESULTS: only weak indication that NR & E have a different uniqueness


## Global distinctiveness by Source II:-----------------------------------------
shapiro.test(residuals(aov(global_di_norm ~ SourceII, data=Diig)))
Diig_fit1 <- kruskal.test(global_di_norm ~ SourceII, data=Diig)
Diig_fit1
pairwise.wilcox.test(Diig$global_di_norm, Diig$SourceII,
                     p.adjust.method = "BH", exact = F)
# RESULTS: 
# 1) Weak indication that E has different distinctivenness with IB & NR
# 2) IA has a different distinctivenness from all the rest


# Av Distinctiveness:-----------------------------------------------------------
dev.off()
boxplot(Di ~ SourceII, data=Dii_v2)
dev.off()
Dii_fit1 <- aov(Di ~ SourceII, data=Dii_v2)
summary(Dii_fit1)
par(mfrow = c(2, 2))
plot(Dii_fit1)
dev.off()
shapiro.test(residuals(Dii_fit1))


Dii_fit2 <- kruskal.test(Di ~ SourceII, data=Dii_v2)
Dii_fit2
pairwise.wilcox.test(Dii_v2$Di, Dii_v2$SourceII,
                     p.adjust.method = "BH", exact=F) 
# Locally extirpated & locally remaining have the same local distinctivenness
# Both IA & IB have a different local distinctivenness in each comparison




#===============================================================================
# Plots:------------------------------------------------------------------------
#===============================================================================

pal3 <- c("#0072B2", "#F0E442", "darkgray")
pal4 <- c("#0072B2","#D55E00","#E69F00","darkgray")
pal4LD <- c("#009E73","#D55E00","#E69F00","#000000")

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
} # base violin plot function


## Uniqueness ------------------------------------------------------------------
uniqueness_violin2 <- violin_p(data=Uii, x=Uii$SourceII, y=Uii$Ui_norm,
                     pal=pal4,
                     legend_title = "Status",
                     labx="", laby="Uniqueness")
ggsave(uniqueness_violin2, filename= paste0(path_plots6, "/S4_uniqueness_violin.jpg"), width = 7, height = 5) 
save(uniqueness_violin2, file="uniqueness_violin2.RData")


## Global distinctiveness-------------------------------------------------------
global_distinctivenness_violin2 <- violin_p(data=Diig, x=Diig$SourceII, y=Diig$global_di_norm,
                      pal=pal4,
                      legend_title = "Status",
                      labx="", laby="Global distinctiveness")

ggsave(global_distinctivenness_violin2, filename= paste0(path_plots6, "/S4_global_distinctivenness_violin.jpg"), width = 7, height = 5) 
save(global_distinctivenness_violin2, file="global_distinctivenness_violin2.RData")


## Local distinctiveness -------------------------------------------------------
Dii_v2$SourceII <- as.factor(Dii_v2$SourceII)
Dii_v2$SourceII <- factor(Dii_v2$SourceII, levels = c("LE", "IA", "IB", "L_NR")) 
local_distinctivenness_violin2 <- ggplot(Dii_v2, aes(x = SourceII, y = Di, fill=SourceII)) +
  geom_point(alpha=0.3)+
  geom_violin(alpha=0.4) +
  geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
  xlab("") +
  ylab("Local distinctiveness") +
  labs(fill="Status")+
  scale_fill_manual(values=pal4LD)+
  theme_bw()
print(local_distinctivenness_violin2)

ggsave(local_distinctivenness_violin2, filename= paste0(path_plots6, "/S4_local_distinctivenness_violin.jpg"), width = 7, height = 5) 
save(local_distinctivenness_violin2, file="local_distinctivenness_violin2.RData")




#===============================================================================
# Main Sps-level figure:--------------------------------------------------------
#===============================================================================

load("trait_space_12_status.RData") # mainc6
load("trait_space_12_status_biplot.RData") # main_biplotc6


(panels_species_level <- ggarrange(main_c6,
                                   uniqueness_violin2,
                                   global_distinctivenness_violin2,
                                   local_distinctivenness_violin2,
                                   common.legend = T, # FALSE (Then panels_speciesII)
                                   legend="bottom",
                                   labels = c("A)", "B)", "C)", "D)"),
                                   align="hv",
                                   font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))
(panels_species_level_legend <- ggarrange(main_c6,
                                   uniqueness_violin2,
                                   global_distinctivenness_violin2,
                                   local_distinctivenness_violin2,
                                   common.legend = F, # FALSE (Then panels_speciesII)
                                   legend="bottom",
                                   labels = c("A)", "B)", "C)", "D)"),
                                   align="hv",
                                   font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))
ggsave(panels_species_level, filename= paste0(path_plots6, "/S4_Figure_SpeciesLevel.jpg"), width = 8, height = 8) # True common legend
ggsave(panels_species_level_legend, filename= paste0(path_plots6, "/S4_Figure_SpeciesLevel_ForLegend.jpg"), width = 8, height = 8) # False common legend




#===============================================================================
# Additional SM plots: ---------------------------------------------------------
#===============================================================================

# Native remaining local distinctiveness: --------------------------------------
#Dii_nat <- subset(Dii, Dii$SourceII == c("NR"))
#Dii_nat$Period <- as.character(Dii_nat$Period)
#Dii_nat$Period[Dii_nat$Period=="B) Contemporary Natives + Exotics"] <- "B) Contemporary (2005)"

#(dviolinNat <- ggplot(Dii_nat, aes(x = Period, y = Di, fill=Period)) +
#    geom_point(alpha=0.3)+
#    geom_violin(alpha=0.4) +
#    geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
#    xlab("Native remaining") +
#    ylab("Local distinctiveness") +
#    scale_fill_manual(values=c("Lightgray", "Darkgray"))+
#    theme_bw())
#ggsave(dviolinNat, filename= paste0(path_plots6, "/S4_LocalDisPanelNR.jpg"), width = 7, height = 5) 

## Uniqueness dimensions:-------------------------------------------------------
Ui_dim <- gather(Ui_dim, key="Trait", value="Uniqueness", -species)
identical(sort(unique(Ui_dim$species)), sort(unique(Uii$species))) # TRUE
Ui_dim$SourceII <- Uii$SourceII[match(Uii$species, Ui_dim$species)]
Ui_dim <- Ui_dim[!Ui_dim$Trait=="Ui_all",]

Ui_dimp <- ggplot(Ui_dim, aes(x=SourceII, y=Uniqueness, fill=SourceII, color=SourceII)) +
  geom_point(alpha=0.3)+
  geom_violin(alpha=0.4)+
  geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
  xlab("") +
  ylab("Uniqueness") +
  labs(fill="Status", color="Status")+
  scale_color_manual(values=pal4)+
  scale_fill_manual(values=pal4)+
  theme_bw()+
  facet_wrap(~Trait, scales = "free")
Ui_dimp

## Global distinctiveness dimensions:-------------------------------------------
Dig_dimgI <- t(do.call(rbind, Dig_dimg))
colnames(Dig_dimgI) <- names(Dig_dimg)
Dig_dimgI <- as.data.frame(Dig_dimgI)
Dig_dimgI$Species <- rownames(Dig_dimgI)

Dig_dimgII <- gather(Dig_dimgI, key="Trait", value="GD", -Species)
Dig_dimgII$SourceII <- Diig$SourceII[match(Diig$species, Dig_dimgII$Species)]
Dig_dimgII <- Dig_dimgII[!Dig_dimgII$Trait=="Di_all",]

Dig_dimp <- ggplot(Dig_dimgII, aes(x=SourceII, y=GD, fill=SourceII, color=SourceII)) +
  geom_point(alpha=0.3)+
  geom_violin(alpha=0.4)+
  geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
  xlab("") +
  ylab("Global Distinctiveness") +
  labs(fill="Status", color="Status")+
  scale_color_manual(values=pal4)+
  scale_fill_manual(values=pal4)+
  theme_bw()+
  facet_wrap(~Trait, scales = "free")
Dig_dimp

## Distinctiveness dimensions:--------------------------------------------------
Di_dim <- lapply(Di_dim, as.data.frame)
traits <- names(Di_dim)
traits
Di_dim <- lapply(Di_dim, function(x) {x$siteperiod <- rownames(x); x})

Di_dim2 <- lapply(Di_dim, function(x) {gather(x, key="Species", value="Distinctiveness", na.rm = TRUE, -c(96))}) # 96 is siteperiod
Di_dim2 <- do.call(rbind, Di_dim2)
Di_dim2$Trait <- str_split_fixed(rownames(Di_dim2), "\\.", 2)[,1]
sum(is.na(Di_dim2$Trait))
Di_dim2$SourceII <- rep(NA, nrow(Di_dim2))


sort(unique(Dii_v2$Period))
Dii_v2$Period[Dii_v2$Period=="Contemporary Natives + Exotics"] <- "ContAll"
Dii_v2$Period[Dii_v2$Period=="Historical"] <- "HNB"
Dii_v2$siteperiodspecies <- paste0(Dii_v2$Site, "_", Dii_v2$Period, "_", Dii_v2$Species)
Di_dim2$siteperiodspecies <- paste0(Di_dim2$siteperiod, "_", Di_dim2$Species)

Di_dim2$SourceII <-  Dii_v2$SourceII[match(Dii_v2$siteperiodspecies,
                                           Di_dim2$siteperiodspecies)]


Di_dim2 <- Di_dim2[!Di_dim2$Trait=="Di_all",]


Di_dim2p <- ggplot(Di_dim2, aes(x=SourceII, y=Distinctiveness, fill=SourceII, color=SourceII)) +
  geom_point(alpha=0.3)+
  geom_violin(alpha=0.4)+
  geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
  xlab("") +
  ylab("Local distinctiveness") +
  labs(fill="Status", color="Status")+
  scale_color_manual(values=pal4LD)+
  scale_fill_manual(values=pal4LD)+
  theme_bw()+
  facet_wrap(~Trait, scales = "free")
Di_dim2p 


ggsave(Ui_dimp, filename= paste0(path_plots6, "/S4_UniquenessDimensions.jpg"), width = 10, height = 8) 
ggsave(Dig_dimp, filename= paste0(path_plots6, "/S4_GlobalDisDimensions.jpg"), width = 10, height = 8)
ggsave(Di_dim2p, filename= paste0(path_plots6, "/S4_LocalDisDimensions.jpg"), width = 10, height = 8) 





# End of script ################################################################
