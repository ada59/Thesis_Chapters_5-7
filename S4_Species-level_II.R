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
library(car)
library(grid)
library(ggpubr)

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

load(paste0(myd, "/dist_mat1.RData"))   # Species x species distance matrix
load(paste0(myd, "/tt2.RData"))         # Species x trait matrix (Status, Source, Family)

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

tt2 <- tt2[order(rownames(tt2)),]
identical(colnames(dt), colnames(dist_mat1))  # TRUE
identical(colnames(dt), rownames(tt2))        # TRUE

tt2$Status[tt2$Status=="Native"] <- "Native Remaining" 

tt2$SourceII[tt2$SourceII=="Extirpated"] <- "E"
tt2$SourceII[tt2$SourceII=="Aquaculture&Sportfishing"] <- "IA"
tt2$SourceII[tt2$SourceII=="Aquarium&Contaminant"] <- "IB"
tt2$SourceII[tt2$SourceII=="Native"] <- "NR"

###########################################################################################
# Trait rarity measures:-------------------------------------------------------------------
Uii <- uniqueness(dt, dist_matrix = dist_mat1) # uniqueness (regional)
Uii$Ui_norm <- (Uii$Ui-min(Uii$Ui))/(max(Uii$Ui)-min(Uii$Ui)) # normalise
Ui_dim <- uniqueness_dimensions(dt, 
                                tt2[,c(1:10)],
                                metric="euclidean") 


Dii <- distinctiveness(dt, dist_matrix = dist_mat1) # distinctiveness (local)
Dii <- as.data.frame(Dii)
warnings()
Di_dim <- distinctiveness_dimensions(dt, 
                                     tt2[,c(1:10)],
                                     metric="euclidean") # warnings when coms have 1 sps only


Diig <- distinctiveness_global(dist_mat1) # global distinctiveness (regional)
Diig$global_di_norm <- (Diig$global_di-min(Diig$global_di))/(max(Diig$global_di)-min(Diig$global_di)) # normalise
dt1 <- as.data.frame(t(dt[1,]))
dt1[1,] <- 1
Dig_dimg <- distinctiveness_dimensions(as.matrix(dt1), 
                                       tt2[,c(1:10)],
                                       metric="euclidean") # warnings when coms have 1 sps only

###########################################################################################
# Period & Status -------------------------------------------------------------------------
# Uniqueness:------------------------------------------------------------------------------
Uii$SourceII <- tt2$SourceII[match(Uii$species, rownames(tt2))]

# Global distinctiveness:------------------------------------------------------------------
Diig$SourceII <- tt2$SourceII[match(Diig$species, rownames(tt2))]

# Distinctiveness:-------------------------------------------------------------------------
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

# Subsetting:
DiiI <- subset(Dii, ! Dii$SourceII == "NR")
DiiII <- subset(Dii, Dii$SourceII == "NR")

DiiIIinsp <- DiiII %>%
  group_by(Site, Species) %>%
  summarise(n=n())

DiiII_1 <- subset(DiiIIinsp, DiiIIinsp$n ==1)
DiiII_2 <- subset(DiiIIinsp, DiiIIinsp$n ==2)

LocE <- merge(x=DiiII_1,y=DiiII,by=c("Site", "Species"), all.x=FALSE, all.y=FALSE)
LocNR <- merge(x=DiiII_2,y=DiiII,by=c("Site", "Species"), all.x=FALSE, all.y=FALSE)

sum(LocE$Period=="Contemporary Natives + Exotics") # Some Translocated
sum(LocE$Period=="Historical")   #142

LocE$SourceII <- ifelse(LocE$Period=="Historical", "E", LocE$SourceII)
LocE$SourceII <- ifelse(LocE$Period=="Contemporary Natives + Exotics", "IB", LocE$SourceII)
LocE$SourceII <- ifelse(LocE$Species=="Chirostoma jordani" & LocE$Site=="Zumpango Lake", "NR", LocE$SourceII)
LocE$SourceII <- ifelse(LocE$Species=="Girardinichthys viviparus" & LocE$Site=="Zempoala Lagoon", "NR", LocE$SourceII)
LocE$SourceII <- ifelse(LocE$Species=="Xenoophorus captivus" & LocE$Site=="Stream at Agua de Enmedio", "NR", LocE$SourceII)

LocE <- within(LocE, rm(n))
LocNR <- within(LocNR, rm(n))

Dii_v2 <- rbind(DiiI, LocNR, LocE) # 507
save(Dii_v2, file="Dii_v2.RData")

length(unique(Dii_v2$Species)) #99
setdiff(unique(colnames(dt)), unique(Dii_v2$Species))
# "Gambusia yucatana" (it means, this species only appears in coms of richness = 1)
# Distinctiveness of this taxon cannot be calculated.

###########################################################################################
# Correlations: ---------------------------------------------------------------------------
count_sps_dis <- Dii_v2 %>%
  group_by(Species) %>%
  summarise(n=n())
range(count_sps_dis$n) # 1 37

sum(count_sps_dis$n==1) #23 species
sum(count_sps_dis$n==2) #28 species

Dis_av <- Dii_v2 %>%
  group_by(Species) %>%
  summarise(av=mean(Di))
Dis_av$av_norm <- (Dis_av$av-min(Dis_av$av))/(max(Dis_av$av)-min(Dis_av$av)) # normalise

Uniqx <- rbind(Uii[!Uii$species == "Gambusia yucatana",], Uii[Uii$species == "Gambusia yucatana",])
Diigx <- rbind(Diig[!Diig$species == "Gambusia yucatana",], Diig[Diig$species == "Gambusia yucatana",])
Gy <- data.frame("Species"="Gambusia yucatana","av"=NA, "av_norm"=NA)
Dis_avx <- rbind(Dis_av, Gy)

identical(Dis_avx$Species, Uniqx$species) #TRUE
identical(Dis_avx$Species, Diigx$species) #TRUE
shapiro.test(Dis_avx$av_norm)
shapiro.test(Uniqx$Ui_norm)
shapiro.test(Diigx$global_di_norm)
cor(data.frame(Dis_avx$av, Uniqx$Ui, Diigx$global_di), method="spearman", use = "complete.obs")

###########################################################################################
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

###########################################################################################
# Tests -----------------------------------------------------------------------------------
# Uniqueness by SourceII:------------------------------------------------------------------
Uii$SourceII <- as.factor(Uii$SourceII)
sum(is.na(Uii))
str(Uii)

Uii_fit1 <- aov(Ui~SourceII, data=Uii)
summary(Uii_fit1)
par(mfrow = c(2, 2))
plot(Uii_fit1)
shapiro.test(residuals(Uii_fit1))

Uii_fit2 <- kruskal.test(Ui~SourceII, data=Uii)
Uii_fit2 # 0.3

# Av Distinctiveness by SourceII:---------------------------------------------------------
boxplot(Di ~ SourceII, data=Dii_v2)
Dii_fit1 <- aov(Di ~ SourceII, data=Dii_v2)
summary(Dii_fit1)
par(mfrow = c(2, 2))
plot(Dii_fit1)
shapiro.test(residuals(Dii_fit1))


Dii_fit2 <- kruskal.test(Di ~ SourceII, data=Dii_v2)
Dii_fit2
pairwise.wilcox.test(Dii_v2$Di, Dii_v2$SourceII,
                     p.adjust.method = "none", exact=FALSE) 

# Global distinctiveness by Source II:-----------------------------------------------------
shapiro.test(residuals(aov(global_di_norm ~ SourceII, data=Diig)))
Diig_fit1 <- kruskal.test(global_di_norm ~ SourceII, data=Diig)
Diig_fit1
pairwise.wilcox.test(Diig$global_di_norm, Diig$SourceII,
                     p.adjust.method = "none")

p.adjust(c(0.3018, 2.667e-09, 2.615e-05), method = "BH")
p1 <- c(1.0e-09, 0.0036, 0.3596, 6.4e-08, 5.3e-07, 0.0475)
p2 <- c(0.00011, 0.01639, 0.03620, 6.3e-06, 4.9e-05, 0.03255)
p.adjust(c(p1, p2), method = "BH")

###########################################################################################
# Plots:-----------------------------------------------------------------------------------
dev.off()
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
}

# Plots Uniqueness ------------------------------------------------------------------------
uviolin2 <- violin_p(data=Uii, x=Uii$SourceII, y=Uii$Ui_norm,
                     pal=pal4,
                     legend_title = "Status",
                     labx="", laby="Uniqueness")
ggsave(uviolin2, filename= paste0(plot_dir, "/SM/Ms/S6_UniquenessPanel.jpg"), width = 7, height = 5) 
save(uviolin2, file="uviolin2.RData")

# Plots Distinctiveness -------------------------------------------------------------------
Dii_v2$SourceII[Dii_v2$SourceII=="NR"] <- "LR"
Dii_v2$SourceII[Dii_v2$SourceII=="E"] <- "LE"
Dii_v2$SourceII <- recode_factor(Dii_v2$SourceII, 
                                 "LE"="LE",
                                 "IA"="IA",
                                 "IB"="IB",
                                 "LR"="LR")
dviolin2 <- ggplot(Dii_v2, aes(x = SourceII, y = Di, fill=SourceII)) +
  geom_point(alpha=0.3)+
  geom_violin(alpha=0.4) +
  geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
  xlab("") +
  ylab("Local distinctiveness") +
  labs(fill="Status")+
  scale_fill_manual(values=pal4LD)+
  theme_bw()
print(dviolin2)

ggsave(dviolin2, filename= paste0(plot_dir, "/SM/Ms/S6_LocalDisPanel.jpg"), width = 7, height = 5) 
save(dviolin2, file="dviolin2.RData")

# Plots Global Distinctiveness-------------------------------------------------------------
dgviolin2 <- violin_p(data=Diig, x=Diig$SourceII, y=Diig$global_di_norm,
                     pal=pal4,
                     legend_title = "Status",
                     labx="", laby="Global distinctiveness")

ggsave(dgviolin2, filename= paste0(plot_dir, "/SM/Ms/S6_GlobDisPanel.jpg"), width = 7, height = 5) 
save(dgviolin2, file="dgviolin2.RData")

###########################################################################################
# Assembly of final plots -----------------------------------------------------------------
load("Main.RData")
(panels_species <- ggarrange(Main,
                             uviolin2,
                             dgviolin2,
                             dviolin2,
                             common.legend = TRUE, # FALSE (Then panels_speciesII)
                             legend="bottom",
                             labels = c("A)", "B)", "C)", "D)"),
                             align="hv",
                             font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))
ggsave(panels_species, filename= paste0(plot_dir, "/S6_Figure2.jpg"), width = 8, height = 8) # True common legend
ggsave(panels_speciesII, filename= paste0(plot_dir, "/S6_Figure2II.jpg"), width = 8, height = 8) 

###########################################################################################
# Additional SM plots: --------------------------------------------------------------------

# Native remaining local distinctiveness: -------------------------------------------------
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

#ggsave(dviolinNat, filename= paste0(plot_dir, "/SM/Ms/S6_LocalDisPanelNR.jpg"), width = 7, height = 5) 

# Uniqueness dimensions:-------------------------------------------------------------------
Ui_dim <- gather(Ui_dim, key="Trait", value="Uniqueness", -species)
identical(sort(unique(Ui_dim$species)), sort(unique(Uii$species))) #TRUE
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

# Global distinctiveness dimensions:------------------------------------------------------
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

# Distinctiveness dimensions:-------------------------------------------------------------
Di_dim <- lapply(Di_dim, as.data.frame)
traits <- names(Di_dim)
Di_dim <- lapply(Di_dim, function(x) {x$siteperiod <- rownames(x); x})

Di_dim2 <- lapply(Di_dim, function(x) {gather(x, key="Species", value="Distinctiveness", na.rm = TRUE, -c(101))})
Di_dim2 <- do.call(rbind, Di_dim2)
Di_dim2$Trait <- rep(traits, each=507)
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


ggsave(Ui_dimp, filename= paste0(plot_dir, "/SM/Ms/S6_UniquenessDimensions.jpg"), width = 10, height = 8) 
ggsave(Dig_dimp, filename= paste0(plot_dir, "/SM/Ms/S6_GlobalDisDimensions.jpg"), width = 10, height = 8)
ggsave(Di_dim2p, filename= paste0(plot_dir, "/SM/Ms/S6_LocalDisDimensions.jpg"), width = 10, height = 8) 

###########################################################################################
# End of script ###########################################################################
