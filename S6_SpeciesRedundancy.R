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
Di_dim <- distinctiveness_dimensions(dt, 
                                     tt2[,c(1:10)],
                                     metric="euclidean") # warnings when coms have 1 sps only
Diig <- distinctiveness_global(dist_mat1) # global distinctiveness (regional)
Diig$global_di_norm <- (Diig$global_di-min(Diig$global_di))/(max(Diig$global_di)-min(Diig$global_di)) # normalise

###########################################################################################
# Period & Status -------------------------------------------------------------------------
# Uniqueness
Uii$Status <- tt2$Status[match(Uii$species, rownames(tt2))]
Uii$SourceII <- tt2$SourceII[match(Uii$species, rownames(tt2))]

# Distinctiveness
Dii$Period <- rep(NA, nrow(Dii))
Dii$Period <- ifelse(rownames(Dii) %like% "HNB", "Historical", Dii$Period)
Dii$Period <- ifelse(rownames(Dii) %like% "ContAll", "Contemporary Natives + Exotics", Dii$Period)
sum(is.na(Dii$Period))

Dii <- gather(Dii, key="Species", value="Di", -Period)
Dii <- Dii[!is.na(Dii$Di),]
length(unique(Dii$Species)) #99
setdiff(unique(colnames(dt)), unique(Dii$Species))
# "Gambusia yucatana" (it means, this species only appears in coms of richness = 1)
# Distinctiveness of this taxon cannot be calculated.
Dii$Period <- as.factor(Dii$Period)
Dii$Period <- recode_factor(Dii$Period, "Historical"="A) Historical", 
                            "Contemporary Natives + Exotics"="B) Contemporary Natives + Exotics")

Dii$Status <- tt2$Status[match(Dii$Species, rownames(tt2))]
Dii$SourceII <- tt2$SourceII[match(Dii$Species, rownames(tt2))]

# Global distinctiveness:
Diig$Status <- tt2$Status[match(Diig$species, rownames(tt2))]
Diig$SourceII <- tt2$SourceII[match(Diig$species, rownames(tt2))]

###########################################################################################
# Correlations: ---------------------------------------------------------------------------
count_sps_dis <- Dii %>%
  group_by(Species, SourceII, Status) %>%
  summarise(n=n())
range(count_sps_dis$n) # 1 37

sum(count_sps_dis$n==1) #23 species
sum(count_sps_dis$n==2) #28 species
sum(count_sps_dis$Status=="Extirpated") #29

sum(count_sps_dis$n[count_sps_dis$SourceII=="E"]==1) #14
sum(count_sps_dis$n[count_sps_dis$SourceII=="E"]==2) #11
14+11 # 25 of 29

Dis_av <- Dii %>%
  group_by(Species, SourceII, Status) %>%
  summarise(av=mean(Di))

Uniqx <- Uii[!Uii$species == "Gambusia yucatana",]
Diigx <- Diig[!Diig$species == "Gambusia yucatana",]

identical(Dis_av$Species, Uniqx$species) #TRUE
identical(Dis_av$Species, Diigx$species) #TRUE
cor(data.frame(Dis_av$av, Uniqx$Ui, Diigx$global_di))

###########################################################################################
# Tests -----------------------------------------------------------------------------------
# Uniqueness by Status:--------------------------------------------------------------------
Uii$Status <- as.factor(Uii$Status)
Uii$SourceII <- as.factor(Uii$SourceII)
sum(is.na(Uii))
str(Uii)

boxplot(Ui~Status, data=Uii)
Uii_fit1 <- aov(Ui ~ Status, data=Uii)
summary(Uii_fit1)

par(mfrow = c(2, 2))
plot(Uii_fit1)
dev.off()

shapiro.test(residuals(Uii_fit1)) # Shapiro-Wilk test of normality
# p-value = 5.507e-10

Uii$logUi <- log(Uii$Ui)
Uii_fit2 <- aov(logUi ~ Status, data=Uii)
summary(Uii_fit2)

par(mfrow = c(2, 2))
plot(Uii_fit2)
shapiro.test(residuals(Uii_fit2))
# p-value = 0.0008383

Uii_fit3 <- kruskal.test(Ui ~ Status, data = Uii)
Uii_fit3

# Uniqueness by SourceII:------------------------------------------------------------------
Uii_fit4 <- aov(Ui~SourceII, data=Uii)
summary(Uii_fit4)
par(mfrow = c(2, 2))
plot(Uii_fit4)

Uii_fit5 <- kruskal.test(Ui~SourceII, data=Uii)
Uii_fit5

# Av Distinctiveness by Status:------------------------------------------------------------
Dis_av$Status <- as.factor(Dis_av$Status)
Dis_av$SourceII <- as.factor(Dis_av$SourceII)
sum(is.na(Dis_av$av)) #0
str(Dis_av)

Dii_fit1 <- aov(av ~ Status, data=Dis_av)
summary(Dii_fit1)

par(mfrow = c(2, 2))
plot(Dii_fit1)
shapiro.test(residuals(Dii_fit1))

Dii_fit2 <- kruskal.test(av ~ Status, data=Dis_av)
Dii_fit2

# Av Distinctiveness by SourceII:---------------------------------------------------------
boxplot(av ~ SourceII, data=Dis_av)
Dii_fit3 <- kruskal.test(av ~ SourceII, data=Dis_av)
Dii_fit3

pairwise.wilcox.test(Dis_av$av, Dis_av$SourceII,
                     p.adjust.method = "BH")

# Global distinctiveness by Status:--------------------------------------------------------
dev.off()
boxplot(global_di_norm ~ Status, data=Diig)
Diig_fit1 <- kruskal.test(global_di_norm ~ Status, data=Diig)
Diig_fit1

# Global distinctiveness by Source II:-----------------------------------------------------
Diig_fit2 <- kruskal.test(global_di_norm ~ SourceII, data=Diig)
Diig_fit2
pairwise.wilcox.test(Diig$global_di_norm, Diig$SourceII,
                     p.adjust.method = "BH")

###########################################################################################
# Plots:-----------------------------------------------------------------------------------
dev.off()
pal3 <- c("#0072B2", "#F0E442", "darkgray")
pal4 <- c("#0072B2","#D55E00","#E69F00","darkgray")

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
uviolin <- violin_p(data=Uii, x=Uii$Status, y=Uii$Ui_norm,
                    pal=pal3,
                    legend_title = "Status",
                    labx="Status", laby="Uniqueness")
uviolin2 <- violin_p(data=Uii, x=Uii$SourceII, y=Uii$Ui_norm,
                     pal=pal4,
                     legend_title = "Status",
                     labx="", laby="Uniqueness")

ggsave(uviolin, filename= paste0(plot_dir, "/SM/Thesis/UniquenessStatus.jpg"), width = 7, height = 5) 
ggsave(uviolin2, filename= paste0(plot_dir, "/SM/Ms/UniquenessPanel.jpg"), width = 7, height = 5) 
save(uviolin2, file="uviolin2.RData")

# Plots Distinctiveness -------------------------------------------------------------------
dviolin <- violin_p(data=Dii, x=Dii$Status, y=Dii$Di,
                    pal=pal3,
                    legend_title = "Status",
                    labx="Status", laby="Local Trait Rarity (Distinctiveness)")
dviolin2 <- violin_p(data=Dii, x=Dii$SourceII, y=Dii$Di,
                     pal=pal4,
                     legend_title = "Status",
                     labx="Status", laby="Distinctiveness")

ggsave(dviolin, filename= paste0(plot_dir, "/SM/Thesis/LocalDisStatus.jpg"), width = 7, height = 5) 
ggsave(dviolin2, filename= paste0(plot_dir, "/SM/Ms/LocalDisSource.jpg"), width = 7, height = 5) 
save(dviolin2, file="dviolin2.RData")

# Plots Global Distinctiveness-------------------------------------------------------------
dgviolin <- violin_p(data=Diig, x=Diig$Status, y=Diig$global_di_norm,
                    pal=pal3,
                    legend_title = "Status",
                    labx="Status", laby="Regional Trait Rarity (Global distinctiveness)")
dgviolin2 <- violin_p(data=Diig, x=Diig$SourceII, y=Diig$global_di_norm,
                     pal=pal4,
                     legend_title = "Status",
                     labx="Status", laby="Global Distinctiveness")


ggsave(dgviolin, filename= paste0(plot_dir, "/SM/Thesis/GlobalDissStatus.jpg"), width = 7, height = 5) 
ggsave(dgviolin2, filename= paste0(plot_dir, "/SM/Ms/GlobalDissSource.jpg"), width = 7, height = 5) 
save(dgviolin2, file="dgviolin2.RData")

###########################################################################################
# Additional SM plots: --------------------------------------------------------------------
Dii_nat <- subset(Dii, Dii$Status == c("Native Remaining"))
Dii_nat$Period <- as.character(Dii_nat$Period)
Dii_nat$Period[Dii_nat$Period=="B) Contemporary Natives + Exotics"] <- "B) Contemporary (2005)"

(dviolinNat <- ggplot(Dii_nat, aes(x = Period, y = Di, fill=Period)) +
    geom_point(alpha=0.3)+
    geom_violin(alpha=0.4) +
    geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+
    xlab("Native remaining") +
    ylab("Distinctiveness") +
    scale_fill_manual(values=c("Lightgray", "Darkgray"))+
    theme_bw())

ggsave(dviolinNat, filename= paste0(plot_dir, "/SM/Ms/DisctinctivenessNativeRemaining.jpg"), width = 8, height = 5) 

###########################################################################################
# Assembly of final plots -----------------------------------------------------------------
load("Main.RData")
(panels_species <- ggarrange(Main,
                             uviolin2,
                             dviolin2,
                             dgviolin2,
                             common.legend = TRUE, 
                             legend="bottom",
                             labels = c("A)", "B)", "C)", "D)"),
                             align="hv",
                             font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))
ggsave(panels_species, filename= paste0(plot_dir, "/Panels_Species.jpg"), width = 8, height = 8) 

###########################################################################################
# End of script ###########################################################################
