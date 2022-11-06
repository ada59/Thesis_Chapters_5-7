###########################################################################################
# Script: Observed vs null trait diversity & redundancy and temporal change
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
library(grid)

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
divT0 <- lapply(ls, function(x) {FD_MLE(x, dist_mat1, min(dist_mat1[dist_mat1>0]), 0)}) # taxonomic diversity
divF0 <- lapply(ls, function(x) {FD_MLE(x, dist_mat1, mean(dist_mat1[dist_mat1>0]), 0)}) # trait diversity

divT0d <- as.data.frame(do.call(rbind, divT0))
divT0d$SiteName <- str_split_fixed(rownames(divT0d), "_", 2)[,1]
divT0d$Period <- str_split_fixed(rownames(divT0d), "_", 2)[,2]

divF0d <- as.data.frame(do.call(rbind, divF0))
divF0d$SiteName <- str_split_fixed(rownames(divF0d), "_", 2)[,1]
divF0d$Period <- str_split_fixed(rownames(divF0d), "_", 2)[,2]

div <- left_join(divT0d, divF0d, by=c("SiteName", "Period"))
names(div) <- c("T0", "SiteName", "Period", "F0")
div <- div %>% relocate(T0, .after = Period)
str(div)

div$SiteName <- as.factor(div$SiteName)
div$Period <- as.factor(div$Period)

div$R0 <- div$T0-div$F0 # trait redundancy

##########################################################################################
# Observed trends: -----------------------------------------------------------------------
# Subsets:
div2 <- div[!div$T0==0,]
save(div2, file = "div2.RData")
hnc <- subset(div2, div2$Period=="HNC")
hnb <- subset(div2, div2$Period=="HNB")
contN <- subset(div2, div2$Period=="ContN")
contAll <- subset(div2, div2$Period=="ContAll")

# Relationship between T0 & Period:-------------------------------------------------------
str(div2)
div2_1T0 <- lm(T0 ~ Period, data=div2)
summary(div2_1T0)

div2B <- subset(div2, div2$Period %in% c("HNB", "ContAll"))
div2_2T0 <- lm(T0 ~ Period, data=div2B)
summary(div2_2T0)

# Relationship between T0 & F0:-----------------------------------------------------------
# Historical conservative:----------------------------------------------------------------
hnc_1 <- lm(F0 ~ T0, data=hnc)
summary(hnc_1)
hnc_2 <- lm(F0 ~ log(T0+1), data=hnc)
summary(hnc_2)

hnc_aic <- as.data.frame(aictab(list("Linear"=hnc_1, "Logarithmic"=hnc_2))) # Curve

# Historical broad:-----------------------------------------------------------------------
hnb_1 <- lm(F0 ~ T0, data=hnb)
summary(hnb_1)
hnb_2 <- lm(F0 ~ log(T0+1), data=hnb)
summary(hnb_2)

hnb_aic <- as.data.frame(aictab(list("Linear H"=hnb_1, "Logarithmic H"=hnb_2))) # Curve

# Current native:-------------------------------------------------------------------------
contN_1 <- lm(F0 ~ T0, data=contN)
summary(contN_1)
contN_2 <- lm(F0 ~ log(T0+1), data=contN)
summary(contN_2)

contN_aic <- as.data.frame(aictab(list("Linear Cont N"=contN_1, "Logarithmic Cont N"=contN_2))) # Curve

# Current All:----------------------------------------------------------------------------
contAll_1 <- lm(F0 ~ T0, data=contAll)
summary(contAll_1)
contAll_2 <- lm(F0 ~ log(T0+1), data=contAll)
summary(contAll_2)

contAll_aic <- as.data.frame(aictab(list("Linear Cont N + E"=contAll_1, "Logarithmic Cont N + E"=contAll_2))) # Curve

# NOTE: A curve is a better fit in all cases.

listaic <- data.frame(rbind(hnb_aic[,c(1,3:4)],contN_aic[,c(1,3:4)], contAll_aic[,c(1,3:4)]))
listModels <- list(hnb_2, hnb_1, contN_2, contN_1, contAll_2, contAll_1)
extract_m <- function(mod=NULL) {
  x <- mod
  df <- x$df.residual
  fstat <- summary(x)$fstatistic[1]
  adjr2 <- summary(x)$adj.r.squared
  pval <- summary(x)$coefficients[2,4]
  df <- as.data.frame(cbind(df, fstat, adjr2, pval))
  return(df)
}                         
l_ext_m <- lapply(listModels, extract_m)
l_ext_m <- do.call(rbind, l_ext_m)

res_m_selection <- data.frame(cbind(listaic, l_ext_m))
write.csv2(res_m_selection, file="res_m_selection.csv", row.names = F)
##########################################################################################
# Assumptions: ---------------------------------------------------------------------------
# Histograms:-----------------------------------------------------------------------------
par(mfrow = c(2, 2))
hist(hnc_2$residuals)
hist(hnb_2$residuals)
hist(contN_2$residuals)
hist(contAll_2$residuals)
dev.off()

# Plots:----------------------------------------------------------------------------------
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
plot(hnc_2)
plot(hnb_2)
plot(contN_2)
plot(contAll_2)
dev.off()

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

div_subset_main <- subset(div, div$Period %in% c("B) Historical broad", "C) Contemporary Natives",
                                                 "D) Contemporary Natives + Exotics"))
div_subset_main$Period <- droplevels(div_subset_main$Period)
div_subset_main$Period <- recode_factor(div_subset_main$Period, 
                                        "B) Historical broad" = "Historical",
                                        "C) Contemporary Natives" = "Cont. Natives",
                                        "D) Contemporary Natives + Exotics" = "Cont. Natives + Exotics")

colorBlind4   <- c("#F0E442", "#E69F00", "#009E73","#0072B2")
colorBlind3   <- c("#E69F00", "#009E73","#0072B2")
scales::show_col(colorBlind4)

(p1 <- ggplot(div, aes(x=T0, y=F0, color=Period))+
  geom_point(aes(color=Period))+
  geom_smooth()+
  theme_classic()+
  scale_color_manual(values=colorBlind4)+
  facet_wrap(~Period)) # free fit

(p2 <- ggplot(div, aes(x=T0, y=F0, color=Period))+
    geom_point(aes(color=Period))+
    geom_smooth(method=lm, formula = y ~ log(x+1))+
    theme_classic()+
    scale_color_manual(values=colorBlind4)+
    facet_wrap(~Period)) # logarithmic fit

ggsave(p1, file= paste0(plot_dir, "/SM/Thesis/F0vsT0FreeFit4.jpg"), width = 8, height = 6) 
ggsave(p2, file= paste0(plot_dir, "/SM/Thesis/F0vsT0Logarithmic4.jpg"), width = 8, height = 6) 

##########################################################################################
# Observed trends (R0 vs T0):-------------------------------------------------------------
# Historical conservative:----------------------------------------------------------------
hnc_1R0 <- lm(R0 ~ T0, data=hnc)
summary(hnc_1R0)
hnc_2R0 <- lm(R0 ~ log(T0+1), data=hnc)

AIC(hnc_1R0, hnc_2R0) # Linear

# Historical broad:-----------------------------------------------------------------------
hnb_1R0 <- lm(R0 ~ T0, data=hnb)
summary(hnb_1R0)
hnb_2R0 <- lm(R0 ~ log(T0+1), data=hnb)

AIC(hnb_1R0, hnb_2R0) # Linear

# Current native:-------------------------------------------------------------------------
contN_1R0 <- lm(R0 ~ T0, data=contN)
summary(contN_1R0)
contN_2R0 <- lm(R0 ~ log(T0+1), data=contN)

AIC(contN_1R0, contN_2R0) # Linear

# Current All:----------------------------------------------------------------------------
contAll_1R0 <- lm(R0 ~ T0, data=contAll)
summary(contAll_1R0)
contAll_2R0 <- lm(R0 ~ log(T0+1), data=contAll)

AIC(contAll_1R0, contAll_2R0) # Linear
# NOTE: A linear regression is a better fit in all cases ( R0 vs T0)

##########################################################################################
# Assumptions: ---------------------------------------------------------------------------
# Histograms:-----------------------------------------------------------------------------
par(mfrow = c(2, 2))
hist(hnc_1R0$residuals)
hist(hnb_1R0$residuals)
hist(contN_1R0$residuals)
hist(contAll_1R0$residuals)
dev.off()

# Plots:----------------------------------------------------------------------------------
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
plot(hnc_1R0)
plot(hnb_1R0)
plot(contN_1R0)
plot(contAll_1R0)
dev.off()

##########################################################################################
# Plot trends:----------------------------------------------------------------------------

(p1R0 <- ggplot(div, aes(x=T0, y=R0, color=Period))+
    geom_point(aes(color=Period))+
    geom_smooth(method=lm, formula = y ~ x)+
    theme_classic()+
    scale_color_manual(values=colorBlind4)+
    facet_wrap(~Period)) 

ggsave(p1R0, file= paste0(plot_dir, "/SM/Thesis/RedundancyLinear4.jpg"), width = 8, height = 6)

##########################################################################################
# Null models: ---------------------------------------------------------------------------
species_list <- unique(c(names(HNB[,-c(1:2)]), names(ContAll[,-c(1:2)]))) 

rand <- list()
rand17 <- list() # 17 sps max observed in any site and period

for (j in 1:999){
  for(i in 1:17){
    df <- as.data.frame(matrix(ncol=100, nrow=1))
    names(df) <- species_list
    samp <- sample(species_list, size = i, replace = FALSE)
    df[as.character(samp)] <- 1
    df[is.na(df)] <- 0
    rand17[[i]] <- df
  } 
  rand[[j]] <- do.call(rbind, rand17)
} 
rand <- do.call(rbind, rand)
sum(rowSums(rand)==17) #999, OK
save(rand, file="rand.RData")

# Compute FD & Redundancy: ---------------------------------------------------------------------
rand <- rand[, order(names(rand))]
dist_mat1 <- dist_mat1[order(rownames(dist_mat1)), order(colnames(dist_mat1))]
identical(colnames(dist_mat1), names(rand)) # TRUE

rand_l <- split(rand, f=rownames(rand))

F0null <- list()
for (i in 1:length(rand_l)){
  assemblage <- as.matrix(t(rand_l[[i]]))
  T0_n <- colSums(assemblage)
  F0_n <- FD_MLE(assemblage, dist_mat1, mean(dist_mat1[dist_mat1>0]), q=0)
  R0_n <- (T0_n - F0_n)
  F0null[[i]] <- cbind(T0_n, F0_n, R0_n)
}
F0null <- as.data.frame(do.call(rbind,  F0null))
save(F0null, file="F0null.RData")

##########################################################################################
# Plots:----------------------------------------------------------------------------------
sum(div$T0==0) # 35
div <- div[!div$T0==0,]
div_subset_main <- div_subset_main[!div_subset_main$T0==0,]

(nullF0 <- ggplot(data=F0null, aes(x=as.factor(T0_n), y=F0_n))+
    geom_boxplot(alpha=0.5, color="gray49")+
    stat_boxplot(geom ='errorbar', color="gray49") + 
    #geom_smooth(aes(x=T0_n, y=F0_n), col="gray49", linetype="dashed")+ # default is gam
    xlab("SR")+
    ylab("T0")+
    theme_minimal()) # null SR vs T0
(nullR0 <- ggplot(data=F0null, aes(x=as.factor(T0_n), y=R0_n))+
    geom_boxplot(alpha=0.5, color="gray49")+
    stat_boxplot(geom ='errorbar', color="gray49") + 
    xlab("SR")+
    ylab("R0")+
    theme_minimal()) # null SR vs R0

(p1nullF0 <- nullF0 + geom_point(data = div, 
                                  aes(x=T0, y=F0, color=Period), size=2, alpha=0.5)+
    geom_smooth(data = div, aes(x=T0, y=F0, color=Period), method=lm, formula = y ~ log(x+1), se=TRUE)+
    scale_color_manual(values=colorBlind4)+
    facet_wrap(~Period)) # SM plot, trait diversity
(p1nullF0_3 <- nullF0 + geom_point(data = div_subset_main, 
                                 aes(x=T0, y=F0, color=Period), size=2, alpha=0.5)+
    geom_smooth(data = div_subset_main, aes(x=T0, y=F0, color=Period), method=lm, formula = y ~ log(x+1), se=TRUE)+
    scale_color_manual(values=colorBlind3)+
    theme(legend.position = "none",
          axis.text.x = element_text(size=8))+
    facet_wrap(~Period)) # Main ms plot, trait diversity
(p1nullR0 <- nullR0 + geom_point(data = div, 
                                 aes(x=T0, y=R0, color=Period), size=2, alpha=0.5)+
    geom_smooth(data = div, aes(x=T0, y=R0, color=Period), method=lm, formula = y ~ x, se=TRUE)+
    scale_color_manual(values=colorBlind4)+
    facet_wrap(~Period)) # SM plot, trait redundancy
(p1nullR0_3 <- nullR0 + geom_point(data = div_subset_main, 
                                   aes(x=T0, y=R0, color=Period), size=2, alpha=0.5)+
    geom_smooth(data = div_subset_main, aes(x=T0, y=R0, color=Period), method=lm, formula = y ~ x, se=TRUE)+
    scale_color_manual(values=colorBlind3)+
    theme(legend.position = "none",
          axis.text.x = element_text(size=8))+
    facet_wrap(~Period)) # Main ms plot, trait redundancy

ggsave(p1nullF0, file= paste0(plot_dir, "/Null-ObservedF0.jpg"), width = 8, height = 6) 
ggsave(p1nullR0, file= paste0(plot_dir, "/Null-ObservedR0.jpg"), width = 8, height = 6) 

###########################################################################################
# SES:-------------------------------------------------------------------------------------
null_l <- split(F0null, f=F0null$T0_n)
null_l_mean <- lapply(null_l, function(x) {c(mean(x$F0_n),
                                           mean(x$R0_n))})
null_l_sd <- lapply(null_l, function(x) {c(sd(x$F0_n),
                                           sd(x$R0_n))})
ll <- list(hnc, hnb, contN, contAll)
ll <- lapply(ll, function(x) {split(x, f=x$T0)})
vecs <- lapply(ll, function(x) {as.vector(names(x))})

means <- lapply(vecs, function(x) {null_l_mean[x]})
sds <- lapply(vecs, function(x) {null_l_sd[x]})

ses_results_hnc <- mapply(function(x,y,z){((x$F0 - y[1])/z[1])}, ll[[1]], means[[1]], sds[[1]], SIMPLIFY = FALSE)
ses_results_hnb <- mapply(function(x,y,z){((x$F0 - y[1])/z[1])}, ll[[2]], means[[2]], sds[[2]], SIMPLIFY = FALSE)
ses_results_cn <- mapply(function(x,y,z){((x$F0 - y[1])/z[1])}, ll[[3]], means[[3]], sds[[3]], SIMPLIFY = FALSE)
ses_results_call <- mapply(function(x,y,z){((x$F0 - y[1])/z[1])}, ll[[4]], means[[4]], sds[[4]], SIMPLIFY = FALSE)

ses_results_hncR0 <- mapply(function(x,y,z){((x$R0 - y[2])/z[2])}, ll[[1]], means[[1]], sds[[1]], SIMPLIFY = FALSE)
ses_results_hnbR0 <- mapply(function(x,y,z){((x$R0 - y[2])/z[2])}, ll[[2]], means[[2]], sds[[2]], SIMPLIFY = FALSE)
ses_results_cnR0 <- mapply(function(x,y,z){((x$R0 - y[2])/z[2])}, ll[[3]], means[[3]], sds[[3]], SIMPLIFY = FALSE)
ses_results_callR0 <- mapply(function(x,y,z){((x$R0 - y[2])/z[2])}, ll[[4]], means[[4]], sds[[4]], SIMPLIFY = FALSE)

as_v <- function(x){as.vector(unlist(x))}
ses_resultsF0 <- c(as_v(ses_results_hnc),as_v(ses_results_hnb),
                   as_v(ses_results_cn),as_v(ses_results_call))
ses_resultsR0 <- c(as_v(ses_results_hncR0),as_v(ses_results_hnbR0),
                   as_v(ses_results_cnR0),as_v(ses_results_callR0))

period <- rep(c("A) Historical Conservative", "B) Historical Broad", "C) Current Native", "D) Current Native + Exotic"), times=c(67,67,46,53))
ses_dt <- data.frame("Period"=period, "SESF0"=ses_resultsF0, "SESR0"=ses_resultsR0)
ses_dt <- ses_dt[!is.na(ses_dt$SESF0),]

summF0 <- summarySE(ses_dt, measurevar=c("SESF0"), groupvars=c("Period"), na.rm = TRUE)
summR0 <- summarySE(ses_dt, measurevar=c("SESR0"), groupvars=c("Period"), na.rm = TRUE)

(summF0p <- ggplot(ses_dt, aes(x=Period, y=SESF0, colour=Period)) + 
    geom_point(alpha=0.5)+
    geom_point(data=summF0, aes(x=Period, y=SESF0, colour=Period), size=2.5)+
    geom_errorbar(data=summF0, aes(ymin=SESF0-ci, ymax=SESF0+ci),
                  width=0.2)+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_manual(values=colorBlind4)+
    labs(y="SES", x="Period")+
    geom_hline(yintercept = 0, linetype="dashed"))

(summR0p <- ggplot(ses_dt, aes(x=Period, y=SESR0, colour=Period)) + 
    geom_point(alpha=0.5)+
    geom_point(data=summR0, aes(x=Period, y=SESR0, colour=Period), size=2.5)+
    geom_errorbar(data=summR0, aes(ymin=SESR0-ci, ymax=SESR0+ci),
                  width=0.2)+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_manual(values=colorBlind4)+
    labs(y="SES", x="Period")+
    geom_hline(yintercept = 0, linetype="dashed"))

ggsave(summF0p, file= paste0(plot_dir, "/SM/Thesis/SES_F0.jpg"), width = 6, height = 4) 
ggsave(summR0p, file= paste0(plot_dir, "/SM/Thesis/SES_R0.jpg"), width = 6, height = 4)

###########################################################################################
# One-tailed t test:---------------------------------------------------------------------------
ses_dt_n <- ses_dt %>%
  group_by(Period) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
ses_dtF0 <- ses_dt[,-c(3)]
ses_dtF0 <- ses_dtF0 %>%
  group_by(Period) %>%
  summarise(res = list(t.test(SESF0, mu=0, alternative="less")))
ses_dtF0$res  

ses_dtR0 <- ses_dt[,-c(2)]
ses_dtR0 <- ses_dtR0 %>%
  group_by(Period) %>%
  summarise(res = list(t.test(SESR0, mu=0, alternative="greater")))
ses_dtR0$res  

###########################################################################################
# Quantile scores:-------------------------------------------------------------------------
divrank <- div[,-c(5)]
divrank <- split(divrank, f=divrank$T0)
vrank <- names(divrank)
nrows <- lapply(divrank, function(x) {nrow(x)})

null_lrank <- lapply(null_l, function(x) {x$F0_n})
null_lrank <- null_lrank[vrank]

null_lrank2 <- mapply(function(x,y) {as.data.frame(t(replicate(x, y)))}, nrows, null_lrank, SIMPLIFY = FALSE)
mat.rank.all <- mapply(function(x,y) {data.frame(cbind(x,y))}, divrank, null_lrank2, SIMPLIFY = FALSE)
mat.rank.all <- do.call(rbind, mat.rank.all)

l <- list()
for (i in 1:nrow(mat.rank.all)){
  ranks <- rank(as.vector(as.numeric(mat.rank.all[i,c(4:1003)])))[1]
  pvals <- ranks/1000
  l[[i]] <- data.frame(cbind(ranks, pvals))
}
l <- do.call(rbind, l)
l <- data.frame(mat.rank.all[,c(1:4)], l) 

l <- l[!l$T0==1,]
(pvalsp <- ggplot(l, aes(x=Period, y=pvals, colour=Period)) + 
    geom_violin(aes(fill=Period), alpha=0.3)+
    geom_boxplot(width=0.1, alpha=0.5)+
    theme_bw()+
    labs(y="P-values")+
    scale_color_manual(values=colorBlind4)+
    scale_fill_manual(values=colorBlind4)+
    theme(legend.position = "none")+
    geom_hline(yintercept = 0.025, linetype="dashed", color="red", alpha=0.5)+
    geom_hline(yintercept = 0.975, linetype="dashed", color="red", alpha=0.5)+
    geom_hline(yintercept = 0.75, linetype="dashed", color="lightgray")+
    geom_hline(yintercept = 0.25, linetype="dashed", color="lightgray")+
    geom_hline(yintercept = 0.5, linetype="dashed", color="lightgray"))

ggsave(pvalsp, file= paste0(plot_dir, "/SM/Thesis/Pvals_F0.jpg"), width = 6, height = 4) 
ggsave(pvalsp, file= paste0(plot_dir, "/SM/Ms/Pvals_F0.jpg"), width = 6, height = 4) 

##########################################################################################
# Figures (Ms) ---------------------------------------------------------------------------

# Main Figure 3: -------------------------------------------------------------------------
(panelsR0T0SR <- ggarrange(p1nullF0_3 + rremove("xlab"),
                           p1nullR0_3 + theme(strip.background = element_blank(),
                                              strip.text.x = element_blank()),
                           align="hv",
                           nrow=2,
                           labels = c("A)", "B)"),
                           font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))
ggsave(panelsR0T0SR, file= paste0(plot_dir, "/panelsR0T0SR.jpg"), width = 8, height = 5) 

# Main Figure 4: -------------------------------------------------------------------------
recfact <- function(x){recode_factor(x, "B) Historical Broad" = "Historical",
                                        "C) Current Native" = "Cont. Natives",
                                        "D) Current Native + Exotic" = "Cont. Natives + Exotics")
  
}
summF03 <- summF0[!summF0$Period=="A) Historical Conservative",]
summR03 <- summR0[!summR0$Period=="A) Historical Conservative",]
ses_dt3 <- ses_dt[!ses_dt$Period=="A) Historical Conservative",]
summF03$Period <- droplevels(recfact(summF03$Period))
summR03$Period <- droplevels(recfact(summR03$Period))
ses_dt3$Period <- droplevels(recfact(ses_dt3$Period))

(summF03p <- ggplot(ses_dt3, aes(x=Period, y=SESF0, colour=Period)) + 
    geom_point(alpha=0.5, size=2)+
    geom_point(data=summF03, aes(x=Period, y=SESF0, colour=Period), size=3.5)+
    geom_errorbar(data=summF03, aes(ymin=SESF0-ci, ymax=SESF0+ci),
                  width=0.2)+
    theme_bw()+
    ylim(c(-3,2.8))+
    theme(legend.position = "none")+
    scale_color_manual(values=colorBlind3)+
    labs(y="SES", x="Period", title="")+
    geom_hline(yintercept = 0, linetype="dashed"))
(summR03p <- ggplot(ses_dt3, aes(x=Period, y=SESR0, colour=Period)) + 
    geom_point(alpha=0.5, size=2)+
    geom_point(data=summR03, aes(x=Period, y=SESR0, colour=Period), size=3.5)+
    geom_errorbar(data=summR03, aes(ymin=SESR0-ci, ymax=SESR0+ci),
                  width=0.2)+
    theme_bw()+
    ylim(c(-3,2.8))+
    theme(legend.position = "none")+
    scale_color_manual(values=colorBlind3)+
    labs(y="SES", x="Period", title ="")+
    geom_hline(yintercept = 0, linetype="dashed"))
(panelsSES <- ggarrange(summF03p + rremove("xlab") + rremove("ylab"),
                        summR03p + rremove("xlab") + rremove("ylab"),
                        align="hv",
                        labels=c("A)", "B)"),
                        nrow=1,
                        font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))

(annotate_figure(panelsSES, left = textGrob("SES", rot = 90, hjust = 0, vjust = 1, gp = gpar(cex = 1)),
                bottom = textGrob("Period", hjust= 0.3, vjust=0.1, gp = gpar(cex = 1))))
ggsave(panelsSES, file= paste0(plot_dir, "/panelsSES.jpg"), width =8, height = 5) 

# Supplementary material, results for HNC (yellow):---------------------------------------
# TBC

##########################################################################################
# End of script ##########################################################################
