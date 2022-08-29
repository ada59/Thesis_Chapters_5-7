###########################################################################################
# Script: Compute redundancies and test for temporal change
# AFE
# August 2022
###########################################################################################

# Libraries:------------------------------------------------------------------------------
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
library(Rmisc)

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
save(ContAll, file="ContAll.RData")

l <- bind_rows(HNC, HNB, ContN, ContAll)
l[is.na(l)] <- 0
sum(is.na(l)) # 0

l$Period <- rep(c("HNC", "HNB", "ContN", "ContAll"), each=67)
All  <- l
save(All, file="All.RData")
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
# Observed trends: -----------------------------------------------------------------------
# Subsets:
hnc <- subset(div, div$Period=="HNC")
hnb <- subset(div, div$Period=="HNB")
contN <- subset(div, div$Period=="ContN")
contAll <- subset(div, div$Period=="ContAll")

str(hnc)

sum(hnc$T0==0)
sum(hnb$T0==0)
sum(contN$T0==0)   #21
sum(contAll$T0==0) #14 (67-53, OK)

contN <- contN[!contN$T0==0,]       # 46
contAll <- contAll[!contAll$T0==0,] # 53

# Historical conservative:
hnc_1 <- lm(F0 ~ T0, data=hnc)
hnc_2 <- lm(F0 ~ log(T0+1), data=hnc)

AIC(hnc_1, hnc_2) # Curve

# Historical broad:
hnb_1 <- lm(F0 ~ T0, data=hnb)
hnb_2 <- lm(F0 ~ log(T0+1), data=hnb)

AIC(hnb_1, hnb_2)

# Current native:
contN_1 <- lm(F0 ~ T0, data=contN)
contN_2 <- lm(F0 ~ log(T0+1), data=contN)

AIC(contN_1, contN_2)

# Current All:
contAll_1 <- lm(F0 ~ T0, data=contAll)
contAll_2 <- lm(F0 ~ log(T0+1), data=contAll)

AIC(contAll_1, contAll_2)

# NOTE:
# A curve is a better fit in all cases.

div$Period <- recode_factor(div$Period, HNC  = "A) Historical conservative", 
                                                HNB = "B) Historical broad",
                                                ContN = "C) Contemporary Natives",
                                                ContAll = "D) Contemporary Natives + Exotics")

# Plot trends:
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

ggsave(p1, file= paste0(plot_dir, "/RedundancyFreeFit.jpg"), width = 12, height = 10) 
ggsave(p2, file= paste0(plot_dir, "/RedundancyLogLinear.jpg"), width = 8, height = 6) 


# Remove locality with richness = 17 in the historical period.

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
# Null TD and FD: ------------------------------------------------------------------------

#Natives <- names(HNB[,-c(1:2)])
All <- unique(c(names(HNB[,-c(1:2)]), names(ContAll[,-c(1:2)]))) 

rand <- list()
rand17 <- list()

# Natives:
#for (j in 1:999){
#  for(i in 1:17){
#    df <- as.data.frame(matrix(ncol=78, nrow=1))
#    names(df) <- Natives
#    samp <- sample(Natives, size = i, replace = FALSE)
#    df[as.character(samp)] <- 1
#    df[is.na(df)] <- 0
#    rand17[[i]] <- df
#  } 
#  rand[[j]] <- do.call(rbind, rand17)
#} 

rand_all <- list()
rand17_all <- list()

#All:
for (j in 1:999){
  for(i in 1:17){
    df <- as.data.frame(matrix(ncol=100, nrow=1))
    names(df) <- All
    samp <- sample(All, size = i, replace = FALSE)
    df[as.character(samp)] <- 1
    df[is.na(df)] <- 0
    rand17_all[[i]] <- df
  } 
  rand_all[[j]] <- do.call(rbind, rand17_all)
} 

#randNatives <- do.call(rbind, rand)
randAll <- do.call(rbind, rand_all)

#save(randNatives, file="randNatives.RData")
save(randAll, file="randAll.RData")

#load(paste0(myd, "/randNatives.RData"))
load(paste0(myd, "/randAll.RData"))

##########################################################################################
# Compute FD: ----------------------------------------------------------------------------

#names(randNatives)[names(randNatives)=="Agonostomus monticola"] <- "Dajaus monticola"
#names(randAll)[names(randAll)=="Agonostomus monticola"] <- "Dajaus monticola"

#randNatives <- randNatives[, order(names(randNatives))]
randAll <- randAll[, order(names(randAll))]

dist_mat1 <- dist_mat1[order(rownames(dist_mat1)), order(colnames(dist_mat1))]
identical(colnames(dist_mat1), names(randAll)) # TRUE

#rem <- setdiff(colnames(dist_mat1), names(randNatives))
#dist_mat2 <- dist_mat1[!rownames(dist_mat1) %in% rem,!colnames(dist_mat1) %in% rem]
#identical(colnames(dist_mat2), names(randNatives)) # TRUE

randAll_l <- split(randAll, f=rownames(randAll))
#randNatives_l <- split(randNatives, f=rownames(randNatives))

FD0All <- list()
for (i in 1:length(randAll_l)){
  assemblage <- as.matrix(t(randAll_l[[i]]))
  TD0_n <- colSums(assemblage)
  FD0_n <- FD_MLE(assemblage, dist_mat1, mean(dist_mat1[dist_mat1>0]), q=0)
  FD0All[[i]] <- cbind(TD0_n, FD0_n)
}

#FD0Natives <- list()
#for (i in 1:length(randNatives_l)){
#  assemblage <- as.matrix(t(randNatives_l[[i]]))
#  TD0_n <- colSums(assemblage)
#  FD0_n <- FD_MLE(assemblage, dist_mat2, mean(dist_mat2[dist_mat2>0]), q=0)
#  FD0Natives[[i]] <- cbind(TD0_n, FD0_n)
#}


RandFD0All <- as.data.frame(do.call(rbind,  FD0All))
#RandFD0Natives <- as.data.frame(do.call(rbind,  FD0Natives))

save(RandFD0All, file="RandFD0All.RData")
#save(RandFD0Natives, file="RandFD0Natives.RData")

load(paste0(myd, "/RandFD0All.RData"))
#load(paste0(myd, "/RandFD0Natives.RData"))

##########################################################################################
# Plots:----------------------------------------------------------------------------------
# https://waterdata.usgs.gov/blog/boxplots/

sum(div$T0==0) # 35
div <- div[!div$T0==0,]


# Null model (all fish)
(nullAll <- ggplot(data=RandFD0All, aes(x=as.factor(TD0_n), y=FD0_n))+
    geom_boxplot(alpha=0.5, color="gray49")+
    stat_boxplot(geom ='errorbar', color="gray49") + 
    xlab("T0")+
    ylab("F0")+
    theme_minimal()) 

(p1 <- nullAll + geom_point(data = div, aes(x=T0,y=F0, color=Period), size=2, alpha=0.5)+
    geom_smooth(data = div, aes(x=T0, y=F0, color=Period), method=lm, formula = y ~ log(x+1), se=TRUE)+
    #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "green"))+
    facet_wrap(~Period))


# Null model (only native fish)
#(nullNatives <- ggplot(data=RandFD0Natives, aes(x=as.factor(TD0_n), y=FD0_n))+
#    geom_boxplot()+
#    stat_boxplot(geom ='errorbar') + 
#    xlab("TD0")+
#    ylab("FD0")+
#    theme_minimal()) 
#(p2 <- nullNatives + geom_point(data = div, aes(x=T0,y=F0, color=Period), size=3, alpha=0.5)+
#    geom_smooth(data = div, aes(x=T0, y=F0, color=Period), method=lm, formula = y ~ log(x+1), se=TRUE)+
#    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "green"))+
#    facet_wrap(~Period))

ggsave(p1, file= paste0(plot_dir, "/Null-ObservedAll.jpg"), width = 10, height = 7) 
#ggsave(p2, file= paste0(plot_dir, "/Null-ObservedNatives.jpg"), width = 12, height = 10) 

###########################################################################################
# SES:-------------------------------------------------------------------------------------
# SES > 0 suggests the prevalence of trait divergence patterns
# SES < 0 suggests the prevalence of trait convergence (REDUNDANCY)

rall_l <- split(RandFD0All, f=RandFD0All$TD0_n)
rall_l_mean <- lapply(rall_l, function(x) {mean(x$FD0_n)})
rall_l_sd <- lapply(rall_l, function(x) {sd(x$FD0_n)})

#rnat_l <- split(RandFD0Natives, f=RandFD0Natives$TD0_n)
#rnat_l_mean <- lapply(rnat_l, function(x) {mean(x$FD0_n)})
#rnat_l_sd <- lapply(rnat_l, function(x) {sd(x$FD0_n)})

ll <- list(hnc, hnb, contN, contAll)
ll <- lapply(ll, function(x) {split(x, f=x$T0)})
vecs <- lapply(ll, function(x) {as.vector(names(x))})

rall_means <- lapply(vecs, function(x) {rall_l_mean[x]})
rall_sds <- lapply(vecs, function(x) {rall_l_sd[x]})

#rnat_means <- lapply(vecs, function(x) {rnat_l_mean[x]})
#rnat_sds <- lapply(vecs, function(x) {rnat_l_sd[x]})

# SES (All species null model)
ses_results_hnc <- mapply(function(x,y,z){((x$F0 - y)/z)}, ll[[1]], rall_means[[1]], rall_sds[[1]], SIMPLIFY = FALSE)
ses_results_hnb <- mapply(function(x,y,z){((x$F0 - y)/z)}, ll[[2]], rall_means[[2]], rall_sds[[2]], SIMPLIFY = FALSE)
ses_results_cn <- mapply(function(x,y,z){((x$F0 - y)/z)}, ll[[3]], rall_means[[3]], rall_sds[[3]], SIMPLIFY = FALSE)
ses_results_call <- mapply(function(x,y,z){((x$F0 - y)/z)}, ll[[4]], rall_means[[4]], rall_sds[[4]], SIMPLIFY = FALSE)


# SES (Native species null model)
#ses_results_hnc_n <- mapply(function(x,y,z){((x$F0 - y)/z)}, ll[[1]], rnat_means[[1]], rnat_sds[[1]], SIMPLIFY = FALSE)
#ses_results_hnb_n <- mapply(function(x,y,z){((x$F0 - y)/z)}, ll[[2]], rnat_means[[2]], rnat_sds[[2]], SIMPLIFY = FALSE)
#ses_results_cn_n <- mapply(function(x,y,z){((x$F0 - y)/z)}, ll[[3]], rnat_means[[3]], rnat_sds[[3]], SIMPLIFY = FALSE)
#ses_results_call_n <- mapply(function(x,y,z){((x$F0 - y)/z)}, ll[[4]], rnat_means[[4]], rnat_sds[[4]], SIMPLIFY = FALSE)

# Plot:
ses_results_all <- c(as.vector(unlist(ses_results_hnc)),
                 as.vector(unlist(ses_results_hnb)),
                 as.vector(unlist(ses_results_cn)),
                 as.vector(unlist(ses_results_call)))
#ses_results_nat <- c(as.vector(unlist(ses_results_hnc_n)),
#                     as.vector(unlist(ses_results_hnb_n)),
#                    as.vector(unlist(ses_results_cn_n)),
#                     as.vector(unlist(ses_results_call_n)))


period <- rep(c("A) Historical Conservative", "B) Historical Broad", "C) Current Native", "D) Current Native + Exotic"), times=c(67,67,46,53))

ses.dt <- data.frame("Period"=period, "SESAll"=ses_results_all)
ses.dt$SESAll[is.nan(ses.dt$SESAll)] <- NA
#ses.dt$SESNat[is.nan(ses.dt$SESNat)] <- NA

sum(is.na(ses.dt$SESAll)) #43
#sum(is.na(ses.dt$SESNat)) #43

# https://www.cyclismo.org/tutorial/R/confidence.html
# 95% confidence interval (normal dis To Be Revised)

summaryAll <- summarySE(ses.dt, measurevar=c("SESAll"), groupvars=c("Period"), na.rm = TRUE)
#summaryNat <- summarySE(ses.dt, measurevar=c("SESNat"), groupvars=c("Period"), na.rm = TRUE)


(psesAll <- ggplot(ses.dt, aes(x=Period, y=SESAll, colour=Period)) + 
  geom_point(alpha=0.5)+
  geom_point(data=summaryAll, aes(x=Period, y=SESAll, colour=Period), size=2.5)+
  geom_errorbar(data=summaryAll, aes(ymin=SESAll-ci, ymax=SESAll+ci),
                width=0.2)+
  theme_bw()+
  labs(y="SES", x="Period")+
  geom_hline(yintercept = 0, linetype="dashed"))

#(psesNat <- ggplot(ses.dt, aes(x=Period, y=SESNat, colour=Period)) + 
#    geom_point(alpha=0.2)+
#   geom_point(data=summaryNat, aes(x=Period, y=SESNat, colour=Period), size=2.5)+
#   geom_errorbar(data=summaryNat, aes(ymin=SESNat-ci, ymax=SESNat+ci),
#                  width=0.2)+
#   theme_bw()+
#   geom_hline(yintercept = 0, linetype="dashed"))
# Removed 43 rows containing missing values (geom_point), OK

ggsave(psesAll, file= paste0(plot_dir, "/SES_All.jpg"), width = 8, height = 5) 
#ggsave(psesNat, file= paste0(plot_dir, "/SES_Natives.jpg"), width = 12, height = 10) 


###########################################################################################
# Quantile scores:-------------------------------------------------------------------------
# The test is two sided:
# P -values less than or equal to 0.025:significantly less F0 than expected
# Greater than or equal to 0.975: significantly more F0 than expected
div$Period <- recode_factor(div$Period,
                            "A) Historical conservative"="A) HNC",
                            "B) Historical broad"="B) HNB",
                            "C) Contemporary Natives"="C) CN",
                            "D) Contemporary Natives + Exotics"="D) CAll")

divII <- split(div, f=div$T0)
vdivII <- names(divII)
nrows <- lapply(divII, function(x) {nrow(x)})

rall_l2 <- lapply(rall_l, function(x) {x$FD0_n})
rall_l2 <- rall_l2[vdivII]

rall_l3 <- mapply(function(x,y) {as.data.frame(t(replicate(x, y)))}, nrows,rall_l2, SIMPLIFY = FALSE)
mat.rank.all <- mapply(function(x,y) {data.frame(cbind(x,y))}, divII, rall_l3, SIMPLIFY = FALSE)

mat.rank.all <- do.call(rbind, mat.rank.all)

l <- list()
for (i in 1:nrow(mat.rank.all)){
  ranks <- rank(c(mat.rank.all[i,4], as.vector(as.numeric(mat.rank.all[i,c(5:1003)]))))[1]
  pvals <- ranks/1000
  l[[i]] <- data.frame(cbind(ranks, pvals))
}
l <- do.call(rbind, l)
l <- data.frame(mat.rank.all[,c(1:4)], l)

l <- l[!l$T0==1,]
(pvalsAll <- ggplot(l, aes(x=Period, y=pvals, colour=Period)) + 
    geom_violin(aes(fill=Period), alpha=0.3)+
    geom_boxplot(width=0.1, alpha=0.5)+
    theme_bw()+
    labs(y="P-values")+
    geom_hline(yintercept = 0.025, linetype="dashed", color="red", alpha=0.5)+
    geom_hline(yintercept = 0.975, linetype="dashed", color="red", alpha=0.5)+
    geom_hline(yintercept = 0.75, linetype="dashed", color="lightgray")+
    geom_hline(yintercept = 0.25, linetype="dashed", color="lightgray")+
    geom_hline(yintercept = 0.5, linetype="dashed", color="lightgray"))

ggsave(pvalsAll, file= paste0(plot_dir, "/pvals_All.jpg"), width = 10, height = 5) 


###########################################################################################
# End of script ###########################################################################
