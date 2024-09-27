#===============================================================================
# Script: Observed vs null trait diversity & redundancy and temporal change
# August/Sept 2022 & August 2024
#===============================================================================




#===============================================================================
# Libraries:--------------------------------------------------------------------
#===============================================================================
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
library(lmtest)

rm(list=ls())
myd <- getwd()

path_lists <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists"   # Raw data lists
path_traits <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Traits" # Dir to trait objects
path_plots6 <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots/Chapter6"   # Dir to save plots Chapter 6




#===============================================================================
# Community data:---------------------------------------------------------------
#===============================================================================
load(paste0(path_lists,"/l67.RData"))     # main list with community data matrices (67 sites)
load(paste0(path_lists, "/l53.RData"))    # main list with community data matrices (53 sites)
load(paste0(path_lists, "/NDT67.RData"))  # List for checks on raw data
load(paste0(path_lists, "/NDT83.RData"))  # List for checks on raw data
load(paste0(path_traits, "/dist_mat1.RData"))


load("div3.RData")
load("l3.RData")
load("pred_redundancy.RData")



#===============================================================================
# read hill functions:----------------------------------------------------------
#===============================================================================
# FD and TD (Hill numbers) use the same units and TD is always greater than or equal to FD
# For q=0, Redundancy is S-FD0
hill_dir <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/Functions/FunD-master"

fskt <- read.table(paste0(hill_dir,"/Plant_Abundance.txt"))
fskt = list(fs=fskt$Fushan, kt=fskt$Kenting)
dij_fskt = read.table(paste0(hill_dir,"/Plant_Distance_Matrix.txt"))
dij_fskt = as.matrix(dij_fskt) # All needed for the code to run

source("C:/Users/afe1/OneDrive - University of St Andrews/PHD/Functions/FunD_Rcode.R") # warning or error here is OK, doesn't affect FD_MLE




#===============================================================================
# Random assemblages: ----------------------------------------------------------
#===============================================================================
#hnb <- subset(l67, l67$Period=="HNB")
#hnb <- within(hnb, rm(SiteNameE, DrainageBasinE, Period))
#all_sps <- names(hnb)         # 95, includes exotics
#hnb <- hnb[,!colSums(hnb)==0] # 77 sps
#sort(unique(names(hnb)))
#all_natives <- names(hnb)
#setdiff(all_sps, all_natives)

#species_list <- all_natives  # or all sps

#rand <- list()
#rand17 <- list() # 17 sps max observed in any site and period

#for (j in 1:999){
#  for(i in 1:17){
#    df <- as.data.frame(matrix(ncol=77, nrow=1)) # or 95, if all_sps
#    names(df) <- species_list
#    samp <- sample(species_list, size = i, replace = FALSE)
#    df[as.character(samp)] <- 1
#    df[is.na(df)] <- 0
#    rand17[[i]] <- df
#  } 
#  rand[[j]] <- do.call(rbind, rand17)
#} 
#rand <- do.call(rbind, rand)
#sum(rowSums(rand)==17) # 999, OK


#randNatives <- rand
#save(randNatives, file="randNatives.RData")


#randAll <- rand
#save(randAll, file="randAll.RData")




#===============================================================================
# Compute Null FD & Redundancy: ------------------------------------------------
#===============================================================================
#load("randNatives.RData")
#load("randAll.RData")

#rand <- randNatives  # or randAll

#rand <- rand[, order(names(rand))]
#dist_mat1 <- dist_mat1[order(rownames(dist_mat1)), order(colnames(dist_mat1))]
#avTau <- mean(dist_mat1[dist_mat1>0])
#dist_mat1 <- dist_mat1[rownames(dist_mat1) %in% all_natives, colnames(dist_mat1) %in% all_natives]
#identical(colnames(dist_mat1), names(rand)) # TRUE

#rand_l <- split(rand, f=rownames(rand))

#F0null <- list()
#for (i in 1:length(rand_l)){
#  assemblage <- as.matrix(t(rand_l[[i]]))
#  T0_n <- colSums(assemblage)
#  F0_n <- FD_MLE(assemblage, dist_mat1, avTau, q=0)
#  R0_n <- (T0_n - F0_n)
#  F0null[[i]] <- cbind(T0_n, F0_n, R0_n)
#}
#F0null <- as.data.frame(do.call(rbind,  F0null))


#F0nullNatives <- F0null
#save(F0nullNatives, file="F0nullNatives.RData")


#F0nullAll <- F0null
#save(F0nullAll, file="F0nullAll.RData")




#===============================================================================
# Plots:------------------------------------------------------------------------
#===============================================================================
load("F0nullNatives.RData")
load("F0nullAll.RData")

F0null <- F0nullNatives

(nullF0 <- ggplot(data=F0null, aes(x=as.factor(T0_n), y=F0_n))+
   geom_boxplot(alpha=0.3, color="gray69")+
   stat_boxplot(geom ='errorbar', color="gray69") + 
   xlab("SR")+
   ylab("T0")+
   theme_bw()) # null SR vs T0 (17)  [HISTORICAL]
(nullF0II <- ggplot(data=F0null[F0null$T0_n<10,], aes(x=as.factor(T0_n), y=F0_n))+
    geom_boxplot(alpha=0.3, color="gray69")+
    stat_boxplot(geom ='errorbar', color="gray69") + 
    xlab("SR")+
    ylab("T0")+
    theme_bw()) # null SR vs T0 (9) [CONTEMPORARY]


div3 <- div3[!div3$Period=="ContE",]   # rm Cont E from rawdata
div3$Period <- recode_factor(div3$Period, 
                             HNC = "Historical Conservative", 
                             HNB = "Historical Broad",
                             ContN = "Contemporary Native",
                             ContAll = "Contemporary Native + Introduced")
div3$Period <- factor(div3$Period, levels = c("Historical Conservative",
                                              "Historical Broad",
                                              "Contemporary Native", 
                                              "Contemporary Native + Introduced"))
unique(pred_redundancy$Period)
str(pred_redundancy)

#### Historical Conservative ####
(p1nullF0_hnc <- nullF0 + 
    geom_point(data = div3[div3$Period=="Historical Conservative",], 
               aes(x=T0, y=F0), size=2, alpha=1, col="#28E2E5")+
    geom_line(data = pred_redundancy[pred_redundancy$Period=="Historical Conservative",], aes(x = T0, y = fit), linetype = 1, size = 1, col="#28E2E5")+
    geom_ribbon(data = pred_redundancy[pred_redundancy$Period=="Historical Conservative",], aes(x = T0, y=fit, ymin = lwr, ymax = upr), alpha=0.03))

#### Historical Broad ####
(p1nullF0_hnb <- nullF0 + 
    geom_point(data = div3[div3$Period=="Historical Broad",], 
               aes(x=T0, y=F0), size=2, alpha=1, col="#61D04F")+
    geom_line(data = pred_redundancy[pred_redundancy$Period=="Historical Broad",], aes(x = T0, y = fit), linetype = 1, size = 1, col="#61D04F")+
    geom_ribbon(data = pred_redundancy[pred_redundancy$Period=="Historical Broad",], aes(x = T0, y=fit, ymin = lwr, ymax = upr), alpha=0.03))

#### Contemporary Native ####
(p1nullF0_contN <- nullF0II + 
   geom_point(data = div3[div3$Period=="Contemporary Native",], 
              aes(x=T0, y=F0), size=2, alpha=1, col="#F5C710")+
   geom_line(data = pred_redundancy[pred_redundancy$Period=="Contemporary Native",], aes(x = T0, y = fit), linetype = 1, size = 1, col="#F5C710")+
   geom_ribbon(data = pred_redundancy[pred_redundancy$Period=="Contemporary Native",], aes(x = T0, y=fit, ymin = lwr, ymax = upr), alpha=0.03))

#### Contemporary Native + Introduced #### 
(p1nullF0_contAll <- nullF0II + 
    geom_point(data = div3[div3$Period=="Contemporary Native + Introduced",], 
               aes(x=T0, y=F0), size=2, alpha=1, col="#CD0BBC")+
    geom_line(data = pred_redundancy[pred_redundancy$Period=="Contemporary Native + Introduced",], aes(x = T0, y = fit), linetype = 1, size = 1, col="#CD0BBC")+
    geom_ribbon(data = pred_redundancy[pred_redundancy$Period=="Contemporary Native + Introduced",], aes(x = T0, y=fit, ymin = lwr, ymax = upr), alpha=0.03))


ggsave(p1nullF0_hnc, file= paste0(path_plots6, "/Null-ObservedHNC.jpg"), width = 8, height = 6) 
ggsave(p1nullF0_hnb, file= paste0(path_plots6, "/Null-ObservedHNB.jpg"), width = 8, height = 6) 
ggsave(p1nullF0_contN, file= paste0(path_plots6, "/Null-ObservedContN.jpg"), width = 8, height = 6) 
ggsave(p1nullF0_contAll, file= paste0(path_plots6, "/Null-ObservedContAll.jpg"), width = 8, height = 6) 




#===============================================================================
# SES:--------------------------------------------------------------------------
#===============================================================================
null_l <- split(F0null, f=F0null$T0_n)
null_l_mean <- lapply(null_l, function(x) {mean(x$F0_n)})  # averages for each species richness level
null_l_sd <- lapply(null_l, function(x) {sd(x$F0_n)})      # sd for each species richness level

#sensi_sites <- as.vector(l3[[3]]$SiteName)  # 33
#l3sensitivity <- lapply(l3, function(x) x[x$SiteName %in% sensi_sites,])
ll <- lapply(l3, function(x) {split(x, f=x$T0)})

vecs <- lapply(ll, function(x) {as.vector(names(x))})
means <- lapply(vecs, function(x) {null_l_mean[x]})
sds <- lapply(vecs, function(x) {null_l_sd[x]})


ses_results_hnc <- mapply(function(x,y,z){((x$F0 - y[1])/z[1])}, ll[[1]], means[[1]], sds[[1]], SIMPLIFY = FALSE)
ses_results_hnb <- mapply(function(x,y,z){((x$F0 - y[1])/z[1])}, ll[[2]], means[[2]], sds[[2]], SIMPLIFY = FALSE)
ses_results_cn <- mapply(function(x,y,z){((x$F0 - y[1])/z[1])}, ll[[3]], means[[3]], sds[[3]], SIMPLIFY = FALSE)
ses_results_call <- mapply(function(x,y,z){((x$F0 - y[1])/z[1])}, ll[[4]], means[[4]], sds[[4]], SIMPLIFY = FALSE)

as_v <- function(x){as.vector(unlist(x))}
ses_resultsF0 <- c(as_v(ses_results_hnc),as_v(ses_results_hnb),
                   as_v(ses_results_cn),as_v(ses_results_call))


period <- rep(c("Historical Conservative", "Historical Broad", "Contemporary Native", "Contemporary Native + Introduced"), times=c(48,60,33,49))
ses_dt <- data.frame("Period"=period, "SESF0"=ses_resultsF0)
ses_dt <- ses_dt[!is.na(ses_dt$SESF0),]
str(ses_dt)
#ses_dt %>% group_by(Period) %>% summarise(count = n(),
#                                          countsAbv = sum(SESF0 > 0, na.rm = TRUE),
#                                          countsUnd = sum(SESF0 < 0, na.rm = TRUE),
#                                          countsAbvPer = (sum(SESF0 > 0, na.rm = TRUE)/count)*100,
#                                          countsUndPer = (sum(SESF0 < 0, na.rm = TRUE)/count)*100)

summF0 <- summarySE(ses_dt, measurevar=c("SESF0"), groupvars=c("Period"), na.rm = TRUE)
summF0$Period <- recode_factor(summF0$Period,
                               "Historical Conservative"="HC",
                               "Historical Broad"="HB",
                               "Contemporary Native"="CN",
                               "Contemporary Native + Introduced"="CN+I")
ses_dt$Period <- recode_factor(ses_dt$Period,
                               "Historical Conservative"="HC",
                               "Historical Broad"="HB",
                               "Contemporary Native"="CN",
                               "Contemporary Native + Introduced"="CN+I")

colorBlind4   <- c("#28E2E5", "#61D04F", "#F5C710", "#CD0BBC")
colorBlind3   <- c("#61D04F", "#F5C710", "#CD0BBC")

(summF0hnc <- ggplot(ses_dt[ses_dt$Period=="HC",], 
                     aes(x=Period, y=SESF0, colour=Period)) +
    geom_point(alpha=0.4, col="#28E2E5") +
    geom_violin(alpha=0.1, col="#28E2E5", fill="#28E2E5") +
    #geom_boxplot(width=0.1, alpha=0.5) +
    geom_point(data=summF0[summF0$Period=="HC",], 
               aes(x=Period, y=SESF0), size=2.5, col="#28E2E5") +
    geom_errorbar(data=summF0[summF0$Period=="HC",],
                  aes(ymin=SESF0-ci, ymax=SESF0+ci),
                  width=0.2, col="#28E2E5")+
    theme_bw()+
    theme(legend.position = "none")+
    labs(y="SES", x="")+
    geom_hline(yintercept = 0, linetype="dashed"))

(summF03p <- ggplot(ses_dt[!ses_dt$Period=="HC",],
                    aes(x=Period, y=SESF0, colour=Period)) + 
    geom_point(alpha=0.5) +
    geom_violin(alpha=0.1, aes(fill=Period)) +
    geom_point(data=summF0[!summF0$Period=="HC",], 
               aes(x=Period, y=SESF0, colour=Period), size=2.5)+
    geom_errorbar(data=summF0[!summF0$Period=="HC",], aes(ymin=SESF0-ci, ymax=SESF0+ci),
                  width=0.2)+
    theme_bw()+
    theme(legend.position = "bottom")+
    scale_color_manual(values=colorBlind3, labels=c("HB (n=60)", "CN (n=33)", "CN+I (n=49)"))+
    scale_fill_manual(values=colorBlind3, labels=c("HB (n=60)", "CN (n=33)", "CN+I (n=49)"))+
    labs(y="SES", x="Period")+
    geom_hline(yintercept = 0, linetype="dashed")) 

ggsave(summF0hnc, file= paste0(path_plots6, "/SES_F0hnc.jpg"), width = 3, height = 4) 
ggsave(summF03p, file= paste0(path_plots6, "/SES_F03.jpg"), width = 6, height = 4)




#===============================================================================
# T-test, different from 0?-----------------------------------------------------
#===============================================================================
ses_dt_n <- ses_dt %>%
  group_by(Period) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
ses_dtF0 <- ses_dt %>%
  group_by(Period) %>%
  summarise(res = list(t.test(SESF0, mu=0, alternative="less")))
ses_dtF0$res  # same p-vals when assessing redundancy




#===============================================================================
# Quantile scores:--------------------------------------------------------------
#===============================================================================
divrank <- div3[,-c(5)] # rm redundancy
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
l$Period <- recode_factor(l$Period, 
                          "Historical Conservative"="HC",
                          "Historical Broad"="HB",
                          "Contemporary Native"="CN",
                          "Contemporary Native + Introduced"="CN+I")
lSigN <- l[l$pvals<0.025,]
lSigP <- l[l$pvals>0.975,]

(pvalsp <- ggplot(l, aes(x=Period, y=pvals, colour=Period)) + 
    geom_violin(aes(fill=Period), alpha=0.1)+
    geom_boxplot(width=0.1, alpha=0.5)+
    theme_bw()+
    labs(y="P-values")+
    scale_color_manual(values=colorBlind4)+
    scale_fill_manual(values=colorBlind4)+
    theme(legend.position = "none")+
    #geom_hline(yintercept = 0.025, linetype="dashed", color="red", alpha=0.5)+
    #geom_hline(yintercept = 0.975, linetype="dashed", color="red", alpha=0.5)+
    geom_hline(yintercept = 0.75, linetype="dashed", color="lightgray")+
    geom_hline(yintercept = 0.25, linetype="dashed", color="lightgray")+
    geom_hline(yintercept = 0.5, linetype="dashed", color="lightgray"))

ggsave(pvalsp, file= paste0(path_plots6, "/QuantileScores_F0.jpg"), width = 6, height = 4) 





#===============================================================================
# Main Figures Assemblage Level ------------------------------------------------
#===============================================================================

## Main Figure: ----------------------------------------------------------------
(MainAssemblageLevel <- ggarrange(p1nullF0_hnb,
                        p1nullF0_contN,
                        p1nullF0_contAll,
                        summF03p,
                        labels=c("A)", "B)", "C)", "D)"),
                        common.legend = T,
                        legend="bottom",
                        ncol=2,
                        nrow=2,
                        font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))

ggsave(MainAssemblageLevel, file= paste0(path_plots6, "/S5_FigureAssemblageLevel.jpg"), width =8, height = 8) 


## Supplementary figure:--------------------------------------------------------
(HNCres <- ggarrange(p1nullF0_hnc,
                     summF0hnc,
                     align="hv",
                     widths = c(1, 0.5),
                     labels=c("A)", "B)"),
                     nrow=1,
                     font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "bottom")))
ggsave(HNCres, file= paste0(path_plots6, "/S5_HNC.jpg"), width =8, height = 5) 





# End of script ################################################################
