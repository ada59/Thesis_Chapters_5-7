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
library(lmtest)

rm(list=ls())
myd <- getwd()
plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots" # Dir to save main plots

# Community data:-------------------------------------------------------------------------
load(paste0(myd, "/HNC.RData"))        # Conservative Historical Native Assemblage
load(paste0(myd, "/HNB.RData"))        # Broad Historical Native Assemblage
load(paste0(myd, "/ContN.RData"))      # Contemporary Native Assemblage
load(paste0(myd, "/ContE.RData"))      # Contemporary Exotic assemblage
load(paste0(myd, "/ContAll.RData"))    # Contemporary Natives + Exotics
load(paste0(myd, "/dist_mat1.RData"))  # Trait distance matrix
load(paste0(myd, "/lsS5.RData"))       # List of assemblages

load("div3.RData")
load("l3.RData")
load("hnc_f.RData")
load("hnb_f.RData")
load("contN_f.RData")
load("contAll_f.RData")

# Functions:------------------------------------------------------------------------------
# FD and TD (Hill numbers) use the same units and TD is always greater than or equal to FD
# For q=0, Redundancy is S-FD0
hill_dir <- "C:/Users/Usuario/Documents/PHD/Functions/FunD-master/"

fskt <- read.table(paste0(hill_dir,"Plant_Abundance.txt"))
fskt = list(fs=fskt$Fushan, kt=fskt$Kenting)
dij_fskt = read.table(paste0(hill_dir,"Plant_Distance_Matrix.txt"))
dij_fskt = as.matrix(dij_fskt) # All needed for the code to run

source("C:/Users/Usuario/Documents/PHD/Functions/FunD_Rcode.R") # Warnings don't affect FD_MLE()

##########################################################################################
# Null models: ---------------------------------------------------------------------------
species_list <- unique(c(names(HNB[,-c(1:2)]), names(ContAll[,-c(1:2)]))) 

colnames(dist_mat1)[colnames(dist_mat1)=="Dajaus monticola"] <-"Agonostomus monticola"
rownames(dist_mat1)[rownames(dist_mat1)=="Dajaus monticola"] <-"Agonostomus monticola"

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
(nullF0 <- ggplot(data=F0null, aes(x=as.factor(T0_n), y=F0_n))+
   geom_boxplot(alpha=0.3, color="gray69")+
   stat_boxplot(geom ='errorbar', color="gray69") + 
   xlab("SR")+
   ylab("T0")+
   theme_bw()) # null SR vs T0 (17)
(nullF0II <- ggplot(data=F0null[F0null$T0_n<10,], aes(x=as.factor(T0_n), y=F0_n))+
    geom_boxplot(alpha=0.3, color="gray69")+
    stat_boxplot(geom ='errorbar', color="gray69") + 
    xlab("SR")+
    ylab("T0")+
    theme_bw()) # null SR vs T0 (9)

(p1nullF0_hnc <- nullF0 + # Historical Conservative
    geom_point(data = div3[div3$Period=="Historical Conservative",], 
               aes(x=T0, y=F0), size=2, alpha=1, col="#F0E442")+
    geom_line(data = hnc_f , aes(x = T0, y = fit), linetype = 1, size = 1, col="#F0E442")+
    geom_ribbon(data = hnc_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3))

(p1nullF0_hnb <- nullF0 + # Historical Broad
    geom_point(data = div3[div3$Period=="Historical Broad",], 
               aes(x=T0, y=F0), size=2, alpha=1, col="#E69F00")+
    geom_line(data = hnb_f , aes(x = T0, y = fit), linetype = 1, size = 1, col="#E69F00")+
    geom_ribbon(data = hnb_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3))

(p1nullF0_contN <- nullF0II + # Contemporary Natives
    geom_point(data = div3[div3$Period=="Contemporary Natives",], 
               aes(x=T0, y=F0), size=2, alpha=1, col="#009E73")+
    geom_line(data = contN_f, aes(x = T0, y = fit), linetype = 1, size = 1, col="#009E73")+
    geom_ribbon(data = contN_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3))

(p1nullF0_contAll <- nullF0II + # Contemporary Natives + Exotics
    geom_point(data = div3[div3$Period=="Contemporary Natives + Exotics",], 
               aes(x=T0, y=F0), size=2, alpha=1, col="#0072B2")+
    geom_line(data = contAll_f , aes(x = T0, y = fit), linetype = 1, size = 1, col="#0072B2")+
    geom_ribbon(data = contAll_f, aes(x = T0, ymin = lwr, ymax = upr), alpha = 0.3))


ggsave(p1nullF0_hnc, file= paste0(plot_dir, "/SM/Ms/Null-ObservedHNC.jpg"), width = 8, height = 6) 
ggsave(p1nullF0_hnb, file= paste0(plot_dir, "/Null-ObservedHNB.jpg"), width = 8, height = 6) 
ggsave(p1nullF0_contN, file= paste0(plot_dir, "/Null-ObservedContN.jpg"), width = 8, height = 6) 
ggsave(p1nullF0_contAll, file= paste0(plot_dir, "/Null-ObservedContAll.jpg"), width = 8, height = 6) 

###########################################################################################
# SES:-------------------------------------------------------------------------------------
null_l <- split(F0null, f=F0null$T0_n)
null_l_mean <- lapply(null_l, function(x) {c(mean(x$F0_n),
                                             mean(x$R0_n))})
null_l_sd <- lapply(null_l, function(x) {c(sd(x$F0_n),
                                           sd(x$R0_n))})
ll <- lapply(l3, function(x) {split(x, f=x$T0)})
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

period <- rep(c("Historical Conservative", "Historical Broad", "Contemporary Natives", "Contemporary Natives + Exotics"), times=c(48,60,33,49))
ses_dt <- data.frame("Period"=period, "SESF0"=ses_resultsF0, "SESR0"=ses_resultsR0)
ses_dt <- ses_dt[!is.na(ses_dt$SESF0),]

ses_dt %>% group_by(Period) %>% summarise(count = n(),
                                          countsAbv = sum(SESF0 > 0, na.rm = TRUE),
                                          countsUnd = sum(SESF0 < 0, na.rm = TRUE),
                                          countsAbvPer = (sum(SESF0 > 0, na.rm = TRUE)/count)*100,
                                          countsUndPer = (sum(SESF0 < 0, na.rm = TRUE)/count)*100)

summF0 <- summarySE(ses_dt, measurevar=c("SESF0"), groupvars=c("Period"), na.rm = TRUE)
summR0 <- summarySE(ses_dt, measurevar=c("SESR0"), groupvars=c("Period"), na.rm = TRUE)
summF0$Period <- recode_factor(summF0$Period,
                               "Historical Conservative"="HNC",
                               "Historical Broad"="HNB",
                               "Contemporary Natives"="CN",
                               "Contemporary Natives + Exotics"="CN+E")
ses_dt$Period <- recode_factor(ses_dt$Period,
                               "Historical Conservative"="HNC",
                               "Historical Broad"="HNB",
                               "Contemporary Natives"="CN",
                               "Contemporary Natives + Exotics"="CN+E")

colorBlind4   <- c("#F0E442", "#E69F00", "#009E73","#0072B2")
colorBlind3   <- c("#E69F00", "#009E73","#0072B2")

(summF0hnc <- ggplot(ses_dt[ses_dt$Period=="HNC",], 
                     aes(x=Period, y=SESF0, colour=Period)) +
    geom_point(alpha=0.4, col="#F0E442") +
    geom_violin(alpha=0.1, col="#F0E442", fill="#F0E442") +
    #geom_boxplot(width=0.1, alpha=0.5) +
    geom_point(data=summF0[summF0$Period=="HNC",], 
               aes(x=Period, y=SESF0), size=2.5, col="#F0E442") +
    geom_errorbar(data=summF0[summF0$Period=="HNC",],
                  aes(ymin=SESF0-ci, ymax=SESF0+ci),
                  width=0.2, col="#F0E442")+
    theme_bw()+
    theme(legend.position = "none")+
    labs(y="SES", x="")+
    geom_hline(yintercept = 0, linetype="dashed"))

(summF03p <- ggplot(ses_dt[!ses_dt$Period=="HNC",],
                    aes(x=Period, y=SESF0, colour=Period)) + 
    geom_point(alpha=0.5) +
    geom_violin(alpha=0.1, aes(fill=Period)) +
    geom_point(data=summF0[!summF0$Period=="HNC",], 
               aes(x=Period, y=SESF0, colour=Period), size=2.5)+
    geom_errorbar(data=summF0[!summF0$Period=="HNC",], aes(ymin=SESF0-ci, ymax=SESF0+ci),
                  width=0.2)+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_manual(values=colorBlind3, labels=c("HNB (n=60)", "CN (n=33)", "CN+E(n=49)"))+
    scale_fill_manual(values=colorBlind3, labels=c("HNB (n=60)", "CN (n=33)", "CN+E(n=49)"))+
    labs(y="SES", x="Period")+
    geom_hline(yintercept = 0, linetype="dashed")) # Same plot for redundancy.

ggsave(summF0hnc, file= paste0(plot_dir, "/SM/Ms/SES_F0hnc.jpg"), width = 6, height = 4) 
ggsave(summF03p, file= paste0(plot_dir, "/SES_F03.jpg"), width = 6, height = 4)

###########################################################################################
# One-tailed t test:-----------------------------------------------------------------------
ses_dt_n <- ses_dt %>%
  group_by(Period) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
ses_dtF0 <- ses_dt %>%
  group_by(Period) %>%
  summarise(res = list(t.test(SESF0, mu=0, alternative="less")))
ses_dtF0$res  # same p-vals when assessing redundancy

###########################################################################################
# Quantile scores:-------------------------------------------------------------------------
divrank <- div3[,-c(5)]
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
                          "HNC"="HNC",
                          "HNB"="HNB",
                          "ContN"="CN",
                          "ContAll"="CN+E")
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
(Figure3 <- ggarrange(p1nullF0_hnb,
                      p1nullF0_contN,
                      p1nullF0_contAll,
                      summF03p,
                      labels=c("A)", "B)", "C)", "D)"),
                      common.legend = T,
                      legend="bottom",
                      ncol=2,
                      nrow=2,
                      font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")))

ggsave(Figure3, file= paste0(plot_dir, "/S5_Figure3.jpg"), width =8, height = 8) 

# Supplementary material, results for HNC (yellow):---------------------------------------
(HNCres <- ggarrange(p1nullF0_hnc,
                     summF0hnc,
                     align="hv",
                     widths = c(1, 0.5),
                     labels=c("A)", "B)"),
                     nrow=1,
                     font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "bottom")))
ggsave(HNCres, file= paste0(plot_dir, "/SM/Ms/S5_HNC.jpg"), width =8, height = 5) 


##########################################################################################
# End of script ##########################################################################
