###########################################################################################
# Script: Sensitivity Analysis/ Observer test
# AFE
# Nov 2022
###########################################################################################

# Libraries:------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpubr)

rm(list=ls())

# Read Data:------------------------------------------------------------------------------
FMorph <- read.csv2("C:/Users/Usuario/Documents/PHD/Datasets/FishMORPH/FISHMORPH_Database.csv", h=T)

setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/Traits/TestObserver/Results")
file_names <- dir() 
mes <- do.call(rbind,lapply(file_names,read.csv))                                  

# Corrected mes:-------------------------------------------------------------------------
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/Traits/TestObserver/Corrected")
file_names_c <- dir() 
mes_c <- do.call(rbind,lapply(file_names_c,read.csv))                                  
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2") 

# Retrieve data: -------------------------------------------------------------------------
FMorph_subset <- subset(FMorph, FMorph$Genus.species %in% c("Astyanax biotae",
                                                            "Astyanax bockmanni",
                                                            "Astyanax lacustris",
                                                            "Astyanax schubarti",
                                                            "Crenicichla britskii",
                                                            "Crenicichla haroldoi",
                                                            "Crenicichla jaguarensis",
                                                            "Crenicichla jupiaensis"))
FMorph_subset <- data.frame("Species"=FMorph_subset$Genus.species,
                            "Observer"=rep("FM", nrow(FMorph_subset)),
                            FMorph_subset[,c(7:15)])

# Retrieve mesurements:
file_names <- str_split_fixed(file_names, "_", 2)[,2]
file_names <- str_split_fixed(file_names, "\\.", 2)[,1]
file_names <- gsub("_", " ", file_names)

mes$S <- rep(file_names, each=11)
mes$M <- rep(c("Bl", "Bd", "Hd", "Ed", "CPd", "CFd", "Eh", "Mo", "Jl", "PFh", "PFl"), times=8)

sum(mes$Length==0)  # 1
sum(mes$Length <0)  # 0
mes$Length[mes$Length==0] <- NA # Mesurements that could not be taken
sum(is.na(mes$Length))

mes <- mes[c("S", "M", "Length")]
names(mes) <- c("Species", "Measure", "Length")

mes_mat <- spread(mes, key="Measure", value="Length")

# Correct some re-taken measurements after data visualisation:-----------------------------
mes_mat$Mo[mes_mat$Species=="A lacustris"] <- mes_c$Length[1]
mes_mat$Jl[mes_mat$Species=="A lacustris"] <- mes_c$Length[2]

mes_mat$Mo[mes_mat$Species=="C jupiaensis"] <- mes_c$Length[3]
mes_mat$Jl[mes_mat$Species=="C jupiaensis"] <- mes_c$Length[4]

# Format: ---------------------------------------------------------------------------------
mes_mat2 <- data.frame(matrix(ncol=11, nrow=8))
names(mes_mat2) <- c("Species", "Observer","BEl", "VEp", "REs", "OGp", "RMl", "BLs", 
                     "PFv", "PFs", "CPt")

mes_mat2$Species <- mes_mat$Species
mes_mat2$Species <- gsub("A ", "Astyanax ", mes_mat2$Species)
mes_mat2$Species <- gsub("C ", "Crenicichla ", mes_mat2$Species)
mes_mat2$Observer <- rep("A", nrow=mes_mat2)

mes_mat2$BEl <- mes_mat$Bl/mes_mat$Bd   #Body elongation
mes_mat2$VEp <- mes_mat$Eh/mes_mat$Bd   #Verical eye position
mes_mat2$REs <- mes_mat$Ed/mes_mat$Hd   #Relative eye size
mes_mat2$OGp <- mes_mat$Mo/mes_mat$Bd   #Oral gape position
mes_mat2$RMl <- mes_mat$Jl/mes_mat$Hd   #Relative maxillary length
mes_mat2$BLs <- mes_mat$Hd/mes_mat$Bd   #Body lateral shape
mes_mat2$PFv <- mes_mat$PFh/mes_mat$Bd  #Pectoral fin vertical position
mes_mat2$PFs <- mes_mat$PFl/mes_mat$Bl  #Pectoral fin size
mes_mat2$CPt <- mes_mat$CFd/mes_mat$CPd #Caudal peduncle throttling

identical(names(mes_mat2), names(FMorph_subset)) # TRUE

dt <- rbind(mes_mat2, FMorph_subset)
dt <- dt %>% hablar::retype()
dt$Observer <- as.factor(dt$Observer)
#dt <- gather(dt, key="Trait", value="Measurement", -c(1:2))

# Plot data: -----------------------------------------------------------------------------
dt_factors <- dt[, c(1:2)]
dt_measures <- dt[, c(3:11)]

dfl <- list()
for(i in 1:ncol(dt_measures)) { 
  mes <- dt_measures[ , i]
  dfl[[i]] <- as.data.frame(cbind(dt_factors, mes))
}
dfl <- lapply(dfl, function(x) {names(x) <- c("Species", "Observer", "Mes");x})
names(dfl) <- names(dt[,c(3:11)])

list_bp <- list()
list_p <- list()
for (i in 1:length(dfl)){
  dataI <- dfl[[i]]
  x <- names(dfl)[[i]]
  bp <- ggplot(dataI, aes(x=Observer, y=Mes)) + geom_violin(fill="lightgray", alpha=0.3) + geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+ theme_bw() + labs(title=paste("Plot of", x), x="Observer", y = "") + theme(plot.title=element_text(size = 10))
  p <- ggplot(dataI, aes(x=Observer, y=Mes)) + geom_point(aes(col=Species)) + theme_bw() + labs(title=paste("Plot of", x), x="Observer", y = "") + theme(plot.title=element_text(size = 10))
  list_bp[[i]] <- bp
  list_p[[i]] <- p
}

boxplotsApp4 <- grid.arrange(grobs = list_bp)
pointsApp4 <- ggarrange(list_p[[1]], list_p[[2]], list_p[[3]], list_p[[4]],
                        list_p[[5]], list_p[[6]], list_p[[7]], list_p[[8]], list_p[[9]],
                        common.legend = TRUE)

plot_dir <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots/SM/Ms" # Supplementary material plots


#ggsave(boxplotsApp4, file= paste0(plot_dir, "/boxplotsApp4.jpg"), width = 8, height = 6) 
#ggsave(pointsApp4, file= paste0(plot_dir, "/pointsApp4.jpg"), width = 8, height = 6) 

###########################################################################################
# Test slightly different measurements of Hd and Eh:---------------------------------------
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/Traits/TestObserver/Hd&Eh")
file_names_hdeh <- dir() 
mes_hdeh <- do.call(rbind,lapply(file_names_hdeh,read.csv))                                  
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2") 

# Correct some re-taken measurements after data visualisation:-----------------------------
mes_mat$Hd[mes_mat$Species=="A biotae"] <- mes_hdeh$Length[1]
mes_mat$Hd[mes_mat$Species=="A bockmanni"] <- mes_hdeh$Length[1]
mes_mat$Hd[mes_mat$Species=="A lacustris"] <- mes_hdeh$Length[1]
mes_mat$Hd[mes_mat$Species=="A schubarti"] <- mes_hdeh$Length[1]

mes_mat$Hd[mes_mat$Species=="C britskii"] <- mes_hdeh$Length[1]
mes_mat$Hd[mes_mat$Species=="C haroldoi"] <- mes_hdeh$Length[1]
mes_mat$Hd[mes_mat$Species=="C jaguarensis"] <- mes_hdeh$Length[1]
mes_mat$Hd[mes_mat$Species=="C jupiaensis"] <- mes_hdeh$Length[1]


mes_mat$Hd[mes_mat$Species=="A biotae"] <- mes_hdeh$Length[2]
mes_mat$Hd[mes_mat$Species=="A bockmanni"] <- mes_hdeh$Length[2]
mes_mat$Hd[mes_mat$Species=="A lacustris"] <- mes_hdeh$Length[2]
mes_mat$Hd[mes_mat$Species=="A schubarti"] <- mes_hdeh$Length[2]

mes_mat$Hd[mes_mat$Species=="C britskii"] <- mes_hdeh$Length[2]
mes_mat$Hd[mes_mat$Species=="C haroldoi"] <- mes_hdeh$Length[2]
mes_mat$Hd[mes_mat$Species=="C jaguarensis"] <- mes_hdeh$Length[2]
mes_mat$Hd[mes_mat$Species=="C jupiaensis"] <- mes_hdeh$Length[2]

mes_mat2B <- data.frame(matrix(ncol=11, nrow=8))
names(mes_mat2B) <- c("Species", "Observer","BEl", "VEp", "REs", "OGp", "RMl", "BLs", 
                     "PFv", "PFs", "CPt")

mes_mat2B$Species <- mes_mat2B$Species
mes_mat2B$Species <- gsub("A ", "Astyanax ", mes_mat2B$Species)
mes_mat2B$Species <- gsub("C ", "Crenicichla ", mes_mat2B$Species)
mes_mat2B$Observer <- rep("A", nrow=mes_mat2B)

mes_mat2B$BEl <- mes_mat$Bl/mes_mat$Bd   #Body elongation
mes_mat2B$VEp <- mes_mat$Eh/mes_mat$Bd   #Verical eye position
mes_mat2B$REs <- mes_mat$Ed/mes_mat$Hd   #Relative eye size
mes_mat2B$OGp <- mes_mat$Mo/mes_mat$Bd   #Oral gape position
mes_mat2B$RMl <- mes_mat$Jl/mes_mat$Hd   #Relative maxillary length
mes_mat2B$BLs <- mes_mat$Hd/mes_mat$Bd   #Body lateral shape
mes_mat2B$PFv <- mes_mat$PFh/mes_mat$Bd  #Pectoral fin vertical position
mes_mat2B$PFs <- mes_mat$PFl/mes_mat$Bl  #Pectoral fin size
mes_mat2B$CPt <- mes_mat$CFd/mes_mat$CPd #Caudal peduncle throttling

dtB <- rbind(mes_mat2B, FMorph_subset)
dtB <- dtB %>% hablar::retype()
dtB$Observer <- as.factor(dtB$Observer)
#dt <- gather(dt, key="Trait", value="Measurement", -c(1:2))

# Plot data: -----------------------------------------------------------------------------
dt_factorsB <- dtB[, c(1:2)]
dt_measuresB <- dtB[, c(3:11)]

dflB <- list()
for(i in 1:ncol(dt_measuresB)) { 
  mes <- dt_measuresB[ , i]
  dflB[[i]] <- as.data.frame(cbind(dt_factorsB, mes))
}
dflB <- lapply(dflB, function(x) {names(x) <- c("Species", "Observer", "Mes");x})
names(dflB) <- names(dtB[,c(3:11)])

list_bp <- list()
list_p <- list()
for (i in 1:length(dflB)){
  dataI <- dflB[[i]]
  x <- names(dflB)[[i]]
  bp <- ggplot(dataI, aes(x=Observer, y=Mes)) + geom_violin(fill="lightgray", alpha=0.3) + geom_boxplot(alpha=0.5, width=0.1, position=position_dodge(1))+ theme_bw() + labs(title=paste("Plot of", x), x="Observer", y = "") + theme(plot.title=element_text(size = 10))
  p <- ggplot(dataI, aes(x=Observer, y=Mes)) + geom_point(aes(col=Species)) + theme_bw() + labs(title=paste("Plot of", x), x="Observer", y = "") + theme(plot.title=element_text(size = 10))
  list_bp[[i]] <- bp
  list_p[[i]] <- p
}
boxplotsApp4B <- grid.arrange(grobs = list_bp)

# Head length measure I correct
# Thus Eh makes more sense centered at the pupil as well.

# End of script###########################################################################
##########################################################################################


