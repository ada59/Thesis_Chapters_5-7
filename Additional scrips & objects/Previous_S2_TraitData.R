
################################################
# previous trait fill code (measurements august):
################################################


setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/Traits/Measurements")
file_names <- dir() # 26
mes <- do.call(rbind,lapply(file_names,read.csv))   # 312 (12*26)                                # All data
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2")   # Go back to main dir

# Retrieve mesurements:
file_names <- str_split_fixed(file_names, "_", 2)[,2]
file_names <- str_split_fixed(file_names, "\\.", 2)[,1]
file_names <- gsub("_", " ", file_names)

mes$S <- rep(file_names, each=12)
mes$M <- rep(c("Bl", "Bd", "Hd", "Ed", "CPd", "CFd", "Eh", "Mo", "Jl", "PFh", "PFi", "PFiII"), times=26)

sum(mes$Length==0)  # 6
sum(mes$Length <0)  # 0
mes$Length[mes$Length==0] <- NA # Mesurements that could not be taken
sum(is.na(mes$Length))

mes <- mes[c("S", "M", "Length")]
names(mes) <- c("Species", "Measure", "Length")

mes_mat <- spread(mes, key="Measure", value="Length")

# Corrected eye measurements (5 sps):
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/Traits/Measurements_Corrected")
file_namesC <- dir() # 5
mesC <- do.call(rbind,lapply(file_namesC,read.csv))
setwd("C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2") 

file_namesC <- gsub("_", " ", str_split_fixed(file_namesC, "_", 2)[,2])
file_namesC <- str_split_fixed(file_namesC, "\\.", 2)[,1]

mesC$S <- rep(file_namesC, each=2)
mesC$M <- rep(c("Hd", "Ed"), times=5)

# Atherinella balsana:
mes_mat$Hd[mes_mat$Species=="Atherinella balsana"] <- mesC$Length[mesC$S=="Atherinella balsana" & mesC$M=="Hd"]
mes_mat$Ed[mes_mat$Species=="Atherinella balsana"] <- mesC$Length[mesC$S=="Atherinella balsana" & mesC$M=="Ed"]

# Chirostoma jordani:
mes_mat$Hd[mes_mat$Species=="Chirostoma jordani"] <- mesC$Length[mesC$S=="Chirostoma jordani" & mesC$M=="Hd"]
mes_mat$Ed[mes_mat$Species=="Chirostoma jordani"] <- mesC$Length[mesC$S=="Chirostoma jordani" & mesC$M=="Ed"]

# Chirostoma labarcae:
mes_mat$Hd[mes_mat$Species=="Chirostoma labarcae"] <- mesC$Length[mesC$S=="Chirostoma labarcae" & mesC$M=="Hd"]
mes_mat$Ed[mes_mat$Species=="Chirostoma labarcae"] <- mesC$Length[mesC$S=="Chirostoma labarcae" & mesC$M=="Ed"]

# Pseudoxiphophorus jonesii:
mes_mat$Hd[mes_mat$Species=="Pseudoxiphophorus jonesii"] <- mesC$Length[mesC$S=="Pseudoxiphophorus jonesii" & mesC$M=="Hd"]
mes_mat$Ed[mes_mat$Species=="Pseudoxiphophorus jonesii"] <- mesC$Length[mesC$S=="Pseudoxiphophorus jonesii" & mesC$M=="Ed"]

# Xiphophorus variatus:
mes_mat$Hd[mes_mat$Species=="Xiphophorus variatus"] <- mesC$Length[mesC$S=="Xiphophorus variatus" & mesC$M=="Hd"]
mes_mat$Ed[mes_mat$Species=="Xiphophorus variatus"] <- mesC$Length[mesC$S=="Xiphophorus variatus" & mesC$M=="Ed"]

# Compute ratios:
fishtrait_fill <- subset(fishtrait, fishtrait$Genus.species %in% file_names) #23, OK
names(fishtrait_fill)

gila_mes_mat <- subset(mes_mat, mes_mat$Species %in% c("Gila conspersa", "Gila minaceae", "Gila pulchra"))
mes_mat <- mes_mat[! mes_mat$Species %in% c("Gila conspersa", "Gila minaceae", "Gila pulchra"),]

fishtrait_fill$BEl <- mes_mat$Bl/mes_mat$Bd #Body elongation
fishtrait_fill$VEp <- mes_mat$Eh/mes_mat$Bd #Verical eye position
fishtrait_fill$REs <- mes_mat$Ed/mes_mat$Hd #Relative eye size
fishtrait_fill$OGp <- mes_mat$Mo/mes_mat$Bd #Oral gape position
fishtrait_fill$RMl <- mes_mat$Jl/mes_mat$Hd #Relative maxillary length
fishtrait_fill$BLs <- mes_mat$Hd/mes_mat$Bd #Body lateral shape
fishtrait_fill$PFv <- mes_mat$PFh/mes_mat$Bd #Pectoral fin vertical position
fishtrait_fill$PFs <- mes_mat$PFiII/mes_mat$Bl #Pectoral fin size
fishtrait_fill$CPt <- mes_mat$CFd/mes_mat$CPd #Caudal peduncle throttling


fishtrait_fill$Reference <- rep("Miller, R.R. (2009) Freshwater Fishes of Mexico (Peces Dulceacuicolas de Mexico). 1st edn.", nrow(fishtrait_fill))
fishtrait_fill$Type.of.illustration <- rep("Picture", nrow(fishtrait_fill))
fishtrait_fill$Type.of.illustration <- ifelse(fishtrait_fill$Genus.species %in% c("Allotoca zacapuensis",
                                                                                  