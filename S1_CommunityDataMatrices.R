
###########################################################################################
# Script: Create community matrices
# AFE
# August 2022
###########################################################################################

# Libreries:-------------------------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(readxl)
library(hablar)
library(data.table)

rm(list=ls())
myd <- getwd()

# Read data:-------------------------------------------------------------------------------
load(paste0(myd, "/NDT53.RData")) # Raw data (53)
load(paste0(myd, "/NDT67.RData")) # Raw data (67)
load(paste0(myd, "/tax.RData"))   # List of species names & taxonomic updates


###########################################################################################
# 1) Create Historical & Contemporary Comm Matrices:---------------------------------------
dim(NDT53)
str(NDT53)

tax$Final <- ifelse(is.na(tax$Update1), tax$Genus_species, tax$Update1)

create_long_df53 <- function(variable=NULL){
  
  # The metadata includes locality and drainage basin info:
  metadata <- NDT53[,c("SiteNameE", "DrainageBasinE")]  
  
  # Elements listed in the variable column get split:
  varI <- data.frame(do.call(rbind, strsplit(variable, ","))) # Will display a warning but this is fixed in the subsequent line
  varI <- as.data.frame(t(apply(varI, 1, function(x) {replace(x, duplicated(x), NA)})))
  
  # Returns data frame in long format
  df <- cbind(metadata, varI)
  df <- gather(df, key="occurrence", value="Species", -c(1:2), na.rm = T)
  df <- within(df, rm(occurrence))
  
  return(df)
} 

varlist53 <- list(NDT53$HistoricalNatCatalog, NDT53$HistoricalNatBroad, NDT53$HistoricalExotics, NDT53$NativeFound2005, NDT53$ExoticsFound2005)
dflist53  <- lapply(varlist53, create_long_df53)
dflist53  <- lapply(dflist53, function(x) {x$Species <- tax$Final[match(x$Species, tax$Short)];x # Updated names
                                           x$Species <- gsub(" cf ", " ", x$Species); x          # Remove " cf " so that e.g. Poecilia cf butleri and Poecilia butleri are regarded as the same sps
                                           x$Abundance <- rep(1, nrow(x));x})                    # Create long data.frame

dflist53 <- lapply(dflist53, function(x) {x <- x %>% distinct()})
dflist53 <- lapply(dflist53, function(x) {x <- spread(x, key="Species", value="Abundance", fill=0)})


# Correct occurrences of A monticola and remove cols that indicate presence of no fish:----
HNC <- dflist53[[1]]    # Historical Catalog Natives
HNB <- dflist53[[2]]    # Historical Broad Natives
HE <- dflist53[[3]]     # Historical Exotics
ContN <- dflist53[[4]]  # Contemporary native
ContE <- dflist53[[5]]  # Contemporary exotics


# Historical Catalog (Natives only):
HNC$SiteName[HNC$`A monticola`==1]
# "De Comala River"/"El Huizcolote River"

names(HNC)[names(HNC)=="A monticola"] <- "Agonostomus monticola"
HNC$`Agonostomus monticola`[HNC$SiteName=="El Huizcolote River"] <- 0

HNC$`Algansea monticola` <- rep(0, nrow(HNC))
HNC$`Algansea monticola`[HNC$SiteName=="El Huizcolote River"] <- 1

# Historical Broad (Natives only):
HNB$SiteName[HNB$`A monticola`==1]
# "De Comala River"/"El Huizcolote River"/"Salado River"       

names(HNB)[names(HNB)=="A monticola"] <- "Agonostomus monticola"
HNB$`Agonostomus monticola`[HNB$SiteName=="El Huizcolote River"] <- 0

HNB$`Algansea monticola` <- rep(0, nrow(HNB))
HNB$`Algansea monticola`[HNB$SiteName=="El Huizcolote River"] <- 1

# Contemporary native:
ContN$SiteName[ContN$`A monticola`==1]
# "De Comala River"/"Salado River"
names(ContN)[names(ContN)=="A monticola"] <- "Agonostomus monticola"
names(ContN)
ContN <- within(ContN, rm("none"))

# Contemporary exotic:
"A monticola" %in% names(ContE) # FALSE
names(ContE)

# Save community presence/absence matrices:
save(HNC, file=paste0(myd, "/HNC.RData"))
save(HNB, file=paste0(myd, "/HNB.RData"))
save(HE, file=paste0(myd, "/HE.RData"))
save(ContN, file=paste0(myd, "/ContN.RData"))
save(ContE, file=paste0(myd, "/ContE.RData"))


# 2) Taxonomic updates in species names: -------------------------------------------------
# These updates are just so that the package fishtree picks up the names.
listData <- list(HNC, HNB, HE, ContN, ContE)    
taxonomicUpdate2 <- function (df=NULL) {
  names(df) <- plyr::revalue(names(df), c("Nosferatu bartoni"="Herichthys bartoni",
                                          "Nosferatu labridens"="Herichthys labridens",
                                          "Mayaheros beani"="Cichlasoma beani",
                                          "Tampichthys dichroma"="Tampichthys dichromus"))
  return(df)
} 
listData <- lapply(listData, taxonomicUpdate2)  # Warnings are fine: sps in the list not present in each comm data

# Check that records of C. jordani include both (C. jordani & C. mezquital):
NDT53$SiteNameE[NDT53$HistoricalNatBroad %like% "mezquital"] #2
NDT53$SiteNameE[NDT53$HistoricalNatBroad %like% "jordani"]   #8

sum(HNB$`Chirostoma jordani`==1) # 10 = 8 + 2, OK

###########################################################################################
# End of script ###########################################################################