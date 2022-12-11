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
dim(NDT67)
str(NDT67)

tax$Final <- ifelse(is.na(tax$Update1), tax$Genus_species, tax$Update1)

create_long_df67 <- function(variable=NULL){
  
  # The metadata includes locality and drainage basin info:
  metadata <- NDT67[,c("SiteNameE", "DrainageBasinE")]  
  
  # Elements listed in the variable column get split:
  varI <- data.frame(do.call(rbind, strsplit(variable, ","))) # Will display a warning but this is fixed in the subsequent line
  varI <- as.data.frame(t(apply(varI, 1, function(x) {replace(x, duplicated(x), NA)})))
  
  # Returns data frame in long format
  df <- cbind(metadata, varI)
  df <- gather(df, key="occurrence", value="Species", -c(1:2), na.rm = T)
  df <- within(df, rm(occurrence))
  
  return(df)
} 

varlist67 <- list(NDT67$HistoricalNatCatalog, NDT67$HistoricalNatBroad, NDT67$HistoricalExotics, NDT67$NativeFound2005, NDT67$ExoticsFound2005)
dflist67  <- lapply(varlist67, create_long_df67)
dflist67  <- lapply(dflist67, function(x) {x$Species <- tax$Final[match(x$Species, tax$Short)];x # Updated names
                                           x$Species <- gsub(" cf ", " ", x$Species); x          # Remove " cf " so that e.g. Poecilia cf butleri and Poecilia butleri are regarded as the same sps
                                           x$Abundance <- rep(1, nrow(x));x})                    # Create long data.frame

dflist67 <- lapply(dflist67, function(x) {x <- x %>% distinct()})
dflist67 <- lapply(dflist67, function(x) {x <- spread(x, key="Species", value="Abundance", fill=0)})


# Correct occurrences of A monticola and remove cols that indicate presence of no fish:----
HNC <- dflist67[[1]]    # Historical Catalog Natives
HNB <- dflist67[[2]]    # Historical Broad Natives
HE <- dflist67[[3]]     # Historical Exotics
ContN <- dflist67[[4]]  # Contemporary native
ContE <- dflist67[[5]]  # Contemporary exotics


# Historical Catalog (Natives only):
HNC$SiteName[HNC$`A monticola`==1]
# "De Comala River"/"El Huizcolote River"/"El Sacristan Spring"

names(HNC)[names(HNC)=="A monticola"] <- "Agonostomus monticola"
HNC$`Agonostomus monticola`[HNC$SiteName=="El Huizcolote River"] <- 0
HNC$`Agonostomus monticola`[HNC$SiteName=="El Sacristan Spring"] <- 0

HNC$`Algansea monticola` <- rep(0, nrow(HNC))
HNC$`Algansea monticola`[HNC$SiteName=="El Huizcolote River"] <- 1
HNC$`Algansea monticola`[HNC$SiteName=="El Sacristan Spring"] <- 1

# Historical Broad (Natives only):
HNB$SiteName[HNB$`A monticola`==1]
# "De Comala River"/"El Huizcolote River"/"El Sacristan Spring"/"Salado River"       

names(HNB)[names(HNB)=="A monticola"] <- "Agonostomus monticola"
HNB$`Agonostomus monticola`[HNB$SiteName=="El Huizcolote River"] <- 0
HNB$`Agonostomus monticola`[HNB$SiteName=="El Sacristan Spring"] <- 0

HNB$`Algansea monticola` <- rep(0, nrow(HNB))
HNB$`Algansea monticola`[HNB$SiteName=="El Huizcolote River"] <- 1
HNB$`Algansea monticola`[HNB$SiteName=="El Sacristan Spring"] <- 1

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
NDT67$SiteNameE[NDT67$HistoricalNatBroad %like% "mezquital"] #2
NDT67$SiteNameE[NDT67$HistoricalNatBroad %like% "jordani"]   #13

sum(HNB$`Chirostoma jordani`==1) # 15= 13 + 2, OK

# 3) Data formatting for assemblage-level analysis:--------------------------------------
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

ls <- split(l, rownames(l)) # 67*4 = 268
save(ls, file="ls.RData")

###########################################################################################
# End of script ###########################################################################