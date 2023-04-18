###########################################################################################
# Script: Create community matrices
# AFE
# April 2023
###########################################################################################

# Libraries:-------------------------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(readxl)
library(hablar)
library(data.table)

rm(list=ls())
myd <- getwd()

# Read data:-------------------------------------------------------------------------------
lists_path <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists"

load(paste0(lists_path, "/NDT53.RData")) # Raw data (53)
load(paste0(lists_path, "/NDT67.RData")) # Raw data (67)
load(paste0(lists_path, "/NDT83.RData")) # Raw data (83)
load(paste0(lists_path, "/tax.RData"))   # List of species names & taxonomic updates


###########################################################################################
# 1) First formatting of data:-------------------------------------------------------------

str(NDT53)
str(NDT67)
str(NDT83)

create_long_df <- function(variable=NULL, n=NULL){
  
  if (n==53){
    ndt <- NDT53
  }
  if (n==67){
    ndt <- NDT67
  }
  if (n==83){
    ndt <- NDT83
  }
  
  # The metadata includes locality and drainage basin info:
  metadata <- ndt[,c("SiteNameE", "DrainageBasinE")]  
  
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
varlist67 <- list(NDT67$HistoricalNatCatalog, NDT67$HistoricalNatBroad, NDT67$HistoricalExotics, NDT67$NativeFound2005, NDT67$ExoticsFound2005)
varlist83 <- list(NDT83$HistoricalNatCatalog, NDT83$HistoricalNatBroad, NDT83$HistoricalExotics, NDT83$NativeFound2005, NDT83$ExoticsFound2005)


##### SUBSET 1:
dflist53  <- lapply(varlist53, function(x) {create_long_df(variable=x, n=53)})
names(dflist53) <- names(varlist53)

dflist53  <- lapply(dflist53, function(x) {x$Species <- tax$Genus_species_ECF[match(x$Species, tax$Short)];x # Updated names
                                           x$Species <- gsub(" cf ", " ", x$Species); x          # Remove " cf " so that e.g. Poecilia cf butleri and Poecilia butleri are regarded as the same sps
                                           x$Abundance <- rep(1, nrow(x));x})                    # Create long data.frame

dflist53 <- lapply(dflist53, function(x) {x <- x %>% distinct()})
dflist53 <- lapply(dflist53, function(x) {x <- spread(x, key="Species", value="Abundance", fill=0)})


##### SUBSET 2:
dflist67  <- lapply(varlist67, function(x) {create_long_df(variable=x, n=67)})
names(dflist67) <- names(varlist67)

dflist67  <- lapply(dflist67, function(x) {x$Species <- tax$Genus_species_ECF[match(x$Species, tax$Short)];x 
                                           x$Species <- gsub(" cf ", " ", x$Species); x          
                                           x$Abundance <- rep(1, nrow(x));x})                    
dflist67 <- lapply(dflist67, function(x) {x <- x %>% distinct()})
dflist67 <- lapply(dflist67, function(x) {x <- spread(x, key="Species", value="Abundance", fill=0)})



##### SUBSET 3:
dflist83  <- lapply(varlist83, function(x) {create_long_df(variable=x, n=83)})
names(dflist83) <- names(varlist83)

dflist83  <- lapply(dflist83, function(x) {x$Species <- tax$Genus_species_ECF[match(x$Species, tax$Short)];x 
                                           x$Species <- gsub(" cf ", " ", x$Species); x          
                                           x$Abundance <- rep(1, nrow(x));x})                    
dflist83 <- lapply(dflist83, function(x) {x <- x %>% distinct()})
dflist83 <- lapply(dflist83, function(x) {x <- spread(x, key="Species", value="Abundance", fill=0)})

# For all three runnings of the function "create_long_df", warnings are fine and fixed within the same code

###########################################################################################
# 2) Correct occurrences of A monticola and remove cols that indicate presence of no fish:-

# ALL SUBSETS (Historical period):

l_historical <- list(dflist53[[1]], dflist53[[2]], 
                 dflist67[[1]], dflist67[[2]], 
                 dflist83[[1]], dflist83[[2]])


l_contemporary <- list(dflist53[[4]], dflist53[[5]], 
                     dflist67[[4]], dflist67[[5]], 
                     dflist83[[4]], dflist83[[5]])


# A monticola is Algansea monticola in "El Huizcolote River" & "El Sacristan Spring"
# A monticola is Agonostomus monticola in "De Comala River", "Salado River"


for (i in 1:length(l_historical)) {
  dt <- l_historical[[i]]
  
  dt$`Dajaus monticola` <- rep(0, nrow(dt))
  dt$`Algansea monticola` <- rep(0, nrow(dt))
  
  dt$`Algansea monticola` <- ifelse(dt$SiteName %in% c("El Huizcolote River", "El Sacristan Spring"), 1, 0)
  dt$`Dajaus monticola` <- ifelse(dt$SiteName %in% c("De Comala River", "Salado River"), 1, 0)
  dt <- dt[,!names(dt)=="A monticola"]
  l_historical[[i]] <- dt
}
for (i in 1:length(l_contemporary)) {
  dt <- l_contemporary[[i]]
  
  dt$`Dajaus monticola` <- rep(0, nrow(dt))
  dt$`Algansea monticola` <- rep(0, nrow(dt))
  
  dt$`Algansea monticola` <- ifelse(dt$SiteName %in% c("El Huizcolote River", "El Sacristan Spring"), 1, 0)
  dt$`Dajaus monticola` <- ifelse(dt$SiteName %in% c("De Comala River", "Salado River"), 1, 0)
  dt <- dt[,!names(dt)=="A monticola"]
  l_contemporary[[i]] <- dt
}

dflist53[[1]] <- l_historical[[1]]
dflist53[[2]] <- l_historical[[2]]
dflist67[[1]] <- l_historical[[3]]
dflist67[[2]] <- l_historical[[4]]
dflist83[[1]] <- l_historical[[5]]
dflist83[[2]] <- l_historical[[6]]

dflist53[[4]] <- l_contemporary[[1]]
dflist53[[5]] <- l_contemporary[[2]]
dflist67[[4]] <- l_contemporary[[3]]
dflist67[[5]] <- l_contemporary[[4]]
dflist83[[4]] <- l_contemporary[[5]]
dflist83[[5]] <- l_contemporary[[6]]


lapply(dflist53, function(x) {unique(names(x))})
lapply(dflist67, function(x) {unique(names(x))})
lapply(dflist83, function(x) {unique(names(x))})

# Remove columns for non-species:
non_species <- c("none", "no aquatic life","Dry", "foul water", 
                 "ceased to exist", "it was not the right locality", 
                 "is the locality", "-", "(?)", "?")
dflist53 <- lapply(dflist53, function(x) {x <- x[,!names(x) %in% non_species]})
dflist67 <- lapply(dflist67, function(x) {x <- x[,!names(x) %in% non_species]})
dflist83 <- lapply(dflist83, function(x) {x <- x[,!names(x) %in% non_species]})

#"Zoogoneticus tequila y Ameca  splendens":
dflist83[[4]]$SiteNameE[dflist83[[4]]$`Zoogoneticus tequila y Ameca  splendens`==1]
"Ameca splendens" %in% names(dflist83[[4]])      # FALSE
"Zoogoneticus tequila" %in% names(dflist83[[4]]) # FALSE

dflist83[[4]]$`Ameca splendens`<- rep(NA, nrow(dflist83[[4]]))
dflist83[[4]]$`Zoogoneticus tequila`<- rep(NA, nrow(dflist83[[4]]))

dflist83[[4]]$`Ameca splendens`[dflist83[[4]]$SiteNameE=="Manantial/balneario El Rincon"] <- 1
dflist83[[4]]$`Zoogoneticus tequila`[dflist83[[4]]$SiteNameE=="Manantial/balneario El Rincon"] <- 1

dflist83[[4]] <- dflist83[[4]][,!names(dflist83[[4]]) %in% "Zoogoneticus tequila y Ameca  splendens"]


# Check column names again:
lapply(dflist53, function(x) {unique(names(x))})
lapply(dflist67, function(x) {unique(names(x))})
lapply(dflist83, function(x) {unique(names(x))})


###########################################################################################
# 3) Generate lists:-----------------------------------------------------------------------
path_list_SM <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists/SM"

# Localities:------------------------------------------------------------------------------
sites <- data.frame(matrix(ncol=3, nrow=83))
names(sites) <- c("Site", "Subset", "Comment")
sites$Site <- sort(unique(NDT83$SiteName))
sites$Subset <- "83"
sites$Subset <- ifelse(sites$Site %in% NDT67$SiteName, "67", sites$Subset)
sites$Subset <- ifelse(sites$Site %in% NDT53$SiteName, "53", sites$Subset)
sites$Comment <- ifelse(sites$Subset=="67", "The historical assemblage could be reconstructed. The contemporary locality didn't sustain life anymore", NA)
sites$Comment <- ifelse(sites$Subset=="83", "Assemblage data for the two periods is not comparable", sites$Comment)
sites$Comment <- ifelse(sites$Site %in% c("rio Marabasco", "rio Coscomate", "laguna El Barril", "rio en 6 de Enero",
                                          "rio Cuarenta", "Los Vergeles", "presa sobre el rio Carrizal"), "NOT IN SM PAPER 2018. Assemblage data for the two periods is not comparable", sites$Comment)

sum(sites$Comment=="NOT IN SM PAPER 2018. Assemblage data for the two periods is not comparable", na.rm=T) # 7, OK
write.csv2(sites, file=paste0(path_list_SM, "/sites.csv"))


# Anecdotics:------------------------------------------------------------------------------

anec <- NDT83$Anectodics
anec <- anec[!is.na(anec)] # % observations
anec

anectax <- c("Sicydium multipunctatum", "Oncorhynchus mykiss", 
             "Rhamdia sp", "Cyprinus carpio", "trompos", "peces limpiadores",
             "Bagre sp", "Micropterus salmoides")

write.csv2(anectax, file=paste0(path_list_SM, "/anectax.csv"))

# Species:
vec_taxa <- sort(unique(c(names(dflist83[[1]]), names(dflist83[[2]]),
                          names(dflist83[[3]]), names(dflist83[[4]]),
                          names(dflist83[[5]]))))
vec_taxa <- vec_taxa[!vec_taxa %in% c("SiteNameE", "DrainageBasinE")]

taxa <- data.frame(matrix(ncol=5, nrow=113))
names(taxa) <- c("Genus_species", "Previous", "Updated", "Regional_Status", "Comment")

taxa$Genus_species <- vec_taxa
taxa$Previous <- tax$Genus_species[match(taxa$Genus_species, tax$Genus_species_ECF)]
# remember that here, for each "taxa$Genus_species" there might be multiple "tax$Genus_species"(i.e. previous)
sum(is.na(taxa$Previous))
taxa$Previous[taxa$Genus_species=="Dajaus monticola"] <- "Agonostomus monticola"
taxa$Previous[taxa$Genus_species=="Algansea monticola"] <- "Algansea monticola"
taxa$Previous[taxa$Genus_species=="Gambusia senilis"] <- "Gambusia senilis"

taxa$Updated <- ifelse(taxa$Previous==taxa$Genus_species, "no", "yes")

# Regional status:
# Introduced:
introduced <- c("Oreochromis aureus", "Oreochromis mossambicus", "Oreochromis niloticus", "Oreochromis sp",
                "Xiphophorus hellerii", "Cyprinus carpio", "Poecilia reticulata", "Pseudoxiphophorus bimaculatus",
                "Pseudoxiphophorus jonesii", "Pseudoxiphophorus sp", "Lepomis macrochirus", "Micropterus salmoides",
                "Xiphophorus maculatus", "Gambusia yucatana", "Poeciliopsis gracilis", "Pomoxis nigromaculatus",
                "Amatitlania nigrofasciata", "Astatotilapia burtoni", "Carassius auratus", "Gambusia affinis", "Xiphophorus variatus")

# Ctenopharyngodon idella (in paper but not in database)

taxa$Regional_Status <- ifelse(taxa$Genus_species %in% introduced, "Introduced", "Native") # OK

# Comments:


# TO FINISH!



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