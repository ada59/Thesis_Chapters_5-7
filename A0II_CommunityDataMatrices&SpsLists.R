#===============================================================================
# Script: Create SM Lists
# AFE
# April 2023 & August 2024
#===============================================================================



#===============================================================================
# Libraries:--------------------------------------------------------------------
#===============================================================================
library(dplyr)
library(tidyverse)
library(readxl)
library(hablar)
library(data.table)

rm(list=ls())
myd <- getwd()



#===============================================================================
# Read data:--------------------------------------------------------------------
#===============================================================================
lists_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists"

load(paste0(lists_path, "/NDT.RData"))   # Raw 
load(paste0(lists_path, "/NDT53.RData")) # Raw data (53 sites)
load(paste0(lists_path, "/NDT67.RData")) # Raw data (67 sites, 53 + 14 with no fish life in 2005)
load(paste0(lists_path, "/NDT83.RData")) # Raw data (83 sites, sampled but including 16 sites with incomplete reconstructions)
load(paste0(lists_path, "/tax.RData"))   # List of species names & taxonomic updates




#===============================================================================
# Deal with unidentified records (addendum): -----------------------------------
#===============================================================================
# Given that the traits for these ocurrences are approoximated by those of other
# species. keeping them in the trait matrix as separate obs would cause problems 
# as indicated by the functions for assessing the quality of trait space in 
# the script Species-level_I.
# therefore, for the purpose of analyses I treat them as the species they've 
# been approximated by.


## Replace Oreochromis sp: -----------------------------------------------------
NDT67$ExoticsFound2005[NDT67$SiteName=="rio de La Pola"] # for check
NDT67$ExoticsFound2005 <- gsub("Oreochromis sp.", "O. aureus", NDT67$ExoticsFound2005)
NDT67$ExoticsFound2005[NDT67$SiteName=="rio de La Pola"] # for check, OK

NDT53$ExoticsFound2005 <- gsub("Oreochromis sp.", "O. aureus", NDT53$ExoticsFound2005)


## Replace Poeciliopsis sp: ----------------------------------------------------
NDT67$ExoticsFound2005[NDT67$SiteName=="canales de Xochimilco"] # for check
NDT67$ExoticsFound2005 <- gsub("Poeciliopsis sp.", "P. gracilis", NDT67$ExoticsFound2005)
NDT67$ExoticsFound2005[NDT67$SiteName=="canales de Xochimilco"] # for check, OK

NDT53$ExoticsFound2005 <- gsub("Poeciliopsis sp.", "P. gracilis", NDT53$ExoticsFound2005)


## Replace Poecilia sp: --------------------------------------------------------
NDT67$ExoticsFound2005[NDT67$SiteName=="rio Turbio cerca de El Nopal"] # for check
NDT67$ExoticsFound2005 <- gsub("Poecilia sp.", "P. sphenops", NDT67$ExoticsFound2005)
NDT67$ExoticsFound2005[NDT67$SiteName=="rio Turbio cerca de El Nopal"] # for check, OK

NDT53$ExoticsFound2005 <- gsub("Poecilia sp.", "P. sphenops", NDT53$ExoticsFound2005)

## Replace Astyanax sp: --------------------------------------------------------
NDT67$HistoricalNatBroad <- gsub("Astyanax sp.", "A. mexicanus", NDT67$HistoricalNatBroad)
NDT67$NativeFound2005 <- gsub("Astyanax sp.", "A. mexicanus", NDT67$NativeFound2005)

NDT53$HistoricalNatBroad <- gsub("Astyanax sp.", "A. mexicanus", NDT53$HistoricalNatBroad)
NDT53$NativeFound2005 <- gsub("Astyanax sp.", "A. mexicanus", NDT53$NativeFound2005)


## Replace Pseudoxiphophorus sp: -----------------------------------------------
NDT67$ExoticsFound2005 <- gsub("Heterandria sp.", "H. bimaculata", NDT67$ExoticsFound2005)
NDT53$ExoticsFound2005 <- gsub("Heterandria sp.", "H. bimaculata", NDT53$ExoticsFound2005)


## Replace Chirostoma sp: ------------------------------------------------------
NDT67$HistoricalNatCatalog <- gsub("Menidia sp.", "M. riojai", NDT67$HistoricalNatCatalog)
NDT53$HistoricalNatCatalog <- gsub("Menidia sp.", "M. riojai", NDT53$HistoricalNatCatalog)

NDT67$HistoricalNatBroad <- gsub("Menidia sp.", "M. riojai", NDT67$HistoricalNatBroad)
NDT53$HistoricalNatBroad <- gsub("Menidia sp.", "M. riojai", NDT53$HistoricalNatBroad)

NDT67$NativeFound2005 <- gsub("Menidia sp.", "M. riojai", NDT67$NativeFound2005)
NDT53$NativeFound2005 <- gsub("Menidia sp.", "M. riojai", NDT53$NativeFound2005)




#===============================================================================
# Initial formatting of data:---------------------------------------------------
#===============================================================================
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


## SUBSET 53:-------------------------------------------------------------------
dflist53  <- lapply(varlist53, function(x) {create_long_df(variable=x, n=53)})
names(dflist53) <- names(varlist53)

dflist53  <- lapply(dflist53, function(x) {x$Species <- tax$Genus_species_ECF[match(x$Species, tax$Short)];x # Updated names
                                           x$Species <- gsub(" cf ", " ", x$Species); x          # Remove " cf " so that e.g. Poecilia cf butleri and Poecilia butleri are regarded as the same sps
                                           x$Abundance <- rep(1, nrow(x));x})                    # Create long data.frame

dflist53 <- lapply(dflist53, function(x) {x <- x %>% distinct()})
dflist53 <- lapply(dflist53, function(x) {x <- spread(x, key="Species", value="Abundance", fill=0)})


## SUBSET 67:-------------------------------------------------------------------
dflist67  <- lapply(varlist67, function(x) {create_long_df(variable=x, n=67)})
names(dflist67) <- names(varlist67)

dflist67  <- lapply(dflist67, function(x) {x$Species <- tax$Genus_species_ECF[match(x$Species, tax$Short)];x 
                                           x$Species <- gsub(" cf ", " ", x$Species); x          
                                           x$Abundance <- rep(1, nrow(x));x})                    
dflist67 <- lapply(dflist67, function(x) {x <- x %>% distinct()})
dflist67 <- lapply(dflist67, function(x) {x <- spread(x, key="Species", value="Abundance", fill=0)})



## SUBSET 83:-------------------------------------------------------------------
dflist83  <- lapply(varlist83, function(x) {create_long_df(variable=x, n=83)})
names(dflist83) <- names(varlist83)

dflist83  <- lapply(dflist83, function(x) {x$Species <- tax$Genus_species_ECF[match(x$Species, tax$Short)];x 
                                           x$Species <- gsub(" cf ", " ", x$Species); x          
                                           x$Abundance <- rep(1, nrow(x));x})                    
dflist83 <- lapply(dflist83, function(x) {x <- x %>% distinct()})
dflist83 <- lapply(dflist83, function(x) {x <- spread(x, key="Species", value="Abundance", fill=0)})



# NOTE ######################################################################### 
# For all three runnings of the function "create_long_df", 
# warnings are fine and fixed within the same code
################################################################################



#===============================================================================
# Correction:-------------------------------------------------------------------
#===============================================================================
# Add occurrences of A monticola and remove cols that indicate presence of no fish
# A monticola is Algansea monticola in "El Huizcolote River" & "El Sacristan Spring"
# A monticola is Agonostomus monticola in "De Comala River", "Salado River"


l_historical <- list(dflist53[[1]],dflist53[[2]],
                 dflist67[[1]],dflist67[[2]],
                 dflist83[[1]],dflist83[[2]])

l_contemporary <- list(dflist53[[4]], 
                     dflist67[[4]],
                     dflist83[[4]])



for (i in 1:length(l_historical)) {
  dt <- l_historical[[i]]
  
  dt$`Dajaus monticola` <- rep(0, nrow(dt))
  dt$`Algansea monticola` <- rep(0, nrow(dt))
  
  dt$`Algansea monticola` <- ifelse(dt$SiteNameE %in% c("El Huizcolote River", "El Sacristan Spring"), 1, 0)
  dt$`Dajaus monticola` <- ifelse(dt$SiteNameE %in% c("De Comala River", "Salado River"), 1, 0) # The occurrence in el Salado should appear in [[1]] so here this is corrected
  dt <- dt[,!names(dt)=="A monticola"]
  l_historical[[i]] <- dt
}

for (i in 1:length(l_contemporary)) {
  dt <- l_contemporary[[i]]
  
  dt$`Dajaus monticola` <- rep(0, nrow(dt))
  #dt$`Algansea monticola` <- rep(0, nrow(dt)) # not found anymore
  
  dt$`Dajaus monticola` <- ifelse(dt$SiteNameE %in% c("De Comala River", "Salado River"), 1, 0) # continues to occur
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
dflist67[[4]] <- l_contemporary[[2]]
dflist83[[4]] <- l_contemporary[[3]]

lapply(dflist53, function(x) {unique(names(x))})
lapply(dflist67, function(x) {unique(names(x))})
lapply(dflist83, function(x) {unique(names(x))})

# Remove columns for non-species:
non_species <- c("none", "no aquatic life","Dry", "foul water", 
                 "ceased to exist", "it was not the right locality", 
                 "is the locality", "-", "(?)", "?", "trompos")
dflist53 <- lapply(dflist53, function(x) {x <- x[,!names(x) %in% non_species]})
dflist67 <- lapply(dflist67, function(x) {x <- x[,!names(x) %in% non_species]})
dflist83 <- lapply(dflist83, function(x) {x <- x[,!names(x) %in% non_species]})

# "Zoogoneticus tequila y Ameca  splendens":
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




#===============================================================================
# Generate lists for SM:--------------------------------------------------------
#===============================================================================
path_list_SM <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists/SM_Chapter5"


## Localities:------------------------------------------------------------------
LocalityList <- NDT83[, names(NDT83) %in% c("SiteNameE", "DrainageBasinE",
                                            "SiteType", "Latitude", "Longitude")]
LocalityList$Subset <- "83"
LocalityList$Subset <- ifelse(LocalityList$SiteNameE %in% NDT67$SiteNameE, "67", LocalityList$Subset)
LocalityList$Subset <- ifelse(LocalityList$SiteNameE %in% NDT53$SiteNameE, "53", LocalityList$Subset)

LocalityList$Comment <- ifelse(LocalityList$Subset=="53", "Assemblage data are comparable for both periods", NA)
LocalityList$Comment <- ifelse(LocalityList$Subset=="67", "The historical assemblage could be reconstructed but in 2005 the locality did not sustain life anymore", LocalityList$Comment)
LocalityList$Comment <- ifelse(LocalityList$Subset=="83", "Assemblage data are not comparable between time periods", LocalityList$Comment)
sort(unique(LocalityList$SiteNameE[LocalityList$Comment=="Assemblage data are not comparable between time periods"]))
sum(LocalityList$Comment != "Assemblage data are not comparable between time periods") # 67

LocalityList <- LocalityList %>% relocate(c(SiteNameE, DrainageBasinE), .before=Latitude)

write.csv(LocalityList, file=paste0(path_list_SM, "/LocalityList.csv"), row.names=F) 
# NOTE: the info in this table was modified manually : 
# excluded any site not in 67 & included more specific comments based on NDT.


## Anecdotics:------------------------------------------------------------------

anec <- NDT83$Anectodics
anec <- anec[!is.na(anec)] # % observations
anec

anectax <- c("Sicydium multipunctatum", "Oncorhynchus mykiss", 
             "Rhamdia sp", "Cyprinus carpio", "trompos", "peces limpiadores",
             "Bagre sp", "Micropterus salmoides")
status <- c("N", "?", 
             "?", "E", "?", "?",
             "?", "E")
anectax <- data.frame("Species"=anectax, "Status"=status)

write.csv(anectax, file=paste0(path_list_SM, "/anectax.csv"), row.names=F)


## Species List: ---------------------------------------------------------------
#sort(unique(dflist83[[4]]$SiteNameE))
vec_taxa83 <- sort(unique(c(names(dflist83[[1]]), names(dflist83[[2]]),
                            names(dflist83[[3]]), names(dflist83[[4]]),
                            names(dflist83[[5]]))))
vec_taxa83 <- vec_taxa83[!vec_taxa83 %in% c("SiteNameE", "DrainageBasinE")]

vec_taxa <- sort(unique(c(names(dflist67[[1]]), names(dflist67[[2]]),
                          names(dflist67[[3]]), names(dflist67[[4]]),
                          names(dflist67[[5]])))) # use 67 site subset (i.e., subset for analyses)
vec_taxa <- vec_taxa[!vec_taxa %in% c("SiteNameE", "DrainageBasinE")] 
sort(unique(vec_taxa)) # 95 species 
# including M. consocia (TBD1) & S. austrinus (TBD3), but not M. patzcuaro (TBD2)
# excludes H cyanogutattus

setdiff(vec_taxa83, vec_taxa) # taxa in 83 list, but not in 67 subset
# "Allotoca meeki"
# "Ameca splendens"
# "Carassius auratus"
# "Chirostoma attenuatum"
# "Chirostoma estor"
# "Cualac tessellatus"
# "Gambusia affinis"
# "Poecilia latipunctata"
# "Skiffia francesae"
# "Chirostoma patzcuaro"
# "Zoogoneticus tequila"
# "Herichthys cyanoguttatus"
# "Zoogoneticus tequila y Ameca splendens"
# PLUS Poecilia sp, Poeciliopsis sp, Oreochromis sp, Chirostoma sp, Astyanax sp, & Pseudoxiphophorus sp

setdiff(vec_taxa, vec_taxa83) # 0, OK

taxa <- data.frame(matrix(ncol=5, nrow=95))
names(taxa) <- c("Genus_species", "Previous", "Updated", "Regional_Status", "Comment")

taxa$Genus_species <- vec_taxa
taxa$Previous <- tax$Genus_species[match(taxa$Genus_species, tax$Genus_species_ECF)]  # re-updated April 2024
# remember that here, for each "taxa$Genus_species" there might be multiple "tax$Genus_species"(i.e. previous)
sum(is.na(taxa$Previous))
taxa$Previous[taxa$Genus_species=="Dajaus monticola"] <- "Agonostomus monticola"
taxa$Previous[taxa$Genus_species=="Algansea monticola"] <- "Algansea monticola"
taxa$Previous[taxa$Genus_species=="Gambusia senilis"] <- "Gambusia senilis"
#taxa$Genus_species[taxa$Previous=="Notropis calientis"] <- "Aztecula calientis"       # April 2024
#taxa$Genus_species[taxa$Previous=="Notropis amecae"] <- "Aztecula amecae"             # April 2024
#taxa$Genus_species[taxa$Previous=="Hybopsis boucardi"] <- "Graodus boucardi"          # April 2024
#taxa$Previous[taxa$Genus_species=="Graodus boucardi"] <- "Hybopsis/Notropis boucardi" # April 2024
taxa$Genus_species[taxa$Previous=="Menidia consocia"] <- "Chirostoma consocium"
taxa$Genus_species[taxa$Previous=="Scartomyzon austrinus"] <- "Moxostoma austrinum"

taxa$Updated <- ifelse(taxa$Previous==taxa$Genus_species, "no", "yes")
sum(is.na(taxa$Updated))
sum(is.na(taxa$Genus_species))


## Add authorities: ------------------------------------------------------------
require(taxize)
listsps <- sort(unique(taxa$Genus_species))
auth <- gnr_resolve(listsps, data_source_ids=11, canonical=FALSE)  
auth <- auth[!auth$matched_name=="Astyanax aeneus (Hensel, 1870)",] 
auth <- auth[!auth$matched_name=="Astyanax Stål, 1867",] 
auth <- auth[!auth$matched_name=="Characodon lateralis Garman, 1895",] 
auth <- auth[!auth$matched_name=="Goodea atripinnis Meek, 1907",] 
auth <- auth[!auth$matched_name=="Oreochromis niloticus (Greenwood, 1960)",] 
auth <- auth[!auth$matched_name=="Oreochromis Carnevale et al., 2003",] 
auth <- auth[!auth$matched_name=="Poecilia mexicana De Filippi, 1940",] 
auth$matched_name[auth$matched_name=="Poecilia"] <- "Poecilia Bloch & Schneider 1801" # ESC Catalog
auth <- auth[!auth$matched_name=="Poecilia Heinemann, 1870",]
auth <- auth[!auth$matched_name=="Poecilia Taczanowski, 1872",]
auth <- auth[!auth$matched_name=="Poecilia Schrank, 1802",]
auth <- auth[!auth$matched_name=="Pseudoxiphophorus bimaculatus (Steindachner, 1863)",]


taxa$Authority <- auth$matched_name[match(taxa$Genus_species, auth$user_supplied_name)]
sum(is.na(taxa$Authority))


## Regional status: -------------------------------------------------------------

introduced <- c("Oreochromis aureus", "Oreochromis mossambicus", "Oreochromis niloticus", 
                "Pseudoxiphophorus bimaculatus", "Pseudoxiphophorus jonesii", 
                "Cyprinus carpio", "Lepomis macrochirus", "Micropterus salmoides", 
                "Pomoxis nigromaculatus", "Carassius auratus", 
                "Xiphophorus maculatus", "Xiphophorus hellerii", "Xiphophorus variatus",
                "Gambusia yucatana", "Gambusia affinis", 
                "Poeciliopsis gracilis", "Poecilia reticulata",
                "Amatitlania nigrofasciata", "Astatotilapia burtoni") 


# NOTE 1: 22 (rm Poecilia sp, bc now that it's P sphenops it is native in some sites)
# rm also Poeciliopsis sp, Oreochromis sp & Pseudoxip. sp bc they've been renamed in addendum


# NOTE 2: Ctenopharyngodon idella (in paper but not in database, considered uncertain)
# additionalSM <- c("Menidia grandocule", "Menidia lucius", "Algansea avia", "Algansea lacustris")
# in Gesundheit & Macias Garcia 2018 (these are updates)

taxa$Regional_Status <- ifelse(taxa$Genus_species %in% introduced, "Introduced", "Native")      # OK
taxa$Regional_Status <- ifelse(taxa$Genus_species %in% c("Chirostoma chapalae"), "Native", taxa$Regional_Status)
#taxa$Comment <- ifelse(taxa$Genus_species %in% c("Chirostoma chapalae"), 
#                      "extirpated in native locs, but translocated in P. Cointzio", taxa$Comment) # this translocation is only in the 83 subset, in the 67 subset the ocurrence of C. chapalae is native

sum(taxa$Regional_Status=="Introduced") #17


## Add families: ---------------------------------------------------------------

# NOTE:
# "sp" except for Gila sp records removed below after addendum

species <- taxa$Genus_species
species <- species[! species %in% c("Gila sp")] 
#fams <- tax_name(species, get = "family")
#save(fams, file=paste0(path_list_SM, "/fams.RData"))
load(paste0(path_list_SM, "/fams.RData"))
taxa$Family <- fams$family[match(taxa$Genus_species, fams$query)]
taxa$Family[taxa$Genus_species %in% c("Pseudoxiphophorus jonesii", "Pseudoxiphophorus bimaculatus")] <- "Poeciliidae"
taxa$Family[taxa$Genus_species %in% c("Astyanax aeneus")] <- "Characidae"
taxa$Family[taxa$Genus_species %in% c("Amatitlania nigrofasciata",
                                      "Amphilophus istlanus",
                                      "Mayaheros beani")] <- "Cichlidae"
taxa$Family[taxa$Genus_species %in% c("Gila sp")] <- "Cyprinidae"
taxa$Family[taxa$Genus_species %in% c("Chirostoma mezquital")] <- "Atherinopsidae"
taxa$Family[taxa$Genus_species %in% c("Dajaus monticola")] <- "Mugilidae"
taxa$Family[taxa$Genus_species %in% c("Neotoca bilineata", "Xenoophorus captivus", 
                                      "Xenotoca variata", "Hubbsina turneri", "Skiffia bilineata")] <- "Goodeidae"
taxa$Family[taxa$Genus_species %in% c("Aztecula sallaei", "Graodus boucardi")] <- "Leuciscidae"


# Orders: ----------------------------------------------------------------------
#ord <- tax_name(species, get = "order")
#save(ord, file=paste0(path_list_SM, "/ord.RData"))
load(paste0(path_list_SM, "/ord.RData"))
taxa$Order <- ord$order[match(taxa$Genus_species, ord$query)]
taxa$Order[taxa$Genus_species %in% c("Pseudoxiphophorus jonesii", "Pseudoxiphophorus bimaculatus")] <- "Cyprinodontiformes"
taxa$Order[taxa$Genus_species %in% c("Astyanax aeneus")] <- "Characiformes"
taxa$Order[taxa$Genus_species %in% c("Amatitlania nigrofasciata",
                                     "Amphilophus istlanus", "Herichthys cyanoguttatus", 
                                     "Mayaheros beani")] <- "Cichliformes"
taxa$Order[taxa$Genus_species %in% c("Gila sp")] <- "Cypriniformes"
taxa$Order[taxa$Genus_species %in% c("Chirostoma mezquital")] <- "Atheriniformes"
taxa$Order[taxa$Genus_species %in% c("Dajaus monticola")] <- "Mugiliformes"
taxa$Order[taxa$Genus_species %in% c("Neotoca bilineata", "Xenoophorus captivus", 
                                     "Xenotoca variata", "Hubbsina turneri", "Skiffia bilineata")] <- "Cyprinodontiformes"
taxa$Order[taxa$Genus_species %in% c("Aztecula sallaei", "Graodus boucardi")] <- "Cypriniformes"


# tbd: -------------------------------------------------------------------------
taxa$Family[taxa$Previous %in% "Menidia consocia"] <- "Atherinopsidae"
#taxa$Family[taxa$Previous %in% "Menidia patzcuaro"] <- "Atherinopsidae"
taxa$Family[taxa$Previous %in% "Scartomyzon austrinus"] <- "Catostomidae"

taxa$Order[taxa$Previous %in% "Menidia consocia"] <- "Atheriniformes"
#taxa$Family[taxa$Previous %in% "Menidia patzcuaro"] <- "Atherinopsidae"
taxa$Order[taxa$Previous %in% "Scartomyzon austrinus"] <- "Cypriniformes"


taxa <- taxa[!names(taxa) %in% "Genus_species"]
taxa <- taxa %>% relocate(c(Order, Family, Authority), .before=Previous)

write.csv(taxa, file=paste0(path_list_SM, "/taxa.csv"), row.names=F) 
# NOTE: the list above reviewed manually for Appendix Table.



#===============================================================================
# Data for assemblage-level analyses:-------------------------------------------
#===============================================================================
path_lists <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Lists"

## SUBSET 67: ------------------------------------------------------------------
ContAll <- full_join(dflist67[[4]], dflist67[[5]], by=c("SiteNameE", "DrainageBasinE")) # contemporary all

ContAll$`Goodea atripinnis`<- ContAll$`Goodea atripinnis.x`| ContAll$`Goodea atripinnis.y`
ContAll$`Poecilia butleri`<- ContAll$`Poecilia butleri.x`| ContAll$`Poecilia butleri.y`
ContAll$`Poecilia sphenops`<- ContAll$`Poecilia sphenops.x`| ContAll$`Poecilia sphenops.y`

ContAll <- ContAll[, ! names(ContAll) %in% c("Goodea atripinnis.x", "Goodea atripinnis.y",
                                             "Poecilia butleri.x", "Poecilia butleri.y",
                                             "Poecilia sphenops.x", "Poecilia sphenops.y")]
check <- sort(unique(names(ContAll)))
sum(is.na(ContAll[,-c(1:2)]))
ContAll[is.na(ContAll)] <- 0
sum(is.na(ContAll[,-c(1:2)]))

# View(ContAll[names(ContAll) %in% c("SiteNameE", "Poecilia butleri", "Goodea atripinnis", "Poecilia sphenops")]) 
# OK P butleri not in Sta Maria Lagoon (historical exotic)
sum(ContAll$`Poecilia butleri`)  # OK (checked visually)
sum(ContAll$`Goodea atripinnis`) # OK
sum(ContAll$`Poecilia sphenops`) # OK
ContAll$Period <- "ContAll"


HNC <- dflist67[[1]]
HNC$Period <- "HNC"
HNB <- dflist67[[2]]
HNB$Period <- "HNB"
# ignore dflist67[[3]] (historical exotics)
ContN <- dflist67[[4]]
ContN$Period <- "ContN"
ContE <- dflist67[[5]]
ContE$Period <- "ContE"
l67 <- bind_rows(HNC, HNB, ContN, ContE, ContAll)
l67 <- l67 %>% relocate(Period, .after=DrainageBasinE)
l67[is.na(l67)] <- 0
sum(is.na(l67)) # 0

save(l67, file=paste0(path_lists, "/l67.RData"))       # data.frame of 305 obs --> 67*4 + 37 (ContE)


## SUBSET 53: ------------------------------------------------------------------
ContAll53 <- full_join(dflist53[[4]], dflist53[[5]], by=c("SiteNameE", "DrainageBasinE")) # contemporary all
names(ContAll53)
ContAll53$`Goodea atripinnis`<- ContAll53$`Goodea atripinnis.x`| ContAll53$`Goodea atripinnis.y`
ContAll53$`Poecilia butleri`<- ContAll53$`Poecilia butleri.x`| ContAll53$`Poecilia butleri.y`
ContAll53$`Poecilia sphenops`<- ContAll53$`Poecilia sphenops.x`| ContAll53$`Poecilia sphenops.y`

ContAll53 <- ContAll53[, ! names(ContAll53) %in% c("Goodea atripinnis.x", "Goodea atripinnis.y",
                                                   "Poecilia butleri.x", "Poecilia butleri.y",
                                                   "Poecilia sphenops.x", "Poecilia sphenops.y")]
sum(is.na(ContAll53[,-c(1:2)]))
ContAll53[is.na(ContAll53)] <- 0
sum(is.na(ContAll53[,-c(1:2)]))
ContAll53$Period <- "ContAll"

HNC <- dflist53[[1]]  #(53)
HNC$Period <- "HNC"
HNB <- dflist53[[2]]  #(53)
HNB$Period <- "HNB"
# ignore dflist67[[3]] (historical exotics)
ContN <- dflist53[[4]]  #(53)
ContN$Period <- "ContN"
ContE <- dflist53[[5]]  #(37)
ContE$Period <- "ContE"
l53 <- bind_rows(HNC, HNB, ContN, ContE, ContAll53)
l53 <- l53 %>% relocate(Period, .after=DrainageBasinE)
l53[is.na(l53)] <- 0
sum(is.na(l53)) # 0

save(l53, file=paste0(path_lists, "/l53.RData"))      # data.frame of 249 obs --> 53*4 + 37 (ContE)


# End of script ################################################################








# NOT INCLUDED #################################################################


#===============================================================================
# IUCN Status: -----------------------------------------------------------------
#===============================================================================
# NOTE: ########################################################################
# To develop for manuscript?
################################################################################

#require(red)
#redlist_api_key <- "YOUR_API_KEY"
#species <- taxa$Genus_species[!]
#status_list <- lapply(, function(s) {
#  rl_search(query = s, token = redlist_api_key)
#})

#for (i in seq_along(species)) {
#  cat(species[i], ": ", status_list[[i]]$result[[1]]$category, "\n")
#}

# https://www.gbif.org/dataset/19491596-35ae-4a91-9a98-85cf505f1bd3
#iucn_tax <- read_delim("C:/Users/afe1/Downloads/iucn-2022-1/taxon.txt", delim = "\t", col_names = F)
#iucn_snap <- read_delim("C:/Users/afe1/Downloads/iucn-2022-1/distribution.txt", delim = "\t", col_names = F)


#species <- taxa$Genus_species
#species <- species[!species %in% c("Astyanax sp", "Chirostoma sp", "Gila sp",
#                                   "Poecilia sp", "TBD1","TBD2","TBD3", "Pseudoxiphophorus sp", "Oreochromis sp",
#                                   "Poeciliopsis sp")]

#iucn_tax$tax <- paste0(str_split_fixed(iucn_tax$X2, " ", 3)[,1], " ", str_split_fixed(iucn_tax$X2, " ", 3)[,2])
#sub_iucn_tax <- iucn_tax[iucn_tax$tax %in% species,]
#sub_iucn_tax <- sub_iucn_tax[!sub_iucn_tax$X1 %in% c("191249_1","134692299_1",
#                                                     "2270", "4529_1", "4530_1",
#                                                     "166052_8", "180896_1", "166066_1",
#                                                     "23116_1", "133768576_4", "133768576_5",
#                                                     "133768576_6", "133768576_7", "82627914_2"),]

#setdiff(species, sub_iucn_tax$tax)
#sub_iucn_tax$status <- iucn_snap$X4[match(sub_iucn_tax$X1, iucn_snap$X1)]

# Aztecula amecae Extinct in the wild (https://www.iucnredlist.org/species/14881/546437)
# Oreochromis aureus Least Concern (https://www.iucnredlist.org/species/166933/6293101)
# Notropis boucardi Endangered (https://www.iucnredlist.org/species/191271/1974646)
# Girardinichthys turneri Critically Endangered (https://www.iucnredlist.org/species/132523146/497499)
# Aztecula sallaei Least Concern (https://www.iucnredlist.org/species/191255/1974392)
# Ictalurus dugesii // Can't find, so Data Deficient (https://www.iucnredlist.org/search?taxonomies=102293&searchType=species)
# Tampichthys dichroma Critically endangered (https://www.iucnredlist.org/species/6624/3135282)


#taxa$IUCN_Status <- sub_iucn_tax$status[match(taxa$Genus_species, sub_iucn_tax$tax)]
#sort(unique(taxa$Genus_species[is.na(taxa$IUCN_Status)]))
#sort(unique(taxa$IUCN_Status))
#taxa$IUCN_Status <- ifelse(taxa$Genus_species %in% c("Aztecula amecae"), "Extinct in the Wild", taxa$IUCN_Status)
#taxa$IUCN_Status <- ifelse(taxa$Genus_species %in% c("Oreochromis aureus",
#                                                     "Aztecula sallaei"), "Least Concern", taxa$IUCN_Status)
#taxa$IUCN_Status <- ifelse(taxa$Genus_species %in% c("Girardinichthys turneri",
#                                                     "Tampichthys dichroma"), "Critically Endangered", taxa$IUCN_Status)
#taxa$IUCN_Status <- ifelse(taxa$Genus_species %in% c("Graodus boucardi"), "Endangered", taxa$IUCN_Status)
#taxa$IUCN_Status <- ifelse(taxa$Genus_species %in% c("Ictalurus dugesii"), "Data Deficient", taxa$IUCN_Status)

