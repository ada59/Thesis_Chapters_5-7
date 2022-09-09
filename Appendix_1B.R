###########################################################################################
# Script: iNext
# AFE
# September 2022
###########################################################################################

# Libraries:------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(iNEXT)
library(gridExtra)


rm(list=ls())
myd <- getwd()
plot_iNext <- "C:/Users/Usuario/Documents/PHD/ThesisChapterMexico_I/TemporalChange_MexicanFish_C2/Plots/iNext" # Dir to save main plots


# Community data:-------------------------------------------------------------------------
# Sampling-unit-based incidence data
load(paste0(myd, "/HNC.RData"))     # Conservative Historical Native Assemblage
load(paste0(myd, "/HNB.RData"))     # Broad Historical Native Assemblage
load(paste0(myd, "/ContN.RData"))   # Contemporary Native Assemblage
load(paste0(myd, "/ContAll.RData")) # Contemporary Native + Exotic assemblage


###########################################################################################
# Example data:----------------------------------------------------------------------------
data(ant)
ant
ant[[1]]

data("ciliates")
ciliates

# list of incidence frequencies.
# the first entry = total number of sampling units.
# followed by the species incidence frequencies (rowsums in a species x sampling unit matrix) as shown below


# STEP 1 (already for ant, not for ciliates)
# Transform incidence raw data to incidence frequencies (iNEXT input format)

step0 <- DataInfo(ciliates, datatype ="incidence_raw") # data info
step0ant <- DataInfo(ant, datatype = "incidence_freq") # data info
step1 <- lapply(ciliates, as.incfreq) # data, correct format

# STEP 2:
# Compute species diversity with a particular level of sample size/coverage

out <- estimateD(ant, datatype = "incidence_freq", base="coverage",
                 level=0.985, conf=0.95) # I guess 0.985 is arbitrary.
# Unlike described in the iNext package, the argument q isn't included

# STEP 3:
# iNterpolation and EXTrapolation of Hill numbers
t <- round(seq(10, 500, length.out=20))
out3 <- iNEXT(ant$h500m, q=1, datatype="incidence_freq", size=t, se=FALSE)
out3$iNextEst

# STEP 4: 
# Plotting method for objects inheriting from class "iNEXT"
# Base R
plot(x=out3, type=1)
plot(x=out3, type=2)
plot(x=out3, type=3)

# ggplot:
ggiNEXT(x=out3, type=1) # sample-size-based rarefaction/extrapolation curve
ggiNEXT(x=out3, type=2) # sample completeness curve
ggiNEXT(x=out3, type=3) # coverage-based rarefaction/extrapolation curve


###########################################################################################
# Central México data (datatype="incidence-raw")
# Each survey is a site following Chao et al., 2014. Each assemblage is a sampling unit.

l <- list(HNC, HNB, ContAll)
sum(is.na(l))

l <- lapply(l, function(x) {rownames(x) <- x$SiteNameE;x})
l <- lapply(l, function(x) {within(x, rm(SiteNameE, DrainageBasinE))})
sum(rowSums(l[[3]])==0) #14, OK
l[[3]] <- l[[3]][rowSums(l[[3]])>0,]
l <- lapply(l, t)
l <- lapply(l, as.data.frame)

l_freq <- lapply(l, as.incfreq)
names(l_freq) <- c("HNC", "HNB", "ContAll")
d_info <- DataInfo(l_freq, datatype = "incidence_freq")
l_freq2 <- lapply(l_freq, function(x) {x[1] <- 2*x[1];x})
d_info2 <- DataInfo(l_freq2, datatype = "incidence_freq")

# size
na <- 67
nb <- 53*2
max(na, nb) #106

# coverage
ca <- max(d_info$SC)
cb <- min(d_info2$SC)
max(ca,cb) #0.9174


divSIZE <- estimateD(l_freq, datatype = "incidence_freq", base="size",
                    level=106, conf=0.95) 
divCOV <- estimateD(l_freq, datatype = "incidence_freq", base="coverage",
                 level=0.9174, conf=0.95) 

inext0 <- lapply(l_freq, function(x) {
  iNEXT(x, q=0, datatype="incidence_freq",se=FALSE)
})

inext0B <- lapply(l_freq, function(x) {
  iNEXT(x, q=0, datatype="incidence_freq", size=106, se=TRUE)
})

# HNC:
hnc1 <- ggiNEXT(x=inext0[[1]], type=1)
hnc2 <- ggiNEXT(x=inext0[[2]], type=2)
hnc3 <- ggiNEXT(x=inext0[[3]], type=3)

hncp <- grid.arrange(hnc1,hnc2, hnc3)

# HNB:
# ...

# ContAll
# ...

# single-assemblage incidence data with three orders q
# data(ant)
# size <- round(seq(10, 500, length.out=20))
# y <- iNEXT(ant$h500m, q=c(0,1,2), datatype="incidence_freq", size=size, se=FALSE)
# ggiNEXT(y, se=FALSE, color.var="Order.q")


###########################################################################################
# End of script ###########################################################################
