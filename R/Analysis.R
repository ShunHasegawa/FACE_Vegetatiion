rm(list=ls(all=TRUE))

library(xlsx)
library(plyr)
library(XLConnect)
library(ggplot2)
library(reshape)
library(lme4)
library(packrat)
library(car)
library(mvabund)

source("R/functions.R")
################
# Process Data #
################
# source("R/FACE.vegetation.2014.management.R")

# Row data for multi variate analysis
load("output/Data//FACE_Vegetation_Raw.RData")

# Data frame with plant functional groups
load("output//Data//FACE_Vegetation_PFG.RData")

# co2 factor
FACE.veg.rslt$co2 <- factor(ifelse(FACE.veg.rslt$ring %in% c(1, 4, 5), "elev", "amb"))

########
# Figs #
########
source("R//Figs.R")

#########
# Stats #
#########
source("R/Stats.R")
