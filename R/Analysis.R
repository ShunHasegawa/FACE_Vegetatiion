rm(list=ls(all=TRUE))

library(xlsx)
library(plyr)
library(XLConnect)
library(ggplot2)
library(reshape)

source("R/functions.R")
################
# Process Data #
################
# source("R/FACE.vegetation.2014.management.R")

# Row data for multi variate analysis
load("output/Data//FACE_Vegetation_Raw.RData")

# Data frame with plant functional groups
load("output//Data//FACE_Vegetation_PFG.RData")

########
# Figs #
########
source("R//Figs.R")


