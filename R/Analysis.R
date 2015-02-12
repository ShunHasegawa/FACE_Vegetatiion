rm(list=ls(all=TRUE))

source("R/Packages.R")
source("R/functions.R")

################
# Process Data #
################
# source("R/FACE.vegetation.2014.management.R")

# Raw data for multi variate analysis
load("output/Data//FACE_Vegetation_Raw.RData")

# Data frame with plant functional groups
load("output//Data//FACE_Vegetation_PFG.RData")

# remove unknown spp
veg <- FACE.veg.rslt[!grepl("Unknown", FACE.veg.rslt$variable), ]

# co2, block and id
veg <- within(veg, {
  block <- recode(ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
  co2 <- factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
  id <- ring:plot
})

########
# Figs #
########
source("R//Figs.R")

#########
# Stats #
#########
source("R/Stats.R")
