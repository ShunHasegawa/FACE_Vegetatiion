rm(list=ls(all=TRUE))

source("R/Packages.R")
source("R/functions.R")

################
# Process Data #
################
source("R/FACE_VegetationDataSheet2015.R")

# Raw data for multi variate analysis
load("output/Data//FACE_Vegetation_Raw_2013_2015.RData")

# Data frame with plant functional groups
load("output//Data//FACE_Vegetation_PFG_2015.RData")

# check unknown spp
all(!grepl("unknown", VegRes15$variable, ignore.case = TRUE))

# co2, block and id
veg <- within(VegRes15, {
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
