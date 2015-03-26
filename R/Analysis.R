rm(list=ls(all=TRUE))

source("R/Packages.R")
source("R/functions.R")

################
# Process Data #
################
# source("R/FACE_VegetationDataSheet2015.R")

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

#######################
# organise data frame #
#######################

# all spp----
veg.face <- within(vdf, {
  co2 = factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
})
SiteName <- c("year", "ring", "co2", "plot", "position", "cell")
SppName <- names(veg.face)[!names(veg.face) %in% SiteName]

# plot sum
PlotSumVeg <- ddply(veg.face, .(year, co2, ring, plot), function(x) colSums(x[, SppName]))

# ring sum
RingSumVeg <- ddply(PlotSumVeg, .(year, co2, ring), function(x) colSums(x[, SppName]))

# PFG matrix----
RingSumPFGMatrix <- dcast(year + co2 + ring ~ PFG, data = subset(veg, !is.na(PFG)), sum)
colSums(RingSumPFGMatrix[4:11])

# remove lichen and wood, also add interaction term
RingSumPFGMatrix <- within(RingSumPFGMatrix, {
  Lichen = NULL
  wood = NULL
  c3_4 = NULL
  yco = year:co2
})
PFGName <- c("c3", "c4", "legume", "moss", "Non_legume")


########
# Figs #
########
source("R//Figs.R")

#########
# Stats #
#########
source("R/Stats.R")
