rm(list=ls(all=TRUE))

## ---- LoadData
source("R/Packages.R")
source("R/functions.R")

################
# Process Data #
################
# source("R/FACE_VegetationDataSheet2015.R")

# Raw data for multi variate analysis
load("output//Data/FACE_DominantVegetation_Raw_2013_2015.RData")

# Data frame with plant functional groups
load("output/Data/FACE_DominantVegetation_PFG_2015.RData")

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
veg.face <- within(DomVdf, {
  co2 <- factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
  block <- recode(ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
  id <- ring:plot
})
SiteName <- c("year", "block", "ring", "co2", "plot", "id", "position", "cell")
SppName <- names(veg.face)[!names(veg.face) %in% SiteName]

# plot sum
PlotSumVeg <- ddply(veg.face, .(year, block, co2, ring, plot, id), function(x) colSums(x[, SppName]))

# ring sum
RingSumVeg <- ddply(PlotSumVeg, .(year, block, co2, ring), function(x) colSums(x[, SppName]))

# PFG matrix----

# plot
PlotSumPFGMatrix <- dcast(year + block + co2 + ring + plot ~ PFG, 
                          data = subset(veg, !is.na(PFG)), sum)
colSums(PlotSumPFGMatrix[,6:11])

# remove lichen, also add interaction term
PlotSumPFGMatrix <- within(PlotSumPFGMatrix, {
  id = ring:plot
  yco = year:co2
})

PFGName <- c("c3", "c4", "legume", "moss", "Non_legume", "wood")

# ring
RingSumPFGMatrix <- ddply(PlotSumPFGMatrix, .(year, block, ring, co2, yco), 
                          function(x) colSums(x[, PFGName]))
  

# Diversity & eveness----
vegDF <- PlotSumVeg[, SppName]
siteDF <- PlotSumVeg[, !names(PlotSumVeg) %in% SppName]

DivDF <- within(siteDF,{
  H <- diversity(vegDF) # Shannon's index
  S <- specnumber(vegDF) # number of spp
  J <- H/log(S)  # Pielou's evenness
})

## ---- CreateFigs
########
# Figs #
########
source("R//Figs.R")

#########
# Stats #
#########
source("R/Stats.R")

# save all objects. This will be used when creating a summary document
save.image(file = "output//Data/AllObj.RData")
