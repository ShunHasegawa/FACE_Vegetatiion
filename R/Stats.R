###########
## Stats ##
###########

# organise data frame
veg.face <- within(veg.face, {
  co2 = factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
})
SiteName <- c("year", "ring", "co2", "plot", "position", "cell")
SppName <- names(veg.face)[!names(veg.face) %in% SiteName]

# plot sum
PlotSumVeg <- ddply(veg.face, .(year, co2, ring, plot), function(x) colSums(x[, SppName]))

# ring sum
RingSumVeg <- ddply(PlotSumVeg, .(year, co2, ring), function(x) colSums(x[, SppName]))

######
# CA #
######
source("R/CA.analysis.R")

#######
# RDA #
#######
source("R/RDA.analysis.R")

########
# CAP  #
########
source("R/CAP.analysis.R")

#######
# GLM #
#######
# source("R/mvabund.analysis.R")

#########
# C3:C4 #
#########
source("R/Stats_C3_4_Grass.R")

#########################
# native vs. introduced #
#########################
source("R/Stat_origin.R")