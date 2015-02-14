###########
## Stats ##
###########

# organise data frame
SiteName <- c("year", "ring", "plot", "position", "cell")
SppName <- names(veg.face)[!names(veg.face) %in% SiteName]

SiteMatrix <- veg.face[, SiteName]
SppMatrix <- veg.face[, SppName]

######
# CA #
######
source("R/CA.analysis.R")

#######
# RDA #
#######


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