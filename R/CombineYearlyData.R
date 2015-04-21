################################
# Combine 2013, 2014 2015 Data #
################################

# 2013 data
# source("R/veg.12.process.R")
load("output/Data/FACE_Vegetation_2013.RData")

# process 2014
# source("R/prcss.veg.2014.R")
load("output/Data/FACE_Veg2014.RData")

veg.face <- rbind.fill(df2013, veg.14)

####################
# organise dataset #
####################
# turn na into 0
veg.face[is.na(veg.face)] <- 0

# sort columns
SiteVec <- c("year", "ring", "plot", "position", "cell")
SppVec <- names(veg.face)[!names(veg.face) %in% SiteVec]
veg.face <- veg.face[, c(SiteVec, sort(SppVec))]

## save
save(veg.face, file = "output/Data/FACE_Vegetation_Raw_2013&2014.RData")
