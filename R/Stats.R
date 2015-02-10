###########
## Stats ##
###########

# organise data frame
FACE.veg.rslt <- within(FACE.veg.rslt, {
  co2 <- factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
  block <- recode(ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
  id <- ring:plot
})

######
# CA #
######
source("R/CA.analysis.R")

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