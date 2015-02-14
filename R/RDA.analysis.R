########################
# Redundancey analysis #
########################

###############
# All species #
###############

# each year separately----
rdaList <- dlply(RingSumVeg, .(year), function(x) rda(log(x[, SppName] + 1) ~ co2, data = x))

par(mfrow = c(1, 2))
l_ply(rdaList, plot)

# two years together----
RingSumVeg$yco <- RingSumVeg$year:RingSumVeg$co2

RDAres <- rda(log(RingSumVeg[, SppName] + 1) ~ yco, data = RingSumVeg[, !names(RingSumVeg) %in% SppName])
RDAres
plot(RDAres)

#######
# PFG #
#######

# PFG matrix
RingSumPFGMatrix <- dcast(year + co2 + ring ~ PFG, data = subset(veg, !is.na(PFG)), sum)

# remove lichen and wood, also add interaction term
RingSumPFGMatrix <- within(RingSumPFGMatrix, {
  Lichen = NULL
  wood = NULL
  yco = year:co2
})

# rda
PFGName <- c("c3", "c3_4", "c4", "legume", "moss", "Non_legume")

# two years together----
PFGrdaRes <- rda(log(RingSumPFGMatrix[, PFGName] + 1) ~ yco, data = RingSumPFGMatrix)
PFGrdaRes
plot(PFGrdaRes)
anova(PFGrdaRes)
