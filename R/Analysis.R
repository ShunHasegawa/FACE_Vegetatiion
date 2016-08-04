rm(list=ls(all=TRUE))

## ---- LoadData
source("R/Packages.R")
source("R/functions.R")

# Process Data ------------------------------------------------------------

# source("R/CombineYearlyData.R")

# Raw data for multi variate analysis
load("output//Data/FACE_FullVegetation_Raw_2013_2016.RData")

# Data frame with plant functional groups
load("output/Data/FACE_FullVegetation_PFG_2016.RData")

# check unknown spp
all(!grepl("unknown", VegRes16$variable, ignore.case = TRUE))

# co2, block and id, combine sedge and grass, wood and shrub
veg <- within(VegRes16, {
  block <- recode(ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
  co2   <- factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
  id    <- ring:plot
  form  <- factor(ifelse(form %in% c("Tree", "Shrub"), "Wood",
                         ifelse(form %in% c("Grass", "Sedge", "Rush"), "Grass",
                                as.character(form)
                         )
  )
  )
}
)

# remove Euc seedlings as it's not reliable
veg <- subsetD(veg, variable != "Euc.seedling")

# remove c3_4 as it's really small number and hard to deal with c3_4..
# (Aristida.warburgii)
unique(veg$variable[veg$PFG == "c3_4"])
sum(veg$value[veg$PFG == "c3_4"])
veg <- subsetD(veg, PFG != "c3_4")

# remove lichen
veg <- subsetD(veg, form != "Lichen")
save(veg, file = "output/Data/EucFACE_understorey_vegetation_2012-2106.RData")
write.csv(veg, file = "output/Data/EucFACE_understorey_vegetation_2012-2106.csv",
          row.names = FALSE)

# organise data frame -----------------------------------------------------


# > all species -----------------------------------------------------------

FullVdf$Euc.seedling <- NULL
FullVdf$Aristida.warburgii <- NULL

veg.face <- within(FullVdf, {
  co2 <- factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
  block <- recode(ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
  id <- ring:plot
})
SiteName <- c("year", "block", "ring", "co2", "plot", "id", "position", "cell")
SppName <- names(veg.face)[!names(veg.face) %in% SiteName]

# plot sum
PlotSumVeg <- ddply(veg.face, .(year, ring, plot, block, co2, id), function(x) colSums(x[, SppName]))

# ring sum
RingSumVeg <- ddply(PlotSumVeg, .(year, ring, block, co2), function(x) colSums(x[, SppName]))


# > PFG -------------------------------------------------------------------

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
  

# diversity indices -------------------------------------------------------

# Diversity & eveness
vegDF <- PlotSumVeg[, SppName]
siteDF <- PlotSumVeg[, !names(PlotSumVeg) %in% SppName]

DivDF <- within(siteDF,{
  H <- diversity(vegDF) # Shannon's index
  S <- specnumber(vegDF) # number of spp
  J <- H/log(S)  # Pielou's evenness
})

# Identify dominant spp
SppSum <- ddply(veg, .(variable), summarise, value = sum(value))
SppSum <- SppSum[order(SppSum$value, decreasing = TRUE),]
SppSum <- within(SppSum, {
  Cov <- round(value * 100/sum(value), 3)
  CumSum <- cumsum(value)
  Dominant <- Cov >= 5 # species with >5 % coverage
})
DmSpp <- droplevels(SppSum$variable[SppSum$Dominant])
DmSpp
sum(SppSum$value[SppSum$Dominant])/sum(SppSum$value)


# figs --------------------------------------------------------------------
source("R//Figs.R")


# stats -------------------------------------------------------------------
source("R/Stats.R")

# save all objects. This will be used when creating a summary document all
# objects. This will be used when creating a summary document
sink("output/session_info.txt", append = TRUE)
paste0("/n/n/n",now(), "/n")
sessionInfo()
sink()

save.image(file = "output//Data/AllObj.RData")
