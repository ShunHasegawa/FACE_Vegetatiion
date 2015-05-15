##########################
# Organise soil variable #
##########################
SoilDf <- read.csv("Data//FACE_Environmtal_Variable.csv", header = TRUE)

#############
# TC and TN #
#############
TcnDF <- read.xlsx2("Data/SoilVariables/Soil.TC.TN.xlsx", 
                    sheetName = "FACE.soil.tc.tn", 
                    header = TRUE, startRow = 1, endRow = 97, 
                    stringsAsFactors = FALSE)
# organise
names(TcnDF)[grepl("tc|tn", names(TcnDF))] <- c("TotalC", "TotalN")
TcnDF[, c("TotalC", "TotalN")] <- apply(TcnDF[, c("TotalC", "TotalN")], 2, as.numeric)
TcnDF$year <- factor(TcnDF$year, labels = c(2013, 2014))
TcnDF <- subsetD(TcnDF, depth == "0-10", select = -depth)

# Ring mean
TcnDF_Ring <- ddply(TcnDF, .(year, ring), summarise,
                    TotalC = mean(TotalC), TotalN = mean(TotalN))

###########################
# Total P measured by Cat #
###########################
TpDF <- read.xlsx2("Data/SoilVariables/EucAFCE P sum 30_05_14.xlsx", 
                   sheetIndex = 1, 
                   header = TRUE, startRow = 1, endRow = 25, 
                   stringsAsFactors = FALSE)
# organise
TpDF <- TpDF[, c("Ring", "Quarter", "Total.CM.mg.kg.1")]
names(TpDF) <- c("ring", "plot", "TotalP_CM")
TpDF <- within(TpDF, {
  TotalP_CM <- as.numeric(TotalP_CM)
  year <- "2013" # not too sure which soil was used... need to check with Cat
})

# Ring mean
TpDF_Ring <- ddply(TpDF, .(year, ring), summarise, TotalP_CM = mean(TotalP_CM))

#######################
# Soil moist and temp #
#######################
load("Data//SoilVariables/soil.var_ring.means.RData")

# 1st year is only from August so just use August-December
SoilMTdf <- subsetD(ring.means, 
                    month(Date) %in% c(8:12) & year(Date) %in% c(2012:2014),
                    select = -co2)
SoilMTdf$year <- factor(year(SoilMTdf$Date), labels = c(2013:2015))
dlply(SoilMTdf, .(year), summary)

# Ring mean
SiteVec <- c("year", "Date", "ring")
Probes <-names(SoilMTdf)[!names(SoilMTdf) %in% SiteVec] 

SoilMTdf_Ring <- ddply(SoilMTdf, .(year, ring), function(x) colMeans(x[, Probes], na.rm = TRUE))

#########
# light #
#########
load("output//Data/FACE_Light_DayMean.RData")
head(FACE_Light_DayMean)
# understorey and canopy PAR is highly correlated so just use understorey
# October 2012 values are little bit weird. It goes up on a sudden.
# Just use Nov-Dec for the time being cause measurements in these months are
# complete in all the three years.
LightNovDec <- subsetD(FACE_Light_DayMean, month(Date) %in% c(11, 12))
LightNovDec$year <- factor(year(LightNovDec$Date), labels = c(2013:2015))

FlorLight_Ring <- ddply(LightNovDec, .(year, ring), summarise, FloorPAR = mean(FloorPAR, na.rm = TRUE))

#######
# IEM #
#######
load("Data/SoilVariables/FACE_IEM.RData")

# use only Nov-Jan
iemNovJan <- subsetD(iem, month(date) == c(1, 11, 12))
# add year; January is counted as an year before
iemNovJan$year <- with(iemNovJan, factor(ifelse(month(date) == 1, year(date), year(date) + 1)))

iem_ring <- ddply(iemNovJan, .(ring, year), summarise,
                  IEM_no = mean(no, na.rm = TRUE), 
                  IEM_nh = mean(nh, na.rm = TRUE), 
                  IEM_p = mean(p, na.rm = TRUE))

#################
# Soil Extracts #
#################
load("Data/SoilVariables/extractable.RData")
tdf <- within(extr, {
  M <- month(date)
  Y <- year(date)
})
xtabs(~ Y + M, data = tdf)

# use december
ExtrDec <- subsetD(extr, month(date) == 12)
ExtrDec$year <- factor(year(ExtrDec$date), labels = c(2013, 2014))
Extract_ring <- ddply(ExtrDec, .(year, ring), summarise, 
                      Ext_no = mean(no, na.rm = TRUE),
                      Ext_nh = mean(nh, na.rm = TRUE),
                      Ext_p = mean(po, na.rm = TRUE))

##################
# Mineralisation #
##################
load("Data/SoilVariables/FACE_mineralisation.RData")
tdf <- within(mine, {
  M <- month(date)
  Y <- year(date)
})
xtabs(~ Y + M, data = tdf)
# use January
MineJan <- subsetD(mine, month(date) == 1)
MineJan$year <- factor(year(MineJan$date))

Mineralisation_ring <- ddply(MineJan, .(year, ring), summarise,
                             n.min = mean(n.min, na.rm = TRUE), 
                             nitrification = mean(nitrification, na.rm = TRUE), 
                             p.min = mean(p.min, na.rm = TRUE))
#############
# Lysimeter #
#############
load("Data/SoilVariables/FACE_lysimeter.Rdata")
tdf <- within(lys, {
  M <- month(date)
  Y <- year(date)
})

tdfS <- subsetD(tdf, depth == "shallow")
xtabs(~ Y + M, data = tdfS)
# no complete months for even two years.. don't use for the time being

###############
# Combine all #
###############
EnvVarDF <- Reduce(function(...) merge(..., by = c("year", "ring"), all.x = TRUE), 
                   list(SoilDf, TcnDF_Ring, TpDF_Ring, SoilMTdf_Ring, FlorLight_Ring, 
                        iem_ring, Extract_ring, Mineralisation_ring))
EnvVarDF <- within(EnvVarDF, {
  year <- factor(year)
  ring <- factor(ring)
  co2 <- factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
})

save(EnvVarDF, file = "output/Data/FACE_EnvironmenVars.RData")
