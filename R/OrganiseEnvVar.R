##########################
# Organise soil variable #
##########################
SoilDf <- read.csv("Data//FACE_Environmtal_Variable.csv", header = TRUE)
summary(SoilDf)

# TC and TN----
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

# Total P measured by Cat----
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

# Soil moist and temp----
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

# light----





Reduce(function(...) merge(..., by = c("year", "ring"), all.x = TRUE), 
       list(SoilDf, TcnDF_Ring, TpDF_Ring, SoilMTdf_Ring))


?Reduce


