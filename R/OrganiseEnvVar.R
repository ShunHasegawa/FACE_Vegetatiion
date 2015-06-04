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
TcnDF$depth <- factor(paste0(TcnDF$depth, "cm"))

#######################
# TC and TN from Hiev #
#######################
HtcnDF <- read.csv("Data/SoilVariables/FACE_P0014_ALL_ SoilCN_June2012-Spet2014_V1.csv")
HtcnDF <- within(HtcnDF, {
  Date <- as.Date(dmy(as.character(Date)))
  ring <- factor(ring)
  year <- year(Date) + 1 # in order to combine with vegetation results from January next year
  Month <- month(Date)
})
HtcnDF <- subsetD(HtcnDF, plot != "Backfill")
names(HtcnDF)[grepl("tc|tn", names(HtcnDF))] <- c("TotalC", "TotalN")
ftable(xtabs( ~ year + Month + depth, data = HtcnDF))
# 2014 is completely missing
# no 20-30 cm for 2013

# compre the above two data frames
llply(list(HtcnDF, TcnDF), function(x) subset(x, ring == 1 & year == 2013))
  # 2013 is the same. so just need 2014 from the above.

# ring mean for each data frames at 0-10 cm from June or May
TcnDF_Ring <- ldply(list(subset(TcnDF, year == 2014), subset(HtcnDF, Month %in% c(5, 6))), 
                    function(x) {
                      ddply(x, .(year, ring, depth), summarise, 
                            TotalC = mean(TotalC), 
                            TotalN = mean(TotalN))
                      })

# check interaction beween ring and depth
par(mfrow = c(2, 3))
d_ply(TcnDF_Ring, .(year), function(x) 
  with(x, interaction.plot(depth, ring, TotalC, main = unique(x$year))))
d_ply(TcnDF_Ring, .(year), function(x) 
  with(x, interaction.plot(depth, ring, TotalN, main = unique(x$year))))
# There doesn't seem to be interaction so just use top 0-10 cm

TcnDF_Ring <- subsetD(TcnDF_Ring, depth == "0-10cm", select = -depth)

###############
# Dry soil pH #
###############
phDF <- read.csv("Data/SoilVariables/FACE_P0014_ALL_ SoilpH_2012to2014_V1.csv")

levels(phDF$Depth)
# unneccessary space. so remove

unique(phDF$Plot)
# some of pots contains ring number so remove

phDF <- within(phDF, {
  Date <- as.Date(dmy(Date))
  Ring <- factor(Ring)
  Plot <- as.character(factor(Plot)) # 1.0->1, 1.1->1.1
  Plot <- factor(gsub(".[.]", "", Plot)) # 1->1, 1.1->1
  Depth <- factor(gsub(" ", "", as.character(Depth)))
  year <- year(Date) + 1  # in order to combine with vegetation results from January next year
  Month <- month(Date)
  })
names(phDF)[c(3, 4, 6)] <- c("ring", "plot", "Drysoil_ph")
# remove NA
phDF <- phDF[!is.na(phDF$Drysoil_ph), ]

# which month and depth to be used
ftable(xtabs( ~ year + Month + Depth, data = phDF))
  # use May/June at 0-10cm

phDF_June <- subsetD(phDF, Depth == "0-10cm" & Month %in% 5:6)

# ring mean
phDF_ring <- ddply(phDF_June, .(year, ring), summarise, Drysoil_ph = mean(Drysoil_ph))

par(mfrow = c(1, 2))

# compare fresh soil pH and dry soil pH
with(phDF_ring, interaction.plot(year, ring, Drysoil_ph, ylim = c(5, 6)))
with(SoilDf, interaction.plot(year, ring, ph, ylim = c(5, 6)))
  # quite different

###############################
# Comprehensive soil analysis #
###############################
SoilChemDF <- read.csv("Data/SoilVariables/FACE_P0014_ALL_ ElementalAnalysis_2012to2014_V2.csv")

# remove unit
names(SoilChemDF) <- gsub("[..].*", "", names(SoilChemDF))

# re-structure
SoilChemDF <- within(SoilChemDF, {
  Date <- as.Date(ymd(Date))
  Ring <- factor(Ring)
  Plot <- factor(Plot)
  Depth <- factor(gsub(" ", "", as.character(Depth)))
  year <- year(Date) + 1
  Month <- month(Date)
  Sulfur <- NULL
})
names(SoilChemDF)[c(3, 4)] <- c("ring", "plot")
summary(SoilChemDF)

# remove rows whose chemical measurements are all NA
chems <- names(SoilChemDF)[!names(SoilChemDF) %in% c("Date", "Sample", "ring", "plot", "Depth", "Month", "year")]
SoilChemDF <- SoilChemDF[apply(SoilChemDF[, chems], 1, function(x) !all(is.na(x))), ]
summary(SoilChemDF)

# which month and depth to be used
ftable(xtabs( ~ year + Month + Depth, data = SoilChemDF))
  # aw what am I gonna do....?

# e.g. phosphorus
par(mfrow = c(2, 3))
l_ply(1:6, function(x)
boxplot(Phosphorus ~ Date, data = subset(SoilChemDF, Depth == "0-10cm" & ring == x), cex.axis = .7, 
        main = paste("Ring", x)))
# there is a monthly variation, so ideally want to use measurment from the same month

# use March at 0-10 cm for time being
SoilChemDF_Mar <- subsetD(SoilChemDF, Depth == "0-10cm" & Month == 3)
summary(SoilChemDF_Mar)

SoilChemDF_ring <- ddply(SoilChemDF_Mar, .(year, ring), function(x) colMeans(x[, chems]))


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
                   list(SoilDf, TcnDF_Ring, phDF_ring, SoilChemDF_ring, 
                        SoilMTdf_Ring, FlorLight_Ring, iem_ring, Extract_ring, 
                        Mineralisation_ring))
EnvVarDF <- within(EnvVarDF, {
  year <- factor(year)
  ring <- factor(ring)
  co2 <- factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
  TotalP <- NULL
})

save(EnvVarDF, file = "output/Data/FACE_EnvironmenVars.RData")
