##########################
# Organise soil variable #
##########################
SoilDf <- read.csv("Data//FACE_Environmtal_Variable.csv", header = TRUE)

#############
# TC and TN #
#############

# 2012-2014
  TcnDF <- read.xlsx2("Data/SoilVariables/Soil.TC.TN.xlsx", 
                      sheetName = "FACE.soil.tc.tn", 
                      header = TRUE, startRow = 1, endRow = 97, 
                      stringsAsFactors = FALSE)
  # organise
  names(TcnDF)[grepl("tc|tn", names(TcnDF))] <- c("TotalC", "TotalN")
  TcnDF[, c("TotalC", "TotalN")] <- apply(TcnDF[, c("TotalC", "TotalN")], 2, as.numeric)
  TcnDF$year <- factor(TcnDF$year, labels = paste0("Year", 0:1))
  TcnDF$depth <- factor(paste0(TcnDF$depth, "cm"))

# 2015
  TcnDF_2015 <- read.xlsx2("Data/SoilVariables/TotalCN_June2015.xlsx",
                           sheetIndex = 1, 
                           colClasses = c("character", rep("numeric", 3)))
  names(TcnDF_2015) <- c("ID", "Mass", "TotalC", "TotalN")
  
  # Subset required rows
  TcnDF_2015$ID <- as.numeric(as.character(TcnDF_2015$ID))
  TcnDF_2015 <- TcnDF_2015[complete.cases(TcnDF_2015), ]
  
  # remove 5 from ID (it was meant to be S but Simmy got it wrong)
  TcnDF_2015$ID <- TcnDF_2015$ID - 500000
  
  # Sample ID
  TcnDF_2015_ID <- read.xlsx2("Data/SoilVariables/TotalCN_June2015_SampleID.xlsx",
                              sheetIndex = 1,
                              startRow = 2,
                              colClasses = c("numeric", rep("character", 2))
                              )
  
  # Merge
  TcnDF_2015dd <- merge(TcnDF_2015, TcnDF_2015_ID, by.x = "ID", by.y = "Sample")
  
  # Organise
  TcnDF_2015dd <- subset(TcnDF_2015dd, 
                         select = c("TotalC", "TotalN", "Ring", "Plot"),
                         !grepl("Backfill", as.character(Ring))
                         )
  names(TcnDF_2015dd)[3:4] <- c("ring", "plot") 
  TcnDF_2015dd <- within(TcnDF_2015dd, {
    year  <- "Year3"
    depth <- "0-10cm"
  })

# Merge all year data
  TcnDF <- rbind.fill(TcnDF, TcnDF_2015dd)
    
#######################
# TC and TN from Hiev #
#######################
HtcnDF <- read.csv("Data/SoilVariables/FACE_P0014_ALL_ SoilCN_June2012-Spet2014_V1.csv")
HtcnDF <- within(HtcnDF, {
  Date <- as.Date(dmy(as.character(Date)))
  ring <- factor(ring)
  year <- factor(year(Date), labels = paste0("Year", c(0, 2))) 
  Month <- month(Date)
})
HtcnDF <- subsetD(HtcnDF, plot != "Backfill")
names(HtcnDF)[grepl("tc|tn", names(HtcnDF))] <- c("TotalC", "TotalN")
ftable(xtabs( ~ year + Month + depth, data = HtcnDF))
# Year1 is completely missing
# no 20-30 cm for Year0

# compre the above two data frames
llply(list(HtcnDF, TcnDF), function(x) subset(x, ring == 1 & year == "Year0"))
  # Year0 is the same. so just need Year1 and Year3 from the above.

# ring mean for each data frames at 0-10 cm from June or May
TcnDF_Ring <- ldply(list(subset(TcnDF, year %in% c("Year1", "Year3")), subset(HtcnDF, Month %in% c(5, 6))), 
                    function(x) {
                      ddply(x, .(year, ring, depth), summarise, 
                            TotalC = mean(TotalC), 
                            TotalN = mean(TotalN))
                      })
# plot mean
dd14 <- subset(TcnDF, year == "Year1" & depth == "0-10cm", 
               select = c("year", "ring", "plot", "TotalC", "TotalN"))
dd13_15 <- subset(HtcnDF, Month %in% c(5, 6) & depth == "0-10cm", 
                  select = c("year", "ring", "plot", "TotalC", "TotalN"))
# fix plot labels
dd13_15$plot <- factor(sub(".[.]", "", as.character(dd13_15$plot)))
TcnDF_Plot <- rbind.fill(list(dd14, dd13_15))
TcnDF_Plot <- TcnDF_Plot[order(TcnDF_Plot$year), ]
save(TcnDF_Plot, file = "output/Data/TcnJune_Plot.RData")

# check interaction beween ring and depth
par(mfrow = c(2, 4))
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
  Date  <- as.Date(dmy(Date))
  Ring  <- factor(Ring)
  Plot  <- as.character(factor(Plot)) # 1.0->1, 1.1->1.1
  Plot  <- factor(gsub(".[.]", "", Plot)) # 1->1, 1.1->1
  Depth <- factor(gsub(" ", "", as.character(Depth)))
  year  <- factor(year(Date), labels = paste0("Year", 0:2)) 
  Month <- month(Date)
  })

names(phDF)[c(3, 4, 6)] <- c("ring", "plot", "Drysoil_ph")
# remove NA
phDF <- phDF[!is.na(phDF$Drysoil_ph), ]

# which month and depth to be used
ftable(xtabs( ~ year + Month + Depth, data = phDF))
  # use May/June at 0-10cm

phDF_June <- subsetD(phDF, Depth == "0-10cm" & Month %in% 5:6)

# pH data from June 2015
ph_2016 <- read.xlsx2(file = "Data/Soil pH June 2015.xlsx", 
                      sheetIndex = 1,
                      colClasses = c("numeric", "character", "character", "numeric"))
names(ph_2016) <- c("Sample.ID", "ring", "plot", "Drysoil_ph") 

ph_2016 <- within(ph_2016, {
  Depth <- "0-10cm"
  Month <- 6
  year  <- "Year3"
})

phDF_June <- rbind.fill(phDF_June, ph_2016)
save(phDF_June, file = "output/Data/FACE_DrysoilPhJune.RData")

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
  year <- factor(year(Date), labels = paste0("Year", 0:2))
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
                    month(Date) %in% c(8:12) & year(Date) %in% c(2012:2015),
                    select = -co2)
SoilMTdf$year <- factor(year(SoilMTdf$Date), labels = paste0("Year", 0:3))
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
  # October 2012 values are little bit weird. It goes up in a sudden.
  # Just use Nov-Dec for the time being cause measurements in these months are
  # complete in all the three years.
LightNovDec <- subsetD(FACE_Light_DayMean, month(Date) %in% c(11, 12))
LightNovDec$year <- factor(year(LightNovDec$Date), labels = paste0("Year", 0:3))

FlorLight_Ring <- ddply(LightNovDec, .(year, ring), summarise, 
                        FloorPAR = mean(FloorPAR, na.rm = TRUE))

#######
# IEM #
#######
load("Data/SoilVariables/FACE_IEM.RData")

# use only Nov-Jan
iemNovJan <- subsetD(iem, time %in% c(6, 7, 13, 14))
# add year; January is counted as an year before
iemNovJan$year <- with(iemNovJan, factor(ifelse(time %in% c(6, 7), "Year0", "Year1")))

iem_ring <- ddply(iemNovJan, .(ring, year), summarise,
                    IEM_no = mean(no, na.rm = TRUE), 
                    IEM_nh = mean(nh, na.rm = TRUE), 
                    IEM_p  = mean(p, na.rm = TRUE)
                  )

#################
# Soil Extracts #
#################
load("Data/SoilVariables/extractable.RData")
tdf <- within(extr, {
  M <- month(date)
  Y <- year(date)
})
xtabs(~ Y + M, data = tdf)

# use December
ExtrDec <- subsetD(extr, month(date) == 12)
ExtrDec$year <- factor(year(ExtrDec$date), labels = paste0("Year", 0:1))
Extract_ring <- ddply(ExtrDec, .(year, ring), summarise, 
                      Ext_no = mean(no, na.rm = TRUE),
                      Ext_nh = mean(nh, na.rm = TRUE),
                      Ext_p  = mean(po, na.rm = TRUE)
                      )

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
MineJan$year <- factor(year(MineJan$date), labels = paste0("Year", 0:1))

Mineralisation_ring <- ddply(MineJan, .(year, ring), summarise,
                               n.min         = mean(n.min, na.rm = TRUE), 
                               nitrification = mean(nitrification, na.rm = TRUE), 
                               p.min         = mean(p.min, na.rm = TRUE)
                             )
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

#######################################################
# Sumamry and simple stats on environmental variables #
#######################################################
# only use totalC, moist, ph, depth of HL, Par and temp
names(EnvVarDF)

envDF <- EnvVarDF[, c("year", "ring", "co2","TotalC", 
                      "moist", "Drysoil_ph", "Depth_HL", 
                      "FloorPAR", "temp")]
envDF$moist <- envDF$moist * 100
# Treatment mean and SE
envDF_mlt <- melt(envDF, id = c("year", "ring", "co2"))
TreatSummary <- ddply(envDF_mlt, .(year, co2, variable), summarise,  
                      Mean = mean(value), 
                      SE = ci(value)[4], 
                      N = sum(!is.na(value)))
TreatSummary$value <- with(TreatSummary, paste0(round(Mean, 2), 
                                                "(",
                                                round(SE, 2),
                                                ")"))
TreatSummary_cst <- dcast(variable ~ year + co2, data = TreatSummary)

#########
# Stats #
#########

# totalC----
bxplts(value = "TotalC", xval = "co2", data = envDF)
m1 <- lmer(TotalC ~ co2 * year + (1|ring), data = envDF)
m2 <- stepLmer(m1)
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))
AnvF_totalC <- Anova(m2, test.statistic = "F")
AnvF_totalC

# moist----
bxplts(value = "moist", xval = "co2", data = envDF)
m1 <- lmer(moist ~ co2 * year + (1|ring), data = envDF)
qqnorm(resid(m2))
qqline(resid(m2))
AnvF_moist <- Anova(m1, test.statistic = "F")
AnvF_moist

# ph----
bxplts(value = "Drysoil_ph", xval = "co2", data = envDF)
m1 <- lmer(Drysoil_ph ~ co2 * year + (1|ring), data = envDF)
m2 <- stepLmer(m1)
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))
AnvF_ph <- Anova(m2, test.statistic = "F")
AnvF_ph

# dpeth HL-----
bxplts(value = "Depth_HL", xval = "co2", data = envDF)
m1 <- lm(Depth_HL ~ co2, data = envDF, subset = year == "Year0")
AnvF_HL <- summary.aov(m1)
AnvF_HL

# PAR----
bxplts(value = "FloorPAR", xval = "co2", data = envDF)
m1 <- lmer(FloorPAR ~ co2 * year + (1|ring), data = envDF)
m2 <- stepLmer(m1)
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))
AnvF_par <- Anova(m2, test.statistic = "F")
AnvF_par

# temp----
bxplts(value = "temp", xval = "co2", data = envDF)
m1 <- lmer(temp ~ co2 * year + (1|ring), data = envDF)
m2 <- stepLmer(m1)
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))
AnvF_temp <- Anova(m2, test.statistic = "F")
AnvF_temp

# Summary----
AnvF_Env <- ldply(list(TotalC = AnvF_totalC, moist = AnvF_moist,
                       temp = AnvF_temp, 
                       Drysoil_ph = AnvF_ph, FloorPAR = AnvF_par),
                  function(x) data.frame(x, term = row.names(x)),
                  .id = "variable")
AnvF_Env_p <- AnvF_Env[, c("variable", "term", "Pr..F.")]
names(AnvF_Env_p)[3] <- "Pr"
AnvF_Env_p <- rbind(AnvF_Env_p, data.frame(variable = "Depth_HL", 
                                           term = "co2",
                                           Pr = AnvF_HL[[1]]$Pr[1]))

AnvF_Env_p$stats <- cut(AnvF_Env_p$Pr, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                        labels = c("***", "**", "*", ".", "ns"))


AnvF_Env_p_cst <- dcast(variable ~ term, value.var = "stats", data = AnvF_Env_p)
AnvF_Env_p_cst[is.na(AnvF_Env_p_cst)] <- "ns"

# combine treatment summary and stats
Env_summary <- merge(TreatSummary_cst, AnvF_Env_p_cst, by = "variable")
Env_summary <- Env_summary[match(c("TotalC", "moist", "Drysoil_ph", 
                                   "Depth_HL", "FloorPAR", "temp"), 
                                 Env_summary$variable), ]
write.csv(Env_summary, file = "output/table/FACE_EnvVarSummary.csv", row.names = FALSE)

##################################
# Extractable NO vs. Dry soil pH #
##################################
extJun <- subset(extr, month(date) == 6)
extJun$year <- factor(year(extJun$date)+1)
phDF_June$year <- factor(phDF_June$year)
llply(list(extJun, phDF_June), nrow)
dd <- merge(extJun, phDF_June, all.x = TRUE, by = c("year", "ring", "plot"))

# get ring mean and remove pseudoreplication
dd_Mean <- ddply(dd, .(ring), function(x) colMeans(x[, c("Drysoil_ph", "no")]))
plot(Drysoil_ph ~ exp(no), data = dd_Mean)
m <- lm(Drysoil_ph ~ exp(no), data = dd_Mean)
summary(m)
anova(m)
abline(m)
plot(Drysoil_ph ~ no, data = dd_Mean)
xval <- seq(.5, 3.5, length.out = 100)
yval <- predict(m, data.frame(no = xval))
lines(xval, yval)
