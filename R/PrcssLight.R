#########
# Light #
#########

# download from HIEv
setToken(tokenfile = "Data/token.txt")

###############
# Understorey #
###############

# download files from HIEv
Light_Hiev <- downloadTOA5("FACE.*AirVars.*dat", 
                           cachefile = "Data/hievdata/tmp.RData",
                           topath = "Data/hievdata/raw_data",
                           maxnfiles = 999)
  # LI190SB_PAR_Den_Avg at 23 m
  # PAR_Den_1-3_Avg at 0.15 m
  # Look at
  # https://sites.google.com/site/hievuws/facilities/eucface/collection-codes for
  # the description of each variables

# save(Light_Hiev, file = "output/Data/FACE_AirVars_Raw.RData")
# load("output/Data/FACE_AirVars_Raw.RData")

# organise data frame. Note that this is really huge so use dplyr----

# Add ring number and remove Source
Light_Hiev <- select(mutate(Light_Hiev, 
                            ring = factor(substr(Light_Hiev$Source, 7, 7))), 
                     -Source)
# remove duplicates
Light_Hiev <- distinct(Light_Hiev)

# save
save(Light_Hiev, file = "output/Data/FACE_AirVars_Processed.RData")

# select only required columns as it's too large
Lightdf <- select(Light_Hiev, DateTime, Date, ring, PAR_Den_1_Avg, PAR_Den_2_Avg, PAR_Den_3_Avg)

# Daily mean
Lightdf_DayMean <- Lightdf %>%
  group_by(Date, ring) %>% 
  summarise(PAR_Den_1_Avg = mean(PAR_Den_1_Avg, na.rm = TRUE), 
            PAR_Den_2_Avg = mean(PAR_Den_2_Avg, na.rm = TRUE), 
            PAR_Den_3_Avg = mean(PAR_Den_3_Avg, na.rm = TRUE) )

# inspection----
Lightdf_DayMean <- data.frame(Lightdf_DayMean)
tdf <- melt(Lightdf_DayMean, id = c("Date", "ring"))
summary(tdf)
p <- ggplot(tdf, aes(x = Date, y = value, col = variable))
p2 <- p + geom_point(size = .05) + facet_grid(ring ~ .)
p2
# some weird values in the beginning
tdf12 <- subset(tdf, Date < as.Date("2012-08-5"))
p <- ggplot(tdf12, aes(x = Date, y = value, col = variable))
p2 <- p + geom_point() + facet_grid(ring ~ .) + 
  scale_x_date(minor_breaks = "days", 
               labels = date_format("%b-%d"))
p2
# don't use before August 2012

tdf <- subset(tdf, Date > as.Date("2012-08-1"))
p <- ggplot(tdf, aes(x = Date, y = value, col = variable))
p2 <- p + geom_point() + facet_grid(ring ~ .) + 
  scale_x_date(breaks = "2 week",
               minor_breaks = "days", 
               labels = date_format("%b-%d")) +
  theme(axis.text.x = element_text(angle = 90))
p2
# Oct2012. Light intensity goes up all of sudden
tdf <- subset(tdf, Date > as.Date("2012-08-1") & Date < as.Date("2012-11-1"))
p <- ggplot(tdf, aes(x = Date, y = value, col = variable))
p2 <- p + geom_point() + facet_grid(ring ~ .) + 
  scale_x_date(breaks = "2 week",
               minor_breaks = "days", 
               labels = date_format("%b-%d")) +
  theme(axis.text.x = element_text(angle = 90))
p2

# what about raw data?
tdf_Oct2012 <- filter(Lightdf, Date > as.Date("2012-10-10") & Date < as.Date("2012-10-20"))
tdf_Oct2012_mlt <- melt(data.frame(tdf_Oct2012), id = c("DateTime", "Date", "ring"))
p <- ggplot(tdf_Oct2012_mlt, aes(x = DateTime, y = value, col = variable))
p2 <- p + geom_point(size = .05) + facet_grid(ring + variable ~ .)
p2
# when I look at the raw data, it's not too weird actually so just stay with
# this. Just remove observations before August 2012

Lightdf <- filter(Lightdf, Date >= as.Date("2012-08-01"))
save(Lightdf, file = "output/Data/FACE_FloorPAR.RData")

# Mean of the all probes for each time point
Lightdf$ProbMean <- rowMeans(Lightdf[, c("PAR_Den_1_Avg", "PAR_Den_2_Avg", "PAR_Den_3_Avg")], na.rm = TRUE)

# Daily Mean
UnderstryPAR <- Lightdf %>% 
  group_by(Date, ring) %>%
  summarise(FloorPAR = mean(ProbMean, na.rm = TRUE))

## anova between ring ##

# load("output/Data/FACE_FloorPAR.RData")
head(Lightdf)
# yearly mean for each ring (use only November and December as these months are
# used in RDA)
ndDf <- subset(Lightdf, month(Date) %in% c(11, 12))
ndDf$year <- factor(year(ndDf$Date))
some(ndDf)

YR_DF <- ndDf %>%
  group_by(year, ring) %>% 
  summarise(PAR_Den_1_Avg = mean(PAR_Den_1_Avg, na.rm = TRUE), 
            PAR_Den_2_Avg = mean(PAR_Den_2_Avg, na.rm = TRUE), 
            PAR_Den_3_Avg = mean(PAR_Den_3_Avg, na.rm = TRUE) )
YR_DF_mlt <- melt(data.frame(YR_DF), id = c("year", "ring"))
boxplot(value ~ ring * year, data = YR_DF_mlt)
dlply(YR_DF_mlt, .(year), function(x) {
  m <- lm(value ~ ring, data = x)
  anova(m)})
# no significant spatial variation

##########################
# canopy light intensity #
##########################

# download files from HIEv
FACE_general_raw <- downloadTOA5("FACE.*general.*dat", 
                                 cachefile = "Data/hievdata/tmp_FACEGeneral.RData",
                                 topath = "Data/hievdata/raw_data//FACE_general",
                                 maxnfiles = 999)
save(FACE_general_raw, file = "output/Data/FACE_general_raw.RData")

# Add ring number and remove duplicates
FACE_general_raw <- select(mutate(FACE_general_raw,
                                  ring = factor(substr(Source, 7, 7))),
               -Source)
FACE_general_RmDup <- distinct(FACE_general_raw)
save(FACE_general_RmDup, file = "output//Data/FACE_general_Processed.RData")

# subset required columns
CnpyLightDF <- select(FACE_general_RmDup, 
                      DateTime, Date, ring, LI190SB_PAR_Den_Avg)
# daily mean
CnpyLight_DayMean <- CnpyLightDF %>% 
  group_by(Date, ring) %>% 
  summarise(CanopyPAR = mean(LI190SB_PAR_Den_Avg, na.rm = TRUE))

tdf <- within(CnpyLight_DayMean, {
  Y = factor(year(Date)) 
  M = factor(month(Date))
}) # for some reasons the following doesn't work,,, so use within
# tdf <- mutate(CnpyLight_DayMean, 
#               Y = factor(year(Date)), 
#               M = factor(month(Date)))

p <- ggplot(data = tdf, aes(x = Date, y = CanopyPAR))
p2 <- p + geom_point(size = .05) + facet_grid(ring  ~ .)
p2
# no obvious problem so just use thsese. 
# BUT NOTE that when I plot raw obsevations I found that there were some dates
# (in Oct, Nov 2013) with no observations.

######################################
# Combine understorey and canopy Par #
######################################
head(UnderstryPAR)
head(CnpyLight_DayMean)
FACE_Light_DayMean <- merge(UnderstryPAR, CnpyLight_DayMean, by = c("Date", "ring"), all = TRUE)
save(FACE_Light_DayMean, file = "output//Data/FACE_Light_DayMean.RData")

load("output//Data/FACE_Light_DayMean.RData")

tdf <- FACE_Light_DayMean
tdf <- within(tdf, {
  M <- month(Date, label = TRUE, abbr = TRUE)
  Y <- year(Date)
})

tdf_mlt <- melt(tdf, id = c("Date", "Y", "M", "ring"))
tdf_mlt_under <- subsetD(tdf_mlt, variable == "FloorPAR")

p <- ggplot(tdf_mlt_under, aes(x = mday(Date), y = value, col = ring))
p2 <- p + 
  geom_line(size = .2, alpha = .7) + 
  facet_grid(Y ~ M, scales = "free_x") +
  scale_x_continuous(breaks = c(10, 20, 30)) +
  scale_color_manual(values = palette()) +
  labs(y = expression(Undersotrey~PAR~(mu*mol~s^"-1"~m^"-2")), x = "Date")
ggsave(p2, filename = "output//figs/FACE_RingDailyUnderstoreyPAR.pdf", width = 10, height = 6)
