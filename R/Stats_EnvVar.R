###############
# Dry soil pH #
###############
load("output/Data/FACE_DrysoilPhJune.RData")
head(phDF_June)
phDF_June <- within(phDF_June, {
                    Year = factor(year(Date)) 
                    Month = factor(month(Date))})
ftable(xtabs(~ Year + ring, data = phDF_June))

# ring variation
ggplot(phDF_June, aes(x = ring, y = Drysoil_ph)) + geom_boxplot() + facet_grid(.~year)
dlply(phDF_June, .(Year), function(x) anova(lm(Drysoil_ph ~ ring, data = x)))
  # significant ring variation

# Year variation
Ph_RingMean <- ddply(phDF_June, .(Year, ring), summarise, Drysoil_ph = mean(Drysoil_ph))
ggplot(Ph_RingMean, aes(x = Year, y = Drysoil_ph)) + geom_boxplot()
Anova(lmer(Drysoil_ph ~ Year + (1|ring), data = Ph_RingMean), test.statistic = "F")
  # substantial year variation

################
# Total Carbon #
################
load("output/Data/TcnJune_Plot.RData")
TcnDF_Plot <- within(TcnDF_Plot, {
  year <- factor(year)
  ring <- factor(ring)
  plot <- factor(plot)
})
ggplot(TcnDF_Plot, aes(x = ring, y = TotalC)) + geom_boxplot() + facet_grid(.~year)

# ring variation
dlply(TcnDF_Plot, .(year), function(x) anova(lm(TotalC ~ ring, data = x)))
# substntaial ring variation

# year variation
tcRingMean <- ddply(TcnDF_Plot, .(year, ring), summarise, TotalC = mean(TotalC))
ggplot(tcRingMean, aes(x = year, y = TotalC)) + geom_boxplot()
Anova(lmer(TotalC ~ year + (1|ring), data = tcRingMean), test.statistic = "F")


#########
# Light #
#########
load("output/Data/FACE_FloorPAR.RData")
head(Lightdf)

# subset November and December
LightNovDec <- subsetD(Lightdf, month(Date) %in% c(11, 12), select = c(-ProbMean, -DateTime))

# Daily Mean
LightNovDec_Daymean <- LightNovDec %>% 
  group_by(Date, ring) %>%
  summarise(Par1 = mean(PAR_Den_1_Avg, na.rm = TRUE), 
            Par2 = mean(PAR_Den_2_Avg, na.rm = TRUE),
            Par3 = mean(PAR_Den_3_Avg, na.rm = TRUE)            
  )
LightNovDec_Daymean_mlt <- melt(data.frame(LightNovDec_Daymean), 
                                id = c("Date", "ring"))
LightNovDec_Daymean_mlt <- within(LightNovDec_Daymean_mlt,{
  Month = factor(month(Date))
  Year = factor(year(Date))})

# inspect probes
str(LightNovDec_Daymean_mlt)
theme_set(theme_bw())
p <- ggplot(LightNovDec_Daymean_mlt, aes(x = Date, y = value, 
                                         col = variable, 
                                         group = variable))
p + geom_line() + facet_grid(ring ~ Year + Month, scale = "free_x")
# they look fine

# plot mean
Light_PltMean <- ddply(LightNovDec_Daymean_mlt, .(Year, ring, variable), summarise, 
                       FloorPAR = mean(value, na.rm = TRUE))
ftable(xtabs(~ Year + ring, data = Light_PltMean))
ggplot(Light_PltMean, aes(x = ring, y = FloorPAR)) + 
  geom_boxplot() + facet_grid(. ~ Year)

# Perform anova and see if there's ring variation
dlply(Light_PltMean, .(Year), function(x) {
  anova(lm(FloorPAR ~ ring, data = x))
})

# year variation
Light_RingMean <- ddply(Light_PltMean, .(Year, ring), summarise,
                        FloorPAR = mean(FloorPAR))
m1 <- lmer(FloorPAR ~ Year + (1|ring), data = Light_RingMean)
Anova(m1, test.statistic = "F")
# strong year variation
