#################
# Dissimilarity #
#################
# Compute dissimiliraity for each plot between 2012 and 2013

# compute dissimilarity for each plot
disDF <- ddply(plt.veg, .(ring, plot), function(x) vegdist(x[, -1:-3], method = "bray"))

# organise data frame
names(disDF)[3] <- "BC"
disDF <- within(disDF, {
  co2 <- factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
  block <- recode(ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
})

boxplot(BC ~ ring, data = disDF)

# perform LMM
m1 <- lmer(BC ~ co2 + (1|block) + (1|ring), data = disDF)
summary(m1)
Anova(m1)
Anova(m1, test.statistic = "F")
plot(m1)
# no co2 effect, but dissimlarity was slightly higher at eCO2

# fit soil moisture data

# load soil moisture data
load("Data/FACE_TDR_ProbeDF.RData")
summary(FACE_TDR_ProbeDF)
soilDF <- subsetD(FACE_TDR_ProbeDF, Sample == "vegetation" & 
                    Date <= as.Date("2013-12-31") &
                    Date >= as.Date("2013-1-1"))
soilDF$plot <- as.factor(soilDF$plot)
MoistDF <- ddply(soilDF, .(ring, plot), summarise, Moist = mean(Moist, na.rm = TRUE))
plot(Moist ~ ring, data = MoistDF)


# merge data frame
vsDF <- merge(disDF, MoistDF, by = c("ring", "plot"))

# fit to the model
plot(BC ~ Moist, vsDF, col = co2, pch = 16)
plot(log(BC) ~ Moist, vsDF, col = co2, pch = 16)

m1 <- lmer(log(BC) ~ co2 * Moist + (1|block), data = vsDF)
m2 <- lmer(log(BC) ~ co2 + Moist + (1|block), data = vsDF)
m3 <- lmer(log(BC) ~ co2 + (1|block), data = vsDF)
m4 <- lmer(log(BC) ~ Moist + (1|block), data = vsDF)
anova(m1, m2, m3, m4)
Anova(m1)
Anova(m1, test.statistic = "F")
plot(m1)
# moisture didn't really fit...
