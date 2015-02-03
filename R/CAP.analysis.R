library(BiodiversityR)

sites$co2 <- factor(ifelse(sites$ring %in% c(1, 4, 5), "elev", "amb"))
sites$YR <- sites$year:sites$ring

transDF <- decostand(vg.data, "log")
head(transDF)

adf <- log(vg.data + 1)

# Jeff's suggestion
# compute dissimilarity for each plot
disDF <- ddply(plt.veg, .(ring, plot), function(x) vegdist(x[, -1:-3], method = "bray"))


disDF <- ddply(plt.veg, .(ring, plot), 
               function(x) vegdist(log(x[, -1:-3] + 1), method = "bray")) # log(n + 1)


# disDF <- ddply(plt.veg, .(ring, plot), 
#                function(x) vegdist(x[, -1:-3], method = "binomial")) # log(n + 1)
# 
# 
# disDF <- ddply(plt.veg, .(ring, plot), 
#                function(x) vegdist(log(x[, -1:-3] + 1), method = "binomial")) # log(n + 1)
# 
# 
# disDF <- ddply(plt.veg, .(ring, plot), 
#                function(x) vegdist(log(x[, -1:-3] + 1), method = "altGower")) # log(n + 1)
# 
# 
# 
# disDF <- ddply(plt.veg, .(ring, plot), 
#                function(x) vegdist(decostand(x[, -1:-3], logbase = 2, "log"), method = "bray")) # log(n + 1)

names(disDF)[3] <- "BC"
disDF <- within(disDF, {
  co2 <- factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
  block <- recode(ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
})
disDF$co2 <- factor(ifelse(disDF$ring %in% c(1, 4, 5), "elev", "amb"))
boxplot(BC ~ ring, data = disDF)
boxplot(log(BC) ~ ring, data = disDF)

m1 <- lmer(log(BC) ~ co2 + (1|block) + (1|ring), data = disDF)
summary(m1)
Anova(m1)
Anova(m1, test.statistic = "F")
plot(m1)
# no co2 effect, but dissimlarity was slightly higher at eCO2, but moisture
# would probably fit quite well

# load soil moisture data
load("Data/FACE_TDR_ProbeDF.RData")
summary(FACE_TDR_ProbeDF)
soilDF <- subsetD(FACE_TDR_ProbeDF, Sample == "vegetation" & 
                    Date <= as.Date("2013-12-31") &
                    Date >= as.Date("2013-1-1"))
soilDF$plot <- as.factor(soilDF$plot)
MoistDF <- ddply(soilDF, .(ring, plot), summarise, Moist = mean(Moist, na.rm = TRUE))
plot(Moist ~ ring, data = MoistDF)

# fit to the model
vsDF <- merge(disDF, MoistDF, by = c("ring", "plot"))
plot(BC ~ Moist, vsDF)
plot(log(BC) ~ Moist, vsDF)

m1 <- lmer(BC ~ co2 * Moist + (1|block) + (1|ring), data = vsDF)
m2 <- lmer(BC ~ co2 + Moist + (1|block) + (1|ring), data = vsDF)
m3 <- lmer(BC ~ Moist + (1|block) + (1|ring), data = vsDF)
Anova(m3)
# moisture didn't really fit...


# transDF <- vegdist(vg.data, method = "altGower") # ln(x + 1)

capDF <- CAPdiscrim(transDF ~ YR, sites, dist = "bray", permutations = 10)

# create a plot
ldDF <- data.frame(sites, ld1 = capDF$x[, 1], ld2 = capDF$x[, 2])
chulDF <- ddply(ldDF, .(year, ring), 
                function(x) {chx <- chull(x[c("ld1", "ld2")]) 
                             chxDF <- data.frame(rbind(x[chx,], x[chx[1], ]))
                             return(chxDF)})

theme_set(theme_bw())
p <- ggplot(ldDF, aes(x = ld1, y = ld2, col = ring, shape = year))
p + geom_point(size = 5) + geom_polygon(data = chulDF, alpha = .1)

# eco2 vegeataion seems to shift towards left along ld1 in 2014
boxplot(ld1 ~ year*ring, data = ldDF)

# vegan
vegDF <- capscale(transDF ~ YR, sites, dist = "bray")
plot(vegDF)
summary(vegDF)
