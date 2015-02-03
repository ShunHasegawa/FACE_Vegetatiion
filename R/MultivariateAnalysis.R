plt.veg <- within(plt.veg, {
  co2 <- factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
  block <- recode(ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
})
Spp <- names(plt.veg)[!grepl("year|ring|plot|block|co2", names(plt.veg))]
Sites <- names(plt.veg)[grepl("year|ring|plot|block|co2", names(plt.veg))]
  

#################
# Dissimilarity #
#################
# Compute dissimiliraity for each plot between 2012 and 2013

# compute dissimilarity for each plot
disDF <- ddply(plt.veg, .(block, ring, co2, plot), 
               function(x) vegdist(x[, Spp], method = "altGower"))
names(disDF)[5] <- "BC"

boxplot(BC ~ ring, data = disDF)

# perform LMM
m1 <- lmer(BC ~ co2 + (1|block), data = disDF)
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
plot(log(BC) ~ log(Moist), vsDF, col = co2, pch = 16)

m1 <- lmer(log(BC) ~ co2 * log(Moist) + (1|block) + (1|ring), data = vsDF)
m2 <- lmer(log(BC) ~ co2 + log(Moist) + (1|block) + (1|ring), data = vsDF)
Anova(m2)
Anova(m2, test.statistic = "F")
summary(m2)
# no variation is explained by random factors
m3 <- lm(log(BC) ~ log(Moist), data = vsDF)
summary.lm(m3)
library(visreg)
visreg(m3, trans = exp)
# negative correlation with moisture

#######################
# Diversity & eveness #
#######################
vegDF <- plt.veg[, Spp]
siteDF <- plt.veg[, Sites]
siteDF$id <- siteDF$ring:siteDF$plot

DivDF <- within(siteDF,{
  H <- diversity(vegDF) # Shannon's index
  S <- specnumber(vegDF) # number of spp
  J <- H/log(S)  # Pielou's evenness
})

###############
## Diversity ##
###############
boxplot(H ~ co2*year, data = DivDF)
boxplot(H ~ year * ring, data = DivDF)
h1 <- lmer(H ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF)
Anova(h1)
Anova(h1, test.statistic = "F")
  # diversity slightly dicreased at eCO2

####################
## number of spp. ##
####################
boxplot(S ~ co2*year, data = DivDF)
boxplot(S ~ year * ring, data = DivDF)
boxplot(log(S) ~ year * ring, data = DivDF)

s1 <- glmer(S ~ co2 * year + (1|block) + (1|ring) + (1|id), family = poisson(), 
            data = DivDF)
s2 <- glmer(S ~ co2 + year + (1|block) + (1|ring) + (1|id), family = poisson(), 
            data = DivDF)
anova(s1, s2)
anova(s1)
anova(s2)
# number of spp decreased in thd 2nd year but no significant co2 effect

s3 <- lmer(log(S) ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF)
anova(s3)

#############
## eveness ##
#############
boxplot(J ~ co2*year, data = DivDF)
boxplot(J ~ year * ring, data = DivDF)

j1 <- lmer(J ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF)
j2 <- lmer(J ~ co2 + year + (1|block) + (1|ring) + (1|id), data = DivDF)
anova(j1, j2)
Anova(j2)
Anova(j2, test.statistic = "F")
plot(j1) # not very good..
# no significant co2 effect, but may decreased eveness slightly

plot(H ~ J, data = DivDF)
