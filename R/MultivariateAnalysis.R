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

###########################
# Plant Functional Groups #
###########################
load("output/Data/FACE_Vegetation_PFG.RData")

# remove rows whicn contain NA in PFG
PFGdf <- FACE.veg.rslt[!is.na(FACE.veg.rslt$PFG),]
library(reshape2)
PFGsum <- dcast(year + ring + plot ~ PFG, sum, data = PFGdf)

# just use c3, c4, legume and non-legume
PFGsum <- PFGsum[, c("year", "ring", "plot","c3", "c4", "legume", "Non_legume")]

# organise
PFGsum <- within(PFGsum, {
  block = recode(ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
  co2 = factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
})

# Compute dissimiliraity for each plot between 2012 and 2013

# compute dissimilarity for each plot
pfgs <- c("c3", "c4", "legume", "Non_legume")
sites <- c("year", "block", "ring", "plot", "co2")

disDF <- ddply(PFGsum, .(block, ring, co2, plot), 
               function(x) vegdist(x[, pfgs], method = "altGower"))
names(disDF)[5] <- "BC"

boxplot(BC ~ ring, data = disDF)
boxplot(BC ~ co2, data = disDF)

# perform LMM
m1 <- lmer(BC ~ co2 + (1|block), data = disDF)
summary(m1)
Anova(m1)
Anova(m1, test.statistic = "F")
plot(m1)
# marginal co2 effect, eCO2 slightly decreased dissimilarity in PFG

# CAP for assembled sample
# ring sum
RngPFGSum <- ddply(PFGsum, .(year, co2, ring), 
                   function(x) colSums(x[, pfgs]))
siteDF <- RngPFGSum[, c("year", "co2", "ring")]
siteDF$YC <- siteDF$year:siteDF$co2
siteDF$YR <- siteDF$year:siteDF$ring

pfgDF <- RngPFGSum[, pfgs]
transDF <- decostand(pfgDF, "log")

capDF <- CAPdiscrim(pfgDF ~ YC, siteDF, dist = "altGower", permutations = 10)

pfg14 <- subsetD(RngPFGSum, year == 2012)
pfg14df <- pfg14[, pfgs]
site14 <- pfg14[, c("year", "co2", "ring")]

capDF <- CAPdiscrim(pfg14df ~ co2, site14, dist = "altGower", permutations = 10)
boxplot(capDF$x[, 1] ~ site14$co2)


# create a plot
boxplot(capDF$x[, 1] ~ siteDF$YC)

ldDF <- data.frame(siteDF, ld1 = capDF$x[, 1], ld2 = capDF$x[, 2])
chulDF <- ddply(ldDF, .(year, co2), 
                function(x) {chx <- chull(x[c("ld1", "ld2")]) 
                             chxDF <- data.frame(rbind(x[chx,], x[chx[1], ]))
                             return(chxDF)})

theme_set(theme_bw())
p <- ggplot(ldDF, aes(x = ld1, y = ld2, col = co2, shape = year))
p + geom_point(size = 5) + geom_polygon(data = chulDF, alpha = .1)

##################################
# assembled sample for each ring #
##################################
vegDF <-  plt.veg[, Spp]
siteDF <-  plt.veg[, Sites]

RngVeg <- ddply(plt.veg, .(year, block, co2, ring), function(x) colSums(x[, Spp]))

# CAP
RngVegdf <- RngVeg[, Spp]
RngVegdf_site <- RngVeg[, c("year", "ring", "block", "co2")]
RngVegdf_site$YC <- RngVegdf_site$year:RngVegdf_site$co2

capDF <- CAPdiscrim(RngVegdf ~ YC, RngVegdf_site, dist = "altGower", permutations = 10)

# create a plot
ldDF <- data.frame(RngVegdf_site, ld1 = capDF$x[, 1], ld2 = capDF$x[, 2])

chulDF <- ddply(ldDF, .(year, co2), 
                function(x) {chx <- chull(x[c("ld1", "ld2")]) 
                             chxDF <- data.frame(rbind(x[chx,], x[chx[1], ]))
                             return(chxDF)})

theme_set(theme_bw())
p <- ggplot(ldDF, aes(x = ld1, y = ld2, col = co2, shape = year))
p + geom_point(size = 5) + geom_polygon(data = chulDF, alpha = .1)

