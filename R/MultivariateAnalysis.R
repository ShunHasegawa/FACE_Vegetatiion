head(PlotSumVeg)

#################
# Dissimilarity #
#################
# Compute dissimiliraity for each plot between 2012 & 2014 and 2014 & 2015

tdf <- PlotSumVeg[1:2, ]

vegdist(log(tdf[, SppName] + 1), method = "bray")

# compute dissimilarity for each plot
disDF <- ddply(PlotSumVeg, .(block, ring, co2, plot), 
               function(x) {
                 df1 <- subset(x, year %in% c(2012, 2014))
                 dis1 <- vegdist(log(df1[, SppName] + 1), method = "bray")
                 df2 <- subset(x, year %in% c(2014, 2015))
                 dis2 <- vegdist(log(df2[, SppName] + 1), method = "bray")
                 dfs <- data.frame(year = c("Year1", "Year2"), BC = c(dis1, dis2))
                 return(dfs)
                 })

# perform LMM----
boxplot(BC ~ ring:year, data = disDF)
bxplts(value = "BC", xval = "co2", data = disDF)
bxplts(value = "BC", xval = "ring", data = disDF)
  # log seems better

disDF$id <- disDF$ring:disDF$plot
m1 <- lmer(log(BC) ~ co2 * year + (1|block) + (1|ring) + (1|id), data = disDF)
m2 <- lmer(log(BC) ~ co2 + year + (1|block) + (1|ring) + (1|id), data = disDF)
anova(m1, m2)
summary(m2)
Anova(m2)
Anova(m2, test.statistic = "F")
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))
# one obvious outlier. what if I remove
which(qqnorm(resid(m2))$x == max(qqnorm(resid(m2))$x))
m3 <- lmer(log(BC) ~ co2 * year + (1|block) + (1|ring) + (1|id), data = disDF[-32, ])
plot(m3)
qqnorm(resid(m3))
qqline(resid(m3))
# improved
summary(m3)
Anova(m3)
# result is not hugely different to the above. so just stay with the above

# fit soil moisture data---

# load soil moisture data
load("Data/FACE_TDR_ProbeDF.RData")
summary(FACE_TDR_ProbeDF)
soilDF <- subsetD(FACE_TDR_ProbeDF, Sample == "vegetation" & 
                    Date >= as.Date("2013-1-1") &
                    Date <= as.Date("2014-12-31"))

soilDF <- within(soilDF, {
  plot <- factor(plot)
  year <- factor(ifelse(year(Date) == "2013", "Year1", "Year2"))
})

MoistDF <- ddply(soilDF, .(year, ring, plot), summarise, Moist = mean(Moist, na.rm = TRUE))
boxplot(Moist ~ year:ring, data = MoistDF)

# merge data frame
vsDF <- merge(disDF, MoistDF, by = c("ring", "plot", "year"))

# fit to the model
theme_set(theme_bw())
p <- ggplot(data = vsDF, aes(x = log(Moist), y = BC, col = co2, shape = year))
p <- ggplot(data = vsDF, aes(x = log(Moist), y = BC, col = co2))
p <- ggplot(data = vsDF, aes(x = Moist, y = BC, col = co2))
p2 <- p + geom_point(size = 4) + geom_smooth(aes(fill = co2), method  = "lm") + 
  facet_grid(.~block)
p2




m1 <- lmer(BC ~ year * co2 * log(Moist) + (1|block) + (1|ring) + (1|id), data = vsDF)
mm2 <- lmer(BC ~ (year + co2 + log(Moist))^2 + (1|block) + (1|ring) + (1|id), data = vsDF)
mm_1 <- lmer(BC ~ year + co2 + Moist + I(Moist^2) + (1|block) + (1|ring) + (1|id), data = vsDF)
mm_2 <- lmer(BC ~ year + co2 + Moist + (1|block) + (1|ring) + (1|id), data = vsDF)
mm_3 <- lmer(BC ~ year + co2 + (1|block) + (1|ring) + (1|id), data = vsDF)
Anova(mm_1)
Anova(mm_2)
Anova(mm4)
AIC(mm_1, mm4, mm_2, mm_3)



AIC(mm4)
drop1(mm2)
mm3 <- update(mm2, ~.-year:co2)
drop1(mm3)
mm4 <- update(mm3, ~.-co2:log(Moist))
drop1(mm4)
Anova(mm4, test.statistic = "F")
plot(allEffects(mm4))
plot(mm4)



anova(m1, mm2)
Anova(mm2)


m2 <- lmer(log(BC) ~ co2 + log(Moist) + (1|block) + (1|ring) + (1|id), data = vsDF)
anova(m1, m2)
Anova(m2)
Anova(m2, test.statistic = "F")
summary(m2)
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))
plot(allEffects(m2))
library(visreg)
visreg(m2)


#############
# PERMANOVA #
#############
# perform permanova for each year separately

# 2012
df12 <- subsetD(plt.veg, year == 2012)
veg12 <- df12[, Spp]
site12 <- df12[, Sites]

a1 <- adonis(veg12 ~ co2, data = site12, strata = site12$ring, 
             perm=100, method = "altGower")
a1
a3 <- adonis(veg12 ~ co2, data = site12, strata = site12$block, 
             perm = 100, method = "altGower")


n1 <- nested.npmanova(veg12 ~ co2 + ring, data = site12, perm=100, method = "altGower")
n2 <- nested.npmanova(veg12 ~ co2 + block, data = site12, perm=100, method = "altGower")
n1
n2

# 2014
df14 <- subsetD(plt.veg, year == 2014)
veg14 <- df14[, Spp]
site14 <- df14[, Sites]

a1 <- adonis(veg14 ~ co2 + ring, data = site14, strata = site14$ring, 
             perm=10, method = "altGower")


a1 <- adonis(veg14 ~ co2, data = site14, perm=100, method = "altGower")
a1
n1 <- nested.npmanova(veg14 ~ co2 + ring, data = site14, perm=5000, method = "altGower")
n1



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



# cap without pooling sample
PFGVeg <- dcast(year + ring + plot ~ PFG, sum, data = PFGdf)
PFGVeg <- PFGVeg[, c("year", "ring", "plot", pfgs)]
PFGVeg$co2 <- factor(ifelse(PFGVeg$ring %in% c(1, 4, 5), "elev", "amb"))
PFGa <- subsetD(PFGVeg, co2 == "amb")
PFGe <- subsetD(PFGVeg, co2 == "elev")
pfgS <- PFGe[,c("year", "ring", "plot")]


summary(pfgS)
summary(df)
df <- data.frame(PFGe[, pfgs])
str(df)
capDF <- CAPdiscrim(df ~ year, pfgS, dist = "altGower", permutations = 10)
?CAPdiscrim

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

