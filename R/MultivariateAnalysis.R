head(PlotSumVeg)

#################
# Dissimilarity #
#################

# compute dissimilarity for each plot
disDF <- ddply(PlotSumVeg, .(block, ring, co2, plot), function(x) YearDssmlrty(x, spp = SppName))

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

MoistTempDF <- ddply(soilDF, .(year, ring, plot), summarise, 
                     Moist = mean(Moist, na.rm = TRUE),
                     Temp = mean(Temp_Mean, na.rm = TRUE))
boxplot(Moist ~ year:ring, data = MoistTempDF)

# merge data frame
vsDF <- merge(disDF, MoistTempDF, by = c("ring", "plot", "year"))
vsDF$id <- with(vsDF, ring:plot)

# fit to the model
theme_set(theme_bw())

# against Moist----
# raw data
p <- ggplot(data = vsDF, aes(x = Moist, y = BC, col = co2, shape = year))
p2 <- p + geom_point(size = 4) + geom_smooth(aes(fill = co2), method  = "lm") 
p2
# log
p <- ggplot(data = vsDF, aes(x = log(Moist), y = log(BC), col = co2, shape = year))
p2 <- p + geom_point(size = 4) + geom_smooth(aes(fill = co2), method  = "lm") 
p2

# against Temp----
# raw
p <- ggplot(data = vsDF, aes(x = Temp, y = BC, col = co2, shape = year))
p2 <- p + geom_point(size = 4) + geom_smooth(aes(fill = co2), method  = "lm") 
p2
# log
p <- ggplot(data = vsDF, aes(x = log(Temp), y = log(BC), col = co2, shape = year))
p2 <- p + geom_point(size = 4) + geom_smooth(aes(fill = co2), method  = "lm") 
p2
# not much difference in temperature between years

# There's small indication of positive correlation in Year2 adn negative
# correlation in Year1.
# Analyse each year separately----
vsDF$logMoist <- log(vsDF$Moist)
mls <- dlply(vsDF, .(year), function(x) {
  tdf <- x
  lmer(log(BC) ~ co2 + logMoist + (1|block) + (1|ring), data = tdf)
})
llply(mls, Anova)
# moisture effects are not significant when analysed separately, so do not worry
# about yearly difference in correlation for the time being

################
## Fit models ##
################
# log----
m1 <- lmer(log(BC) ~ co2 * logMoist + (1|block) + (1|ring) + (1|id), data = vsDF)
Anova(m1)
m2 <- update(m1, ~. - co2:logMoist)
Anova(m2)
Anova(m2, test.statistic = "F")

# 2nd polynomial----
m3 <- lmer(BC ~ co2 * (Moist + I(Moist^2)) + (1|block) + (1|ring) + (1|id), data = vsDF)
Anova(m3)
m4 <- lmer(BC ~ co2 + Moist + I(Moist^2) + (1|block) + (1|ring) + (1|id), data = vsDF)
Anova(m4)
Anova(m4, test.statistic = "F")

# compare two models
ListMl <- list(m2, m4)
names(ListMl) <- c("log", "2nd Polynomial")

# model diagnosis
plot(m2)
plot(m4)
par(mfrow = c(1, 2))
l_ply(names(ListMl), function(x) {
  l <- ListMl[[x]]
  qqnorm(resid(l), main = x)
  qqline(resid(l))
})
ldply(ListMl, r.squared)

visreg(m2, xvar = "logMoist", xtrans = exp, trans = exp, by = "co2", overlay = TRUE)
visreg(m4, xvar = "Moist", by = "co2", overlay = TRUE)
# no huge difference but log seems slightly better so use log


###########################
# Plant Functional Groups #
###########################
head(PlotSumVeg)
PFGName

# compute dissimilarity for each plot

disDF <- ddply(PlotSumPFGMatrix, .(block, ring, co2, plot), 
               function(x) YearDssmlrty(x, spp = PFGName))
disDF$id <- with(disDF, ring:plot)

# combine with moisture data
summary(MoistTempDF)
PFG_vsDF <- merge(disDF, MoistTempDF, by = c("ring", "plot", "year"))

theme_set(theme_bw())
p <- ggplot(data = PFG_vsDF, aes(x = log(Moist), y = log(BC), col = co2, shape = year))
p2 <- p + geom_point(size = 4) + geom_smooth(method = "lm")
p2
# it seems that BC goes up till moisture of 0.1, then goes down. may be
# secon-plynomial might be appropriate
p <- ggplot(data = PFG_vsDF, aes(x = Moist, y = BC, col = co2))
p2 <- p + geom_point(size = 4) + 
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE))
p2

############################
## what if I include year ##
############################
PFG_vsDF$logMoist <- log(PFG_vsDF$Moist)
m1 <- lmer(log(BC) ~ co2 * year * logMoist + (1|block) + (1|ring) + (1|id), data = PFG_vsDF)
m2 <- stepLmer(m1)
Anova(m2)
Anova(m2, test.statistic = "F")
plot(allEffects(m2))
# no moisture effect. Large proportion of variability is explained by moisture

# each year----
Mlist <- dlply(PFG_vsDF, .(year), function(x) {
  dfs <- x
  lmer(log(BC) ~ co2 + logMoist + (1|block) + (1|ring), data = x)
  })
llply(Mlist, Anova)
llply(Mlist, summary)
# negative moisture effect only in the first year... this is contradictive to
# the overall trend (i.e., it was really dry in the Year2 and dissimilarity was
# quite low compared to Year1, so overall treand of dissimilarity against
# moisture is positive)--> need more inspection

#################
## perform LMM ##
#################

# fit log and 2nd polynomial

# log-log----
Pfg_ll1 <- lmer(log(BC) ~ co2 * log(Moist) + (1|block) + (1|ring) + (1|id), data = PFG_vsDF)
Anova(Pfg_ll1)
Anova(Pfg_ll1, test.statistic = "F")

# 2nd polynomial----
Pfg_pl1 <- lmer(BC ~ co2 * (Moist + I(Moist^2)) + (1|block) + (1|ring) + (1|id), data = PFG_vsDF)
Pfg_pl2 <- update(Pfg_pl1, ~. - co2:I(Moist^2))
anova(Pfg_pl1, Pfg_pl2)
Anova(Pfg_pl2)

# model diagnosis----
ListMl <- list(Pfg_ll1, Pfg_pl2)
names(ListMl) <- c("log-log", "2nd polynomial")

plot(Pfg_ll1)
plot(Pfg_pl2)
  # indication of curvature..
par(mfrow = c(1, 2))
l_ply(names(ListMl), function(x) {
  l <- ListMl[[x]]
  qqnorm(resid(l), main = x)
  qqline(resid(l))
})
# squared R
ldply(ListMl, r.squared)
llply(ListMl, summary)

par(mfrow = c(1, 2))
visreg()

# although r.squared is slightly higher for 2nd plynomial, the number of 
# parameters is also bigger (5 vs. 4). Also, qqplot looks slightly better (?)
# for log-log. So use log.
visreg(Pfg_ll1, xvar = "Moist", trans = exp, by = "co2", overlay = TRUE)
visreg(Pfg_pl2, xvar = "Moist", by = "co2", overlay = TRUE)

# create a plot with confidense intervals
range(PFG_vsDF$Moist)
expDF <- expand.grid(co2 = c("amb", "elev"), Moist = seq(0.03, 0.2, length.out = 50))
bb <- bootMer(Pfg_ll1, FUN = function(x) predict(x, expDF, re.form = NA), nsim=500)
lci <- exp(apply(bb$t, 2, quantile, 0.025))
uci <- exp(apply(bb$t, 2, quantile, 0.975))
BC <- exp(bb$t0)
predDF <- cbind(lci, uci, BC, expDF)

p <- ggplot(data = PFG_vsDF, aes(x = log(Moist), y = log(BC), col = co2, fill = co2))
p2 <- p + geom_point(size = 3) +
  geom_line(aes(x = log(Moist), y = log(BC)), data = predDF) +
  geom_ribbon(aes(ymin = log(lci), ymax = log(uci), x = log(Moist)), 
              alpha = .2, color = NA, data = predDF) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"))
p2
ggsavePP(filename = "output//figs/PFG_DissimPredVal_Moist", plot = p2, width = 6, height = 4)

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

