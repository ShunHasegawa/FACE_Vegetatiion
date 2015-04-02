head(PlotSumVeg)

#################
# Dissimilarity #
#################
# Compute dissimiliraity for each plot between 2012 & 2014 and 2014 & 2015


YearDssmlrty <- function(x, spp) {
  df1 <- subset(x, year %in% c(2012, 2014))
  dis1 <- vegdist(log(df1[, spp] + 1), method = "bray")
  df2 <- subset(x, year %in% c(2014, 2015))
  dis2 <- vegdist(log(df2[, spp] + 1), method = "bray")
  dfs <- data.frame(year = c("Year1", "Year2"), BC = c(dis1, dis2))
  return(dfs)
}


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

MoistDF <- ddply(soilDF, .(year, ring, plot), summarise, Moist = mean(Moist, na.rm = TRUE))
boxplot(Moist ~ year:ring, data = MoistDF)

# merge data frame
vsDF <- merge(disDF, MoistDF, by = c("ring", "plot", "year"))
vsDF$id <- with(vsDF, ring:plot)

# fit to the model
theme_set(theme_bw())
p <- ggplot(data = vsDF, aes(x = Moist, y = BC, col = co2, shape = year))
p2 <- p + geom_point(size = 4) + geom_smooth(aes(fill = co2), method  = "lm") 
p2
# log
p <- ggplot(data = vsDF, aes(x = log(Moist), y = log(BC), col = co2, shape = year))
p2 <- p + geom_point(size = 4) + geom_smooth(aes(fill = co2), method  = "lm") 
p2
# log or 2nd polynomial

# log
Sp1 <- lmer(log(BC) ~ co2 * log(Moist) + (1|block) + (1|ring) + (1|id), data = vsDF)
Sp2 <- update(Sp1, ~. -co2:log(Moist))
Anova(Sp2)
plot(Sp2)

# 2nd polynomial
Sp_pl1 <- lmer(BC ~ co2 * (Moist + I(Moist^2)) + (1|block) + (1|ring) + (1|id), data = vsDF)
Sp_pl2 <- stepLmer(Sp_pl1)
summary(Sp_pl2)
Anova(Sp_pl2)
plot(Sp_pl2)

# model diagnosis
MList <- list(Sp2, Sp_pl2)
names(MList) <- c("log-log", "2nd polynomial")

par(mfrow = c(1, 2))
l_ply(names(MList), function(x) {
  ls <- MList[[x]]
  qqnorm(resid(ls), main = x)
  qqline(resid(ls))
})

# squared r
ldply(MList, r.squared.merMod)

# log-log is better
visreg(Sp2, xvar = "Moist", by = "co2", overlay = TRUE, 
       line = list(col = c("blue", "red")), points = list(col = c("blue", "red")))

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
summary(MoistDF)
PFG_vsDF <- merge(disDF, MoistDF, by = c("ring", "plot", "year"))

theme_set(theme_bw())
p <- ggplot(data = PFG_vsDF, aes(x = log(Moist), y = log(BC), col = co2, shape = year))
p2 <- p + geom_point(size = 4) + geom_smooth(method = "lm") + facet_grid(.~ block)
p2
# it seems that BC goes up till moisture of 0.1, then goes down. may be
# secon-plynomial might be appropriate
p <- ggplot(data = PFG_vsDF, aes(x = Moist, y = BC, col = co2))
p2 <- p + geom_point(size = 4) + 
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE))

#################
## perform LMM ##
#################

# log-log----
Pfg_ll1 <- lmer(log(BC) ~ co2 * log(Moist) + (1|block) + (1|ring) + (1|id), data = PFG_vsDF)
plot(Pfg_ll1)
Anova(Pfg_ll1)
Anova(Pfg_ll1, test.statistic = "F")

# 2nd polynomial----
Pfg_pl1 <- lmer(BC ~ co2 * (Moist + I(Moist^2)) + (1|block) + (1|ring) + (1|id), data = PFG_vsDF)
Pfg_pl2 <- stepLmer(Pfg_pl1)
summary(Pfg_pl2)
Anova(Pfg_pl2)
plot(Pfg_pl2)

# model diagnosis
Pfg_pl2 <- lmer(BC ~ co2 * I(Moist^2) + Moist + (1|block) + (1|ring) + (1|id), data = PFG_vsDF)
pfgMlList <- list(Pfg_ll1, Pfg_pl2)
names(pfgMlList) <- c("log-log", "2nd polynomial")

par(mfrow = c(1, 2))
l_ply(names(pfgMlList), function(x) {
  ls <- MList[[x]]
  qqnorm(resid(ls), main = x)
  qqline(resid(ls))
})

# what about squared r
ldply(pfgMlList, r.squared.merMod)
# 2n polynomial is very slightly better, it also uses one more parameter than
# the log-log

# not much difference betwween the above two models so just stay with log-log to
# be consisitent with species dissimilarity

# predicted values
par(mfrow = c(1, 2))
l_ply(names(pfgMlList), function(x) {
  visreg(pfgMlList[[x]], xvar = "Moist", by = "co2", overlay = TRUE, 
       line = list(col = c("blue", "red")), 
       points = list(col = c("blue", "red")),
       main = x,
       legend = FALSE)
  legend("bottomright", legend = c("amb", "co2"), col = c("blue", "red"), 
         lty = 1, bty = "n")
  })

# create a plot with confidense intervals
range(PFG_vsDF$Moist)
expDF <- expand.grid(co2 = c("amb", "elev"), Moist = seq(0.03, 0.2, length.out = 100))
bb <- bootMer(Pfg_ll1, FUN = function(x) predict(x, expDF, re.form = NA), nsim=500)
lci <- apply(bb$t, 2, quantile, 0.025)
uci <- apply(bb$t, 2, quantile, 0.975)
BC <- exp(bb$t0)
predDF <- cbind(lci, uci, BC, expDF)

p <- ggplot(data = PFG_vsDF, aes(x = log(Moist), y = log(BC), col = co2, fill = co2))
p2 <- p + geom_point(size = 4) +
  geom_ribbon(aes(ymin = lci, x = log(Moist), ymax = uci), alpha = .4, data = predDF)
p2


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

