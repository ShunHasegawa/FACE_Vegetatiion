summary(DivDF)

# orgnaise data frame
DivDF <- within(DivDF, {
  co2 <- factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
  block <- recode(ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
  id <- ring:plot
})

###########
# Eveness #
###########
bxplts(value = "J", xval = "co2", data = DivDF)
bxplts(value = "J", xval = "ring", data = DivDF)
par(mfrow = c(2,2))
boxplot(J ~ ring:year, data= DivDF, main = "raw")

Eml1 <- lmer(J ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF)
Eml2 <- stepLmer(Eml1)
summary(Eml2)
plot(Eml2)
qqnorm(resid(Eml2))
qqline(resid(Eml2))

#############
# Diversity #
#############
bxplts(value = "H", xval = "co2", data = DivDF)
bxplts(value = "H", xval = "ring", data = DivDF)

Dml1 <- lmer(H ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF)
Dml2 <- stepLmer(Dml1)
Anova(Dml2)
Anova(Dml2, test.statistic = "F")

####################
# Species richness #
####################
bxplts(value = "S", xval = "co2", data = DivDF)
bxplts(value = "S", xval = "ring", data = DivDF)

Sml1 <- glmer(S ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF,
              family = "poisson")
summary(Sml1)
# Devience >> df; highly overdispersed

Sml2 <- lmer(log(S) ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF)
Sml3 <- stepLmer(Sml2)
Anova(Sml3)
plot(Sml3)
qqnorm(resid(Sml3))
qqline(resid(Sml3))

##########
# Ancova #
##########
# reorganise data frame
DivDF_mlt <- melt(DivDF, id = c("year", "block", "co2", "ring", "plot", "id"))
DivDF_cst <- dcast(DivDF_mlt, block + co2 + ring + plot + id  + variable ~ year)
names(DivDF_cst)[7:8] <- c("Init", "Fin")

p <- ggplot(DivDF_cst, aes(x = Init, y = Fin, col = co2))
p2 <- p + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", aes(group = co2)) +
  facet_wrap(~ variable, scales =  "free")
p2

Hdf <- subsetD(DivDF_cst, variable == "H")
Jdf <- subsetD(DivDF_cst, variable == "J")
Sdf <- subsetD(DivDF_cst, variable == "S")

# Diveristy ----
h1 <- lmer(Fin ~ co2 * Init + (1|block) + (1|ring), data = Hdf)
summary(h1)
# no variation is explained by ring and block
h2 <- lm(Fin ~ co2 * Init, data = Hdf)
summary(h2)
par(mfrow = c(2, 2))
plot(h2)
# no co2 effect

# Eveness----
j1 <- lmer(Fin ~ co2 * Init + (1|block) + (1|ring), data = Jdf)
summary(j1)
# no variation is explained by ring and block
j2 <- lm(Fin ~ co2 * Init, data = Jdf)
summary(j2)
par(mfrow = c(2, 2))
plot(j2)
# no co2 effect

# Species----
s1 <- glmer(Fin ~ co2 * Init + (1|block) + (1|ring), data = Sdf, family = "poisson")
summary(s1)
# highly overdispersed
plot(s1)

# log transformation
par(mfrow = c(1, 2))
plot(Fin ~ Init, pch = 16, col = co2, data = Sdf)
plot(log(Fin) ~ Init, pch = 16, col = co2, data = Sdf)
s2 <- lmer(log(Fin) ~ co2 * Init + (1|block) + (1|ring), data = Sdf)
summary(s2)
# no co2 effect
