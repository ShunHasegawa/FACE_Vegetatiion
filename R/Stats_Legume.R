head(FACE.veg.rslt)

# Create legume, non_legume summary DF

# Only forbs (remove shrubs)

# PFG_plot sum for legume and non_legume forbes wihtought unknown spp
PFGPltSum <- ddply(subsetD(FACE.veg.rslt, form == "Forb" & PFG %in% c("legume", "Non_legume")),
                   .(year, co2, block, ring, plot, id, PFG), summarise, 
                   value = sum(value, na.rm = TRUE))

# cast and make data frame for anlysis
LgmDF <- cast(PFGPltSum, year + co2 + block +  ring + plot + id ~ PFG)

# y value for glm
LgmDF$yv <- cbind(LgmDF$legume, LgmDF$Non_legume)

##################################
# try glm without rondom factors #
##################################
glm1 <- glm(yv ~ year * co2, data = LgmDF, family = binomial)
summary(glm1)
# highly overdispersed
glm2 <- glm(yv ~ year * co2, data = LgmDF, family = quasibinomial)
summary(glm2)
glm3 <- update(glm2, ~. - year:co2)
anova(glm2, glm3, test = "Chisq")
anova(glm2, glm3, test = "F")
par(mfrow = c(2, 2))
plot(glm3)
# no co2 effect

##########################
# include random factors #
##########################
# add Residual id for each observation (i.e., each observation has individual
# ids)
LgmDF$ResID <- factor(1:nrow(LgmDF))

m1 <- glmer(yv ~ year * co2 +  (1|block) +(1|ring) + (1|id), family = "binomial", 
            data = LgmDF)
summary(m1)
plot(m1)
# deviance >> df.resid, so it's highly overdispersed..
# also unknown warnings

m2 <- glmer(yv ~ year * co2
            + (1|block) +  (1|ring) + (1|id) + (1|ResID), 
            family = "binomial", 
            data = LgmDF)
summary(m2)
plot(m2)
# slightly improved, but still overdispersed

#######################################################
# no good solution now so, just simply use proportion #
#######################################################

LgmDF$LegPr <- with(LgmDF, legume/(legume + Non_legume))

par(mfrow = c(2, 2))

boxplot(LegPr ~ year:ring, data = LgmDF, main = "raw")
boxplot(asin(LegPr) ~ year:ring, data = LgmDF, main = "arcsin")
boxplot(logit(LegPr) ~ year:ring, data = LgmDF, main = "logit")

pm1 <- lmer(LegPr ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = LgmDF)
summary(pm1)
Anova(pm1)
Anova(pm1, test.statistic = "F")
plot(pm1)
# little bit wedge patter

# arc sin transformation
pm2 <- lmer(asin(LegPr) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = LgmDF)
summary(pm2)
Anova(pm2)
Anova(pm2, test.statistic = "F")
plot(pm2)
qqnorm(resid(pm2))
qqline(resid(pm2))
# slightly improved...

# logit
pm3 <- lmer(logit(LegPr, adjust = 0.001) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = LgmDF)
summary(pm3)
Anova(pm3)
plot(pm3)
qqnorm(resid(pm3))
qqline(resid(pm3))
# slightly better than above
Anova(pm3)
# there don't seem to be co2 effect anyway

# model simplification
pm4 <- stepLmer(pm3)
summary(pm4)
plot(pm4)
qqnorm(resid(pm4))
qqline(resid(pm4))

# anothe glm which handles overdispersion(?)
library(MASS)
t1 <- glmmPQL(yv ~ co2 * year, random = ~1|block/ring/id, family = binomial, 
              data = LgmDF)
summary(t1)
plot(t1)
qqnorm(resid(t1))
qqline(resid(t1))
