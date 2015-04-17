head(PlotSumPFGMatrix)

# y value for glm
LgmDF <- within(PlotSumPFGMatrix,{yv <- cbind(legume, Non_legume)})

##################################
# try glm without rondom factors #
##################################
glm1 <- glm(yv ~ year * co2, data = LgmDF, family = binomial)
summary(glm1)
# highly overdispersed
glm2 <- glm(yv ~ year * co2, data = LgmDF, family = quasibinomial)
summary(glm2)
plot(glm2)
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

# plot var against mean
tdf <- ddply(LgmDF, .(year, ring), summarise, 
             Mean = mean(LegPr), Var = var(LegPr))
plot(Var ~ Mean, data = tdf) # not looking like binormial pattern

# boxplot
par(mfrow = c(2, 2))
boxplot(LegPr ~ year:ring, data = LgmDF, main = "raw")
boxplot(asin(LegPr) ~ year:ring, data = LgmDF, main = "arcsin")
boxplot(logit(LegPr) ~ year:ring, data = LgmDF, main = "logit")
# no distictive difference, try one by one

# raw
PropFm1 <- lmer(LegPr ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = LgmDF)
plot(PropFm1)
qqnorm(resid(PropFm1))
qqline(resid(PropFm1))
# wedged pattern

# arc sin transformation
PropFm2 <- lmer(asin(LegPr) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = LgmDF)
plot(PropFm2)
qqnorm(resid(PropFm2))
qqline(resid(PropFm2))
# no difference from the above

# logit
PropFm3 <- lmer(logit(LegPr, adjust = 0.001) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = LgmDF)
plot(PropFm3)
qqnorm(resid(PropFm3))
qqline(resid(PropFm3))
# improved a lot so use this
summary(PropFm3)
Anova(PropFm3)
# no co2 effect

# what if I remove one value off the q-q line
which(qqnorm(resid(PropFm3))$y == max(qqnorm(resid(PropFm3))$y))
tdf <- LgmDF[-4, ]
PropFm3_2 <- lmer(logit(LegPr, adjust = 0.001) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = tdf)
plot(PropFm3_2)
qqnorm(resid(PropFm3_2))
qqline(resid(PropFm3_2))
Anova(PropFm3_2)
# no difference

# different adjustment value
TryAdj <- function(adj) {
  m <- lmer(logit(LegPr, adjust = adj) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = LgmDF)
  qqnorm(resid(m), main = paste("Adjust = ", adj))
  qqline(resid(m))
}
par(mfrow = c(3, 3))
l_ply(seq(0.001, 0.002, length.out = 9), TryAdj)
# not much difference around .001 so just sty with .001

summary(PropFm3)
Anova(PropFm3)
# model simplification
PropFm4 <- stepLmer(PropFm3)
summary(PropFm4)
plot(PropFm4)
qqnorm(resid(PropFm4))
qqline(resid(PropFm4))

## ---- Stats_ForbPropSmmry
# The model
PropFm4@call
  # nothing was significant

# model diagnosis
plot(PropFm4)
qqnorm(resid(PropFm4))
qqline(resid(PropFm4))
