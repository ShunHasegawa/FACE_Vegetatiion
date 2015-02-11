################
# plant origin #
################
# create data frame for the analysis (also remove moss and lichen)
OrgnDat <- cast(subset(FACE.veg.rslt, !form %in% c("Moss", "Lichen")), year + co2 + block + ring + plot + id ~ origin, function(x) sum(x, na.rm = TRUE))

colSums(OrgnDat[, c("native", "naturalised", "NA")])
# some NA observations, but only a few so ignore for the time being

OrgnDat[names(OrgnDat) == "NA"] <- NULL

########
# GLMM #
########
OrgnDat$yv <- cbind(OrgnDat$native, OrgnDat$naturalised)

# add Residual id for each observation (i.e., each observation has individual
# ids)
OrgnDat$ResID <- factor(1:nrow(OrgnDat))

glm1 <- glmer(yv ~ year * co2 + (1|block) +(1|ring) + (1|id), 
              family = "binomial", 
              data = OrgnDat)
summary(glm1)
plot(glm1)
# deviance >> df.resid, so it's highly overdispersed..

glm2 <- glmer(yv ~ year * co2 + (1|block) +(1|ring) + (1|id) + (1|ResID), 
              family = "binomial", 
              data = OrgnDat)
summary(glm2)
# slightly improved, but still overdispersed

# ring donesn't seeme to be needed so remove
glm3 <- glmer(yv ~ year * co2 + (1|block) + (1|id) + (1|ResID), 
              family = "binomial", 
              data = OrgnDat)
summary(glm3)
# still overdispersed

#######################################################
# no good solution now so, just simply use proportion #
#######################################################

OrgnDat$NatPr <- with(OrgnDat, native/(native + naturalised))

par(mfrow = c(1, 2))
boxplot(NatPr ~ year:ring, data = OrgnDat)
boxplot(asin(NatPr) ~ year:ring, data = OrgnDat)

pm1 <- lmer(NatPr ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = OrgnDat)
summary(pm1)
Anova(pm1)
Anova(pm1, test.statistic = "F")
plot(pm1)

# arc sin transformation
pm2 <- lmer(asin(NatPr) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = OrgnDat)
summary(pm2)
Anova(pm2)
Anova(pm2, test.statistic = "F")
plot(pm2)
# slightly improved...
qqnorm(resid(pm2))
qqline(resid(pm2))
# no co2 effect
