head(FACE.veg.rslt)

# Create C3, C4 grass summary DF

# PFG_plot sum for c3 and c4 grass wihtought unknown spp
PFGPltSum <- ddply(subsetD(FACE.veg.rslt, form == "Grass" & PFG %in% c("c3", "c4")),
                  .(year, co2, block, ring, plot, id, PFG), summarise, 
                  value = sum(value, na.rm = TRUE))

# cast and make data frame for anlysis
C34grassDF <- cast(PFGPltSum, year + co2 + block +  ring + plot + id ~ PFG)

# glm
C34grassDF$yv <- cbind(C34grassDF$c3, C34grassDF$c4)

##################################
# try glm without rondom factors #
##################################
glm1 <- glm(yv ~ year * co2, data = C34grassDF, family = binomial)
summary(glm1)
  # highly overdispersed
glm2 <- glm(yv ~ year * co2, data = C34grassDF, family = quasibinomial)
summary(glm2)
glm3 <- update(glm2, ~. - year:co2)
anova(glm2, glm3, test = "Chisq")
anova(glm2, glm3, test = "F")
# no co2 effect

##########################
# include random factors #
##########################
# add Residual id for each observation (i.e., each observation has individual
# ids)
C34grassDF$ResID <- factor(1:nrow(C34grassDF))

m1 <- glmer(yv ~ year * co2
            +  (1|block) +(1|ring) + (1|id), 
            family = "binomial", 
            data = C34grassDF)
summary(m1)
plot(m1)
# deviance >> df.resid, so it's highly overdispersed..

m2 <- glmer(yv ~ year * co2
            + (1|block) +  (1|ring) + (1|id) + (1|ResID), 
            family = "binomial", 
            data = C34grassDF)
summary(m2)
plot(m2)
# slightly improved, but still overdispersed

#######################################################
# no good solution now so, just simply use proportion #
#######################################################

C34grassDF$C3Pr <- with(C34grassDF, c3/(c3 + c4))

par(mfrow = c(1, 2))
boxplot(C3Pr ~ year:ring, data = C34grassDF)
boxplot(asin(C3Pr) ~ year:ring, data = C34grassDF)

pm1 <- lmer(C3Pr ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = C34grassDF)
summary(pm1)
Anova(pm1)
Anova(pm1, test.statistic = "F")
plot(pm1)
# little bit wedge patter

# arc sin transformation
pm2 <- lmer(asin(C3Pr) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = C34grassDF)
summary(pm2)
Anova(pm2)
Anova(pm2, test.statistic = "F")
plot(pm2)
# slightly improved...
qqnorm(resid(pm2))
qqline(resid(pm2))

plot(allEffects(pm2))

# contrast

# contrast doesn't work with lmer. so use lme
lmeMod <- lme(asin(C3Pr) ~ year * co2, random = ~1|block/ring/id, data = C34grassDF)

cntrst<- contrast(lmeMod, 
                  a = list(year = levels(C34grassDF$year), co2 = "amb"),
                  b = list(year = levels(C34grassDF$year), co2 = "elev"))
cntrst
# no significant differece between treatments, but still c3 proportion decreased
# at amb

