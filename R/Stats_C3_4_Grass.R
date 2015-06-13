head(veg)

# Create C3, C4 grass summary DF

# PFG_plot sum for c3 and c4 grass wihtought unknown spp
PFGPltSum <- ddply(subsetD(veg, form == "Grass" & PFG %in% c("c3", "c4")),
                  .(year, co2, block, ring, plot, id, PFG), summarise, 
                  value = sum(value, na.rm = TRUE))

# cast and make data frame for anlysis
C34grassDF <- dcast(PFGPltSum, year + co2 + block + ring + plot + id ~ PFG)

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
# no co2 effect

##########################
# include random factors #
##########################
# add Residual id for each observation (i.e., each observation has individual
# ids)
C34grassDF$ResID <- factor(1:nrow(C34grassDF))

m1 <- glmer(yv ~ year * co2 +  (1|block) +(1|ring) + (1|id), 
            family = "binomial", data = C34grassDF)
summary(m1)
Anova(m1)
plot(m1)
# deviance >> df.resid, so it's highly overdispersed.. but P value is really 
# really low so it would probably be still significant even if I take
# overdispersion into account

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

par(mfrow = c(1, 3))
boxplot(C3Pr ~ year:ring, data = C34grassDF, main = "raw")
boxplot(asin(C3Pr) ~ year:ring, data = C34grassDF, main = "arcsin")
boxplot(logit(C3Pr) ~ year:ring, data = C34grassDF, main = "logit")
  # logit might be slightly better

pm1 <- lmer(logit(C3Pr) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = C34grassDF)
summary(pm1)
Anova(pm1)
Anova(pm1, test.statistic = "F")
plot(pm1)
# little bit wedge patter
qqnorm(resid(pm1))
qqline(resid(pm1))
# there one off the line. what if I remove
which(qqnorm(resid(pm1))$y == min(qqnorm(resid(pm1))$y))
pm2 <- lmer(logit(C3Pr) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = C34grassDF[-1, ])
plot(pm2)
qqnorm(resid(pm2))
qqline(resid(pm2))
AnvF_PropC3 <- Anova(pm2, test.statistic = "F")
AnvF_PropC3
# improved a lot. but needs to inspect more

# try different adjust values for logit transformation
logFun <- function(x){
  pm3 <- lmer(logit(C3Pr, adjust = x) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = C34grassDF)
  qqnorm(resid(pm3), main = x)
  qqline(resid(pm3))
}
par(mfrow = c(3, 3))
sapply(seq(0.02, 0.06, length = 9), logFun)
# not much difference

# contrast----

# contrast doesn't work with lmer. so use lme
tdf <- C34grassDF[-1, ]
lmeMod <- lme(logit(C3Pr) ~ year * co2, random = ~1|block/ring/id, data = tdf)
cntrst<- contrast(lmeMod, 
                  a = list(year = levels(tdf$year), co2 = "amb"),
                  b = list(year = levels(tdf$year), co2 = "elev"))
PropC3_CntrstRes <- cntrstTbl(cntrst, data = tdf, variable = "C3Prop")
PropC3_CntrstRes
# no significant differece between treatments, but still c3 proportion decreased
# at amb

########################################
# C3 proportion in the whole community #
########################################
head(veg)

# remove c3_4
veg2 <- subsetD(veg, PFG != "c3_4")

# get C3 proportion fo each plot
PropDF <- ddply(veg2, .(year, co2, block, ring, plot, id), function(x){
  C3prop = 1-sum(x$value[x$PFG == "c4"])/sum(x$value)
  Total = sum(x$value)
  data.frame(C3prop, Total)
  })

# glm
m1 <- glmer(C3prop ~ year * co2 +  (1|block) +(1|ring) + (1|id), 
            family = "binomial", data = PropDF, weight = Total)
summary(m1)
plot(m1)
# overdispersed so try transformation with LMM

# LMM
bxplts(value = "C3prop", xval = "co2", data = PropDF)
bxplts(value = "C3prop", xval = "ring", data = PropDF)

par(mfrow = c(1, 3))
boxplot(C3prop ~ year:ring, data = PropDF, main = "raw")
boxplot(asin(C3prop) ~ year:ring, data = PropDF, main = "arcsin")
boxplot(logit(C3prop) ~ year:ring, data = PropDF, main = "logit")
# try logit

rm1 <- lmer(logit(C3prop) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = PropDF)
summary(rm1)
Anova(rm1)
Anova(rm1, test.statistic = "F")
plot(rm1)
# slightly wedged pattern
qqnorm(resid(rm1))
qqline(resid(rm1))

# one complete outlier, so remove them
rmv <- which(qqnorm(resid(rm1))$y %in% min(qqnorm(resid(rm1))$y))
rm2 <- update(rm1, subset = -rmv)
plot(rm2)
qqnorm(resid(rm2))
qqline(resid(rm2))

# improved a lot
AnvF_PropC3_total <- Anova(rm2, test.statistic = "F")
AnvF_PropC3_total

# contrast
# contrast doesn't work with lmer. so use lme
tdf <- PropDF[-rmv, ]
lmeMod <- lme(logit(C3prop) ~ year * co2, random = ~1|block/ring/id, data = tdf)
cntrst<- contrast(lmeMod, 
                  a = list(year = levels(tdf$year), co2 = "amb"),
                  b = list(year = levels(tdf$year), co2 = "elev"))
PropC3_Total_CntrstRes <- cntrstTbl(cntrst, data = tdf, variable = "C3Prop")
PropC3_Total_CntrstRes
# no significant difference between treatment for each year

## ---- Stats_C34PropSmmry
# The model
pm2@call

# Chisq
Anova(pm2)

# F test
AnvF_PropC3

# Contrast
PropC3_CntrstRes

# Model diagnosis
plot(pm2)
qqnorm(resid(pm2))
qqline(resid(pm2))

## ---- Stats_C3_TotalPropSmmry
# The model
rm2@call

# Chisq
Anova(rm2)

# F test
AnvF_PropC3_total

# Contrast
PropC3_Total_CntrstRes

# Model diagnosis
plot(rm2)
qqnorm(resid(rm2))
qqline(resid(rm2))
