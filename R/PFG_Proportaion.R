head(veg)

# remove c3_4
veg2 <- subsetD(veg, PFG != "c3_4")

# get C3, grass, forb etc. proportion fo each plot
PFGpropDF <- dlply(veg, .(form), function(x) {
  dd <- ddply(x, .(year, co2, block, ring, plot, id), )
  
})

rm(dd)
PropDF <- ddply(veg2, .(year, co2, block, ring, plot, id), function(x){
  Total     = sum(x$value)
  C3prop    = 1-sum(x$value[x$PFG == "c4"])/Total
  Forbprop  = sum(x$value[x$form == "Forb"])/Total
  Grassprop = sum(x$value[x$form == "Grass"])/Total
  Woodprop  = sum(x$value[x$form == "Tree"])/Total
  Mossprop  = sum(x$value[x$form == "Moss"])/Total
  data.frame(Total, C3prop, Forbprop, Grassprop, Woodprop, Mossprop)
  })

########################################
# C3 proportion in the whole community #
########################################

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

#########
# Grass #
#########
bxplts(value = "Grassprop", xval = "co2", data = PropDF)
bxplts(value = "Grassprop", xval = "ring", data = PropDF)

par(mfrow = c(1, 3))
boxplot(Grassprop ~ year:ring, data = PropDF, main = "raw")
boxplot(asin(Grassprop) ~ year:ring, data = PropDF, main = "arcsin")
boxplot(logit(Grassprop) ~ year:ring, data = PropDF, main = "logit")
# try logit
gm1 <- lmer(logit(Grassprop) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = PropDF)

Anova(gm1)
Anova(gm1, test.statistic = "F")
# model simplification
gm2 <- stepLmer(gm1)
AnvF_grass <- Anova(gm2, test.statistic = "F")
AnvF_grass
plot(allEffects(gm2))

# model diagnosis
plot(gm2)
qqnorm(resid(gm2))
qqline(resid(gm2))

########
# Forb #
########
bxplts(value = "Forbprop", xval = "co2", data = PropDF)
bxplts(value = "Forbprop", xval = "ring", data = PropDF)

par(mfrow = c(1, 3))
boxplot(Forbprop ~ year:ring, data = PropDF, main = "raw")
boxplot(asin(Forbprop) ~ year:ring, data = PropDF, main = "arcsin")
boxplot(logit(Forbprop) ~ year:ring, data = PropDF, main = "logit")
# try logit
fm1 <- lmer(logit(Forbprop) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = PropDF)
Anova(fm1)
Anova(fm1, test.statistic = "F")

# model simplification
fm2 <- stepLmer(fm1)
AnvF_forb <- Anova(fm2, test.statistic = "F")
AnvF_forb

# model diagnosis
plot(fm2)
qqnorm(resid(fm2))
qqline(resid(fm2))

# what if I remove two outliers
rmv <- which(qqnorm(resid(fm2), plot.it = FALSE)$y %in% 
               sort(qqnorm(resid(fm2), plot.it = FALSE)$y)[1:2])

fm3 <- update(fm1, subset = -rmv)
plot(fm3)
qqnorm(resid(fm3))
qqline(resid(fm3))
Anova(fm3) # no difference so just present fm2

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


