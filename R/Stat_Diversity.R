###########
# Eveness #
###########
bxplts(value = "J", xval = "co2", data = DivDF)
bxplts(value = "J", xval = "ring", data = DivDF)

Eml1 <- lmer(J ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF)
Anova(Eml1)
AnvF_Eml <- Anova(Eml1, test.statistic = "F")
AnvF_Eml
summary(Eml1)
plot(Eml1)
qqnorm(resid(Eml1))
qqline(resid(Eml1))

#############
# Diversity #
#############
bxplts(value = "H", xval = "co2", data = DivDF)
bxplts(value = "H", xval = "ring", data = DivDF)

Dml1 <- lmer(H ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF)
Anova(Dml1)
Anova(Dml1, test.statistic = "F")
plot(Dml1)
qqnorm(resid(Dml1))
qqline(resid(Dml1))
# One on the left bottom is quite out of the line. what if I remove it.
which(qqnorm(resid(Dml1))$y == min(qqnorm(resid(Dml1))$y))
Dml2 <- lmer(H ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF[-6, ])
plot(Dml2)
qqnorm(resid(Dml2))
qqline(resid(Dml2))
# looks a lot better. Need to inspect more about this point.
Anova(Dml2)
AnvF_Dml <- Anova(Dml2, test.statistic = "F")
AnvF_Dml
plot(allEffects(Dml2))
 # Diversity decreased in eCO2

# contrast----
# contrast doesn't work with lmer so rewrite the model with lme
tdf <- DivDF[-6, ]
Dml_lme <- lme(H ~ co2 * year, random = ~1|block/ring/id, data = tdf)

cntrst <- contrast(Dml_lme,
                   a = list(year = levels(tdf$year), co2 = "amb"),
                   b = list(year = levels(tdf$year), co2 = "elev"))
H_CntrstRes <- cntrstTbl(cntrst, data = tdf, variable = "H", digit = 2)
H_CntrstRes

####################
# Species richness #
####################
bxplts(value = "S", xval = "co2", data = DivDF)
bxplts(value = "S", xval = "ring", data = DivDF)
Sml1 <- glmer(S ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF,
              family = "poisson")
summary(Sml1)
plot(Sml1)
# Devience >> df; highly overdispersed

Sml2 <- lmer(log(S) ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF)
Anova(Sml2, test.statistic = "F")
plot(Sml2)
qqnorm(resid(Sml2))
qqline(resid(Sml2))
# left bottom is off the line. what if I remove
which(qqnorm(resid(Sml2))$y == min(qqnorm(resid(Sml2))$y))
Sml3 <- lmer(log(S) ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF[-6,])
AnvF_Sml <- Anova(Sml3, test.statistic = "F")
AnvF_Sml
plot(Sml3)
qqnorm(resid(Sml3))
qqline(resid(Sml3))
# improved a lot. but need to inspect more

# contrast----
tdf <- DivDF[-6, ]
Sml_lme <- lme(log(S) ~ co2 * year, random = ~1|block/ring/id, data = tdf)

cntrst <- contrast(Sml_lme,
                   a = list(year = levels(tdf$year), co2 = "amb"),
                   b = list(year = levels(tdf$year), co2 = "elev"))
S_CntrstRes <- cntrstTbl(cntrst, data = tdf, variable = "S", digit = 2)
S_CntrstRes

## ---- Stats_DiversityInx_Evenness
# Evenness
# The initial model
Eml1@call
Anova(Eml1)
# F
AnvF_Eml

## ---- Stats_DiversityInx_Diversity
# Shannon's index
# The model
Dml2@call
# Chisq
Anova(Dml2)
# F
AnvF_Dml
# Contrast
H_CntrstRes

## ---- Stats_DiversityInx_SpRichness
Sml3@call
# Chisq
Anova(Sml3)
# F test
AnvF_Sml
# Contrast
S_CntrstRes

# Model diagnosis
plot(Sml3)
qqnorm(resid(Sml3))
qqline(resid(Sml3))
