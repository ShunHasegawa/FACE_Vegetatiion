###########
## Stats ##
###########

######
# CA #
######
source("R/CA.analysis.R")

########
# GLMs #
########
source("functions/mvabund.analysis.R")

#########
# C3:C4 #
#########
head(FACE.veg.rslt)

# Create C3, C4 grass summary DF

# PFG_plot sum for c3 and c4 grass wihtought unknown spp
PFGPltSum<- ddply(subsetD(FACE.veg.rslt, form == "Grass" & PFG %in% c("c3", "c4")),
                          .(year, co2, ring, plot, PFG), summarise, value = sum(value, na.rm = TRUE))

# cast and make data frame for anlysis
C34grassDF <- cast(PFGPltSum, year + co2 + ring + plot ~ PFG)
C34grassDF$id <- C34grassDF$ring:C34grassDF$plot

## linear model ##
# proportion data 
df <- ddply(C34grassDF, .(ring, year), summarise, p = mean(c3/(c3 + c4)), V = var(c3/(c3 + c4)))

plot(V ~ p, data = df, xlim = c(0,1))

C34grassDF$prp <- C34grassDF$c3/(C34grassDF$c3 + C34grassDF$c4)
boxplot(prp ~ co2 * year, data = C34grassDF)



m1 <- lme(prp ~ co2 * year, random = ~1|ring/plot, data = C34grassDF)
m1 <- lme(asin(sqrt(prp)) ~ co2 * year, random = ~1|ring/plot, data = C34grassDF)
Anova(m1)
plot(m1)
qqnorm(m1, ~resid(.) | ring)
qqnorm(residuals.lm(m1))
qqline(residuals.lm(m1))


# glm
C34grassDF$yv <- cbind(C34grassDF$c3, C34grassDF$c4)

# different rondom factor structure
m1 <- glmer(yv ~ year * co2 + (1|ring/plot/id), family = "binomial", data = C34grassDF)
m2 <- glmer(yv ~ year * co2 + (1|ring/plot), family = "binomial", data = C34grassDF)
m3 <- glmer(yv ~ year * co2 + (1|ring), family = "binomial", data = C34grassDF)
m4 <- glmer(yv ~ year * co2 + (1|id), family = "binomial", data = C34grassDF)
anova(m1, m2, m3, m4)
#use m2

m2_2 <- glmer(yv ~ year + co2 + (1|ring/plot), family = "binomial", data = C34grassDF)
anova(m2, m2_2, test = "Chi")

summary(m2)
#df.resid <<< deviance... overdispersion..?

# model diagnosis
plot(m2)
qqnorm(resid(m2))
qqline(resid(m1))
# hmm, not too bad though,, not sure if it's correct..

