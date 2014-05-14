################
# plant origin #
################
# glm

OrgnDat <- cast(FACE.veg.rslt, year + co2 + ring + plot ~ origin, function(x) sum(x, na.rm = TRUE))
OrgnDat <- OrgnDat[, c("year", "co2", "ring", "plot", "native", "naturalised"), drop = TRUE]
OrgnDat$yv <- cbind(OrgnDat$native, OrgnDat$naturalised*2)
OrgnDat$id <- OrgnDat$ring:OrgnDat$plot

# different rondom factor structure
m1 <- glmer(yv ~ year * co2 + (1|ring/plot/id), family = "binomial", data = OrgnDat, 
            control=glmerControl(optimizer="bobyqa"))
m2 <- glmer(yv ~ year * co2 + (1|ring/plot), family = "binomial", data = OrgnDat, 
            control=glmerControl(optimizer="bobyqa"))
m3 <- glmer(yv ~ year * co2 + (1|ring), family = "binomial", data = OrgnDat, 
            control=glmerControl(optimizer="bobyqa"))
m4 <- glmer(yv ~ year * co2 + (1|id), family = "binomial", data = OrgnDat, 
            control=glmerControl(optimizer="bobyqa"))
anova(m1, m2, m3, m4)
#use m4

m4_2 <- glmer(yv ~ year + co2 + (1|id), family = "binomial", data = OrgnDat, 
              control=glmerControl(optimizer="bobyqa"))
anova(m4, m4_2, test = "Chi")

summary(m4)
#df.resid <<< deviance... overdispersion..?

# model diagnosis
plot(m4)
qqnorm(resid(m4))
qqline(resid(m4))
