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
