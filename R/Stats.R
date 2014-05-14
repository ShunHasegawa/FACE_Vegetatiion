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


## linear model ##

boxplot()


# cast and make data frame for anlysis
C34grassDF <- cast(PFGPltSum, year + co2 + ring + plot ~ PFG)




# glm
C34grassDF$yv <- cbind(C34grassDF$c3, C34grassDF$c4)
m1 <- glmer(yv ~ year * co2 + (1|ring/plot), family = "binomial", data = C34grassDF)
m2 <- glmer(yv ~ year + co2 + (1|ring/plot), family = "binomial", data = C34grassDF)
anova(m1, m2, test = "F")
# unable to remove year:co2 interaction

anova(m1)
summary(m1)


plot(m1)
summary(m1)
AIC(m1, m2)
anova(m1, test = "F")
plot(m1)

ds <- ddply(plt.pfg.cst, .(year, co2), summarise, ms = mean(p))
ddply(ds, .(year), summarise, ms[co2 == "elev"]/ms[co2 == "amb"])


library(effects)
plot(allEffects(m1))
plt.pfg.cst$p <- plt.pfg.cst$C3/(plt.pfg.cst$C3 + plt.pfg.cst$C4)
boxplot(p ~ co2 * year, data = plt.pfg.cst)




print(model.tables(aov.ex2,"means"),digits=3)  
ddply(plt.pfg.cst, .(year, co2), summarise, mean(p))

library(effects)
??plotallEffect

?model.tables
plot(m1)
?glmer

m <-  ddply(plt.pfg.cst, .(year, ring), summarise, means = mean(p), ss = var(p))
plot(ss ~ means, data = m, xlim = c(0, 1))
