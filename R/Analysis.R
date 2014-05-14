rm(list=ls(all=TRUE))

library(xlsx)
library(plyr)
library(XLConnect)
library(ggplot2)
library(reshape)

source("R/functions.R")
################
# Process Data #
################
# source("R/FACE.vegetation.2014.management.R")

# Row data for multi variate analysis
load("output/Data//FACE_Vegetation_Raw.RData")

# Data frame with plant functional groups
load("output//Data//FACE_Vegetation_PFG.RData")

########
# Figs #
########
source("R//Figs.R")


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

# plot, PFG mean for grass
plt.pfg <- ddply(subset(FACE.veg.rslt, form == "Grass"), .(year, ring, plot, fgs), summarise, frq = sum(value))

# remove Grass:C3.4, not known
plt.pfg <- subset(plt.pfg, !(fgs %in% c("Grass:C3.4", "Grass:notknown")))
plt.pfg <- droplevels(plt.pfg)

# cast
names(plt.pfg)[5]  <- "value"
plt.pfg.cst <- cast(plt.pfg, year + ring + plot ~ fgs)
names(plt.pfg.cst)[4:5] <- c("C3", "C4")

plt.pfg.cst$yv <- cbind(plt.pfg.cst$C3, plt.pfg.cst$C4)
plt.pfg.cst$co2 <- factor(ifelse(plt.pfg.cst$ring %in% c("1", "4", "5"), "elev", "amb"))

# glm
library(lme4)
m1 <- glmer(yv ~ year * co2 + (1|ring/plot), family = "binomial", data = plt.pfg.cst)
m2 <- glmer(yv ~ year * co2 + (1|year/ring/plot), family = binomial, data = plt.pfg.cst)
plot(m1)

anova(m1, m2, test = "F")
summary(m2)

m2 <- glmer(yv ~ year + co2 + (1|ring/plot), family = binomial, data = plt.pfg.cst)
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
