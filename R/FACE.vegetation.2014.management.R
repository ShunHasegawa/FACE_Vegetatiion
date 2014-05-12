rm(list=ls(all=TRUE))

source("functions/list_library.R")
source("functions//functions.R")
library(XLConnect)

################
# process 2014 #
################
# source("functions//prcss.veg.2014.R")

#############
# 2012 data #
#############
load("output/veg.12.Rdata")

# remove plant type
sp.12 <- names(veg.12)[-1:-4]
splt.12 <- strsplit(sp.12, "[.]")
new.sp.12 <- paste(sapply(splt.12, "[", 2),
                sapply(splt.12, "[", 3), sep = ".")
names(veg.12)[-1:-4] <- new.sp.12

# year
veg.12$year <- factor("2012")


############### 
# 2012 & 2014 #
###############
load("output/veg.14.Rdata")

veg.face <- rbind.fill(veg.12, veg.14)

# turn na into 0
veg.face[is.na(veg.face)] <- 0
save(veg.face, file = "output/FACE.Vegetation.2012&2013.Raw.Rdata")

veg.face.mlt <- melt(veg.face, id = c("year", "ring", "plot", "position", "cell"))

spp <- data.frame(sp = sort(levels(veg.face.mlt$variable)))
write.csv(spp, file = "output/spp.csv",row.names = FALSE)

# plant properties
spList <- read.csv("Data//FACE_Vegetation_sp.list.csv")

FACE.veg.rslt <- merge(veg.face.mlt, spList, by.x = "variable", by.y = "sp", all = TRUE)

save(FACE.veg.rslt, file = "output/FACE_Vegetation_2012&2014.Rdata")
load("output/FACE_Vegetation_2012&2014.Rdata")
########
# Figs #
########
source("functions//crt_brgrph.R")
dev.off()

###########
## Stats ##
###########

######
# CA #
######
source("functions/CA.analysis.R")

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

