# Environmental variables
load("output//Data/FACE_EnvironmenVars.RData")

# correlations
CorMatrix <- cor(Rm_ymc(EnvVarDF), use = "pairwise.complete.obs")
corrplot(CorMatrix)

# some values are really redundant so remove. Also some of them are hard to
# interpret so remove. e.g. EC, Theta75
RmVar <- c("sand", "silt", "clay", "TotalP", "OrganicC", "TotalP_CM", "TotalC",
           "TotalN", "Theta5", "Theta30", "Theta75", "EC", 
           "T5", "T10", "T20", "T30", "T50", "T100")
EnvVarDF <- EnvVarDF[, !names(EnvVarDF) %in% RmVar]
CorMatrix <- cor(Rm_ymc(EnvVarDF), use = "pairwise.complete.obs")
corrplot(CorMatrix)

# some of the variables are not complete for three years

# Variables measured only in the 1st year---
naCol_1stYear <-cbind(EnvVarDF[, c("year", "ring", "co2")], 
                      EnvVarDF[, apply(subset(EnvVarDF, year == 2014),
                                       2, function(x) any(is.na(x)))])
names(naCol_1stYear)[4] <- "OrganicMatter"

# these are all soil variable and shouldn't change dramatically year by year so
# just repeat the 1st year values
NewnaCol <- ldply(1:3, function(...) subset(naCol_1stYear, year == 2013))

# replace na columns with new columns
EnvVarDF[, names(NewnaCol)[c(-1:-3)]] <- Rm_ymc(NewnaCol)

# Variables measured only in the 1 and 2nd years----
naCol_2ndYear <-cbind(EnvVarDF[, c("year", "ring")], 
                      EnvVarDF[, apply(EnvVarDF, 2, function(x) any(is.na(x)))])
# these variables may vary year by year

###########################
# Analyse 3-year data set #
###########################
naCol_2ndYear
# these variables may vary year by year. so remove them for the time being
EnvDF_3df <- EnvVarDF[, apply(EnvVarDF, 2, function(x) !any(is.na(x)))]


#########
## RDA ##
#########

rda1 <- rda(log(RingSumVeg[, SppName] + 1) ~ year + co2 + OrganicMatter + Depth_HL + 
              ThetaHL + moist + FloorPAR + temp, EnvDF_3df)
rda1   # Constrained (canonical) axes explains 81%
summary(rda1)
anova(rda1)  # Significant association between species and environmental variables
anova(rda1, by = "axis")  # RDA1-5 are significant
anova(rda1, by = "margin") # not all explanatory variables is significant

# model simplification
rdaRes <- llply(list("backward", "forward", "both"), 
                function(x) ordistep(rda1, direction = x, trace = 0))
SimplModAnv <- llply(rdaRes, function(x) anova(x, by = "margin"))
SimplModAnv
# co2, OrganicMatter, Depth_HL, ThetaHL, moist and temp are significant
rda2 <- rdaRes[[3]]

anova(rda2, by = "axis")

# plot
p <- TriPlot(MultValRes = rda2, env = EnvDF_3df, yaxis = "RDA axis", axispos = c(1, 2, 3))
p2 <- p + scale_color_hue(labels = paste0(1:6, c("a", "e", "a", "e", "e", "a")))
ggsavePP(filename = "output//figs/RDA_EnvVar", plot = p2, width = 6, height = 6)

# Envrironmental variables have significant association with plant community, 
# but masks CO2 effect. Environmental variables and associated plant community 
# seems really similar between adjacent FACE rings within the same block (1&2, 
# 3&4 and 5&6). So conduct Partial RDA (fit block first to remove initial
# difference and fit co2 and year)

# Partial RDA
EnvDF_3df$block <- recode(EnvDF_3df$ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
prda1 <- rda(log(RingSumVeg[, SppName] + 1) ~ year + co2 + Condition(block), EnvDF_3df)
prda1 
# as expected condition accounts for a large propotion of variation (48 %)
# Semipartial R2 = 25 %
anova(prda1) # significant association
plot(prda1)
anova(prda1, by = "axis") # only RDA1 is significant
anova(prda1, by = "margin")
summary(prda1)

p <- TriPlot(MultValRes = prda1, env = EnvDF_3df, yaxis = "RDA axis", axispos = c(1, 2, 4), EnvNumeric = FALSE)
p2 <- p + scale_color_hue(labels = paste0(1:6, c("a", "e", "a", "e", "e", "a")))
ggsavePP(filename = "output//figs/PartialRDA_EnvVar", plot = p2, width = 6, height = 6)

# Plant community is different between treatments, but simply becuase of initial
# difference

#########
## CAP ##
#########
# bray-curtis is used to compute dissimilarity
Cap1 <- capscale(log(RingSumVeg[, SppName] + 1) ~ year + co2 + OrganicMatter + Depth_HL + 
                   moist + temp + FloorPAR, EnvDF_3df, dist = "bray")
anova(Cap1) # significant association between plant community and environmental variables
anova(Cap1, by = "axis") # CAP1-4 are significant
Cap1
# 76 % is explained by constrained axes
TriPlot(Cap1, env = EnvDF_3df, yaxis = "CAP", axispos = c(1:3), biplcons = 1,
        lowx = .2, lowy = .2)

# model simplification
anova(Cap1, by = "margin")
capRes <- llply(list("backward", "forward", "both"), 
                function(x) ordistep(Cap1, direction = x, trace = 0))
SimplModAnv <- llply(capRes, function(x) anova(x, by = "margin"))
SimplModAnv
# slightly different results depending on the order of deletion. year and temp
# probably have similar contribution. keep year for time being
cap2 <- capRes[[3]]
p <- TriPlot(cap2, env = EnvDF_3df, yaxis = "CAP", axispos = c(1:3), biplcons = 1,
        lowx = .2, lowy = .2)
p2 <- p 
ggsavePP(filename = "output//figs/FACE_CAP_EnvVar", plot = p2, width = 6, height = 6)

# Partial CAP
pcap1 <- capscale(log(RingSumVeg[, SppName] + 1) ~ year + co2 + Condition(block), EnvDF_3df, dist = "bray")
pcap1
# block explains 55 % of variation
anova(pcap1)
anova(pcap1, by = "margin")
p <- TriPlot(MultValRes = pcap1, env = EnvDF_3df, yaxis = "CAP", axispos = c(1:3), biplcons = 1, 
             lowx = .1, lowy = .1, EnvNumeric = FALSE)
ggsavePP(filename = "output//figs/FACE_PartialCAP_EnvVar", plot = p, width = 6, height = 6)





EnvDF_3df$block <- recode(EnvDF_3df$ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")

a1 <- anova(CapRes1, by = "term", permu = 999)
a1
a2 <- anova(CapRes1, by = "margin", permu = 999)
a2

CapRes2 <- capscale(log(RingSumVeg[, SppName] + 1) ~ year + co2 + OrganicMatter + Depth_HL +
                      moist + temp, EnvDF_3df, dist = "bray")
a3 <- anova(CapRes2, by = "margin", permu = 999)
a3
CapRes3 <- capscale(log(RingSumVeg[, SppName] + 1) ~ co2 + OrganicMatter + 
                      Depth_HL + moist + year, EnvDF_3df , dist = "bray")

a4 <- anova(CapRes3, by = "margin", permu = 999)
a4
plot(CapRes3)
plot(CapRes3, choices = c(1, 3))

CapRes4 <- capscale(log(RingSumVeg[, SppName] + 1) ~ co2 + OrganicMatter + 
                      Depth_HL + moist + year, Condition = EnvDF_3df$ring,
                    EnvDF_3df, dist = "bray", )
plot(CapRes4)



CapRes3 <- capscale(log(RingSumVeg[, SppName] + 1) ~ co2 + OrganicMatter + Depth_HL + moist
                    + temp, EnvDF_3df , dist = "bray")
?varpart
a3 <- anova(CapRes3, by = "margin", permu = 999)
a4 <- anova(CapRes3, by = "axis", permu = 999)
a3
a4
plot(CapRes3, choices = a2$Variance/a2$Variance[length(a2$Variance)]

str(a2)


a <- goodness(rda1, summ = TRUE)


str(CapRes1)

fit <- envfit(CapRes1, EnvDF_3df[, c("year", "FloorPAR", "OrganicMatter", "Depth_HL", 
                                  "ThetaHL", "moist", "temp")],
              perm = 999, display = "lc")



CapRes2 <- capscale(log(RingSumVeg[, SppName] + 1) ~ co2 + year + OrganicMatter + Depth_HL + ThetaHL +
                   moist + temp + FloorPAR, EnvDF_3df, dist = "bray",add = TRUE)
ordistep(CapRes2)
?prc



data(pyrifos)
week <- gl(11, 12, labels=c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
dose <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))
ditch <- gl(12, 1, length=132)

summary(data.frame(week, dose, ditch))
head(data.frame(week, dose, ditch), n = 30)
ctrl <- how(plots = Plots(strata = ditch,type = "free"),
            within = Within(type = "series"), nperm = 99)


summary(a2)
a3
anova(rda1, by = "margin")

?ordistep


rda3 <- rda(log(RingSumVeg[, SppName] + 1) ~ co2 * year + Condition(year + block), EnvDF_3df)

rda4 <- rda(log(RingSumVeg[, SppName] + 1) ~ yc + Condition(block), EnvDF_3df)
anova(rda4)
anova(rda4, by = "axis")
summary(rda4)
plot(rda4, scaling = 1)

CapRes2 <- capscale(log(RingSumVeg[, SppName] + 1) ~ yc + Condition(block), EnvDF_3df, dist = "bray",add = TRUE)
CapRes2 <- capscale(log(RingSumVeg[, SppName] + 1) ~ year + co2 + Condition(block), EnvDF_3df, dist = "bray",add = TRUE)
anova(CapRes2, by = "axis")
plot(CapRes2, scaling = 1)

CapRes3 <- capscale(log(RingSumVeg[, SppName] + 1) ~ OrganicMatter + Depth_HL + ThetaHL + 
                      moist + temp + FloorPAR, EnvDF_3df, dist = "bray",add = TRUE)



anova(rda3, by = "term", strata = EnvDF_3df$year)

plot(rda3, scaling = 1)

EnvDF_3df$block <- recode(EnvDF_3df$ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
EnvDF_3df$yc <- with(EnvDF_3df, year:co2)

rda2 <- rda(log(RingSumVeg[, SppName] + 1) ~ co2 + year + Condition(block), EnvDF_3df)
plot(rda2, scaling = 1)
anova(rda2, by = "axis", permutation = 999)
anova(rda2, by = "margin", permutation = 999)


spdf <- data.frame(
  sp1 = c(1,0,0,11,11,9,9,7,7,5),
  sp2 = c(0,0,1,4,5,6,7,8,9,10),
  sp3 = c(0,0,0,0,17,0,13,0,10,0),
  sp4 = c(0,0,0,0,7,0,10,0,13,0),
  sp5 = c(0,0,0,8,0,6,0,4,0,2),
  sp6 = c(0,0,0,1,0,2,0,3,0,4)
  )
xdf <- data.frame(Depth = 1:10, 
                  Substrate = c("sand", "sand", "sand", 
                                "other","coral", "other", 
                                "coral", "other", "coral", 
                                "other"))

rda.ts <- rda(spdf ~ Depth + Condition(Substrate), xdf)
rda.ts2 <- rda(spdf ~Substrate  + Condition(Depth), xdf)
rda.ts
rda.ts2
summary(rda.ts)
plot(rda.ts, scaling = 1)


?rda
vegan::scores(rda.ts, scaling = 1)

rda
smrdsa <- summary(rda.ts)
smrdsa$sites
smrdsa$constraints
for (i in 1:3) print(cor(smrdsa$sites[, i], smrdsa$constraints[, i]))
spenvcor (rda.ts)

str(smrdsa)

r2 <- 0.95971
ar2 <- 1 - (1 - r2) * (10-1)/(10-3-1)
anova(rda.ts, permutation = 999)
aa <- anova(rda.ts, permutation = 999, by = "axis")
aa2 <- anova(rda.ts, permutation = 999, by = "terms", type = "reduced")
aa2
aa2$Variance/sum(aa2$Variance)
str(aa)
?anova.ccv


rda3 <- rda(log(RingSumVeg[, SppName] + 1) ~ co2 + year + OrganicMatter + Depth_HL + 
              ThetaHL + moist + temp + FloorPAR, EnvDF_3df)
summary(rda3)
spenvcor(rda3)
anova(rda3, by = "axis", permutation = 999)
plot(rda3, scaling = 1)
anova(rda3, by = "margin", permutation = 999)
rda2 <- rda(log(RingSumVeg[, SppName] + 1))


CapRes1
str(CapRes1)
CapRes1$CCA$rank

Smmry_CAP <- summary(CapRes1)
str(Smmry_CAP)
summary(CapRes1)
SiteDF <- cbind(EnvDF_3df, Smmry_CAP$sites[, c("CAP1", "CAP2", "CAP3")])
envCap <- data.frame(predictor = row.names(Smmry_CAP$biplot), 
                     Smmry_CAP$biplot[, c("CAP1", "CAP2", "CAP3")], 
                     row.names = NULL, 
                     ring = "1", year = "2013")
envCap <- envCap[-(1:2), ]  
TrtCap <- data.frame(treatment = row.names(Smmry_CAP$centroids), 
                     Smmry_CAP$centroids[, c("CAP1", "CAP2", "CAP3")], 
                     ring = "1",
                     year = "2013",
                     row.names = NULL)
TrtCap$treatment <- factor(TrtCap$treatment, labels = paste0("Year", 1:3))
SpDf <- data.frame(species = gsub("[.]", "\n", row.names(Smmry_CAP$species)),
                   Smmry_CAP$species[, c("CAP1", "CAP2", "CAP3")], 
                   ring = "1",
                   year = "2013",
                   row.names = NULL)

listdf <- llply(list(SiteDF = SiteDF, envCap = envCap,  TrtCap = TrtCap, SpDf=  SpDf), 
                function(x) { NonCAPCol <- names(x)[!names(x) %in% c("CAP2", "CAP3")]
                              x_mlt <- melt(x, id = NonCAPCol)
                              
                              return(x_mlt)
                              }
                )

theme_set(theme_bw())
# CAP1 vs. CAP2
p <- ggplot(data = listdf$SiteDF, aes(x = CAP1, y = value, col = ring, shape = year))
p2 <- p + geom_point(size = 3)  
p3 <- p2 + 
  geom_segment(data = listdf$envCap,
               aes(x = 0, y = 0, xend = CAP1, yend = value), 
               arrow = arrow(length = unit(.1, "cm")), 
               alpha = .6,
               color = "black") +
  geom_text(data = listdf$envCap, 
            aes(x = CAP1 * 1.1 , y = value * 1.1, label = predictor), 
            alpha = .6, size = 2, lineheight = .7, color = "black")
p4 <- p3 + 
  geom_segment(data = listdf$TrtCap,
                        aes(x = 0, y = 0, xend = CAP1, yend = value), 
                        arrow = arrow(length = unit(.1, "cm")), 
                        alpha = .6,
                        color = "blue") +
  geom_text(data = listdf$TrtCap, 
            aes(x = CAP1 * 1.1 , y = value * 1.1, label = treatment), 
            alpha = .6, lineheight = .7, 
            color = "blue", size = 2, 
            fontface = "bold") +
  facet_grid(variable ~ .)
tdf <- subset(listdf$SpDf, abs(CAP1) > .22 | abs(value) > .22)
tdf <- within(tdf, {
  CAP1 <- CAP1 * 1.5
  value <- value * 1.5
})
p5 <- p4 + 
  geom_segment(data = tdf, 
               aes(x = 0, y = 0, xend = CAP1, yend = value), 
               arrow = arrow(length = unit(.1, "cm")),
               alpha = .6) +
  geom_text(data = tdf, 
            aes(x = CAP1 * 1.3, y = value * 1.3, label = species),
            size = 2, fontface = "bold.italic", 
            colour = "red", alpha = .6, lineheight = .7) +
  labs(x = "CAP1", y = "CAP axis")
p5
ggsavePP(filename = "output/figs/FACE_CAP_EnvVar", plot = p5, width = 6.5, height = 7)

# PERMANOVA
RngSppDF <- log(RingSumVeg[, SppName] + 1)
PermRes <- TypeIIpermanova(c("OrganicMatter", "Depth_HL", "ThetaHL", "moist", "temp",
                             "FloorPAR", "year"), EnvironmentalDF = EnvDF_3df)
PermRes
# only four variables seem to be important. 












# RDA
rda1 <- rda(log(RingSumVeg[, SppName] + 1) ~ year + co2 + OrganicMatter + 
              Depth_HL + ThetaHL + moist + temp + FloorPAR, EnvDF_3df)
Smmry_rda <- summary(rda1)

plot(rda1)

SiteDF <- cbind(EnvDF_3df, Smmry_rda$sites)
envCap <- data.frame(predictor = row.names(Smmry_rda$biplot), Smmry_rda$biplot * 2, 
                     row.names = NULL, 
                     ring = "1", year = "2013")
envCap <- envCap[-(1:3), ]  
TrtCap <- data.frame(treatment = row.names(Smmry_rda$centroids), 
                     Smmry_rda$centroids, 
                     ring = "1",
                     year = "2013",
                     row.names = NULL)
SpDf <- data.frame(species = row.names(Smmry_rda$species),
                   Smmry_rda$species, 
                   ring = "1",
                   year = "2013",
                   row.names = NULL)

# RDA1 vs. RDA2
p <- ggplot(data = SiteDF, aes(x = RDA1, y = RDA2, col = ring, shape = year))
p2 <- p + geom_point(size = 4)  
p3 <- p2 + 
  geom_segment(data = envCap,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(.1, "cm")), 
               alpha = .7,color = "black") +
  geom_text(data = envCap, 
            aes(x = RDA1 * 1.1 , y = RDA2 * 1.1, label = predictor), 
            alpha = .7, lineheight = .7, color = "black")
p4 <- p3 + geom_segment(data = TrtCap,
                        aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                        arrow = arrow(length = unit(.1, "cm")), 
                        alpha = .7,color = "red") +
  geom_text(data = TrtCap, 
            aes(x = RDA1 * 1.1 , y = RDA2 * 1.1, label = treatment), 
            alpha = .7, lineheight = .7, color = "red")
p4





# principal reponse curve

prc1 <- prc(log(RingSumVeg[, SppName] + 1), EnvDF_3df$co2, EnvDF_3df$year)
prc1
summary(prc1)

splogsum <- colSums(log(RingSumVeg[, SppName] + 1))
plot(prc1, type = "b", pch = 19)


 ?rda
theme(plot.margin=unit(c(0, 0.5, 0, 0), "lines"),
        strip.background = element_blank(),
        strip.text.x = element_blank(), 
        axis.title.y = element_text(lineheight = 1.2))
anova(prc1, strata = EnvDF_3df$by, first = TRUE)


EnvDF_3df$block <- factor(ifelse(EnvDF_3df$ring %in% c(1, 2), "A", 
                                 ifelse(EnvDF_3df$ring %in% c(3, 4), "B", "C")))

EnvDF_3df$by <- with(EnvDF_3df, block:year)



###########################
# Analyse 2-year data set #
###########################
# Anlyse only first 2-year data set so environmental variables are more complete
EnvDF_2df <- subsetD(EnvVarDF, year != 2015)
names(EnvDF_2df)

# Performe PCoA----
RngVegdf <- ddply(veg.face, .(year, ring, co2), function(x) colSums(x[, SppName]))
RngVegdf <- subsetD(RngVegdf, year != 2015)
RngSppDF <- RngVegdf[, SppName]
RngSiteDF <- RngVegdf[, !names(RngVegdf) %in% SppName]

PCoA <- cmdscale(d = vegdist(log(RngSppDF + 1), method = "bray"), eig = TRUE, k = 3)
PCoA_SiteScoreDF <- cbind(RngSiteDF, PCoA$points)
names(PCoA_SiteScoreDF)[4:6] <- paste0("PCoA", 1:3)

# Fit environmental variables----
# 1st two axes
par(mfrow = c(1, 2))
efLst <- llply(2:3, function(x) envfit(PCoA$points[, c(1, x)], 
                                       EnvDF_2df[, c("OrganicMatter", "moist", "temp")], permu = 999))

plot(PCoA2 ~ PCoA1, data = PCoA_SiteScoreDF, type = "n")
d_ply(PCoA_SiteScoreDF, .(year, ring), function(x) {
  points(PCoA2 ~ PCoA1, col = ring, pch = as.numeric(year), data = x, cex = 2)
})
plot(efLst[[1]], p.max = .05)

plot(PCoA3 ~ PCoA1, data = PCoA_SiteScoreDF, type = "n")
d_ply(PCoA_SiteScoreDF, .(year, ring), function(x) {
  points(PCoA3 ~ PCoA1, col = ring, pch = as.numeric(year), data = x, cex = 2)
})
plot(efLst[[2]], p.max = .05)

# perform permanova----

# too many variables for only 12 sites so remove some of the variables

PermRes_Year2 <- TypeIIpermanova(c("OrganicMatter", "moist", "temp"), EnvironmentalDF = EnvDF_2df)
PermRes_Year2

##





EnvDF_2df <- subsetD(EnvVarDF, year == 2013)
names(EnvDF_2df)

# Performe PCoA----
RngVegdf <- ddply(veg.face, .(year, ring, co2), function(x) colSums(x[, SppName]))
RngVegdf <- subsetD(RngVegdf, year == 2013)
RngSppDF <- RngVegdf[, SppName]
RngSiteDF <- RngVegdf[, !names(RngVegdf) %in% SppName]

PCoA <- cmdscale(d = vegdist(log(RngSppDF + 1), method = "bray"), eig = TRUE, k = 3)
PCoA_SiteScoreDF <- cbind(RngSiteDF, PCoA$points)
names(PCoA_SiteScoreDF)[4:6] <- paste0("PCoA", 1:3)

# Fit environmental variables----
# 1st two axes
par(mfrow = c(1, 2))
efLst <- llply(2:3, function(x) envfit(PCoA$points[, c(1, x)], EnvDF_2df[c("moist", "temp", "OrganicMatter")],
                                       permu = 720))

plot(PCoA2 ~ PCoA1, data = PCoA_SiteScoreDF, type = "n")
d_ply(PCoA_SiteScoreDF, .(year, ring), function(x) {
  points(PCoA2 ~ PCoA1, col = ring, pch = as.numeric(year), data = x, cex = 2)
})
plot(efLst[[1]], p.max = .05)

plot(PCoA3 ~ PCoA1, data = PCoA_SiteScoreDF, type = "n")
d_ply(PCoA_SiteScoreDF, .(year, ring), function(x) {
  points(PCoA3 ~ PCoA1, col = ring, pch = as.numeric(year), data = x, cex = 2)
})
plot(efLst[[2]], p.max = .05)

adonis(log(RngSppDF + 1) ~ moist + temp, 
       data = EnvDF_2df, method = "bray", permutation = 720)



# Peform permanova
TypeIIpermanova <- function(terms, EnvironmentalDF) {
  # reorder terms. the tem of interest is placed at the end
  TermList <- llply(1:length(terms), function(x) c(terms[-x], terms[x]))
  names(TermList) <- terms
  # create formulas
  formulaList <- llply(TermList, function(x) {
    form1 <- paste(x, collapse = "+")
    form <- paste("log(RngSppDF + 1) ~", form1)
    return(as.formula(form))
  })
  # perform permanova for each formula
  PermDF <- ldply(names(formulaList), function(x) {
    perm <- adonis(formulaList[[x]],  
                   method = "bray", 
                   data = EnvironmentalDF,
                   strata = EnvironmentalDF$year, 
                   permutations = 999)
    res <- data.frame(tem = x, perm$aov.tab[row.names(perm$aov.tab) == x, ])
    return(res)
  })
  return(PermDF)
}
# SS is TypeI SS so variance is allocated sequentially; thereby significant 
# levels depend on the order of terms. Hence obtain F and associated P values by
# fitting the term of interest at the end. Repeat this for all the variables. 

PermRes <- TypeIIpermanova(names(EnvDF_3df)[c(-1, -2)], EnvironmentalDF = EnvDF_3df)
PermRes
