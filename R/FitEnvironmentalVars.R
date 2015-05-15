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
              moist + FloorPAR + temp, EnvDF_3df)
rda1   # Constrained (canonical) axes explains 76%
anova(rda1)  # Significant association between species and environmental variables
anova(rda1, by = "axis")  # RDA1-4 are significant
anova(rda1, by = "margin") # not all explanatory variables is significant

# model simplification
rdaRes <- llply(list("backward", "forward", "both"), 
                function(x) ordistep(rda1, direction = x, trace = 0))
SimplModAnv <- llply(rdaRes, function(x) anova(x, by = "margin"))
SimplModAnv
# co2, OrganicMatter, Depth_HL, moist and temp are significant
rda2 <- rdaRes[[3]]
anova(rda2, by = "axis")

# plot
p <- TriPlot(MultValRes = rda2, env = EnvDF_3df, yaxis = "RDA axis", axispos = c(1, 2, 3))
ggsavePP(filename = "output//figs/FACE_RDA_EnvVar", plot = p, width = 6, height = 6)
# moisture is very well representing high moisutre in ring5 followed by ring6
spenvcor(rda2)


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
ggsavePP(filename = "output//figs/PartialRDA_EnvVar", plot = p, width = 6, height = 6)

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
spenvcor(cap2)

p <- TriPlot(cap2, env = EnvDF_3df, yaxis = "CAP", axispos = c(1:3), biplcons = 1,
        lowx = .2, lowy = .2)
ggsavePP(filename = "output//figs/FACE_CAP_EnvVar", plot = p, width = 6, height = 6)

# Partial CAP
pcap1 <- capscale(log(RingSumVeg[, SppName] + 1) ~ year + co2 + Condition(block), EnvDF_3df, dist = "bray")
pcap1
# block explains 55 % of variation
anova(pcap1)
anova(pcap1, by = "axis") # only 1st axis is significant
anova(pcap1, by = "margin")
p <- TriPlot(MultValRes = pcap1, env = EnvDF_3df, yaxis = "CAP", axispos = c(1, 2, 4), biplcons = 1, 
             lowx = .1, lowy = .1, EnvNumeric = FALSE)
ggsavePP(filename = "output//figs/FACE_PartialCAP_EnvVar", plot = p, width = 6, height = 6)


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
