# Environmental variables
load("output//Data/FACE_EnvironmenVars.RData")

summary(EnvVarDF)

# correlations
CorMatrix <- cor(EnvVarDF[c(-1, -2)], use = "pairwise.complete.obs")
corrplot(CorMatrix)
# some values are really redundant so remove. Also some of them are hard to
# interpret so remove. e.g. EC, Theta75
RmVar <- c("silt", "clay", "OrganicC", "TotalP_CM",
           "Theta5", "Theta30", "Theta75", "EC", 
           "T5", "T10", "T20", "T30", "T50", "T100")
EnvVarDF <- EnvVarDF[, !names(EnvVarDF) %in% RmVar]
CorMatrix <- cor(EnvVarDF[c(-1, -2)], use = "pairwise.complete.obs")
corrplot(CorMatrix)

# some of the variables are not complete for three years

# Variables measured only in the 1st year---
naCol_1stYear <-cbind(EnvVarDF[, c("year", "ring")], 
                      EnvVarDF[, apply(subset(EnvVarDF, year == 2014),
                                       2, function(x) any(is.na(x)))])
# these are all soil variable and shouldn't change dramatically year by year so
# just repeat the 1st year values
NewnaCol <- ldply(1:3, function(...) subset(naCol_1stYear, year == 2013))

# replace na columns with new columns
EnvVarDF[, names(NewnaCol)[c(-1, -2)]] <- NewnaCol[, c(-1, -2)]

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

# Performe PCoA----
RngVegdf <- ddply(veg.face, .(year, ring, co2), function(x) colSums(x[, SppName]))
RngSppDF <- RngVegdf[, SppName]
RngSiteDF <- RngVegdf[, !names(RngVegdf) %in% SppName]

PCoA <- cmdscale(d = vegdist(log(RngSppDF + 1), method = "bray"), eig = TRUE, k = 3)
PCoA_SiteScoreDF <- cbind(RngSiteDF, PCoA$points)
names(PCoA_SiteScoreDF)[4:6] <- paste0("PCoA", 1:3)

# Fit environmental variables----
# 1st two axes
par(mfrow = c(1, 2))
efLst <- llply(2:3, function(x) envfit(PCoA$points[, c(1, x)], EnvDF_3df[, c(-1, -2)], permu = 999))

plot(PCoA2 ~ PCoA1, data = PCoA_SiteScoreDF, type = "n")
d_ply(PCoA_SiteScoreDF, .(year, ring), function(x) {
  points(PCoA2 ~ PCoA1, col = ring, pch = as.numeric(year), data = x, cex = 2)
})
plot(efLst[[1]], p.max = .05)
legend("bottomleft", 
       col = c(1:6, 1, 1, 1), pch = c(rep(19, 6), 1:3), 
       legend = c(paste("Ring", 1:6), paste("Year", 1:3)), 
       bty = "n", cex = .7)

plot(PCoA3 ~ PCoA1, data = PCoA_SiteScoreDF, type = "n")
d_ply(PCoA_SiteScoreDF, .(year, ring), function(x) {
  points(PCoA3 ~ PCoA1, col = ring, pch = as.numeric(year), data = x, cex = 2)
})
plot(efLst[[2]], p.max = .05)

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
#                    strata = EnvironmentalDF$ring, 
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
