## ---- OrganiseEnvVars
# Environmental variables
load("output//Data/FACE_EnvironmenVars.RData")

names(EnvVarDF)
# remove all minor metals
rmMetal <- c("Silver", "Arsenic", "Lead", "Cadmium", "Chromium", "Copper", "Manganese", 
             "Nickel", "Selenium", "Zinc", "Mercury", "Iron", "Aluminium", "Boron", 
             "Silicon", "Vanadium", "Cobalt", "Molybdenum", "Barium", "Calcium", "Magnesium",
             "Potassium", "Sodium")

EnvVarDF <- EnvVarDF[, !names(EnvVarDF) %in% rmMetal]

# correlations
CorMatrix <- cor(Rm_ymc(EnvVarDF), use = "pairwise.complete.obs")
p = CorMatrix
p[is.na(CorMatrix)]=0.2 
p[is.na(CorMatrix)==F]=0
CorMatrix[is.na(CorMatrix)]=0
corrplot.mixed(CorMatrix, p.mat=p, tl.cex = .5)

# some values are really redundant so remove. Also some of them are hard to
# interpret so remove. e.g. EC, Theta75
RmVar <- c("sand", "silt", "clay", "OrganicMatter",
           "OrganicC", "OrganicC", "TotalP_CM", 
           "Theta5", "Theta30", "Theta75", "EC", 
           "T5", "T10", "T20", "T30", "T50", "T100", 
           "TotalN", "Phosphorus", "nitrification")
EnvVarDF <- EnvVarDF[, !names(EnvVarDF) %in% RmVar]
CorMatrix <- cor(Rm_ymc(EnvVarDF), use = "pairwise.complete.obs")
corrplot.mixed(CorMatrix, tl.cex = .5, mar = par()$mar)

# some of the variables are not complete for three years

# Variables measured only in the 1 and 2nd years----
naCol_2ndYear <-cbind(EnvVarDF[, c("year", "ring")], 
                      EnvVarDF[, apply(EnvVarDF, 2, function(x) any(is.na(x)))])

###########################
# Analyse 3-year data set #
###########################
naCol_2ndYear
# these variables may vary year by year. so remove them for the time being
EnvDF_3df <- EnvVarDF[, apply(EnvVarDF, 2, function(x) !any(is.na(x)))]

EnvDF_3df$block <- recode(EnvDF_3df$ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
## ---- RunRDA

#########
## RDA ##
#########
# combine environment and spp df
seDF <- merge(RingSumVeg, EnvDF_3df, by = c("year", "ring", "block", "co2"))

# possible explanatory variables
expl <- c("co2",  "TotalC", "moist", "Drysoil_ph", "Depth_HL", "FloorPAR", "temp")

###############
# Single term #
###############
# R2adj for each single term

adjR_singl_Lst <- list()
for (i in 1:3) {
  dd <- subsetD(seDF, year == levels(seDF$year)[i])
  spdd <- dd[ , SppName]
  # formula for each variable
  singl_fmls <- llply(paste("log(spdd + 1) ~", expl), as.formula)
  names(singl_fmls) <- expl
  adjR_singl <- ldply(singl_fmls, function(y) 
    RsquareAdj(rda(y, data = dd))$adj.r.squared, 
    .id = "variable")
  adjR_singl_Lst[[i]] <- adjR_singl
  rm(dd, spdd, singl_fmls, adjR_singl)
}
names(adjR_singl_Lst) <- paste0("Year", 1:3)

# Get variables with positive R2adj
PosAdjR <- llply(adjR_singl_Lst, function(x) as.character(x$variable[x$V1 > 0]))

# Formula for full models
FullFormula <- llply(PosAdjR, function(x) {
  comb_exp <- combn(x, 4) # combination of four
  expl_fml <-apply(comb_exp, 2, function(x) paste(x, collapse = "+"))
  return(expl_fml)
  })

# Y matric (lefthand part)
LH <- list("log(df2013[ , SppName] + 1) ~",
           "log(df2014[ , SppName] + 1) ~",
           "log(df2015[ , SppName] + 1) ~")
fmls <- llply(list(Year1 = 1, Year2 = 2, Year3 = 3), 
              function(x) llply(paste(LH[[x]], FullFormula[[x]]), as.formula))

##############
## 1st year ##
##############
df2013 <- subsetD(seDF, year == 2013)

# There are too many environmental variables to fit. so choose four which showed
# highest R2adj 
# adjusted R2
adjR <- ldply(fmls$Year1, function(x) RsquareAdj(rda(x, data = df2013))$adj.r.squared)

# highest R2
rr <- rda(fmls$Year1[[which(max(adjR) == adjR)]], df2013)

# check multicollinearity
vif.cca(rr)
anova(rr, permutations = allPerms(6))
rr2 <- rda(log(df2013[ , SppName] + 1) ~ 1, df2013)
rr3 <- ordiR2step(rr2, rr, permutations = allPerms(6), direction = "forward", Pin = .1)

# summary result
rda2013 <- list(IniRda = rr, FinRda = rr3)

##############
## 2nd year ##
##############
df2014 <- subsetD(seDF, year == 2014)

# adjusted R2
adjR <- laply(fmls$Year2, function(x) RsquareAdj(rda(x, data = df2014))$adj.r.squared)

# highest R2
rr <- rda(fmls$Year2[[which(max(adjR) == adjR)]], df2014)

# check multicolliniarity using vif
vif.cca(rr)
anova(rr, permutations = allPerms(6))
  # not significant

# choose only three terms
comb_exp <- combn(PosAdjR$Year2, 3)
expl_fml <-apply(comb_exp, 2, function(x) paste(x, collapse = "+"))
fmls_3 <- llply(paste("log(df2014[ , SppName] + 1) ~", expl_fml), as.formula)

adjR <- laply(fmls_3, function(x) RsquareAdj(rda(x, data = df2014))$adj.r.squared)
rr <- rda(fmls_3[[which(max(adjR) == adjR)]], df2014)
vif.cca(rr)
anova(rr, permutations = allPerms(6))
# good
rr2 <- rda(log(df2014[ , SppName] + 1) ~ 1, df2014)
rr3 <- ordiR2step(rr2, rr, permutations = allPerms(6), direction = "forward", Pin = .1)

# summary result
rda2014 <- list(IniRda = rr, FinRda = rr3)

##############
## 3rd year ##
##############
df2015 <- subsetD(seDF, year == 2015)

# adjusted R2
adjR <- laply(fmls$Year3, function(x) RsquareAdj(rda(x, data = df2015))$adj.r.squared)

# highest R2
rr <- rda(fmls$Year3[[which(max(adjR) == adjR)]], df2015)

# check multicollinearity
vif.cca(rr)
# TotalC, Depth_HL >10. 

# Second highset R2adj
adjR[which(max(adjR) == adjR)] <- NA
rr <- rda(fmls$Year3[[which(max(adjR, na.rm = TRUE) == adjR)]], df2015)
vif.cca(rr)
anova(rr, permutations = allPerms(6))
# good

rr2 <- rda(log(df2015[ , SppName] + 1) ~ 1, df2015)
rr3 <- ordiR2step(rr2, rr, permutations = allPerms(6), direction = "forward", Pin = .1)

# summary result
rda2015 <- list(IniRda = rr, FinRda = rr3)

#############
## Summary ##
#############
# R2adj for initial full model
FuladjR_pv <- ldply(RdaLst, function(x){ 
  data.frame(
  Full = RsquareAdj(x$IniRda)$adj.r.squared, 
  Pr = anova(x$IniRda, permutations = allPerms(6))$Pr[1])
  }, 
  .id = "year")
FuladjR_pv <- dcast(variable ~ year, data = melt(FuladjR_pv, id = "year"))

# Adjusted R2 for each term
AdjTbl <- dcast(variable ~ .id, data = ldply(adjR_singl_Lst), value.var = "V1")
AdjTbl <- rbind(AdjTbl, FuladjR_pv)
AdjTbl[, 2:4] <- round(AdjTbl[, 2:4], 3)
# replace negative values with <0
AdjTbl[AdjTbl < 0] <- "<0" 
# save
write.csv(AdjTbl, file = "output/table/RDA_AdjR_Speices.csv")

RdaLst <- list(Year1 = rda2013, Year2 = rda2014, Year3 = rda2015)
llply(RdaLst, function(x) x$IniRda)

AnovaRes <- ldply(RdaLst, 
                  function(x) {
                    dd <- anova(x$FinRda, permutations = allPerms(6), by = "margin")
                    nr <- row.names(dd)
                    data.frame(variable = nr, dd)
                    }, 
                  .id = "year")

AdjR2 <- llply(RdaLst, function(x) RsquareAdj(x$FinRda)$adj.r.squared)
AdjR2 <- laply(RdaLst, function(x) RsquareAdj(x$IniRda)$adj.r.squared)
round(AdjR2, 3)
llply(RdaLst, function(x) anova(x$IniRda, permutations = allPerms(6)))
llply(RdaLst, function(x) anova(x$FinRda, permutations = allPerms(6)))

#######
# PFG #
#######
peDF <- merge(RingSumPFGMatrix, EnvDF_3df, by = c("year", "ring", "block", "co2")) 

###############
# Single term #
###############
# R2adj for each single term

adjR_singl_Lst <- list()
for (i in 1:3) {
  dd <- subsetD(peDF, year == levels(peDF$year)[i])
  spdd <- dd[ , PFGName]
  # formula for each variable
  singl_fmls <- llply(paste("log(spdd + 1) ~", expl), as.formula)
  names(singl_fmls) <- expl
  adjR_singl <- ldply(singl_fmls, function(y) 
    RsquareAdj(rda(y, data = dd))$adj.r.squared, 
    .id = "variable")
  adjR_singl_Lst[[i]] <- adjR_singl
  rm(dd, spdd, singl_fmls, adjR_singl)
}
names(adjR_singl_Lst) <- paste0("Year", 1:3)

# Get variables with positive R2adj
PosAdjR <- llply(adjR_singl_Lst, function(x) as.character(x$variable[x$V1 > 0]))
 # less than 4 terms left

# Formula for full models
FullFormula <- llply(PosAdjR, function(x) paste(x, collapse = "+"))

# Y matric (lefthand part)
LH <- list("log(df2013[ , PFGName] + 1) ~",
           "log(df2014[ , PFGName] + 1) ~",
           "log(df2015[ , PFGName] + 1) ~")
fmls <- llply(list(Year1 = 1, Year2 = 2, Year3 = 3), 
              function(x) llply(paste(LH[[x]], FullFormula[[x]]), as.formula))

##############
## 1st year ##
##############
df2013 <- subsetD(peDF, year == 2013)

# adjusted R2
adjR <- laply(fmls$Year1, function(x) RsquareAdj(rda(x, data = df2013))$adj.r.squared)

# highest R2
rr <- rda(fmls$Year1[[which(max(adjR) == adjR)]], df2013)
# check multicollinearity
vif.cca(rr)
anova(rr, permutations = allPerms(6))
rr2 <- rda(log(df2013[ , PFGName] + 1) ~ 1, df2013)
rr3 <- ordiR2step(rr2, rr, permutations = allPerms(6), direction = "forward", Pin = .1)

# summary result
rda2013 <- list(IniRda = rr, FinRda = rr3)

##############
## 2nd year ##
##############
df2014 <- subsetD(peDF, year == 2014)

# adjusted R2
adjR <- ldply(fmls$Year2, function(x) RsquareAdj(rda(x, data = df2014))$adj.r.squared)

# highest R2
rr <- rda(fmls$Year2[[which(max(adjR) == adjR)]], df2014)
anova(rr, permutations = allPerms(6))
rr2 <- rda(log(df2014[ , PFGName] + 1) ~ 1, df2014)
rr3 <- ordiR2step(rr2, rr, permutations = allPerms(6), direction = "forward", Pin = .1)

# summary result
rda2014 <- list(IniRda = rr, FinRda = rr3)

##############
## 3rd year ##
##############
df2015 <- subsetD(peDF, year == 2015)

# adjusted R2
adjR <- laply(fmls$Year3, function(x) RsquareAdj(rda(x, data = df2015))$adj.r.squared)

# highest R2
rr <- rda(fmls$Year3[[which(max(adjR) == adjR)]], df2015)
# check multicollinearity
vif.cca(rr)
anova(rr, permutations = allPerms(6))
# not significant so use two terms instead

# try only two variables
comb_exp <- combn(PosAdjR$Year3, 2) 
expl_fml <-apply(comb_exp, 2, function(x) paste(x, collapse = "+"))
fmls_2 <- llply(paste("log(df2015[ , PFGName] + 1) ~", expl_fml), as.formula)
adjR <- laply(fmls_2, function(x) RsquareAdj(rda(x, data = df2015))$adj.r.squared)
rr <- rda(fmls_2[[which(max(adjR) == adjR)]], df2015)
vif.cca(rr)
anova(rr, permutations = allPerms(6))
# good

rr2 <- rda(log(df2015[ , PFGName] + 1) ~ 1, df2015)
rr3 <- ordiR2step(rr2, rr, permutations = allPerms(6), direction = "forward", Pin = .1)

# summary result
rda2015 <- list(IniRda = rr, FinRda = rr3)

#############
## Summary ##
#############
RdaLst_pfg <- list(Year1 = rda2013, Year2 = rda2014, Year3 = rda2015)

# R2adj for initial full model
FuladjR_pv_pfg <- ldply(RdaLst_pfg, function(x){ 
  data.frame(
    Full = RsquareAdj(x$IniRda)$adj.r.squared, 
    Pr = anova(x$IniRda, permutations = allPerms(6))$Pr[1])
}, 
.id = "year")
FuladjR_pv_pfg <- dcast(variable ~ year, data = melt(FuladjR_pv_pfg, id = "year"))

# Adjusted R2 for each term
AdjTbl_pfg <- dcast(variable ~ .id, data = ldply(adjR_singl_Lst), value.var = "V1")
AdjTbl_pfg <- rbind(AdjTbl_pfg, FuladjR_pv_pfg)
AdjTbl_pfg[, 2:4] <- round(AdjTbl_pfg[, 2:4], 3)
# replace negative values with <0
AdjTbl_pfg[AdjTbl_pfg < 0] <- "<0" 
AdjTbl_pfg
write.csv(AdjTbl_pfg, file = "output/table/RDA_AdjR_PFG.csv")

AnovaRes_pfg <- ldply(RdaLst_pfg, function(x) {
                        aa <- anova(x$FinRda, permutations = allPerms(6), by = "margin")
                        rn <- row.names(aa)
                        data.frame(variable = rn, aa)},
                  .id = "year")
AdjR2 <- laply(RdaLst_pfg, function(x) RsquareAdj(x$FinRda)$adj.r.squared)
round(AdjR2, 3)
llply(RdaLst_pfg, function(x) anova(x$FinRda, permutations = allPerms(6)))


#################
# CO2 treatment #
#################

#################
## All speices ##
#################

#############
## Ambient ##
#############
ambDF <- subsetD(seDF, co2 == "amb")

# permutaiton Possible number of permutation is (3!)^3 = 216. it's small so need
# to use exact permutation test

allPerms(9, 
         how(within = Within(type = "free"), plot = Plots(strata = ambDF2$ring)))
# Null hypothesis is no difference between years. Ring is no exchangeable but 
# subplots are. NOTE # First three colmuns of the above indices are for the 1st 
# factor of ring, and four to six are for the 2nd and seven to nine are for 3rd.
# So that data frame needs to be reordered to match this. 

# reorder by ring
ambDF <- ambDF[order(ambDF$ring), ]

# Run RDA
rda_amb <- rda(log(ambDF[, SppName] + 1) ~ year + moist + TotalC + 
                 Drysoil_ph + FloorPAR + Condition(ring) , ambDF)

hh <- allPerms(9, how(within = Within(type = "free"), plot = Plots(strata = ambDF$ring)))

# model simplificaiton
rdaRes <- llply(list("backward", "forward", "both"), 
                function(x) ordistep(rda_amb, direction = x, trace = 0, 
                                     permutations = hh))
SimplModAnv <- llply(rdaRes, function(x) anova(x, by = "margin", permutations = hh))
SimplModAnv
# only year is left
rda_amb2 <- rdaRes[[1]]
anova(rda_amb2, permutations = hh)

# Axis
anova(rda_amb2, permutations = hh, by = "axis")
# RDA1 is significant

# replace time with temp
rda_amb <- rda(log(ambDF[, SppName] + 1) ~ temp + moist + TotalC + 
                 Drysoil_ph + FloorPAR + Condition(ring) , ambDF)
rdaRes <- llply(list("backward", "forward", "both"), 
                function(x) ordistep(rda_amb, direction = x, trace = 0, 
                                     permutations = hh))
SimplModAnv <- llply(rdaRes, function(x) anova(x, by = "margin", permutations = hh))
SimplModAnv

rda_amb_temp <- rdaRes[[1]]
anova(rda_amb_temp, permutations = hh)

# plot
p <- TriPlot(MultValRes = rda_amb_temp, env = ambDF, yaxis = "RDA axis", axispos = c(1, 2, 3))
ggsavePP(filename = "output//figs/FACE_RDA_EnvVar_Ambient", plot = p, width = 6, height = 6)

#########
## co2 ##
#########

elvDF <- subsetD(seDF, co2 == "elev")
# reorder by ring
elvDF <- elvDF[order(elvDF$ring), ]

# Run RDA
rda_elev <- rda(log(elvDF[, SppName] + 1) ~ year + moist + TotalC + 
                  Drysoil_ph + Condition(ring) , elvDF)

hh <- allPerms(9, how(within = Within(type = "free"), plot = Plots(strata = elvDF$ring)))

# model simplificaiton
rdaRes <- llply(list("backward", "forward", "both"), 
                function(x) ordistep(rda_elev, direction = x, trace = 0, 
                                     permutations = hh))
SimplModAnv <- llply(rdaRes, function(x) anova(x, by = "margin", permutations = hh))
SimplModAnv

# year and moist
rda_elev2 <- rdaRes[[3]]
anova(rda_elev2, permutations = hh)
anova(rda_elev2, permutations = hh, by = "margin")

# year maybe replaceable with temp
summary(lm(temp ~ year, data = seDF))

# use temp itead of time
rda_elev <- rda(log(elvDF[, SppName] + 1) ~ temp + moist + TotalC + 
                  Drysoil_ph + Condition(ring) , elvDF)
rdaRes <- llply(list("backward", "forward", "both"), 
                function(x) ordistep(rda_elev, direction = x, trace = 0, 
                                     permutations = hh))
SimplModAnv <- llply(rdaRes, function(x) anova(x, by = "margin", permutations = hh))
SimplModAnv

rda_elev_temp <- rdaRes[[1]]
anova(rda_elev_temp, permutations = hh)
anova(rda_elev_temp, permutations = hh, by = "margin")

# plot
p <- TriPlot(MultValRes = rda_elev_temp, env = ambDF, yaxis = "RDA axis", axispos = c(1, 2, 3))
ggsavePP(filename = "output//figs/FACE_RDA_EnvVar_eCO2", plot = p, width = 6, height = 6)

########
# Plot #
########
p2 <- PlotRDA_Year(rdaResLst = list(rda_amb2, rda_elev2))
ggsavePP(filename = "output/figs/FACE_RDAvsYearbyCO2", plot = p2, 
         width = 6.65, height = 4)


#####################
## 3-year data set ## 
#####################
# From the above
#analysis, moist, temp, TotalC and Dry_soilph are determied to be imporatnt
#driver

rda_all <- rda(log(seDF[, SppName] + 1) ~ TotalC + moist + Drysoil_ph + temp + co2, 
               seDF)
# can't run anova as it is. cause different perumutation units need to be
# defined for year and co2. so anyway create a triplot and see the pattern.
rda_all

# plot
p <- TriPlot(MultValRes = rda_all, env = seDF, yaxis = "RDA axis", axispos = c(1, 2, 3), centcons = 2)
ggsavePP(filename = "output//figs/FACE_RDA_EnvVar_Year1_3", plot = p, width = 6, height = 6)


#########
## PFG ##
#########

#########
## CO2 ##
#########
summary(lm(FloorPAR ~ year, data = peDF)) 
summary(lm(temp ~ year, data = peDF)) 
cor(peDF$FloorPAR, peDF$temp)
# significant temporal change in FloorPAR and temperature
plot(FloorPAR ~ year, data = peDF)
plot(temp ~ year, data = peDF)

#############
## Ambient ##
#############
ambDF_pfg <- subsetD(peDF, co2 == "amb")

# reorder by ring
ambDF_pfg <- ambDF_pfg[order(ambDF_pfg$ring), ]

# Run RDA
rda_amb <- rda(log(ambDF_pfg[, PFGName] + 1) ~ year + Condition(ring) , ambDF_pfg)

# Define permutation
hh <- allPerms(9, how(within = Within(type = "free"), plot = Plots(strata = ambDF_pfg$ring)))

anova(rda_amb, permutations = hh)
# no year effect

# Use temp and FloorPAR intstead
rda_amb2 <- rda(log(ambDF_pfg[, PFGName] + 1) ~ temp + FloorPAR + Condition(ring) , ambDF_pfg)
# model simplificaiton
rda_amb3 <- ordistep(rda_amb2, direction = "both", trace = 0, permutations = hh)
anova(rda_amb3, permutations = hh, by = "terms")
# significant tempearture effect

##########
## eCO2 ##
##########

elevDF_pfg <- subsetD(peDF, co2 == "elev")
elevDF_pfg <- elevDF_pfg[order(elevDF_pfg$ring), ]

rda_elev <- rda(log(elevDF_pfg[, PFGName] + 1) ~ year + Condition(ring) , elevDF_pfg)
hh <- allPerms(9, how(within = Within(type = "free"), plot = Plots(strata = elevDF_pfg$ring)))
anova(rda_elev, permutations = hh)
# no year effect

# fit FloorPAR and temp
rda_elev2 <- rda(log(elevDF_pfg[, PFGName] + 1) ~ FloorPAR + temp + Condition(ring) , elevDF_pfg)
anova(rda_elev2, permutations = hh, by = "terms")

# remove FloorPAR or temp
rda_elev3 <- rda(log(elevDF_pfg[, PFGName] + 1) ~ temp + Condition(ring) , elevDF_pfg)
rda_elev4 <- rda(log(elevDF_pfg[, PFGName] + 1) ~ FloorPAR + Condition(ring) , elevDF_pfg)
anova(rda_elev3, permutations = hh, by = "terms")
anova(rda_elev4, permutations = hh, by = "terms")
# no temp effect or FloorPAR effect

# Plot against year by co2
p2 <- PlotRDA_Year(rdaResLst = list(rda_amb, rda_elev), env = list(ambDF_pfg, elevDF_pfg), spscore = 0)
ggsavePP(filename = "output/figs/FACE_RDAvsYearbyCO2_PFG", plot = p2, 
         width = 6.65, height = 4)

#####################
## 3-year data set ## 
#####################
# From the above
#analysis, moist, temp, TotalC and Dry_soilph are determied to be imporatnt
#driver
rda_pfg_all <- rda(log(peDF[, PFGName] + 1) ~ TotalC + moist + Drysoil_ph + Depth_HL + temp, peDF)

# plot
p <- TriPlot(MultValRes = rda_pfg_all, env = peDF, yaxis = "RDA axis", axispos = c(1, 2, 3), 
             centcons = 2, spcons = .5, biplcons = 1, lowx = 0, lowy = 0)
ggsavePP(filename = "output//figs/FACE_RDA_EnvVar_PFG_Year1_3", plot = p, width = 6, height = 6)

##################
# Summary result #
##################
Res_list <- list(Res_Year1, Res_Year2, Res_Year3)
names(Res_list) <- c("Year1", "Year2", "Year3")
str(Res_Year1)

SummaryRDAanova <- function(x){
  data.frame(source = row.names(x), 
             DF = x$Df,
             Var = round(x$Variance*100/sum(x$Variance), 2),  
             F = round(x$F, 2), 
             P =  round(x$Pr, 3))
}

ResDF <- ldply(Res_list, SummaryRDAanova, .id = "Year")
write.csv(ResDF, file = "output/table/Result_RDA_Anova.csv", row.names = FALSE)
PfgResDF <- ldply(RdaPfgRes, SummaryRDAanova, .id = "Year")
write.csv(PfgResDF, file = "output/table/Result_RDA_Anova_PFG.csv", row.names = FALSE)


