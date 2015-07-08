###############
# All species #
###############

########################
# Each year separately #
########################

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
names(adjR_singl_Lst) <- paste0("Year", 0:2)

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
fmls <- llply(list(Year0 = 1, Year1 = 2, Year2 = 3), 
              function(x) llply(paste(LH[[x]], FullFormula[[x]]), as.formula))

##############
## 1st year ##
##############
df2013 <- subsetD(seDF, year == "Year0")

# There are too many environmental variables to fit. so choose four which showed
# highest R2adj 
# adjusted R2
adjR <- ldply(fmls$Year0, function(x) RsquareAdj(rda(x, data = df2013))$adj.r.squared)

# highest R2
rr <- rda(fmls$Year0[[which(max(adjR) == adjR)]], df2013)

# check multicollinearity
vif.cca(rr)
anova(rr, permutations = allPerms(6))
rr2 <- rda(log(df2013[ , SppName] + 1) ~ 1, df2013)
rr3 <- ordiR2step(rr2, rr, permutations = allPerms(6), direction = "forward", Pin = .1)
summary(rr3)

# summary result
rda2013 <- list(IniRda = rr, FinRda = rr3)

##############
## 2nd year ##
##############
df2014 <- subsetD(seDF, year == "Year1")

# adjusted R2
adjR <- laply(fmls$Year1, function(x) RsquareAdj(rda(x, data = df2014))$adj.r.squared)

# highest R2
rr <- rda(fmls$Year1[[which(max(adjR) == adjR)]], df2014)

# check multicolliniarity using vif
vif.cca(rr)
anova(rr, permutations = allPerms(6))
# not significant

# choose only three terms
comb_exp <- combn(PosAdjR$Year1, 3)
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
df2015 <- subsetD(seDF, year == "Year2")

# adjusted R2
adjR <- laply(fmls$Year2, function(x) RsquareAdj(rda(x, data = df2015))$adj.r.squared)

# highest R2
rr <- rda(fmls$Year2[[which(max(adjR) == adjR)]], df2015)

# check multicollinearity
vif.cca(rr)
# TotalC, Depth_HL >10. 

# Second highset R2adj
adjR[which(max(adjR) == adjR)] <- NA
rr <- rda(fmls$Year2[[which(max(adjR, na.rm = TRUE) == adjR)]], df2015)
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
RdaLst <- list(Year0 = rda2013, Year1 = rda2014, Year2 = rda2015)
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
llply(RdaLst, function(x) anova(x$IniRda, permutations = allPerms(6), by = "term"))
llply(RdaLst, function(x) anova(x$FinRda, permutations = allPerms(6)))
llply(RdaLst, function(x) anova(x$FinRda, permutations = allPerms(6), by = "margin"))

#################
# CO2 treatment #
#################

#############
## Ambient ##
#############
ambDF <- subsetD(seDF, co2 == "amb")

# permutaiton Possible number of permutation is (3!)^3 = 216. it's small so need
# to use exact permutation test

allPerms(9, 
         how(within = Within(type = "free"), plot = Plots(strata = ambDF$ring)))
# Null hypothesis is no difference between years. Ring is no exchangeable but 
# subplots are. NOTE # First three colmuns of the above indices are for the 1st 
# factor of ring, and four to six are for the 2nd and seven to nine are for 3rd.
# So that data frame needs to be reordered to match this. 

# reorder by ring
ambDF <- ambDF[order(ambDF$ring), ]

# Run RDA
rda_amb <- rda(log(ambDF[, SppName] + 1) ~ year + Condition(ring) , ambDF)

hh <- allPerms(9, how(within = Within(type = "free"), plot = Plots(strata = ambDF$ring)))

anova(rda_amb, permutations = hh)
# significant year effect

# Axis
anova(rda_amb, permutations = hh, by = "axis")
# RDA1 is significant

summary(rda_amb)

#########
## co2 ##
#########

elvDF <- subsetD(seDF, co2 == "elev")
# reorder by ring
elvDF <- elvDF[order(elvDF$ring), ]

# Run RDA
rda_elev <- rda(log(elvDF[, SppName] + 1) ~ year + Condition(ring) , elvDF)

hh <- allPerms(9, how(within = Within(type = "free"), plot = Plots(strata = elvDF$ring)))

anova(rda_elev, permutations = hh)
# significant year effect
anova(rda_elev, permutations = hh, by = "axis")
# only RDA1 is significant

############
# Summary  #
############
RdaResLst <- list(amb = rda_amb, elev = rda_elev)

# Driving spp
SpScores <- ldply(RdaResLst, function(x) {
  sp <- vegan::scores(x)$species
  data.frame(sp, species = row.names(sp))
}, .id = "co2")

drvSpp <- unique(SpScores$specie[abs(SpScores$RDA1) > .4])
drvSppDF <- subset(SpScores, species %in% drvSpp)

# add covrage column and PFG
drvSppDF <- Reduce(function(x, y) {
  merge(x, y, by.x = "species", by.y = "variable", all.x = TRUE)},
  list(drvSppDF, SppSum[, c("variable", "Cov")], unique(veg[, c("variable", "PFG")])))
  
drvSppDF$species <- with(drvSppDF, 
                         paste0(gsub("[.]", " ", as.character(species)), 
                                "(", format(Cov, digit = 1, nsmall = 1), ")"))
drvSppDF_cst <- dcast(species + PFG ~ co2, value.var = "RDA1", data = drvSppDF)
drvSppDF_cst <- drvSppDF_cst[order(drvSppDF_cst$amb, decreasing = TRUE), ]
drvSppDF_cst <- within(drvSppDF_cst, {
  amb <- format(amb, digit = 2, nsmall = 2)
  elev <- format(elev, digit = 2, nsmall = 2)
})
write.csv(drvSppDF_cst, file = "output/table/RDA_AllSpp_SpScore.csv", row.names = FALSE)

########
# Plot #
########
p2 <- PlotRDA_Year(rdaResLst = list(rda_amb, rda_elev), env = list(ambDF, elvDF), 
                   spscore = .4)
ggsavePP(filename = "output/figs/FACE_RDAvsYearbyCO2", plot = p2, 
         width = 6.65, height = 4)

## Plot for thesis ##
RdaResLst
envDFLst <- list(amb = ambDF, elev = elvDF)

SummaryRda <- llply(RdaResLst, summary)

# Site score
siteDD <- ldply(c("amb", "elev"), function(x) data.frame(SummaryRda[[x]]$site, envDFLst[[x]]))
siteDD$year <- factor(siteDD$year, labels = paste0("Year", 0:2))

# create a plot
p <- ggplot(siteDD, aes(x = year, y = RDA1))
SiteScorePlot <- p + 
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~co2, scale = "free_y") +
  labs(x = NULL, y = "RDA1") +
  science_theme
SiteScorePlot
# change facet_label

# % variance for RDA1
Rda1Prop <- laply(SummaryRda, function(x) 
  format(x$cont$importance["Proportion Explained", "RDA1"] * 100, 
         digits = 2, nsmall = 2)
)

labs <- c(paste0("Ambient (", Rda1Prop[1], "%)"), 
          parse(text = paste("eCO[2]~(", Rda1Prop[2], "*'%')")))
Rda_Year_AllSp <- facet_wrap_labeller(SiteScorePlot, labels = labs)

ggsavePP(Rda_Year_AllSp, filename = "output/figs/Fig_Thesis/RDAvsYearbyCO2_AllSpp", 
         width = 4, height = 3)

#####################
## 3-year data set ## 
#####################
# From the above
#analysis, moist, temp, TotalC and Dry_soilph are determied to be imporatnt
#driver

rda_all <- rda(log(seDF[, SppName] + 1) ~ TotalC + moist +year, seDF)

# can't run anova as it is. cause different perumutation units need to be
# defined for year and co2. so anyway create a triplot and see the pattern.
rda_all

# plot
p <- TriPlot(MultValRes = rda_all, env = seDF, yaxis = "RDA axis", axispos = c(1, 2, 3), centcons = 2)
ggsavePP(filename = "output//figs/FACE_RDA_EnvVar_Year0_3", plot = p, width = 6, height = 6)

##################
# Fig for thesis #
##################

#############
## All spp ##
#############
RdaAllRes <- summary(rda_all)
seDF$year <- factor(seDF$year, labels = paste0("Year", 0:2))
sitedd <- data.frame(RdaAllRes$site, seDF)

centdd <- data.frame(RdaAllRes$centroids, co2 = "amb", year = "Year0", 
                     variable = row.names(RdaAllRes$centroids))
centdd$variable <- factor(centdd$variable, labels = paste0("Year", 0:2))

bipldd <- data.frame(RdaAllRes$biplot, co2 = "amb", year = "Year0", 
                     variable = row.names(RdaAllRes$biplot))
bipldd <- subsetD(bipldd, variable %in% c("TotalC", "moist"))
bipldd$variable <- factor(bipldd$variable, labels = c("Moist", "Total C"))

VarProp <- RdaAllRes$cont$importance["Eigenvalue",]/RdaAllRes$tot.chi
axislabs <- paste0(c("RDA1", "RDA2"), "(", round(VarProp[c(1, 2)] * 100, 2), "%)")

# make a plot
p <- ggplot(data = sitedd, aes(x = RDA1, y = RDA2, shape = year))
p2 <- p + 
  geom_path(aes(group = ring), col = "black") +
  geom_point(aes(fill = co2), size = 4) + 
  #   geom_text(aes(x = RDA1 + 0.1,label = ring)) +
  scale_fill_manual(values = c("black", "white"), 
                    labels = c("Ambient", expression(eCO[2])),
                    guide = guide_legend(override.aes = list(shape = 21))) +
  # need to add overrisde here to make white circle in the legend with shape of
  # 21
  scale_shape_manual(values = c(21, 22, 24)) + 
  geom_segment(data = bipldd,
               aes(x = 0, y = 0, xend = RDA1 * 2, yend = RDA2 * 2), 
               arrow = arrow(length = unit(.2, "cm")), 
               color = "red") +
  geom_text(data = bipldd, 
            aes(x = RDA1 * 2.3 , y = RDA2 * 2.3, label = variable), 
            alpha = .6, lineheight = .7, 
            color = "red", size = 4, 
            fontface = "bold") +
  geom_segment(data = centdd,
               aes(x = 0, y = 0, xend = RDA1 * 2, yend = RDA2 * 2), 
               arrow = arrow(length = unit(.2, "cm")), 
               color = "blue") +
  geom_text(data = centdd, 
            aes(x = RDA1 * 2.3 , y = RDA2 * 2.3, label = variable), 
            alpha = .6, lineheight = .7, 
            color = "blue", size = 4, 
            fontface = "bold") +
  geom_hline(xintercept = 0, linetype = "dashed") +
  geom_vline(yintercept = 0, linetype = "dashed") +
  science_theme +
  theme(legend.position = c(.17, .15), 
        legend.box = "horizontal", 
        legend.box.just = "top") +
  labs(x = axislabs[1], y = axislabs[2])
RDA_Plot_AllSpp <- p2
RDA_Plot_AllSpp
ggsavePP(plot = RDA_Plot_AllSpp, filename = "output/figs/Fig_Thesis/RDA_3yr_AllSpp", width = 6, height = 4)
