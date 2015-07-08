#######
# PFG #
#######

########################
# Each year separately #
########################

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
names(adjR_singl_Lst) <- paste0("Year", 0:2)

# Get variables with positive R2adj
PosAdjR <- llply(adjR_singl_Lst, function(x) as.character(x$variable[x$V1 > 0]))
# less than 4

# Formula for full models
FullFormula <- llply(PosAdjR, function(x) paste(x, collapse = "+"))

# Y matric (lefthand part)
LH <- list("log(df2013[ , PFGName] + 1) ~",
           "log(df2014[ , PFGName] + 1) ~",
           "log(df2015[ , PFGName] + 1) ~")
fmls <- llply(list(Year0 = 1, Year1 = 2, Year2 = 3), 
              function(x) llply(paste(LH[[x]], FullFormula[[x]]), as.formula))

##############
## 1st year ##
##############
df2013 <- subsetD(peDF, year == "Year0")

# adjusted R2
adjR <- laply(fmls$Year0, function(x) RsquareAdj(rda(x, data = df2013))$adj.r.squared)

# highest R2
rr <- rda(fmls$Year0[[which(max(adjR) == adjR)]], df2013)
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
df2014 <- subsetD(peDF, year == "Year1")

# adjusted R2
adjR <- ldply(fmls$Year1, function(x) RsquareAdj(rda(x, data = df2014))$adj.r.squared)

# highest R2
rr <- rda(fmls$Year1[[which(max(adjR) == adjR)]], df2014)
anova(rr, permutations = allPerms(6))
rr2 <- rda(log(df2014[ , PFGName] + 1) ~ 1, df2014)
rr3 <- ordiR2step(rr2, rr, permutations = allPerms(6), direction = "forward", Pin = .1)

# summary result
rda2014 <- list(IniRda = rr, FinRda = rr3)

##############
## 3rd year ##
##############
df2015 <- subsetD(peDF, year == "Year2")

# adjusted R2
adjR <- laply(fmls$Year2, function(x) RsquareAdj(rda(x, data = df2015))$adj.r.squared)

# highest R2
rr <- rda(fmls$Year2[[which(max(adjR) == adjR)]], df2015)
# check multicollinearity
vif.cca(rr)
anova(rr, permutations = allPerms(6))
# not significant so use two terms instead

# try only two variables
comb_exp <- combn(PosAdjR$Year2, 2) 
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
RdaLst_pfg <- list(Year0 = rda2013, Year1 = rda2014, Year2 = rda2015)

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
hh <- allPerms(9, how(within = Within(type = "free"), 
                      plot = Plots(strata = ambDF_pfg$ring)))
anova(rda_amb, permutations = hh)
# no year effect

# Use temp and FloorPAR intstead
rda_amb2 <- rda(log(ambDF_pfg[, PFGName] + 1) ~ temp + Condition(ring) , ambDF_pfg)
rda_amb3 <- rda(log(ambDF_pfg[, PFGName] + 1) ~ FloorPAR + Condition(ring) , ambDF_pfg)
rda_amb4 <- rda(log(ambDF_pfg[, PFGName] + 1) ~ FloorPAR + temp + Condition(ring) , 
                ambDF_pfg)
vif.cca(rda_amb2)
vif.cca(rda_amb3)
vif.cca(rda_amb4)

RsquareAdj(rda_amb2)$adj.r.squared
RsquareAdj(rda_amb3)$adj.r.squared
RsquareAdj(rda_amb4)$adj.r.squared
anova(rda_amb2, permutations = hh, by = "margin")
anova(rda_amb3, permutations = hh, by = "margin")
anova(rda_amb4, permutations = hh, by = "margin")
# fitting two variables at the same time ends up no significant association so
# don't use
# floorpar has slightly higher R2
anova(rda_amb3, permutations = hh, by = "margin")

# significant tempearture or FloorPAR effect

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
rda_elev3 <- rda(log(elevDF_pfg[, PFGName] + 1) ~ temp + Condition(ring) , elevDF_pfg)
rda_elev4 <- rda(log(elevDF_pfg[, PFGName] + 1) ~ FloorPAR + Condition(ring) , elevDF_pfg)
anova(rda_elev2, permutations = hh, by = "terms")
anova(rda_elev3, permutations = hh, by = "terms")
anova(rda_elev4, permutations = hh, by = "terms")
# nothing is significant

# no temp effect or FloorPAR effect

# Plot against year by co2
p2 <- PlotRDA_Year(rdaResLst = list(rda_amb, rda_elev), env = list(ambDF_pfg, elevDF_pfg), spscore = 0)
ggsavePP(filename = "output/figs/FACE_RDAvsYearbyCO2_PFG", plot = p2, 
         width = 6.65, height = 4)

## Plot for thesis ##
RdaResPfgLst <- list(amb = rda_amb, elev = rda_elev)
envDFLst <- list(amb = ambDF_pfg, elev = elevDF_pfg)

SummaryRdaPfg <- llply(RdaResPfgLst, summary)

# Site score
siteDD <- ldply(c("amb", "elev"), function(x) data.frame(SummaryRdaPfg[[x]]$site, envDFLst[[x]]))
siteDD$year <- factor(siteDD$year, labels = paste0("Year", 0:2))

# create a plot

# Site score plot
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
Rda1Prop <- laply(SummaryRdaPfg, function(x) 
  format(x$cont$importance["Proportion Explained", "RDA1"] * 100, 
         digits = 2, nsmall = 2)
  )

labs <- c(paste0("Ambient (", Rda1Prop[1], "%)"), 
          parse(text = paste("eCO[2]~(", Rda1Prop[2], "*'%')")))
Rda_Year_PFG <- facet_wrap_labeller(SiteScorePlot, labels = labs)
ggsavePP(Rda_Year_PFG, filename = "output/figs/Fig_Thesis/RDAvsYearbyCO2_PFG", 
         width = 4.5, height = 2.5)

# Species score Table
sppdd <- ldply(SummaryRdaPfg, function(x) 
  data.frame(x$species, variable = row.names(x$species)), 
  .id = "co2")
sppdd$RDA1 <- round(sppdd$RDA1, 2)
sppdd_cst <- dcast(variable ~ co2, value.var = "RDA1", data = sppdd)
sppdd_cst <- sppdd_cst[order(sppdd_cst$amb, decreasing = TRUE), ]
write.csv(sppdd_cst, file = "output/table/RDA_PFG_SpScore.csv", row.names = FALSE)

#####################
## 3-year data set ## 
#####################
# From the above
#analysis, moist, temp, TotalC and Dry_soilph are determied to be imporatnt
#driver
rda_pfg_all <- rda(log(peDF[, PFGName] + 1) ~ Drysoil_ph, peDF)

# plot
p <- TriPlot(MultValRes = rda_pfg_all, env = peDF, yaxis = "RDA axis", axispos = c(1, 2, 3), 
             centcons = 2, spcons = .5, biplcons = 1, lowx = 0, lowy = 0)
ggsavePP(filename = "output//figs/FACE_RDA_EnvVar_PFG_Year0_3", plot = p, width = 6, height = 6)

# Fig for thesis

RdaAllRes <- summary(rda_pfg_all)
peDF$year <- factor(peDF$year, labels = paste0("Year", 0:2))
sitedd <- data.frame(RdaAllRes$site, peDF)

sppdd <- data.frame(RdaAllRes$species, year = "Year0", 
                    PFG = factor(row.names(RdaAllRes$species)))
sppdd$PFG <- factor(sppdd$PFG, 
                    labels = c("Non~legume", 
                               "C[3*'\u005F'*grass]","C[4*'\u005F'*grass]",
                               "Legume", "Moss",  "Woody~plants"))

bipldd <- data.frame(RdaAllRes$biplot, co2 = "amb", year = "Year0", 
                     variable = row.names(RdaAllRes$biplot))
bipldd$variable <- factor(bipldd$variable, labels = c("Soil pH"))
VarProp <- RdaAllRes$cont$importance["Eigenvalue",]/RdaAllRes$tot.chi
axislabs <- paste0(c("RDA1", "PC1"), "(", round(VarProp[c(1, 2)] * 100, 2), "%)")

# make a plot
p <- ggplot(data = sitedd, aes(x = RDA1, y = PC1, shape = year))
p2 <- p + 
  geom_path(aes(group = ring), col = "black") +
  geom_point(aes(fill = co2), size = 4) + 
  scale_fill_manual(values = c("black", "white"), 
                    labels = c("Ambient", expression(eCO[2])),
                    guide = guide_legend(override.aes = list(shape = 21))) +
  # need to add overrisde here to make white circle in the legend with shape of
  # 21
  scale_shape_manual(values = c(21, 22, 24)) + 
  # species score
  geom_segment(data = sppdd,
               aes(x = 0, y = 0, xend = RDA1 * .9, yend = PC1 * .9), 
               arrow = arrow(length = unit(.2, "cm")), 
               color = "darkgreen") +
  geom_text(data = sppdd, 
            aes(x = RDA1 * 1, y = PC1 * 1, label = PFG), 
            alpha = .6, lineheight = .7, 
            color = "darkgreen", size = 4, 
            fontface = "bold", 
            parse = TRUE) +
  # environmental variable
  geom_segment(data = bipldd,
               aes(x = 0, y = 0, xend = RDA1 * 1.3, yend = PC1 * 1.3), 
               arrow = arrow(length = unit(.2, "cm")), 
               color = "red") +
  geom_text(data = bipldd, 
            aes(x = RDA1 * 1.5 , y = PC1 * 1.5, label = variable), 
            alpha = .6, lineheight = .7, 
            color = "red", size = 4, 
            fontface = "bold") +
  geom_hline(xintercept = 0, linetype = "dashed") +
  geom_vline(yintercept = 0, linetype = "dashed") +
  science_theme+ 
  theme( legend.key = element_blank(),
         legend.position = c(.8, .15), 
         legend.box = "horizontal", 
         legend.box.just = "top",
         legend.key.height = unit(1, "lines")) +
  labs(x = axislabs[1], y = axislabs[2])
RDA_Plot_PFG <- p2
RDA_Plot_PFG

ggsavePP(plot = RDA_Plot_PFG, filename = "output/figs/Fig_Thesis/RDA_3yr_PFG", width = 6, height = 4)

