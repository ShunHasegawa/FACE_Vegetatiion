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

##############
## 1st year ##
##############

# determine environmental variables driving initial difference in species
# composition using the 1st year data

names(EnvDF_3df)
df2013 <- subsetD(seDF, year == 2013)

rda1 <- rda(log(df2013[ , SppName] + 1) ~ TotalC + moist + Drysoil_ph + Depth_HL + FloorPAR, df2013)

# model simplification
rdaRes <- llply(list("backward", "forward", "both"), 
                function(x) ordistep(rda1, direction = x, trace = 0, 
                                     permutations = allPerms(6)))
  # possible permutation is 6! = 720 (small), so carry out exact permutation
  # test it says number of permutations 198 but actually 720 (I think.. and
  # don't know why it says 198)

SimplModAnv <- llply(rdaRes, function(x) anova(x, by = "margin", 
                                               permutations = allPerms(6)))
SimplModAnv
# TotalC and moist
rda2 <- rdaRes[[1]] # constrained axis explains 71% of variation
anova(rda2, permutations = allPerms(6))
Res_Year1 <- anova(rda2, permutations = allPerms(6), by = "terms")

# include co2
rda3 <- rda(log(df2013[ , SppName] + 1) ~ co2 + Condition(TotalC + moist), df2013)
anova(rda3, permutations = allPerms(6))
# no co2 effect

# plot
p <- TriPlot(MultValRes = rda2, env = subset(seDF, year == 2013), 
             yaxis = "RDA axis", axispos = c(1, 2, 3))
ggsavePP(filename = "output//figs/FACE_RDA_EnvVar_Year1", plot = p, width = 6, height = 6)
spenvcor(rda2)

##############
## 2nd year ##
##############
df2014 <- subsetD(seDF, year == 2014)
rda1 <- rda(log(df2014[ , SppName] + 1) ~ TotalC + moist + Drysoil_ph + Depth_HL, df2014)
rda1   # Constrained (canonical) axes explains 84%
anova(rda1, permutations = allPerms(6))  # 6! = 720
anova(rda1, permutations = allPerms(6), by = "terms")  # 6! = 720
# no significant association

# model simplification
rdaRes <- llply(list("backward", "forward", "both"), 
                function(x) ordistep(rda1, direction = x, trace = 0, 
                                     permutations = allPerms(6), 
                                     Pin = .1, Pout = .11))
SimplModAnv <- llply(rdaRes, function(x) anova(x, by = "margin", 
                                               permutations = allPerms(6)))
SimplModAnv
# moist and Drysoil_ph is marginally significant

rda2 <- rdaRes[[1]]
anova(rda2, permutations =0, allPerms(6))
Res_Year2 <- anova(rda2, permutations = allPerms(6), by = "margin")

# include co2
rda3 <- rda(log(df2014[ , SppName] + 1) ~ co2 + Condition(moist + Drysoil_ph) , df2014)
anova(rda3, permutations = allPerms(6))
# no

# plot
p <- TriPlot(MultValRes = rda2, env = df2014, yaxis = "RDA axis", axispos = c(1, 2, 3))
ggsavePP(filename = "output//figs/FACE_RDA_EnvVar_Year2", plot = p, 
         width = 6, height = 6)

##############
## 3rd year ##
##############
df2015 <- subsetD(seDF, year == 2015)
rda1 <- rda(log(df2015[ , SppName] + 1) ~ TotalC + moist + Drysoil_ph + Depth_HL, df2015)
rda1   # Constrained (canonical) axes explains 89%
anova(rda1, permutations = allPerms(6))  # marginaly significant
anova(rda1, permutations = allPerms(6), by = "terms")  # 6! = 720
# no significant association

# model simplification
rdaRes <- llply(list("backward", "forward", "both"), 
                function(x) ordistep(rda1, direction = x, trace = 0))
SimplModAnv <- llply(rdaRes, function(x) anova(x, by = "margin", 
                                               permutations = allPerms(6)))
SimplModAnv
# TotalC, moist and Drysoil_ph
rda2 <- rdaRes[[1]] 
rda2 # 78%
Res_Year3 <- anova(rda2, permutations = allPerms(6), by = "terms")  # 6! = 720

# include co2
rda3 <- rda(log(df2015[ , SppName] + 1) ~ co2 + Condition(TotalC + moist + Drysoil_ph), df2015)
rda3
anova(rda3, permutations = allPerms(6))  
# no co2 effect

# plot
p <- TriPlot(MultValRes = rda2, env = df2015, yaxis = "RDA axis", axispos = c(1, 2, 3))
ggsavePP(filename = "output//figs/FACE_RDA_EnvVar_Year3", plot = p, width = 6, height = 6)

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

#######
# PFG #
#######
peDF <- merge(RingSumPFGMatrix, EnvDF_3df, by = c("year", "ring", "block", "co2")) 

simpleRDAFun <- function(x) {
  dd <- x
  rda1 <- rda(log(dd[ , PFGName] + 1) ~ TotalC + moist + Drysoil_ph + Depth_HL, data = dd)
  rda2 <- ordistep(rda1, direction = "backward", trace = 0, permutations = allPerms(6))
  return(rda2)
}

dd1 <- subsetD(peDF, year == 2013)
dd2 <- subset(peDF, year == 2014)
dd3 <- subset(peDF, year == 2015)

rda1 <- rda(log(dd1[ , PFGName] + 1) ~ TotalC + moist + Drysoil_ph + Depth_HL, data = dd1)
rda2 <- rda(log(dd2[ , PFGName] + 1) ~ TotalC + moist + Drysoil_ph + Depth_HL, data = dd2)
rda3 <- rda(log(dd3[ , PFGName] + 1) ~ TotalC + moist + Drysoil_ph + Depth_HL, data = dd3)

rda1_2 <- ordistep(rda1, direction = "both", trace = 0, permutations = allPerms(6), Pin = .1, Pout = .11)
rda2_2 <- ordistep(rda2, direction = "both", trace = 0, permutations = allPerms(6), Pin = .1, Pout = .11)
rda3_2 <- ordistep(rda3, direction = "both", trace = 0, permutations = allPerms(6), Pin = .1, Pout = .11)

peDF_res <- list(rda1_2, rda2_2, rda3_2)
llply(peDF_res, function(x) anova(x, permutations = allPerms(6)))
RdaPfgRes <- llply(peDF_res, function(x) anova(x, permutations = allPerms(6), by = "margin"))
RdaPfgRes

# TotalC, moist, Drysoil_ph and Depth_HL

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


#######################################
# investigate environmental variables #
#######################################

# define permutation
seDF <- seDF[order(seDF$ring, seDF$year), ] # reorder by ring
ctrl <- allPerms(nrow(seDF), 
                 how(within = Within(type = "series", constant = TRUE),
                     # keep order (Year1-3) to take autocorrelation into account
                     plot = Plots(strata = seDF$ring, type = "free")) )
                          # ring is exchangeable.

# check if permutation is proparly designed
seDF[, c("year", "ring")]
some(ctrl)
# correct as each ring is exchanged and the order 1-3, 4-6, 7-9 etc. are kept


# use temp as condition as temp explains yearly variation
rda3 <- rda(log(seDF[, SppName] + 1) ~ TotalC + moist + Drysoil_ph + Depth_HL + Condition(temp), seDF)

anova(rda3, permutations = ctrl)
# marginally significant
anova(rda3, by = "margin", permutations = ctrl)

# model simplification
rdaRes <- llply(list("backward", "forward", "both"), 
                function(x) ordistep(rda3, direction = x, trace = 0, 
                                     permutations = ctrl))
SimplModAnv <- llply(rdaRes, function(x) anova(x, by = "margin", permutations = ctrl))
SimplModAnv

# moist, totalC and Drysoil_ph
rda4 <- rdaRes[[1]]
anova(rda4, permutations = ctrl)
anova(rda4, by = "margin", permutations = ctrl)

# plot
p <- TriPlot(MultValRes = rda4, env = seDF, yaxis = "RDA axis", axispos = c(1, 2, 3))
ggsavePP(filename = "output//figs/FACE_RDA_EnvVar_AllYear", plot = p, width = 6, height = 6)



## ---- RunCap

#########
## CAP ##
#########
# bray-curtis is used to compute dissimilarity
Cap1 <- capscale(log(RingSumVeg[, SppName] + 1) ~ TotalC +  Depth_HL + 
                   moist + temp + FloorPAR, EnvDF_3df, dist = "bray")

Cap1 <- capscale(log(RingSumVeg[, SppName] + 1) ~ TotalC +  moist, 
                 EnvDF_3df, dist = "bray")
Cap1
# 66 % is explained by constrained axes

anova(Cap1, permutations = ctrl) # significant association between plant community and environmental variables
anova(rda1, permutations = ctrl) # significant association between plant community and environmental variables
anova(Cap1, by = "axis") # CAP1-4 are significant

theme_set(theme_bw())
TriPlot(Cap1, env = EnvDF_3df, yaxis = "CAP", axispos = c(1:3), biplcons = 1, lowx = .2, lowy = .2)

# model simplification; remove each term one by one backwards, forwards and both
anova(Cap1, by = "margin")
capRes <- llply(list("backward", "forward", "both"), 
                function(x) ordistep(Cap1, direction = x, trace = 0))
SimplModAnv <- llply(capRes, function(x) anova(x, by = "margin"))
SimplModAnv
# consistent results

cap2 <- capRes[[3]]

p <- TriPlot(cap2, env = EnvDF_3df, yaxis = "CAP", axispos = c(1:3), biplcons = 1,
        lowx = .2, lowy = .2)
p
ggsavePP(filename = "output//figs/FACE_CAP_EnvVar", plot = p, width = 6, height = 6)

## ---- RunPartialCap

# Partial CAP Environmental variables have significant association with plant 
# community. Adjacent FACE rings (Block) seems to have similar plant community.
# So in order to take between-block variation and run partial CAP

pcap1 <- capscale(log(RingSumVeg[, SppName] + 1) ~ year + co2 + Condition(block), EnvDF_3df, dist = "bray")
pcap1
# block explains 55 % of variation
anova(pcap1)
anova(pcap1, by = "axis") # only 1st axis is significant
anova(pcap1, by = "margin")
p <- TriPlot(MultValRes = pcap1, env = EnvDF_3df, yaxis = "CAP", axispos = c(1, 2, 4), biplcons = 1, 
             lowx = .1, lowy = .1, EnvNumeric = FALSE)
p
ggsavePP(filename = "output//figs/FACE_PartialCAP_EnvVar", plot = p, width = 6, height = 6)

## ---- RunCAPbyCO2_AllSpp
CapByCo2 <- llply(list(amb = "amb", elev = "elev"), function(x) {
  capscale(log(RingSumVeg[, SppName] + 1) ~ year + Condition(ring), 
           EnvDF_3df, dist = "bray", subset = EnvDF_3df$co2 == x)
})
CapByCo2
llply(CapByCo2, anova)
# There was significant association between year and plant community at both CO2
# treatments

## ---- RunCAPbyCO2_AllSpp_plot
SmmryCapByCO2 <- llply(CapByCo2, summary)

# score for spp, site and centroids
CapScoreList <- llply(CapByCo2, vegan::scores)

# combine the above scores for each of co2 treatments
CapScores <- llply(names(CapScoreList$amb), function(x) {
  tdf <- ldply(list(CapScoreList$amb[[x]], CapScoreList$elev[[x]]), 
               function(y) data.frame(y, treatment = row.names(y)))
  tdf$co2 <- rep(c("amb", "elev"), each = nrow(tdf)/2)
  return(tdf)
  })
names(CapScores) <- names(CapScoreList$amb)

CapScores$species$year <- factor("2013")
CapScores$centroids$year <- factor("2013")
CapScores$sites <- cbind(CapScores$sites, 
                         ldply(c("amb", "elev"), function(x) subset(EnvDF_3df, co2 == x)))

p <- ggplot(data = CapScores$sites, aes(x = CAP1, y = CAP2, shape = year))
p2 <- p + geom_point(size = 3) + facet_grid(co2 ~.)
p3 <- p2 + 
  geom_segment(data = CapScores$centroids, 
                   aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
                   arrow = arrow(length = unit(.1, "cm")), 
                   alpha = .6,
                   color = "blue") + 
      geom_text(data = CapScores$centroids, 
                aes(x = CAP1 * 1.1 , y = CAP2 * 1.1, label = treatment), 
                alpha = .6, lineheight = .7, 
                color = "blue", size = 2, 
                fontface = "bold")
spdf <- subset(CapScores$species, abs(CAP1) > .08| abs(CAP2) > .08)
spdf <- within(spdf, {
  treatment <- gsub("[.]", "\n", as.character(treatment))
  CAP1 <- CAP1 * 2.5
  CAP2 <- CAP2 * 2.5
  })

p4 <- p3 + 
  geom_segment(data = spdf, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(.1, "cm")),
               alpha = .6,
               colour = "red") +
    geom_text(data = spdf, 
              aes(x = CAP1 * 1.2, y = CAP2 * 1.2, label = treatment),
              size = 2, fontface = "bold.italic", 
              colour = "red", alpha = .6, lineheight = .7) +
    labs(x = "CAP1", y = "CAP2") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  geom_vline(aes(xintercept = 0), linetype = "dotted")
p4
ggsavePP(filename = "output/figs/FACE_CAPvsYear_byCO2", plot = p4, width = 6, height = 6)

## ---- RunCAPbyCO2_PFG

#######
# PFG #
#######
cap_pfg1 <- capscale(log(RingSumPFGMatrix[, PFGName] + 1) ~ OrganicMatter + moist + temp + 
                       FloorPAR + Depth_HL, EnvDF_3df, dist = "bray")

anova(cap_pfg1)
anova(cap_pfg1, by = "axis")
anova(cap_pfg1, by = "margin")

capRes_pfg <- llply(list("backward", "forward", "both"), 
                    function(x) ordistep(cap_pfg1, direction = x, trace = 0))

SimplModAnv_pfg <- llply(capRes_pfg, function(x) anova(x, by = "margin"))
cap_pfg2 <- capRes_pfg[[3]]

p <- TriPlot(cap_pfg2, yaxis = "CAP", env = EnvDF_3df, axispos = c(1:3),
             lowx = 0, lowy = 0, spcons = .7, biplcons = .7)
p
ggsavePP(filename = "output/figs/FACE_CAP_EnvVar_PFG", plot = p, width = 6, height = 6)


## ---- others
# Partial CAP with PFG
cap_pfg2 <- capscale(log(RingSumPFGMatrix[, PFGName] + 1) ~ co2 + year + Condition(block),
                     EnvDF_3df, dist = "bray")
cap_pfg3 <- capscale(log(RingSumPFGMatrix[, PFGName] + 1) ~ co2 + year,
                     EnvDF_3df, dist = "bray")
cap_pfg4 <- capscale(log(RingSumPFGMatrix[, PFGName] + 1) ~ year + Condition(ring),
                     EnvDF_3df, dist = "bray")
anova(cap_pfg2)
anova(cap_pfg3)
anova(cap_pfg4)
anova(cap_pfg4, by = "axis")
plot(cap_pfg4)


CapByCo2_pfg <- llply(list(amb = "amb", elev = "elev"), function(x) {
  capscale(log(RingSumPFGMatrix[, PFGName] + 1) ~ year + Condition(ring), 
           EnvDF_3df, dist = "bray", subset = EnvDF_3df$co2 == x)
  })

llply(CapByCo2_pfg, anova)
par(mfrow = c(1, 2))
llply(CapByCo2_pfg, function(x) plot(x, scaling = 1))
# no evidence of association


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
