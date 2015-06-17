head(veg)

#########################
# co2 effect (% change) #
#########################
# plant forms

# get C3, grass, forb etc. abundance fo each plot
PfgAbundDF <- ddply(veg, .(year, co2, block, ring), function(x){
  Total     = sum(x$value)
  C3    = sum(x$value[!x$PFG %in% c("c4", "moss")])
  C4grass   = sum(x$value[x$PFG == "c4"])
  Forb  = sum(x$value[x$form == "Forb"])
  Grass = sum(x$value[x$form == "Grass"])
  Wood  = sum(x$value[x$form == "Wood"])
  Moss  = sum(x$value[x$form == "Moss"])
  data.frame(Total, C3, C4grass, Forb, Grass, Wood, Moss)
})
PfgAbundDF_mlt <- melt(PfgAbundDF, id = c("year", "co2", "block", "ring", "Total"))
PfgAbundDF_mlt_II <- melt(PfgAbundDF_mlt, id = c("year", "co2", "block", "ring", "variable"), 
                          variable_name = "count")

# cast to pair each block
PfgAbundDF_cst <- dcast(PfgAbundDF_mlt_II, variable + year + block ~ co2 + count)

ratio <- function(d, w) {
  AmbM  <- with(d, sum(amb_value * w) / sum(amb_Total * w))
  elevM <- with(d, sum(elev_value * w)/ sum(elev_Total * w))
  elevM/AmbM- 1
}

RatioSE <- ddply(PfgAbundDF_cst, .(year, variable), function(x) {
  b <- boot::boot(x, ratio, R = 999, stype = "w")
  summary(b)
})
RatioSE$co2R <- RatioSE$original
RatioSE[RatioSE$variable == "C3", ]
# plot

# block mean
blockMean <- ddply(PfgAbundDF_mlt, .(year, variable, block), 
                   function(x) {
                     am <- with(x, value[co2 == "amb"] / Total[co2 == "amb"])
                     em <- with(x, value[co2 == "elev"]/ Total[co2 == "elev"])
                     co2R <- em/am - 1
                     data.frame(co2R) 
                   })

p <- ggplot(RatioSE, aes(x = year, y = co2R))
p2 <- p + 
  geom_point(aes(x = year, y = co2R), size = 4)+
  geom_errorbar(aes(x = year, ymin = co2R - bootSE, ymax = co2R + bootSE), width = 0) + 
  geom_point(data = blockMean, size = 2, col = "red", alpha = .7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "CO2 response ratio", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 7)) +
  facet_wrap( ~ variable, scale = "free_y")
ggsavePP(plot = p2, filename = "output/figs/FACE_CO2ResponseRatio_PFGfraction", width = 6, height = 4)

# organise it to export as a table
RatioSE$co2R_se <- with(RatioSE, paste0(round(co2R, 2), "(", round(bootSE, 2), ")"))

RatioSE_cst <- dcast(RatioSE[c("variable", "year", "co2R_se")], variable ~ year, value.var = "co2R_se")
write.csv(RatioSE_cst, file = "output/table/CO2ResponseRatio_PFGfraction.csv", row.names = FALSE)

############
# Analysis #
############

# get C3, grass, forb etc. proportion fo each plot
PropDF <- ddply(veg, .(year, co2, block, ring, plot, id), function(x){
  Total     = sum(x$value)
  C3prop    = 1-sum(x$value[x$PFG %in% c("c4", "moss")])/Total
  C4grass   = sum(x$value[x$PFG == "c4"])/Total
  Forbprop  = sum(x$value[x$form == "Forb"])/Total
  Grassprop = sum(x$value[x$form %in% c("Grass", "Sedge")])/Total
  Woodprop  = sum(x$value[x$form %in% c("Tree", "Shrub")])/Total
  Mossprop  = sum(x$value[x$form == "Moss"])/Total
  data.frame(Total, C3prop, C4grass, Forbprop, Grassprop, Woodprop, Mossprop)
  })
PropDF$obs <- 1:nrow(PropDF)

########################################
# C3 proportion in the whole community #
########################################

# glm
m1 <- glmer(C3prop ~ year * co2 +  (1|block) +(1|ring) + (1|id), 
            family = "binomial", data = PropDF, weight = Total)
overdisp.glmer(m1)
# overdispersed
# add obs
m2 <- update(m1, ~. + (1|obs))
overdisp.glmer(m2)
Anova(m2)
plot(m2)

# highly wedged..
qqnorm(resid(m2))
qqline(resid(m2))
# compare AIC
c3PropComAIC <- CompAIC(m2)
c3PropComAIC

# LMM
bxplts(value = "C3prop", xval = "co2", data = PropDF)
bxplts(value = "C3prop", xval = "ring", data = PropDF)

par(mfrow = c(1, 3))
boxplot(C3prop ~ year:ring, data = PropDF, main = "raw")
boxplot(asin(C3prop) ~ year:ring, data = PropDF, main = "arcsin")
boxplot(logit(C3prop) ~ year:ring, data = PropDF, main = "logit")

# try logit
rm1 <- lmer(logit(C3prop) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = PropDF)
summary(rm1)
Anova(rm1)
Anova(rm1, test.statistic = "F")
plot(rm1)
# slightly wedged pattern
qqnorm(resid(rm1))
qqline(resid(rm1))

# one complete outlier, so remove them
rmv <- which(qqnorm(resid(rm1))$y %in% min(qqnorm(resid(rm1))$y))
rm2 <- update(rm1, subset = -rmv)
Anova(rm2)
qqnorm(resid(rm2))
qqline(resid(rm2))

# improved a lot
AnvF_PropC3_total <- Anova(rm2, test.statistic = "F")
AnvF_PropC3_total

# contrast
# contrast doesn't work with lmer. so use lme
tdf <- PropDF[-rmv, ]
lmeMod <- lme(logit(C3prop) ~ year * co2, random = ~1|block/ring/id, data = tdf)
cntrst<- contrast(lmeMod, 
                  a = list(year = levels(tdf$year), co2 = "amb"),
                  b = list(year = levels(tdf$year), co2 = "elev"))
PropC3_Total_CntrstRes <- cntrstTbl(cntrst, data = tdf, variable = "C3Prop")
PropC3_Total_CntrstRes
# no significant difference between treatment for each year

# improve the outlier and try glmm again
m3 <- update(m1, subset = -rmv)
overdisp.glmer(m3)
m4 <- update(m3, ~. + (1|obs))
overdisp.glmer(m4)
Anova(m4)
# compare AIC
CompAIC(m4)

plot(m4)
qqnorm(resid(m4))
qqline(resid(m4))
summary(m4)

############
# C4 grass #
############
# glm
m1 <- glmer(C4grass ~ year * co2 +  (1|block) +(1|ring) + (1|id), 
            family = "binomial", data = PropDF, weight = Total)
overdisp.glmer(m1)
# overdispersed
# add obs
summary(m1)
m2 <- update(m1, ~. + (1|obs))
overdisp.glmer(m2)
Anova(m2)
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))
# Compare AIC
C4grassCompAic <- CompAIC(m2)
C4grassCompAic

# indication of co2xyear interaction

# LMM
boxplot(logit(C4grass) ~ year:ring, data = PropDF, main = "logit")
m2Lmer <- lmer(logit(C4grass) ~ year * co2 +  (1|block) +(1|ring) + (1|id), data = PropDF, weight = Total)
Anova(m2Lmer, test.statistic = "F")

plot(m2Lmer)
qqnorm(resid(m2Lmer))
qqline(resid(m2Lmer))
# one obvious outlier
rmv <- which(qqnorm(resid(m2Lmer))$y == max(qqnorm(resid(m2Lmer))$y))
m2Lmer2 <- update(m2Lmer, subset = -rmv)
plot(m2Lmer2)
qqnorm(resid(m2Lmer2))
qqline(resid(m2Lmer2))
AnvF_C4grass <- Anova(m2Lmer2, test.statistic = "F")
AnvF_C4grass

# update glm, remove the outlier
m3 <- update(m1, subset = -rmv)
Anova(m3)
summary(m3)
overdisp.glmer(m3)
# overdispersed
# add obs
m4 <- update(m3, ~. + (1|obs))
plot(m4)
qqnorm(resid(m4))
qqline(resid(m4))
C4grassCompAic <- CompAIC(m4)
C4grassCompAic

#################
# Grass + Sedge #
#################
m1 <- glmer(Grassprop ~ year * co2 +  (1|block) +(1|ring) + (1|id), 
            family = "binomial", data = PropDF, weight = Total)
overdisp.glmer(m1)
# overdispersed
# add obs
m2 <- update(m1, ~. + (1|obs))
overdisp.glmer(m2)
Anova(m2)
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))

# compare AIC
GrassPropCompAIC <- CompAIC(m2)
GrassPropCompAIC

# LMM
par(mfrow = c(1, 3))
boxplot(Grassprop ~ year:ring, data = PropDF, main = "raw")
boxplot(asin(Grassprop) ~ year:ring, data = PropDF, main = "arcsin")
boxplot(logit(Grassprop) ~ year:ring, data = PropDF, main = "logit")

# try logit
gm1 <- lmer(logit(Grassprop) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = PropDF)
plot(gm1)
qqnorm(resid(gm1))
qqline(resid(gm1))

Anova(gm1)
AnvF_Grass <- Anova(gm1, test.statistic = "F")
AnvF_Grass
# maybe not

########
# Forb #
########
m1 <- glmer(Forbprop ~ year * co2 +  (1|block) +(1|ring) + (1|id), 
            family = "binomial", data = PropDF, weight = Total)
overdisp.glmer(m1)
# overdispersed
# add obs
m2 <- update(m1, ~. + (1|obs))
overdisp.glmer(m2)
Anova(m2)
# compare AIC
ForbPropCompAIC <- CompAIC(m2)
ForbPropCompAIC
# year looks important

# LMM
bxplts(value = "Forbprop", xval = "co2", data = PropDF)
bxplts(value = "Forbprop", xval = "ring", data = PropDF)

par(mfrow = c(1, 3))
boxplot(Forbprop ~ year:ring, data = PropDF, main = "raw")
boxplot(asin(Forbprop) ~ year:ring, data = PropDF, main = "arcsin")
boxplot(logit(Forbprop) ~ year:ring, data = PropDF, main = "logit")
# try logit
fm1 <- lmer(logit(Forbprop) ~ year * co2 + (1|block) +  (1|ring) + (1|id), data = PropDF)
Anova(fm1)
Anova(fm1, test.statistic = "F")

# model simplification
fm2 <- stepLmer(fm1)
AnvF_forb <- Anova(fm1, test.statistic = "F")
AnvF_forb

# model diagnosis
plot(fm2)
qqnorm(resid(fm2))
qqline(resid(fm2))

# what if I remove two outliers
rmv <- which(qqnorm(resid(fm2), plot.it = FALSE)$y %in% 
               sort(qqnorm(resid(fm2), plot.it = FALSE)$y)[1:2])

fm3 <- update(fm1, subset = -rmv)
plot(fm3)
qqnorm(resid(fm3))
qqline(resid(fm3))
Anova(fm3) # no difference so just present fm2

###########
# Summary #
###########

# Model comparison from GLMM
SummaryCompAIC <- ldply(list(C3total = c3PropComAIC, 
                             C4grass = C4grassCompAic, 
                             Grass = GrassPropCompAIC, 
                             Forb = ForbPropCompAIC), 
                        function(x)data.frame(x, terms = row.names(x)),
                        .id = "variable")
SummaryCompAIC

# F from LMM
SummaryAnvF_PFG <- ldply(list(Forb = AnvF_forb,
                              Grass = AnvF_grass, 
                              C3total = AnvF_PropC3_total, 
                              C4grass = AnvF_C4grass), 
                         function(x) data.frame(x, terms = row.names(x)), .id = "variable")
SumamryStat <- merge(SummaryCompAIC, SummaryAnvF_PFG, by = c("variable", "terms"), all = TRUE)
SumamryStat$terms <- factor(SumamryStat$terms, levels = c("Full", "co2", "year", "year:co2"))
SumamryStat <- SumamryStat[order(SumamryStat$variable, SumamryStat$terms), ]
SumamryStat$Pr <- round(SumamryStat$Pr..F., 3)
SumamryStat_PvalAic <- subset(SumamryStat, terms != "Full", select = c("variable", "terms", "dAIC", "Pr"))
SumamryStat_PvalAic_mlt <- melt(SumamryStat_PvalAic, id = c("variable", "terms"), variable_name = "stats")
SumamryStat_PvalAic_mlt$stats <- factor(SumamryStat_PvalAic_mlt$stats, levels = c("Pr", "dAIC"))
PFG_Sumamry <- dcast(SumamryStat_PvalAic_mlt, variable ~ terms + stats)
write.csv(PFG_Sumamry, file = "output/table/FACE_PFG_Prop_LMMGlmm.csv", row.names = FALSE)

## ---- Stats_C3_TotalPropSmmry
# The model
rm2@call

# Chisq
Anova(rm2)

# F test
AnvF_PropC3_total

# Contrast
PropC3_Total_CntrstRes

# Model diagnosis
plot(rm2)
qqnorm(resid(rm2))
qqline(resid(rm2))

