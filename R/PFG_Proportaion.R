head(veg)


# co2 effect (% change) ---------------------------------------------------

#########################
# co2 effect (% change) #
#########################
# plant forms
summary(veg)
# get C3, grass, forb etc. abundance fo each plot
PfgAbundDF <- ddply(veg, .(year, co2, block, ring), function(x){
  Total   = sum(x$value)
  C3      = sum(x$value[!x$PFG %in% c("c4", "moss", "fern")])
  C4grass = sum(x$value[x$PFG  == "c4"])
  Forb    = sum(x$value[x$form == "Forb"])
  Grass   = sum(x$value[x$form == "Grass"])
  Wood    = sum(x$value[x$form == "Wood"])
  Moss    = sum(x$value[x$form == "Moss"])
  data.frame(Total, C3, C4grass, Forb, Grass, Wood, Moss)
})
PfgAbundDF_mlt <- melt(PfgAbundDF, 
                       id = c("year", "co2", "block", "ring", "Total"))
PfgAbundDF_mlt_II <- melt(PfgAbundDF_mlt, 
                          id            = c("year", "co2", "block", "ring", 
                                            "variable"), 
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
  geom_point(aes(x = year, y = co2R), size = 4) +
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

RatioSE_cst <- dcast(RatioSE[c("variable", "year", "co2R_se")], 
                     variable  ~ year, 
                     value.var = "co2R_se")
write.csv(RatioSE_cst, 
          file      = "output/table/CO2ResponseRatio_PFGfraction.csv", 
          row.names = FALSE)


# Analysis ----------------------------------------------------------------

############
# Analysis #
############

# get C3, grass, forb etc. proportion for each plot
PropDF <- ddply(veg, .(year, co2, block, ring, plot, id), function(x){
  Total     = sum(x$value)
  C3prop    = 1-sum(x$value[x$PFG %in% c("c4", "moss", "fern")])/Total
  C4grass   = sum(x$value[x$PFG == "c4"])/Total
  Forbprop  = sum(x$value[x$form == "Forb"])/Total
  Grassprop = sum(x$value[x$form %in% c("Grass", "Sedge")])/Total
  Woodprop  = sum(x$value[x$form %in% "Wood"])/Total
  Mossprop  = sum(x$value[x$form == "Moss"])/Total
  data.frame(Total, C3prop, C4grass, Forbprop, Grassprop, Woodprop, Mossprop)
  })

# Move Year0 value to a new column to be used as covariate in the analyssis
# below

  PropDF_mlt <- melt(PropDF, id = c("year", "co2", "block", "ring", "plot", "id", "Total"))
  
  # Year0
  year0_dd <- PropDF_mlt %>%                   
    filter(year == "Year0") %>%
    select(variable, id, value) %>%
    rename(value0 = value)
  
  # subsequent years
  subyear_dd <- filter(PropDF_mlt, year != "Year0") 
  
  # merge
  PropDF_year0 <- merge(subyear_dd, year0_dd, by = c("id", "variable")) 
  PropDF_year0$obs <- 1:nrow(PropDF_year0) 

  unique(PropDF_year0$variable)
  par(mfrow = c(2, 3), mar = c(5, 4, 2, 1))
  d_ply(PropDF_year0, .(variable), function(x) {
    plot(logit(value) ~ value0, pch = 19, col = year, data = x, main = unique(x$variable))
    })


# C3 proportion in the whole community ------------------------------------

########################################
# C3 proportion in the whole community #
########################################

c3prop_dd <- filter(PropDF_year0, variable == "C3prop")
par(mfrow = c(1, 3))
plot(value ~ value0, data = c3prop_dd, pch = 19, col = year)
plot(asin(value) ~ value0, data = c3prop_dd, pch = 19, col = year)
plot(logit(value) ~ value0, data = c3prop_dd, pch = 19, col = year)

rm1 <- lmer(value ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id),
            data = c3prop_dd)
rm2 <- lmer(asin(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id),
            data = c3prop_dd)
rm3 <- lmer(logit(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id),
            data = c3prop_dd)

ldply(list(rm1, rm2, rm3), r.squared)

Anova(rm1, test.statistic = "F")
plot(rm1)
qqnorm(resid(rm1))
qqline(resid(rm1))

AnvF_PropC3_total <- Anova(rm1, test.statistic = "F")
AnvF_PropC3_total


# C4 grass ----------------------------------------------------------------

############
# C4 grass #
############
c4gprop_dd <- filter(PropDF_year0, variable == "C4grass")

par(mfrow = c(2, 2))
boxplot(logit(value) ~ year:ring, data = c4gprop_dd, main = "logit")
plot(logit(value) ~ value0, data = c4gprop_dd, pch = 19, col = year)
plot(value ~ value0, data = c4gprop_dd, pch = 19, col = year)

m2Lmer <- lmer(logit(value) ~ year * co2 + value0 + 
                 (1|block) +(1|ring) + (1|id), data = c4gprop_dd)
m3Lmer <- lmer(value ~ year * co2 + value0 + 
                 (1|block) +(1|ring) + (1|id), data = c4gprop_dd)
ldply(list(m2Lmer, m3Lmer), r.squared)

plot(m3Lmer)
qqnorm(resid(m3Lmer))
qqline(resid(m3Lmer))
AnvF_C4grass <- Anova(m3Lmer, test.statistic = "F")
AnvF_C4grass


# . post-hoc test ---------------------------------------------------------
plot(lmerTest::lsmeans(m3Lmer))
lsmeans::lsmeans(m3Lmer, pairwise ~ co2 | year)


# Grass + Sedge + rush ----------------------------------------------------

########################
# Grass + Sedge + rush #
########################
gProp_dd <- filter(PropDF_year0, variable == "Grassprop")

par(mfrow = c(1, 3))
plot(value ~ value0, pch = 19, col = year,  data = gProp_dd)
plot(asin(value) ~ value0, pch = 19, col = year, data = gProp_dd)
plot(logit(value) ~ value0, pch = 19, col = year, data = gProp_dd)

gm1 <- lmer(value ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
            data = gProp_dd)
gm2 <- lmer(logit(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
            data = gProp_dd)
gm3 <- lmer(asin(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
            data = gProp_dd)
ldply(list(gm1, gm2, gm3), r.squared)

plot(gm2)
qqnorm(resid(gm2))
qqline(resid(gm2))
AnvF_Grass <- Anova(gm2, test.statistic = "F")
AnvF_Grass

# . post-hoc test ---------------------------------------------------------

plot(lmerTest::lsmeans(gm2))
lsmeans::lsmeans(gm2, pairwise ~ co2 | year)


# Forb --------------------------------------------------------------------


########
# Forb #
########

forbProp_dd <- filter(PropDF_year0, variable == "Forbprop")

par(mfrow = c(1, 3))
plot(value ~ value0, pch = 19, col = year, data = forbProp_dd)
plot(asin(value) ~ value0,  pch = 19, col = year, data = forbProp_dd)
plot(logit(value) ~ value0, pch = 19, col = year, data = forbProp_dd)

fm1 <- lmer(value ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
            data = forbProp_dd)
fm2 <- lmer(logit(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
            data = forbProp_dd)
fm3 <- lmer(asin(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
            data = forbProp_dd)
ldply(list(fm1, fm2, fm3), r.squared)

Anova(fm2)
AnvF_forb <- Anova(fm2, test.statistic = "F")
AnvF_forb

# model diagnosis
plot(fm2)
qqnorm(resid(fm2))
qqline(resid(fm2))


# Summary -----------------------------------------------------------------

###########
# Summary #
###########

# F from LMM
SummaryAnvF_PFG <- ldply(list(Forb    = AnvF_forb,
                              Grass   = AnvF_Grass, 
                              C3total = AnvF_PropC3_total, 
                              C4grass = AnvF_C4grass), 
                         function(x) data.frame(x, terms = row.names(x)), 
                         .id = "variable")

# summary of anova
PFGResAnvF <- SummaryAnvF_PFG
names(PFGResAnvF)[5] <- "Pr"
PFGResAnvF <- within(PFGResAnvF, {
              F      <- round(F, 2)
              Df.res <- round(Df.res, 0)
              Pr     <- round(Pr, 3)
            })
PFGResAnvF$terms <- factor(PFGResAnvF$terms, 
                           levels = c("value0", "co2", "year", "year:co2"))
PFGResAnvF <- PFGResAnvF[order(PFGResAnvF$variable, PFGResAnvF$terms), ]
write.csv(PFGResAnvF, "output/table/FACE_PFG_AnvF.csv", row.names = FALSE)

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
