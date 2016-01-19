########################
# Redundancey analysis #
########################

###############
# All species #
###############

# each year separately----
rdaList <- dlply(RingSumVeg, .(year), function(x) rda(log(x[, SppName] + 1) ~ co2, data = x))

par(mfrow = c(1, 2))
l_ply(names(rdaList), function(x) plot(rdaList[[x]], main = x))

# each co2 treatment separately----
rdaList_year <- dlply(RingSumVeg, .(co2), function(x) rda(log(x[, SppName] + 1) ~ year, data = x))

par(mfrow = c(1, 2))
l_ply(names(rdaList_year), function(x) plot(rdaList_year[[x]], main = x))

# highest and lowest three spp scores
Smmry_rdaList_year <- llply(rdaList_year, summary)

llply(Smmry_rdaList_year, function(x){
  spsc <- x$species[, 1]
  sort(spsc)[c(1:3, 70:72)]
})

# two years together----
RingSumVeg$yco <- RingSumVeg$year:RingSumVeg$co2

RDAres <- rda(log(RingSumVeg[, SppName] + 1) ~ yco, data = RingSumVeg[, !names(RingSumVeg) %in% SppName])
RDAres <- rda(RingSumVeg[, SppName] ~ yco, data = RingSumVeg[, !names(RingSumVeg) %in% SppName])
RDAres

#######
# PFG #
#######

# rda

# each year separately----
pfg_co2 <- dlply(RingSumPFGMatrix, .(year), function(x) rda(log(x[, PFGName] + 1) ~ co2, data = x))
par(mfrow = c(1, 2))
l_ply(names(pfg_co2), function(x) plot(pfg_co2[[x]], main = x))

# each co2 treatment separately----
pfg_year <- dlply(RingSumPFGMatrix, .(co2), function(x) rda(log(x[, PFGName] + 1) ~ year, data = x))
par(mfrow = c(1, 2))
l_ply(names(pfg_year), function(x) plot(pfg_year[[x]], main = x))

SmmryPfgYear <- llply(pfg_year, summary)
PfgYearSitedf <- ldply(names(SmmryPfgYear), function(x) {
  data.frame(rda1 = SmmryPfgYear[[x]]$sites[, 1],
             subset(RingSumPFGMatrix, co2 == x)[, c("year", "ring", "co2")])})
PfgYearSitedf <- within(PfgYearSitedf, {
  co2 <- factor(co2, labels = c("Ambient RDA(7.9%)", 
                                "eCO2 RDA(4.1%)"))
  year <- factor(year, labels = c("2013\nPre-CO2", "2014"))
})

PfgYearSppdf <- ldply(names(SmmryPfgYear), function(x) {
  data.frame(rda1 = SmmryPfgYear[[x]]$species[, 1], 
             PFG = row.names(SmmryPfgYear[[x]]$species),
             co2 = x, 
             year = "Species score")})
PfgYearSppdf$co2 <- factor(PfgYearSppdf$co2, labels = c("Ambient RDA(7.9%)", "eCO2 RDA(4.1%)"))

theme_set(theme_bw())
p <- ggplot(PfgYearSitedf, aes(x = year, y = rda1))
p2 <- p + geom_point(aes(shape = year), size = 5, alpha = .7) + 
  geom_line(aes(group = ring), linetype = 3) + 
  geom_point(data = PfgYearSppdf, aes(x = year, y = rda1, col = PFG), 
             size = 5, alpha = .7, position = position_dodge(.5)) +
  facet_grid(~co2) +
  labs(x = NULL, y = "RDA1") +
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = .7)
p2

# two years together----
PFGrdaRes <- rda(log(RingSumPFGMatrix[, PFGName] + 1) ~ yco, data = RingSumPFGMatrix)
PFGrdaRes
plot(PFGrdaRes)
