head(PlotSumVeg)

#################
# Dissimilarity #
#################

# compute dissimilarity for each plot
disDF <- ddply(PlotSumVeg, .(block, ring, co2, plot), function(x) YearDssmlrty(x, spp = SppName))
disDF$id <- disDF$ring:disDF$plot

################
## CO2 x Year ##
################

# perform LMM----
bxplts(value = "EU", xval = "co2", data = disDF)
bxplts(value = "EU", xval = "ring", data = disDF)
  
# raw data
m1 <- lmer(EU ~ co2 * year + (1|block) + (1|ring) + (1|id), data = disDF)
DisAllSp_AnvF <- Anova(m1, test.statistic = "F")
plot(m1)
qqnorm(resid(m1))
qqline(resid(m1))

###########################
# Plant Functional Groups #
###########################
head(PlotSumVeg)
PFGName

# compute dissimilarity for each plot

Pfg_disDF <- ddply(PlotSumPFGMatrix, .(block, ring, co2, plot), 
               function(x) YearDssmlrty(x, spp = PFGName))
Pfg_disDF$id <- with(Pfg_disDF, ring:plot)

bxplts(value = "EU", xval = "ring", data = Pfg_disDF)
bxplts(value = "EU", xval = "co2", data = Pfg_disDF)

#################
## perform LMM ##
#################
# raw data
m1 <- lmer(EU ~ co2 * year + (1|block) + (1|ring) + (1|id), data = Pfg_disDF)
Anova(m1)
m2 <- stepLmer(m1)
DisPfg_AnvF <- Anova(m2, test.statistic = "F")
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))

##########
# Figure #
##########
# summary df---
disDF$type <- "All species"
PFG_vsDF$type <- "PFG"

disDF_Ring <- ldply(list(disDF, PFG_vsDF), function(x) {
  ddply(x, .(type, year, block, co2, ring), summarise, EU = mean(EU))
})

disDF_CO2 <- ddply(disDF_Ring, .(type, year, co2), summarise, 
                   Mean = mean(EU), 
                   SE = ci(EU)[4], 
                   N = sum(!is.na(EU)))

# plot----
PsDdg <- .3
p <- ggplot(data = disDF_CO2, aes(x = year, y = Mean, fill = co2))
p2 <- p +  
  geom_errorbar(aes(x = year, ymax = Mean + SE, ymin = Mean - SE), 
                position = position_dodge(PsDdg), 
                width = 0) +
  geom_point(size = 3, position = position_dodge(PsDdg), shape = 21) +
  scale_fill_manual(values = c("black", "white"), labels = c("Ambient", expression(eCO[2]))) + 
  labs(x = "Year", y = "Euclidean dissimilarity") + 
  facet_wrap(~ type, scale = "free_y") + 
  science_theme + 
  theme(legend.position = c(.85, .85))
ggsavePP(filename = "output/figs/FACE_EU_Dissimilarity_CO2", plot = p2,
         width = 4, height = 2.5)


#######################
# CO2 response ratios #
#######################
disDF_Ring

# cast to get each block ratio
disDF_Ring_cst <- dcast(disDF_Ring, type + year + block ~ co2, value.var = "EU")
disDF_Ring_cst$co2R <- disDF_Ring_cst$elev/disDF_Ring_cst$amb - 1

# compute SE for ratio using non-parametric bootstrap
ratio <- function(d, w) sum(d$elev * w)/sum(d$amb * w) - 1
RatioSE <- ddply(disDF_Ring_cst, .(type, year), function(x) {
  b <- boot::boot(x, ratio, R = 999, stype = "w")
  summary(b)
})

# plot
tdf <- within(RatioSE, {co2R = original})

theme_set(theme_bw())
p <- ggplot(tdf, aes(x = year, y = co2R))
p2 <- p + geom_point(size = 3) + 
  geom_errorbar(aes(x = year, ymin = co2R - bootSE, ymax = co2R + bootSE), width = 0) +
  geom_point(data = disDF_Ring_cst, aes(x = year, y = co2R), col = "red", alpha = .7, size = 2)+
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap( ~ type, scale = "free_y", ncol = 5) +
  labs(y = "CO2 response ratio", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 7))
p2
ggsavePP(plot = p2, filename = "output/figs/FACE_CO2ResponseRatio_Dissimilarity", width = 4, height = 3)

# organise it to export as a table
RatioSE$co2R <- with(RatioSE, paste0(round(original, 2), "(", round(bootSE, 2), ")"))

# reorder
RatioSE_cst <- dcast(RatioSE, type ~ year, value.var = "co2R")
write.csv(RatioSE_domspp_cst, file = "output/table/CO2ResponseRatio_Dissimilarity.csv", row.names = FALSE)

###########
# Summary #
###########
DisAllSp_AnvF
DisPfg_AnvF
