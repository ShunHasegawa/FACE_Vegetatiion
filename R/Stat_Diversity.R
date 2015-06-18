###########
# Eveness #
###########
bxplts(value = "J", xval = "co2", data = DivDF)
bxplts(value = "J", xval = "ring", data = DivDF)

Eml1 <- lmer(J ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF)
Anova(Eml1)
AnvF_Eml <- Anova(Eml1, test.statistic = "F")
AnvF_Eml
summary(Eml1)
plot(Eml1)
qqnorm(resid(Eml1))
qqline(resid(Eml1))

#############
# Diversity #
#############
bxplts(value = "H", xval = "co2", data = DivDF)
bxplts(value = "H", xval = "ring", data = DivDF)

Dml1 <- lmer(H ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF)
Anova(Dml1)
AnvF_Dml <- Anova(Dml1, test.statistic = "F")
AnvF_Dml
plot(Dml1)
qqnorm(resid(Dml1))
qqline(resid(Dml1))
plot(allEffects(Dml1))
 # Diversity decreased in eCO2

# contrast----
# contrast doesn't work with lmer so rewrite the model with lme
Dml_lme <- lme(H ~ co2 * year, random = ~1|block/ring/id, data = DivDF)

cntrst <- contrast(Dml_lme,
                   a = list(year = levels(DivDF$year), co2 = "amb"),
                   b = list(year = levels(DivDF$year), co2 = "elev"))
H_CntrstRes <- cntrstTbl(cntrst, data = DivDF, variable = "H", digit = 2)
H_CntrstRes

####################
# Species richness #
####################
bxplts(value = "S", xval = "co2", data = DivDF)
bxplts(value = "S", xval = "ring", data = DivDF)

Sml1 <- glmer(S ~ co2 * year + (1|block) + (1|ring) + (1|id), 
              data = DivDF, family = poisson)
Sml2 <- glmer(S ~ co2 * year + (1|block) + (1|ring) + (1|id), 
              data = DivDF, family = poisson(link = sqrt))
Sml3 <- glmer(S ~ co2 * year + (1|block) + (1|ring) + (1|id), 
              data = DivDF, family = poisson(link = power(1/3)))
l_ply(list(Sml1, Sml2, Sml3), overdisp.glmer)
# use Sml1
summary(Sml1)
plot(Sml1)
qqnorm(resid(Sml1))
qqline(resid(Sml1))
Anova(Sml1)
S_CompAic <- CompAIC(Sml1)
S_CompAic
# AIC decrease when co2:year removed but

# LMM
SmlLmm <- lmer(log(S) ~ co2 * year + (1|block) + (1|ring) + (1|id), data = DivDF)
summary(SmlLmm)
plot(SmlLmm)
qqnorm(resid(SmlLmm))
qqline(resid(SmlLmm))
AnvF_Sml <- Anova(SmlLmm, test.statistic = "F")
AnvF_Sml

###########
# Summary #
###########
Anv_lst <- list('Species richness' = AnvF_Sml, 
                'Diveristy' = AnvF_Dml,
                'Evenness' = AnvF_Eml)

Anv_df <- ldply(Anv_lst, function(x) {
  predictors = factor(row.names(x), levels = c("co2", "year", "co2:year"))
  pval = x$Pr
  stat = as.character(cut(pval, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                          labels = c("***", "**", "*", ".", "ns")))
  stat[stat == "ns"] <- round(pval[which(stat == "ns")], 2)
  data.frame(predictors, stat)
  }, .id = "variable")
Anv_cst <- dcast(variable ~ predictors, data = Anv_df, value.var = "stat")
write.csv(Anv_cst, file = "output/table/FACE_Diversity_StatSummary.csv")

######################
# Co2 response ratio #
######################
DivDF_mlt <- melt(DivDF, id = c("year", "block", "co2", "ring", "plot", "id"))
DivDF_RingMean <- ddply(DivDF_mlt, .(year, block, co2, ring, variable), 
                        summarise, value = mean(value))
DivDF_RingMean_cst <- dcast(DivDF_RingMean, year + block + variable ~ co2)

# compute SE for ratio using non-parametric bootstrap
ratio <- function(d, w) sum(d$elev * w) / sum(d$amb * w) -1 

# not the following function is wrong.. (don't know why...) for some reasons it
# estimates smaller SE ratio2 <- function(d, w) with(d, sum(value[co2 == "elev"]
# * w)/sum(value[co2 == "amb"] * w)) -1

RatioSE <- ddply(DivDF_RingMean_cst, .(year, variable), function(x) {
  b <- boot::boot(x, ratio, R = 999, stype = "w")
  summary(b)
})
RatioSE$co2R <- RatioSE$original

# organise it to export as a table
RatioSE$co2R_se <- with(RatioSE, paste0(round(co2R, 2), "(", round(bootSE, 2), ")"))

RatioSE_cst <- dcast(RatioSE[c("variable", "year", "co2R_se")], variable ~ year, value.var = "co2R_se")
write.csv(RatioSE_cst, file = "output/table/CO2ResponseRatio_DivIndx.csv", row.names = FALSE)

# plot

# block mean
subset(DivDF_RingMean_cst, year == 2013 & variable == "H" & block == "A")
DivDF_RingMean_cst$co2R <- with(DivDF_RingMean_cst, elev/amb -1 )

p <- ggplot(RatioSE, aes(x = year, y = co2R))
p2 <- p + 
  geom_point(aes(x = year, y = co2R), size = 4)+
  geom_errorbar(aes(x = year, ymin = co2R - bootSE, ymax = co2R + bootSE), width = 0) + 
  geom_point(data = DivDF_RingMean_cst, size = 2, col = "red", alpha = .7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "CO2 response ratio", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 7)) +
  facet_wrap( ~ variable, scale = "free_y")
p2
ggsavePP(plot = p2, filename = "output/figs/FACE_CO2ResponseRatio_DivIndx", width = 6, height = 3)

## ---- Stats_DiversityInx_Evenness
# Evenness
# The initial model
Eml1@call
Anova(Eml1)
# F
AnvF_Eml

## ---- Stats_DiversityInx_Diversity
# Shannon's index
# The model
Dml2@call
# Chisq
Anova(Dml2)
# F
AnvF_Dml
# Contrast
H_CntrstRes

## ---- Stats_DiversityInx_SpRichness
Sml3@call
# Chisq
Anova(Sml3)
# F test
AnvF_Sml
# Contrast
S_CntrstRes

# Model diagnosis
plot(Sml3)
qqnorm(resid(Sml3))
qqline(resid(Sml3))
