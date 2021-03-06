
# Year0 = covariate -------------------------------------------------------

# Move Year0 value to a new column to be used as covariate for the analysis

DivDF_year0_list <- llply(DivDF_list, function(x){
  
  DivDF_mlt <- melt(x, id = c("year", "ring", "plot", "block", "co2", "id"))
  
  # year0
  year0_dd <- DivDF_mlt %>%
    filter(year == "Year0") %>%
    select(id, value, variable) %>%
    rename(value0 = value)
  
  # subseqent years
  subyear_dd <- filter(DivDF_mlt, year != "Year0")
  
  # merge
  DivDF_year0 <- merge(subyear_dd, year0_dd, by = c("id", "variable")) 
  
  return(DivDF_year0)
  })

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

l_ply(names(DivDF_year0_list), function(x){
  d_ply(DivDF_year0_list[[x]], .(variable), function(y){
    figtitle <- paste(x, unique(y$variable), sep = "-")
    plot(value ~ value0, pch = 19, col = year, data = y, main = figtitle)
  })
})


# all species -----------------------------------------------------------
DivDF_year0 <- DivDF_year0_list[["all_spp"]]

# . Evenness ----------------------------------------------------------------

############
# Evenness #
############

j_dd <- filter(DivDF_year0, variable == "J")
plot(value ~ value0, pch = 19, col = year, data = j_dd)

Eml1 <- lmer(value ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), data = j_dd)

AnvF_Eml <- Anova(Eml1, test.statistic = "F")
AnvF_Eml
plot(Eml1)
qqnorm(resid(Eml1))
qqline(resid(Eml1))


# . Diversity ---------------------------------------------------------------

#############
# Diversity #
#############
h_dd <- filter(DivDF_year0, variable == "H")
plot(value ~ value0, pch = 19, col = year, data = h_dd)

Dml1 <- lmer(value ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), 
             data = h_dd)

AnvF_Dml <- Anova(Dml1, test.statistic = "F")
AnvF_Dml

# CO2 is not quite significant, but F is reasonably hight (>5). worth to more
# explore
options(na.action = "na.fail")
tm <- update(Dml1, REML = FALSE)
mm <- dredge(tm)
mm
# removing co2 increase more than 2 (>4) AICc units, indicating that co2 is
# animportant term
tm2 <- get.models(mm, 1)[[1]]
tm3 <- update(tm2, REML = TRUE)
Anova(tm3, test.statistic = "F")

plot(Dml1)
qqnorm(resid(Dml1))
qqline(resid(Dml1))


# . Species richness --------------------------------------------------------

####################
# Species richness #
####################
s_dd <- filter(DivDF_year0, variable == "S")
par(mfrow = c(1, 3))
plot(value ~ value0, pch = 19, col = year, data = s_dd)
plot(log(value) ~ value0, pch = 19, col = year, data = s_dd)
plot(log(value) ~ log(value0), pch = 19, col = year, data = s_dd)

SmlLmm1 <- lmer(value ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), 
                data = s_dd)
SmlLmm2 <- lmer(log(value) ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), 
                data = s_dd)
SmlLmm3 <- lmer(log(value) ~ co2 * year + log(value0) + (1|block) + (1|ring) + (1|id), 
                data = s_dd)
ldply(list(SmlLmm1, SmlLmm2, SmlLmm3), r.squared)

plot(SmlLmm3)
qqnorm(resid(SmlLmm3))
qqline(resid(SmlLmm3))
AnvF_Sml <- Anova(SmlLmm3, test.statistic = "F")
AnvF_Sml


# Grass -------------------------------------------------------------------

DivDF_year0_grass <- DivDF_year0_list[["grass_spp"]]

# . Evenness ----------------------------------------------------------------

j_dd <- filter(DivDF_year0_grass, variable == "J" & !is.na(value))
plot(value ~ value0, pch = 19, col = year, data = j_dd)

Eml1 <- lmer(value ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), 
             data = j_dd)

plot(j_dd$value)

AnvF_Eml_grass <- Anova(Eml1, test.statistic = "F")
AnvF_Eml_grass
plot(Eml1)
qqnorm(resid(Eml1))
qqline(resid(Eml1))

# one obvious outlier. what if we remove
rm_val <- which.min(resid(Eml1))
Eml2 <- update(Eml1, subset = -rm_val)
summary(Eml2)
plot(Eml2)
qqnorm(resid(Eml2))
qqline(resid(Eml2))
AnvF_Eml_grass <- Anova(Eml2, test.statistic = "F")
AnvF_Eml_grass 

# . Diversity -------------------------------------------------------------
h_dd <- filter(DivDF_year0_grass, variable == "H")
plot(value ~ value0, pch = 19, col = year, data = h_dd)

Dml1 <- lmer(value ~ year* co2 + value0 + (1|block) + (1|ring) + (1|id), 
             data = h_dd)

AnvF_Dml_grass <- Anova(Dml1, test.statistic = "F")
AnvF_Dml_grass

plot(Dml1)
qqnorm(resid(Dml1))
qqline(resid(Dml1))

plot(lmerTest::lsmeans(Dml1))
lsmeans::lsmeans(Dml1, pairwise ~ co2 | year)

# . Species richness --------------------------------------------------------

s_dd <- filter(DivDF_year0_grass, variable == "S")
plot(value ~ value0, pch = 19, col = year, data = s_dd)
plot(log(value) ~ value0, pch = 19, col = year, data = s_dd)
plot(log(value) ~ log(value0), pch = 19, col = year, data = s_dd)

SmlLmm1 <- lmer(value ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), 
                data = s_dd)
SmlLmm2 <- lmer(log(value) ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), 
                data = s_dd)
SmlLmm3 <- lmer(log(value) ~ co2 * year + log(value0) + (1|block) + (1|ring) + (1|id), 
                data = s_dd)
ldply(list(SmlLmm1, SmlLmm2, SmlLmm3), r.squared)

plot(SmlLmm1)
qqnorm(resid(SmlLmm1))
qqline(resid(SmlLmm1))
AnvF_Sml_grass <- Anova(SmlLmm1, test.statistic = "F")
AnvF_Sml_grass

plot(lmerTest::lsmeans(SmlLmm1))


# Forb --------------------------------------------------------------------
DivDF_year0_forb <- DivDF_year0_list[["forb_spp"]]

# . Evenness ----------------------------------------------------------------

j_dd <- filter(DivDF_year0_forb, variable == "J" & !is.na(value))
plot(value ~ value0, pch = 19, col = year, data = j_dd)

Eml1 <- lmer(value ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), 
             data = j_dd)

AnvF_Eml_forb <- Anova(Eml1, test.statistic = "F")
AnvF_Eml_forb
# COVARIATE IS NOT SIGNIFICANT!!! RECONSIDER COVARIATE

plot(Eml1)
qqnorm(resid(Eml1))
qqline(resid(Eml1))

# . Diversity -------------------------------------------------------------
h_dd <- filter(DivDF_year0_forb, variable == "H")
plot(value ~ value0, pch = 19, col = year, data = h_dd)

Dml1 <- lmer(value ~ year* co2 + value0 + (1|block) + (1|ring) + (1|id), 
             data = h_dd)

AnvF_Dml_forb <- Anova(Dml1, test.statistic = "F")
AnvF_Dml_forb
# COVARIATE IS NOT SIGNIFICANT!!! RECONSIDER COVARIATE

plot(Dml1)
qqnorm(resid(Dml1))
qqline(resid(Dml1))


# . Species richness ------------------------------------------------------

s_dd <- filter(DivDF_year0_forb, variable == "S")
plot(value ~ value0, pch = 19, col = year, data = s_dd)
plot(log(value) ~ value0, pch = 19, col = year, data = s_dd)
plot(log(value) ~ log(value0), pch = 19, col = year, data = s_dd)

SmlLmm1 <- lmer(value ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), 
                data = s_dd)
SmlLmm2 <- lmer(log(value) ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), 
                data = s_dd)
SmlLmm3 <- lmer(log(value) ~ co2 * year + log(value0) + (1|block) + (1|ring) + (1|id), 
                data = s_dd)
ldply(list(SmlLmm1, SmlLmm2, SmlLmm3), r.squared)

plot(SmlLmm3)
qqnorm(resid(SmlLmm3))
qqline(resid(SmlLmm3))
AnvF_Sml_forb <- Anova(SmlLmm3, test.statistic = "F")
AnvF_Sml_forb

# number of species -----------------------------------------------------

totalSum <- ddply(veg, .(variable), summarise, value = sum(value, na.rm = TRUE))
summary(totalSum)
length(unique(totalSum$variable))


treatSum <- ddply(veg, .(variable, co2), summarise, value = sum(value, na.rm = TRUE))
nrow(treatSum)
# remove 0
treatSum <- subset(treatSum, value != 0)
nrow(treatSum)

vs <- dlply(treatSum, .(co2), function(x) droplevels(unique(x$variable)))

sapply(vs, length)
intersect(vs[[1]], vs[[2]])


# Summary -----------------------------------------------------------------


###########
# Summary #
###########
Anv_lst <- list('Species richness' = AnvF_Sml, 
                'Diveristy' = AnvF_Dml,
                'Evenness' = AnvF_Eml)

DivAnvF <- ldply(Anv_lst, function(x) {
  x$terms <- factor(row.names(x), levels = c("value0", "co2", "year", "co2:year"))
  return(x)}, .id = "Response")
names(DivAnvF)[5] <- "Pr"

DivAnvF <- within(DivAnvF, {
  F  <- round(F, 2)
  Df.res <- round(Df.res, 2)
  Pr <- round(Pr, 3)
  })
DivAnvF <- arrange(DivAnvF, Response, terms)

write.csv(DivAnvF, "output/table/SummaryResultDiversity.csv", row.names = FALSE)

Anv_df <- ldply(Anv_lst, function(x) {
  predictors = factor(row.names(x), 
                      levels = c("value0", "co2", "year", "co2:year"))
  pval = x$Pr
  stat = as.character(cut(pval, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                          labels = c("***", "**", "*", ".", "ns")))
  stat[stat == "ns"] <- round(pval[which(stat == "ns")], 2)
  data.frame(predictors, stat)
  }, .id = "variable")
Anv_cst <- dcast(variable ~ predictors, data = Anv_df, value.var = "stat")
write.csv(Anv_cst, file = "output/table/FACE_Diversity_StatSummary.csv")


# CO2 response ratio ------------------------------------------------------

######################
# CO2 response ratio #
######################
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
Dml1@call
# Chisq
Anova(Dml1)
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

# Model diagnosis
plot(Sml3)
qqnorm(resid(Sml3))
qqline(resid(Sml3))
