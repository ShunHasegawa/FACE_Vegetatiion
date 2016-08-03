head(veg)

# Identify dominant spp
DmSpp

# plot sum of dominant species
SppPlotSum <- ddply(subset(veg, variable %in% DmSpp), 
                    .(variable, year, co2, block, ring, plot, id), 
                    summarise, value = sum(value))
SppPlotSum$obs <- with(SppPlotSum, year:id)

# Move Year0 value to a new column to be used as covariate for the analysis

  # year0
  year0_dd <-SppPlotSum %>%
    filter(year == "Year0") %>%
    select(id, value, variable) %>%
    rename(value0 = value)
  
  # subseqent years
  subyear_dd <- filter(SppPlotSum, year != "Year0")
  
  # merge
  SppPlotSum_year0 <- merge(subyear_dd, year0_dd, by = c("id", "variable")) 
  
  unique(SppPlotSum_year0$variable)
  par(mfrow = c(2, 2), mar = c(5, 4, 2, 1))
  d_ply(SppPlotSum_year0, .(variable),  function(x) 
    plot(value ~ value0, pch = 19, col = year, data = x, 
         main = unique(x$variable)))
  
# co2 effect (% change) ---------------------------------------------------

#########################
# co2 effect (% change) #
#########################

RingSumSpp <- ddply(veg, .(year, block, co2, variable), summarise, value = sum(value))
RingSumSpp
# there're lots of 0s which cause trboule with calculating ratios so add 1
RingSumSpp$value <- RingSumSpp$value + 1

# cast for to get each block ratio
RingSumSpp_cst <- dcast(RingSumSpp, variable + year + block ~ co2)
RingSumSpp_cst$co2R <- RingSumSpp_cst$elev/RingSumSpp_cst$amb - 1

# compute SE for ratio using non-parametric bootstrap
ratio <- function(d, w) sum(d$elev * w)/sum(d$amb * w) - 1
RatioSE <- ddply(RingSumSpp_cst, .(variable, year), function(x) {
  b <- boot::boot(x, ratio, R = 999, stype = "w")
  summary(b)
  })

# plot
tdf <- within(RatioSE, {co2R = original
                        variable = gsub("[.]", "\n", as.character(variable))
                        })
tdf2 <- within(RingSumSpp_cst, {variable = gsub("[.]", "\n", as.character(variable))
                                })

theme_set(theme_bw())
p <- ggplot(tdf, aes(x = year, y = co2R))
p2 <- p + geom_point(size = 3) + 
  geom_errorbar(aes(x = year, ymin = co2R - bootSE, ymax = co2R + bootSE), width = 0) +
  geom_point(data = tdf2, aes(x = year, y = co2R), col = "red", alpha = .7, size = 2)+
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap( ~ variable, scale = "free_y", ncol = 5) +
  labs(y = "CO2 response ratio", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 7))
p2
ggsavePP(plot = p2, filename = "output/figs/FACE_CO2ResponseRatio_Spp", width = 7.5, height = 20)

# organise it to export as a table
RatioSE$co2R <- with(RatioSE, paste0(round(original, 2), "(", round(bootSE, 2), ")"))
RatioSE_domspp <- subset(RatioSE, variable %in% DmSpp, select = c("variable", "year", "co2R"))
RatioSE_domspp$variable <- factor(RatioSE_domspp$variable, levels = DmSpp)
# reorder
RatioSE_domspp <- RatioSE_domspp[order(as.numeric(RatioSE_domspp$variable)), ]
RatioSE_domspp_cst <- dcast(RatioSE_domspp, variable ~ year, value.var = "co2R")
write.csv(RatioSE_domspp_cst, file = "output/table/CO2ResponseRatio_DominantSpp.csv", row.names = FALSE)


# Microlaena.stipoides ----------------------------------------------------

DmSpp[[1]]
msDF <- subset(SppPlotSum_year0, variable == "Microlaena.stipoides")
bxplts(value = "value", xval = "co2", data = msDF)
bxplts(value = "value", xval = "ring", data = msDF)
par(mfrow = c(1, 2))
boxplot(value ~ year:ring, data = msDF, main = "raw")
boxplot(logit(value) ~ year:ring, data = msDF, main = "logit")
plot(value ~ value0, data = msDF, pch = 19, col = year)
plot(sqrt(value) ~ value0, data = msDF, pch = 19, col = year)
plot(sqrt(value) ~ sqrt(value0), data = msDF, pch = 19, col = year)

msDF$sqrt_value0 <- sqrt(msDF$value0)
m2lmr <- lmer(value ~ year * co2 + value0 + 
                (1|block) + (1|ring)  + (1|id), data = msDF)
m3lmr <- lmer(sqrt(value) ~ year * co2 + value0 + 
                (1|block) + (1|ring)  + (1|id), data = msDF)
m4lmr <- lmer(sqrt(value) ~ year * co2 + sqrt_value0 + 
                (1|block) + (1|ring)  + (1|id), data = msDF)
ldply(list(m2lmr, m3lmr, m4lmr), r.squared)
AICc(m2lmr, m3lmr, m4lmr)

plot(m4lmr)
qqnorm(resid(m4lmr))
qqline(resid(m4lmr))
AnvF_ms <- Anova(m4lmr, test.statistic = "F")
AnvF_ms


# Pratia.purpurascens -----------------------------------------------------

DmSpp[[2]]
ppDF <- subset(SppPlotSum_year0, variable == "Pratia.purpurascens")
bxplts(value = "value", xval = "co2", ofst = 1, data = ppDF)
par(mfrow = c(2, 2))
plot(value ~ value0, pch = 19, col = year, data = ppDF)
plot(sqrt(value) ~ value0, pch = 19, col = year, data = ppDF)
plot(sqrt(value) ~ sqrt(value0), pch = 19, col = year, data = ppDF)

ppDF$sqrt_value0 <- sqrt(ppDF$value0)

m5lm <- lmer(value ~ year * co2 + value0 + (1|block) + (1|ring)  + (1|id), 
             data = ppDF)
m6lm <- lmer(sqrt(value) ~ year * co2 + value0 + (1|block) + (1|ring)  + (1|id), 
             data = ppDF)
m7lm <- lmer(sqrt(value) ~ year * co2 + sqrt_value0 + (1|block) + (1|ring)  + (1|id), 
             data = ppDF)
ldply(list(m5lm, m6lm, m7lm), r.squared)
AICc(m5lm, m6lm, m7lm)

AnvF_pp <- Anova(m7lm, test.statistic = "F")
AnvF_pp 
plot(m7lm)
qqnorm(resid(m7lm))
qqline(resid(m7lm))


# Cynodon.dactylon --------------------------------------------------------

DmSpp[[3]]
cdDF <- subsetD(SppPlotSum_year0, variable == "Cynodon.dactylon")
bxplts(value = "value", xval = "co2", ofst = 1, data = cdDF)
par(mfrow = c(2, 2))
plot(value ~ value0, pch = 19, col = year, data = cdDF)
plot(sqrt(value) ~ value0, pch = 19, col = year, data = cdDF)
plot(sqrt(value) ~ sqrt(value0), pch = 19, col = year, data = cdDF)

cdDF$sqrt_valu0 <- sqrt(cdDF$value0)
m2lmr <- lmer(value ~ year * co2 + value0 + 
                (1|block) + (1|ring)  + (1|id), data = cdDF)
m3lmr <- lmer(sqrt(value) ~ year * co2 + value0 + 
                (1|block) + (1|ring)  + (1|id), data = cdDF)
m4lmr <- lmer(sqrt(value) ~ year * co2 + sqrt_valu0 + 
                (1|block) + (1|ring)  + (1|id), data = cdDF)
ldply(list(m2lmr, m3lmr, m4lmr), r.squared)
AICc(m2lmr, m3lmr, m4lmr)

plot(m4lmr)
qqnorm(resid(m4lmr))
qqline(resid(m4lmr))

AnvF_cd <- Anova(m3lmr, test.statistic = "F")
AnvF_cd

# . post-hoc test ---------------------------------------------------------

plot(lmerTest::lsmeans(m4lmr))
lsmeans::lsmeans(m4lmr, pairwise ~ co2 | year)


# Commelina.cyanea --------------------------------------------------------
DmSpp[4]
ccDF <- subsetD(SppPlotSum_year0, variable == "Commelina.cyanea")
bxplts(value = "value", xval = "co2", ofst = 1, data = ccDF)
par(mfrow = c(2, 2))
plot(value ~ value0, pch = 19, col = year, data = ccDF)
plot(sqrt(value) ~ value0, pch = 19, col = year, data = ccDF)
plot(sqrt(value) ~ sqrt(value0), pch = 19, col = year, data = ccDF)

ccDF$sqrt_value0 <- sqrt(ccDF$value0)
m5lmer <- lmer(value ~ year * co2 + value0 + 
                 (1|block) + (1|ring)  + (1|id), data = ccDF)
m6lmer <- lmer(sqrt(value) ~ year * co2 + value0 + 
                 (1|block) + (1|ring)  + (1|id), data = ccDF)
m7lmer <- lmer(sqrt(value) ~ year * co2 + sqrt_value0 + 
                 (1|block) + (1|ring)  + (1|id), data = ccDF)
ldply(list(m5lmer, m6lmer, m7lmer), r.squared)
AICc(m5lmer, m6lmer, m7lmer)

plot(m7lmer)
qqnorm(resid(m7lmer))
qqline(resid(m7lmer))
AnvF_cc <- Anova(m7lmer, test.statistic = "F")
AnvF_cc


# Summary -----------------------------------------------------------------

###########
# Summary #
###########

# Combine Species names and object names of the results
a <- llply(strsplit(as.character(DmSpp), split = "[.]"))
names(a) <- DmSpp
AbrSpp <- ldply(a, 
                function(x) paste(tolower(substring(x, 1, 1)), collapse = ""), 
                .id = "variable")
AbrSpp$AmvF <- as.character(paste0("AnvF_", AbrSpp$V1)) # Anova results

# combine all results

# Anova F
DomSppAnvF <- ddply(AbrSpp, .(variable), function(x) {
  d <- get(x$AmvF)
  d <- within(d, {
          F = round(F, 3)
          Df.res = round(Df.res, 0)
          Pr = round(d$Pr, 3)
          'Pr(>F)' = NULL
          terms = row.names(d)
        })
  return(d)
  })

# organise
DomSppAnvF <- arrange(DomSppAnvF, variable, terms) 
write.csv(DomSppAnvF, file = "output/table/DominantSpp_Stats.csv", row.names = FALSE)
