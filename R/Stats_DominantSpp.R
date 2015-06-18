head(veg)

# Identify dominant spp
DmSpp

# plot sum of dominant species
SppPlotSum <- ddply(subset(veg, variable %in% DmSpp), 
                    .(variable, year, co2, block, ring, plot, id), 
                    summarise, value = sum(value))
SppPlotSum$obs <- with(SppPlotSum, year:id)

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
                        year = factor(year, labels = paste0("Year", 1:3))})
tdf2 <- within(RingSumSpp_cst, {variable = gsub("[.]", "\n", as.character(variable))
                                year = factor(year, labels = paste0("Year", 1:3))})

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

#######
# GLM #
#######

# Microlaena.stipoides----
DmSpp[[1]]
msDF <- subset(SppPlotSum, variable == "Microlaena.stipoides")
bxplts(value = "value", xval = "co2", data = msDF)

m1 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson, data = msDF)
m2 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson(link = sqrt), data = msDF)
m3 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson(link = power(1/3)), data = msDF)
l_ply(list(m1, m2, m3), overdisp.glmer)

# m2
Anova(m2)
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))

# compare AIC
CompAIC_ms <- CompAIC(m2)
CompAIC_ms

# LMM
m2lmr <- lmer(sqrt(value + 1) ~ year * co2 + (1|block) + (1|ring)  + (1|id), data = msDF)
plot(m2lmr)
qqnorm(resid(m2lmr))
qqline(resid(m2lmr))
AnvF_ms <- Anova(m2lmr, test.statistic = "F")
AnvF_ms

# Pratia.purpurascens----
DmSpp[[2]]
ppDF <- subset(SppPlotSum, variable == "Pratia.purpurascens")
bxplts(value = "value", xval = "co2", ofst = 1, data = ppDF)

m1 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson, data = ppDF)
m2 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson(link = sqrt), data = ppDF)
m3 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson(link = power(1/3)), data = ppDF)
l_ply(list(m1, m2, m3), overdisp.glmer)
# overdispersed, include obs

m4 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id) + (1|obs), 
            family = poisson, data = ppDF)
m5 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id)+ (1|obs), 
            family = poisson(link = sqrt), data = ppDF)
m6 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id)+ (1|obs), 
            family = poisson(link = power(1/3)), data = ppDF)
l_ply(list(m4, m5, m6), overdisp.glmer)
# use m5

Anova(m5)
plot(m5)
qqnorm(resid(m5))
qqline(resid(m5))

# compare AIC
CompAIC_pp <- CompAIC(m5)
CompAIC_pp
  # year shouldn't be removed. co2 don't seemt to be imporatant as no change in
  # Conditiaonl R2

# tiny indication of co2 and year effect
m5lm <- lmer(sqrt(value + 1) ~ year * co2 + (1|block) + (1|ring)  + (1|id), data = ppDF)
AnvF_pp <- Anova(m5lm, test.statistic = "F")
AnvF_pp
# probably year but not co2

#Cynodon.dactylon----
DmSpp[[3]]
cdDF <- subsetD(SppPlotSum, variable == "Cynodon.dactylon")
bxplts(value = "value", xval = "co2", ofst = 1, data = cdDF)

m1 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson, data = cdDF)
m2 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson(link = sqrt), data = cdDF)
m3 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson(link = power(1/3)), data = cdDF)
l_ply(list(m1, m2, m3), overdisp.glmer)
# overdispersed
mls <- llply(list(m1, m2, m3), function(x) update(x, ~ . + (1|obs)))
l_ply(mls, overdisp.glmer)
# use the 2nd one
m2 <- mls[[2]]
Anova(m2)
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))
# compare AIC
CompAIC_cd <- CompAIC(m2)
CompAIC_cd

# small indication of year
m2lmr <- lmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), data = cdDF)
plot(m2lmr)
qqnorm(resid(m2lmr))
qqline(resid(m2lmr))

AnvF_cd <- Anova(m2lmr, test.statistic = "F")
AnvF_cd
# maybe not..

#Commelina.cyanea----
DmSpp[4]
ccDF <- subsetD(SppPlotSum, variable == "Commelina.cyanea")
bxplts(value = "value", xval = "co2", ofst = 1, data = ccDF)

m1 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson, data = ccDF)
m2 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson(link = sqrt), data = ccDF)
m3 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson(link = power(1/3)), data = ccDF)
l_ply(list(m1, m2), overdisp.glmer)
# log is better, but overdispersed
mls <- llply(list(m1, m2), function(x) update(x, ~. + (1|obs), 
                                              control = glmerControl(optimizer = "bobyqa")))
l_ply(mls, overdisp.glmer)
# use the 2nd model
m4 <- mls[[2]]
Anova(m4)
plot(m4)
qqnorm(resid(m4))
qqline(resid(m4))
# compare AIC
CompAIC_cc <- CompAIC(m4)
CompAIC_cc

# indication of massive year effect
m4lmer <- lmer(sqrt(value + 1) ~ year * co2 + (1|block) + (1|ring)  + (1|id), data = ccDF)
plot(m4lmer)
qqnorm(resid(m4lmer))
qqline(resid(m4lmer))
AnvF_cc <- Anova(m4lmer, test.statistic = "F")
AnvF_cc

#Hydrocotyle.peduncularis----
DmSpp[5]
hpDF <- subsetD(SppPlotSum, variable == "Hydrocotyle.peduncularis")
bxplts(value = "value", xval = "co2", ofst = 1, data = hpDF)

m1 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson, data = hpDF)
m2 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id), 
            family = poisson(link = sqrt), data = hpDF)
m3 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id) + (1|obs), 
            family = poisson(link = power(1/3)), data = hpDF)
# convergent problem in m1 and m3
l_ply(list(m1, m2, m3), overdisp.glmer)

# slightly overdespersed. m3 is probably not reliable
m4 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id) + (1|obs), 
            family = poisson, data = hpDF)
m5 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id) + (1|obs), 
            family = poisson(link = sqrt), data = hpDF)
m6 <- glmer(value ~ year * co2 + (1|block) + (1|ring)  + (1|id) + (1|obs), 
            family = poisson(link = power(1/3)), data = hpDF)
l_ply(list(m4, m5, m6), overdisp.glmer)
# use m5
Anova(m5)
plot(m5)
qqnorm(resid(m5))
qqline(resid(m5))
# compare AIC
CompAIC_hp <- CompAIC(m5)
CompAIC_hp

# LMM
m5lmr <- lmer(sqrt(value + 1) ~ year * co2 + (1|block) + (1|ring)  + (1|id), data = hpDF)
plot(m5lmr)
qqnorm(resid(m5lmr))
qqline(resid(m5lmr))
AnvF_hp <- Anova(m5lmr, test.statistic = "F")
AnvF_hp

###########
# Summary #
###########

# Combine Species names and object names of the results
a <- llply(strsplit(as.character(DmSpp), split = "[.]"))
names(a) <- DmSpp
AbrSpp <- ldply(a, function(x) paste(tolower(substring(x, 1, 1)), collapse = ""), .id = "variable")
AbrSpp$AmvF <- as.character(paste0("AnvF_", AbrSpp$V1)) # Anova results
AbrSpp$CmpAic <- as.character(paste0("CompAIC_", AbrSpp$V1)) # AIC results

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

# Coimpare AIC
DomSppCompAic <- ddply(AbrSpp, .(variable), function(x) {
  d <- get(x$CmpAic)
  d$terms <- row.names(d)
  return(d)
})

# merge and organise
DomSppSummary <- merge(DomSppAnvF, DomSppCompAic, by = c("variable", "terms"), all = TRUE)
DomSppSummary$terms <- factor(DomSppSummary$terms, levels = c("Full", "co2", "year", "year:co2"))
DomSppSummary <- DomSppSummary[order(DomSppSummary$variable, DomSppSummary$terms), ]
write.csv(DomSppSummary, file = "output/table/DominantSpp_Stats.csv", row.names = FALSE)

