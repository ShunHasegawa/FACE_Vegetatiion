
# fit polinomial time -----------------------------------------------------


library(lattice)
tdf <- PropDF
tdf$timeN <-  as.numeric(tdf$year)

xyplot(logit(C4grass) ~ timeN | block , data=tdf, group=co2)
print(xyplot(logit(C4grass) ~ timeN | ring + plot, group=co2,  tdf, type = c("r", "p")))


m2Lmer <- lmer(logit(C4grass) ~ co2 + 
                 (1|block) + (1|ring) + (timeN|id), 
               data = tdf)
m3Lmer <- lmer(logit(C4grass) ~ co2 + 
                 (1|block) + (1|ring) + (1|id), 
               data = tdf)
m4Lmer <- lmer(logit(C4grass) ~ co2 * timeN + 
                 (1|block) + (1|ring) + (1|id), 
               data = tdf)


m5Lmer <- lmer(logit(C4grass) ~ co2 * poly(timeN, 2) + 
                 (1|block) + (1|ring) + (1|id), 
               data = tdf)


lplist <- llply(1:3, 
                function(x) {
                fl <- lmer(logit(C3prop) ~ co2 * poly(timeN, x) + 
                             (1|block) + (1|ring) + (1|id), 
                           data = tdf)
                return(fl)
                })
anova(lplist[[1]], lplist[[3]])
ldply(lplist, r.squared)

anova(m2Lmer, m3Lmer, m4Lmer, m5Lmer)
Anova(m4Lmer, test.statistic = "F")


ml3 <- lmer(logit(C3prop) ~ co2 * poly(timeN, 1) + 
             (1|block) + (1|ring) + (1|id), 
           data = tdf)
visreg(ml3, xvar = "timeN", by = "co2", overlay = "TRUE")


?visreg
Anova(m2Lmer, test.statistic = "F")
plot(m2Lmer)

expDF <- expand.grid(co2  = c("amb", "elev"),
                    timeN = 1:4)

bb <- bootMer(ml3,
              FUN  = function(x) predict(x, expDF, re.form = NA),
              nsim = 500)
lci <- apply(bb$t, 2, quantile, 0.025)
uci <- apply(bb$t, 2, quantile, 0.975)
PredVal <- bb$t0
ddf <- cbind(lci, uci, PredVal, expDF)

p <- ggplot(ddf, aes(x = timeN, y = PredVal, col = co2, fill = co2))
p2 <- p + 
      geom_line() +
      geom_ribbon(aes(x = timeN, ymin = lci, ymax = uci), 
                  alpha = .5)

p2
#





m1 <- glmer(C4grass ~ co2 +  (timeN|block) +(1|ring) + (1|id), 
            family = "binomial", data = tdf, weight = Total)
overdisp.glmer(m1)

anova(m1)

summary(m1)



cdDF <- subsetD(SppPlotSum, variable == "Cynodon.dactylon")

tdf2 <- cdDF
tdf2$timeN <- as.numeric(tdf2$year)

mll <- llply(1:3, 
              function(x) lmer(sqrt(value + 1) ~ poly(timeN, x) * co2 + 
                                 (1|block) + (1|ring)  + (1|id), data = tdf2))
ldply(mll, r.squared)

m2 <- lmer(sqrt(value + 1) ~ poly(timeN, 1) * co2 + 
             (1|block) + (1|ring)  + (1|id), data = tdf2)
Anova(m2, test.statistic = "F")
visreg(m2, xvar = "timeN", by = "co2", overlay = "TRUE")
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))


# scatter plot for PFG composition ----------------------------------------


head(PFG_Fraction)

p <- ggplot(PFG_Fraction, 
            aes(x = year, y = original, col = co2, group = co2))
p2 <- p + 
  geom_point(size = 3, position = position_dodge(.2)) +
  geom_line() +
  geom_errorbar(aes(ymin = original - bootSE, ymax = original + bootSE),
                width = 0.2, position = position_dodge(.2)) +
  facet_grid(. ~ PFG)
p2
?geom_errorbar

# ANCOVA with initial data as a covariate ---------------------------------

PfgAbundDF_cst


ttd <- subsetD(PfgAbundDF_cst, variable == "C4grass")
ttd$R <- with(ttd, (elev_value/elev_Total)/(amb_value/amb_Total)) 
plot(log(R) ~  as.numeric(year), col = block, data = ttd)
head(ttd)

ttd$yearN <- as.numeric(ttd$year)


m1 <- lmer(log(R) ~ yearN + (1|block) , data = ttd)
visreg(m1, xvar = "yearN")

ttd2 <- cbind(ttd, R1 = subset(ttd, year == "Year0")[, "R"])
ttd2 <- subsetD(ttd2, year != "Year0")

plot(R ~ R1, col = year, data = ttd2)

m2 <- aov(log(R) ~ yearN + block, data = ttd)


anova(m2)
visreg(m2, xvar = "yearN", by = "block", overlay = TRUE)

?visreg
Anova(m1, test.statistic = "F")
plot(m1)
qqnorm(residuals(m1))
qqline(residuals(m1))




subset(RatioSE, variable == "C4grass")



init <- subset(PropDF, year == "Year0")[, c("id", "C4grass")]
names(init)[2] <- "C4_init"
ttd3 <- merge(PropDF, init, by = "id", all.x = TRUE)
ttd3 <- subsetD(ttd3, year != "Year0")
plot(log(C4grass+1) ~ log(C4_init+1), data = ttd3, col = co2, pch = 19)

ttd3$yearN <- as.numeric(ttd3$year)
m1 <- lmer(C4grass ~ co2 * yearN + C4_init + (1|block) + (1|ring) + (1|id), 
           data = ttd3)
m2 <- lmer(C4grass ~ co2 * yearN + (C4_init|block) + (1|ring) + (1|id), 
           data = ttd3)
m3 <- lmer(C4grass ~ co2 * yearN + (1|block) + (C4_init|ring) + (1|id), 
           data = ttd3)
m4 <- lmer(C4grass ~ co2 * yearN + (C4_init|block) + (C4_init|ring) + (1|id), 
           data = ttd3)
anova(m1, m2, m3, m4)
Anova(m1, test.statistic = "F")
visreg(m1, xvar = "yearN", by = "co2", overlay = TRUE)
visreg(m1, xvar = "C4_init", by = "co2", overlay = TRUE)


summary(m3)

m5 <- lme(C4grass ~ co2, random = ~as.numeric(year)|id, data = PropDF)
anova(m5)



ldply(list(m2), r.squared)

Anova(m1, test.statistic = "F")
visreg(m1, xvar = "yearN", by = "co2", overlay = TRUE)


