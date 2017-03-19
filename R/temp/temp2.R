c34sum$logmoist <- log(c34sum$totalmoist)
c4d_m6 <- lmer(scale(c4_ddiff) ~ co2 * (scale(logmoist) + scale(annual_temp2m) + scale(c3)) + (1|ring) + (1|RY) + (1|id), data = c34sum)
c4d_m7 <- lmer(scale(c4_ddiff) ~ co2 * (scale(logmoist) + scale(annual_temp2m) + scale(log(c3 + 1))) + (1|ring) + (1|RY) + (1|id), data = c34sum)



c34sum2 <- c34sum
c34sum2$smoist <- scale(c34sum2$logmoist)[, 1]
c34sum2$sc3    <- scale(log(c34sum2$c3 + 1))[, 1]
c34sum2$stemp  <- scale(c34sum2$annual_temp2m)[, 1]


c4d_m8 <- lmer(c4_ddiff ~ co2 * (smoist + stemp + sc3) +
                 (1+smoist|ring) + (1|RY) + (1|id), data = c34sum2)
c4d_m8 <- lmer(c4_ddiff ~ co2 * (smoist + sc3) +
                 (1|ring) + (1|RY) + (1+smoist|id), data = c34sum2)

c4d_m8 <- lmer(c4_ddiff ~ co2 + smoist + sc3 +
                 (1|ring) + (1+smoist|RY) + (1|id), data = c34sum2)
Anova(c4d_m8, test.statistic = "F")

# # summary(c34sum2)
# summary(c34sum)
summary(c4d_m8)
c4d_m6full <- dredge(c4d_m8, REML = F)
c4d_m6full
subset(c4d_m6full, !nested(.))
newmod <- model.avg(get.models(subset(c4d_m6full, !nested(.)), subset = delta < 4))
summary(newmod)
importance(newmod)
confint(newmod, level = .95)

bm <- get.models(c4d_m6full, subset = 1)[[1]]
confint(bm, level = .95, method = "boot", nsim = 500, boot.type = "basic")


