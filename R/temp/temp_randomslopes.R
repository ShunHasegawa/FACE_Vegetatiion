names(c34sum)
# c4d_m0_slp_9   <- lmer(c4_ddiff ~ co2 * (log(c4moist)+ s_c3 +s_temp) + (0+s_logmiost+s_c3|ring/id) + (1+s_logmiost + s_c3|year:ring), data = c34sum)
c4d_m0_slp_9   <- lmer(c4_ddiff ~ co2 * (log(c4moist)+ log(c3)) + (1+s_logmiost+s_c3|ring) + (0+s_logmiost+s_c3|id:ring) + (0+s_logmiost + s_c3|year:ring), data = c34sum)
c4d_m0_slp_9   <- lmer(c4_ddiff ~ co2 * (log(c4moist)+ log(c3)) + (1+s_logmiost+s_c3|ring) + (0+s_logmiost+s_c3|id:ring) + (0 + s_c3|year:ring), data = c34sum)
# c4d_m0_slp_10   <- lmer(c4_ddiff ~ co2 * (log(c4moist)+ log(c3) +annual_temp2m) + (1+s_logmiost+s_c3|ring/id) + (1+s_logmiost + s_c3|year:ring), data = c34sum)

c4d_m0_slp_9   <- lmer(c4_ddiff ~ co2 * (log(totalmoist)+ log(c3)) + 
                         (0+log(totalmoist)+log(c3)|ring) + (0+log(totalmoist)+log(c3)|id:ring) + (0+log(totalmoist) + log(c3)|year:ring), data = c34sum)

?scale
plot(s_logmiost ~ log(totalmoist), data = c34sum)
aaaa <- (log(c34sum$c3) - mean(log(c34sum$c3)))/sd(c34sum$c3)
aaaa <- scale(log(c34sum$c3))[, 1]
c34sum$s_c3
plot(s_c3 ~ log(c3), data = c34sum)
plot(aaaa ~ log(c3), data = c34sum)

c4d_m0_slp_9   <- lmer(c4_ddiff ~ co2*(s_logmiost+s_c3)+
                         (1           +s_c3|year) +
                         (0+s_logmiost     |id:ring) + 
                         (0           +s_c3|year:ring), data = c34sum)

# c4d_m0_slp_9   <- lmer(c4_ddiff ~ 1+
#                          (0+s_logmiost+s_c3|ring) +
#                          (1+s_logmiost+s_c3|year) +
#                          (0+s_logmiost+s_c3|id:ring) + 
#                          (1+s_logmiost+s_c3|year:ring), data = c34sum)
summary(c4d_m0_slp_9)
# c4d_m0_slp_9   <- lmer(c4_ddiff ~ co2*(s_logmiost+s_c3)+
#                          (0+s_logmiost+s_c3|ring) +
#                          (0+s_logmiost|id:ring) + 
#                          (0+s_c3|year:ring), data = c34sum)

c34sum$s_c4_ddiff <- scale(c34sum$c4_ddiff)[,1]
c4d_m0_slp_9   <- lmer(s_c4_ddiff ~ co2*(s_logmiost+s_c3)+
                         # (1|ring) +
                         (0+s_logmiost|id:ring)+
                         (1+s_c3|year:ring),
                       data = c34sum)
# c4d_m0_slp_92   <- lmer(c4_ddiff ~ co2*(s_logmiost+s_c3)+
#                          (1|ring) +
#                          (1|id:ring) +
#                          (1|year:ring), data = c34sum)

# library(arm)
# std_m9 <- standardize(c4d_m0_slp_9, standardize.y = TRUE)
# std_m92 <- standardize(c4d_m0_slp_92, standardize.y = FALSE)
summary(c4d_m0_slp_9)
r.squaredGLMM(c4d_m0_slp_9)
# summary(c4d_m0_slp_9)
std_m9full <- dredge(c4d_m0_slp_9, REML = F)
mbs <- get.models(std_m9full, subset = 1)[[1]]
# std_m9nest <- subset(std_m9full, !nested(.))
?nested
# std_m92full <- dredge(std_m92, REML = F)
plot(std_m9)
qqnorm(resid(std_m9))
qqline(resid(std_m9))

m10_avg <- model.avg(get.models(std_m9full, subset = delta <= 2))
# m10_avg1 <- model.avg(get.models(std_m92full, subset = delta <= 2))
summary(m10_avg)
# summary(m10_avg1)
confint(m10_avg, level = .95)
summary(mbs)

ggplot(c34sum, aes(x = s_c3, y = s_c4_ddiff, col = ring))+
  geom_point()+
  geom_smooth(method = "lm", se = F)+
  facet_wrap(~ year)



summary(c4d_m0_slp_9)
plot(c4d_m0_slp_9)
qqnorm(resid(c4d_m0_slp_9))
qqline(resid(c4d_m0_slp_9))
which.min(resid(c4d_m0_slp_9))
c4d_m0_slp_10 <- update(c4d_m0_slp_9, subset = -40)
qqnorm(resid(c4d_m0_slp_10))
qqline(resid(c4d_m0_slp_10))
summary(c4d_m0_slp_10)
m10full <- dredge(c4d_m0_slp_9, REML = F)
m10_avg <- model.avg(get.models(m10full, subset = delta <= 4))
summary(m10_avg)
m10bs <- get.models(m10full, subset = 1)[[1]]
confint(m10_avg)
confint(m10_avg)
coef(m10_avg)
predict(m10_avg)
summary(m10bs)
r.squaredGLMM(m10bs)
r.squaredGLMM(c4d_m0_slp_10)

# c4d_m0_slp_10  <- lmer(c4_ddiff ~ co2 + log(c4moist) + log(c3)  + (0+s_logmiost+s_c3|ring/id) + (1+s_logmiost + s_c3|year:ring), data = c34sum)
# c4d_m0_slp_11  <- lmer(c4_ddiff ~ co2                + log(c3)  + (0+s_logmiost+s_c3|ring/id) + (1+s_logmiost + s_c3|year:ring), data = c34sum)
# c4d_m0_slp_12  <- lmer(c4_ddiff ~ co2 + log(c4moist)            + (0+s_logmiost+s_c3|ring/id) + (1+s_logmiost + s_c3|year:ring), data = c34sum)
# c4d_m0_slp_13  <- lmer(c4_ddiff ~ co2 +                         + (0+s_logmiost+s_c3|ring/id) + (1+s_logmiost + s_c3|year:ring), data = c34sum)
# summary(c4d_m0_slp_11)
# Anova(c4d_m0_slp_11, test.statistic = "F")
# AICc(c4d_m0_slp_11)
# r.squaredGLMM(c4d_m0_slp_11)

plot(c4d_m0_slp_9)
qqnorm(resid(c4d_m0_slp_9))
qqline(resid(c4d_m0_slp_9))

anova(c4d_m0_slp_9, c4d_m0_slp_10, c4d_m0_slp_11, c4d_m0_slp_12, c4d_m0_slp_13)
aful <- dredge(c4d_m0_slp_9, REML = F, extra = "r.squaredGLMM")

a <- pdredge(c4d_m0_slp_9, REML = F)
# model0, 2, and 10 caused convergence issues

mb <- get.models(aful, subset = 1)[[1]]
summary(mb)

m0 <- get.models(a, row.names(.) == 0)[[1]]
summary(m0)
m10 <- get.models(a, row.names(.) == 10)[[1]]

summary(m10)
m0_2 <- lmer(c4_ddiff ~ 1 + (0+s_logmiost+s_c3|ring/id) + (1+s_logmiost+s_c3|year:ring), data = c34sum, REML = F)
summary(m10_2)
anova(m10, m10_2)

mbb   <- lmer(c4_ddiff ~ co2 + log(c4moist) + log(c3)  + (0+s_logmiost+s_c3|ring/id) + (0+s_logmiost+s_c3|year:ring), data = c34sum, REML = F)
model.sel(m10_2, mbb)



summary(m2)
m2_2   <- lmer(c4_ddiff ~ s_c3 + (0+s_logmiost+s_c3|id) + (1+s_logmiost + s_c3|year:ring), data = c34sum)
summary(m2_2)

m10_2  <- update(m10, control = lmerControl(optimizer ="Nelder_Mead"))

summary(m0)
summary(m2)
summary(m10)



a

print(a, warnings = F)
dredge(c4d_m0_slp_10, REML = F)
