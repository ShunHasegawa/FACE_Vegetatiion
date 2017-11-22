newsoil2 <- c43_ratio_iem %>% 
  filter(year != "Year0") %>% 
  select(year, ring, s_n, s_p) %>% 
  left_join(grass_DS_acr) %>% 
  left_join(c43_ratio_year0)

head(c43_ratio_year0)


head(newsoil2)

names(newsoil2)


nm1 <- lmer(c43_r ~ s_logmoist + s_temp + s_logpar + s_n + s_p + (1|ring) + (1|year), data = newsoil2)
summary(nm1)
Anova(nm1, test.statistic = "F")

nm1 <- lmer(s_sc4_ddiff ~ s_logmoist + s_temp + s_logpar + s_n + s_p + (1|ring) + (1|year), data = newsoil2)
nm2 <- lmer(s_dc4_ddiff ~ s_logmoist + s_temp + s_logpar + s_n + s_p + (1|ring) + (1|year), data = newsoil2)



nm2 <- lmer(s_dc4_ddiff ~ s_logmoist + s_temp + s_logpar + s_n + s_p + (1|ring) + (1|year), data = newsoil2)
summary(nm2)



Anova(nm2, test.statistic = "F")

nm3 <- lmer(s_sc3_ddiff ~ s_logmoist + s_temp + s_logpar + s_n + s_p + (1|ring) + (1|year), data = newsoil2)
nm4 <- lmer(s_dc3_ddiff ~ s_logmoist + s_temp + s_logpar + s_n + s_p + (1|ring) + (1|year), data = newsoil2)
Anova(nm1, test.statistic = "F")
Anova(nm3, test.statistic = "F")
Anova(nm4, test.statistic = "F")





newsoil <- c43_ratio_iem %>% 
  filter(year != "Year0") %>% 
  left_join(c34growth_moist) %>% 
  left_join(grass_DS_dd_ed) %>% 
  left_join(select(c43_ratio_year0, -c43_r)) %>% 
  mutate(s_moist = scale(totalmoist)[, 1],
         s_logmoist = scale(log(totalmoist))[, 1],
         s_temp  = scale(annual_temp2m)[, 1],
         s_par   = scale(PAR)[, 1],
         s_logpar   = scale(log(PAR))[, 1],
         s_y0 = scale(ratios0)[, 1])




c43_soil_m21 <- lmer(I(scale(c4)[, 1]) ~ s_n + s_p + s_moist + s_temp + s_par  + (1|ring) + (1|year), data = newsoil)
c43_soil_m22 <- lmer(I(scale(c4)[, 1]) ~ s_n + s_p + s_logmoist + s_temp + s_logpar + (1|ring) + (1|year), data = newsoil)
r.squared(c43_soil_m21)
r.squared(c43_soil_m22)
summary(c43_soil_m22)
plot(c43_soil_m22)
qqPlot(resid(c43_soil_m22))







names(c34growth_moist)
c43_soil_m31 <- lmer(s_c43_r ~ s_n + s_p + s_moist + s_temp + s_par  + (1|ring) + (1|year), data = newsoil)
c43_soil_m32 <- lmer(s_c43_r ~ s_n + s_p + s_logmoist + s_temp + s_logpar + (1|ring) + (1|year), data = newsoil)
r.squared(c43_soil_m31)
r.squared(c43_soil_m32)

summary(c43_soil_m31)
summary(c43_soil_m32)

r.squared(c43_soil_m3)

summary(c43_soil_m3)
c43_soil_m3_full <- dredge(c43_soil_m3, REML  = F)
importance(c43_soil_m3_full)

vif(c43_soil_m3)


names(newsoil)
c43_soil_m41 <- lmer(S_c3 ~ s_n + s_p + s_moist + s_temp + s_par + (1|ring) + (1|year), data = newsoil)
c43_soil_m42 <- lmer(S_c3 ~ s_n + s_p + s_logmoist + s_temp + s_logpar + (1|ring) + (1|year), data = newsoil)
r.squared(c43_soil_m41)
r.squared(c43_soil_m42)
summary(c43_soil_m41)
summary(c43_soil_m42)

plot(c43_soil_m41)
plot(c43_soil_m42)

qqPlot(resid(c43_soil_m41))
qqPlot(resid(c43_soil_m42))

mcp.fnc(c43_soil_m42)
oldf <- romr.fnc(c43_soil_m41, data.frame(newsoil))
c43_soil_m43 <- update(c43_soil_m41, data = oldf$data)
mcp.fnc(c43_soil_m43)


plot(c43_soil_m43)
qqPlot(resid(c43_soil_m43))
summary(c43_soil_m43)
Anova(c43_soil_m43, test.statistic = "F")
visreg(c43_soil_m6, xvar = "s_moist")
visreg(c43_soil_m6, xvar = "s_par")




c43_soil_m51 <- lmer(D_c3 ~ s_n + s_p + s_moist + s_temp + s_par + (1|ring) + (1|year), data = newsoil)
c43_soil_m52 <- lmer(D_c3 ~ s_n + s_p + s_logmoist + s_temp + s_logpar + (1|ring) + (1|year), data = newsoil)
r.squared(c43_soil_m51)
r.squared(c43_soil_m52)
plot(c43_soil_m5)
qqPlot(resid(c43_soil_m5))
summary(c43_soil_m51)
summary(c43_soil_m52)
visreg(c43_soil_m52, xvar = "s_n")


c43_soil_m61 <- lmer(S_c4 ~ s_n + s_p + s_moist + s_temp + s_par + (1|ring) + (1|year), data = newsoil)
c43_soil_m62<- lmer(S_c4 ~ s_n + s_p + s_logmoist + s_temp + s_logpar + (1|ring) + (1|year), data = newsoil)
r.squared(c43_soil_m61)
r.squared(c43_soil_m62)

plot(c43_soil_m6)
qqPlot(resid(c43_soil_m6))
summary(c43_soil_m61)
summary(c43_soil_m62)
visreg(c43_soil_m62, xvar = "s_n")


c43_soil_m71 <- lmer(D_c4 ~ s_n + s_p + s_moist + s_temp + s_par + (1|ring) + (1|year), data = newsoil)
c43_soil_m72 <- lmer(D_c4 ~ s_n + s_p + s_logmoist + s_temp + s_logpar + (1|ring) + (1|year), data = newsoil)
r.squared(c43_soil_m71)
r.squared(c43_soil_m72)

plot(c43_soil_m72)
qqPlot(resid(c43_soil_m72))
summary(c43_soil_m71)
summary(c43_soil_m72)
visreg(c43_soil_m8, xvar = "s_logmoist")
visreg(c43_soil_m8, xvar = "s_temp")






summary(c43_soil_m5)
summary(c43_soil_m6)
summary(c43_soil_m7)

Anova(c43_soil_m4, test.statistic = "F")
Anova(c43_soil_m5, test.statistic = "F")
Anova(c43_soil_m6, test.statistic = "F")
Anova(c43_soil_m7, test.statistic = "F")





visreg(c43_soil_m3, xvar = "s_n")
visreg(c43_soil_m6, xvar = "s_n")




plot(c43_soil_m3)
qqPlot(resid(c43_soil_m3))
summary(c43_soil_m3)
dredge(c43_soil_m3, REML = F)

confint(c43_soil_m3, method = "boot", nsim = 99)


Anova(c43_soil_m2, test.statistic = "F")

grass_DS_dd_ed <- grass_DS_dd %>% 
  select(-PFG, -type, -value0, -co2) %>% 
  spread(pfg_type, value)




gh_m1 <- lmer(value ~ co2 * year * value0 + (1|ring), data = gh_d)
gj_m1 <- lmer(value ~ co2 * year * value0 + (1|ring), data = gj_d)
Anova(gj_m1, test.statistic = "F")
