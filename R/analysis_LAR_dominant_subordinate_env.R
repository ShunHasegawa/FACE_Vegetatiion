
# Here, I analyse LAR of dominant/subordinate C3/C4 species against
# environmental variables



# prepare data frame ------------------------------------------------------


# compute annual change rate (acr) for each group
grass_DS_acr <- grass_DS %>%
  mutate(type_pfg = paste(type, PFG, sep = "_")) %>% 
  group_by(year, block, ring, plot, co2, id, RY, type_pfg, variable) %>% 
  summarise(value = sum(value)) %>%
  mutate(value = log(value + 1)) %>% 
  group_by(year, block, ring, plot, co2, id, RY, type_pfg) %>% 
  summarise(value = sum(value)) %>% 
  spread(key = type_pfg, value = value) %>% 
  ungroup() %>% 
  arrange(id, year) %>%
  group_by(id) %>%
  mutate_each(funs(ddiff = . - lag(., 1)),
              D_c3, D_c4, S_c3, S_c4) %>%
  filter(year != "Year0") %>%
  left_join(c34growth_moist) %>% 
  ungroup() %>% 
  mutate(s_dc4_ddiff = scale(D_c4_ddiff)[, 1],
         s_dc3_ddiff = scale(D_c3_ddiff)[, 1],
         s_sc4_ddiff = scale(S_c4_ddiff)[, 1],
         s_sc3_ddiff = scale(S_c3_ddiff)[, 1],
         s_logmoist  = scale(log(totalmoist))[, 1],
         s_temp      = scale(annual_temp2m)[, 1],
         s_logpar    = scale(log(PAR))[, 1],
         s_dc3       = scale(log(D_c3 + 1))[, 1],
         s_dc4       = scale(log(D_c4 + 1))[, 1],
         s_sc3       = scale(log(S_c3 + 1))[, 1],
         s_sc4       = scale(log(S_c4 + 1))[, 1])

plot(D_c4_ddiff ~ D_c3, data = grass_DS_acr)
plot(D_c4_ddiff ~ S_c3, data = grass_DS_acr)


# m0 <- lmer(s_dc4_ddiff ~ s_dc3 + s_sc4  + s_sc3 + s_logmoist+s_temp + s_logpar +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
m0 <- lmer(s_sc4_ddiff ~ s_dc3 + s_dc4  + s_sc3 + s_logmoist+s_temp + s_logpar +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
m0 <- lmer(s_sc4_ddiff ~ s_logmoist+s_temp + s_logpar +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
# m0 <- lmer(s_sc3_ddiff ~ s_dc3 + s_dc4  + s_sc4 + s_logmoist+s_temp + s_logpar +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr[-32, ])
# m0 <- lmer(s_dc3_ddiff ~ co2 + s_sc3 + s_dc4  + s_sc4 + s_logmoist+s_temp + s_logpar +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr[-25, ])
xvars <- c("s_dc3", "s_dc4", "s_sc3", "s_sc4", "s_logmoist", "s_temp", "s_logpar")
par(mfrow = c(2, 4))
l_ply(xvars[-1:-4], function(x) visreg(m0, xvar = x))
m1 <- lmer(I(s_logmoist+s_temp+s_logpar) ~ co2 * year + (1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
Anova(m1, test.statistic = 'F')

# m0 <- lmer(s_dc3_ddiff ~ co2 + s_logmoist+s_temp+s_logpar+(1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
# m0 <- lmer(s_sc4_ddiff ~ co2 + s_logmoist+s_temp+s_logpar + s_dc3 + s_dc4 + (1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
m0_full <- dredge(m0)



m0 <- lmer(D_c4_ddiff ~ co2 * (s_logmoist + s_temp + s_logpar) + (1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
m0 <- lmer(S_c4_ddiff ~ co2 * (s_logmoist + s_temp + s_logpar) + (1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
m0 <- lmer(D_c3_ddiff ~ co2 * (s_logmoist + s_temp + s_logpar) + (1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
m0 <- lmer(S_c3_ddiff ~ co2 + s_logmoist + s_temp + s_logpar + (1|ring)+(1|RY)+(1|id), data = grass_DS_acr)


m0 <- lmer(S_c4_ddiff ~ totalmoist + annual_temp2m + PAR + (1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
# m0 <- lmer(S_c3_ddiff ~ co2 * (s_logmoist + s_temp + s_logpar) + (1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
# m0 <- lmer(D_c3_ddiff ~ co2 * (s_logmoist + s_temp + s_logpar) + (1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
Anova(m0, test.statistic = "F")
summary(m0)
mfull <- dredge(m0)
mbest <- get.models(mfull, subset = 1)[[1]]
Anova(mbest, test.statistic = "F")
visreg(mbest, xvar = "annual_temp2m")
visreg(mbest, xvar = "PAR")

r.squared(m0)
rmean <- grass_DS_acr %>% 
  group_by(year, co2) %>% 
  summarise_each(funs(mean), totalmoist, annual_temp2m, PAR)
# rmean <- cbind(rmean, grass_DS_acr[1, c("ring", "RY", "id")])

bb <- bootMer(m0, 
              FUN = function(x) predict(x, rmean, re.form = NA),
              nsim = 499)
lci     <- apply(bb$t, 2, quantile, 0.05)
uci     <- apply(bb$t, 2, quantile, 0.95)
predval <- bb$t0
rmean %>% 
  bind_cols(data.frame(lci, uci, predval))


BtsCI <- function(model, MoistVal, TempVal, variable){
  expDF <- data.frame(co2 = c("amb", "elev"),
                      Moist = rep(MoistVal, 2),
                      Temp_Mean = rep(TempVal, 2))
  bb <- bootMer(model,
                FUN=function(x) predict(x, expDF, re.form = NA),
                nsim=500)
  lci <- apply(bb$t, 2, quantile, 0.025)
  uci <- apply(bb$t, 2, quantile, 0.975)
  PredVal <- bb$t0
  df <- cbind(lci, uci, PredVal, expDF, variable)
  return(df)
} 












?confint
pd <- predict(m0, rmean)
?predict
?predict
rmean$pred <- pd

Anova(m0, test.statistic = "F")
qqPlot(resid(m0))
summary(m0)
which.min(qqnorm(resid(m0))$y)

