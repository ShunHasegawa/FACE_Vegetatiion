
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


# subordinate c4  -------------------------------------------------------------

lar_sc4_m1 <- lmer(s_sc4_ddiff ~ co2 * (s_logmoist+s_temp + s_logpar) +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
lar_sc4_m1_full <- dredge(lar_sc4_m1, fixed = c("co2", "s_logmoist", "s_temp",  "s_logpar"))
plot(lar_sc4_m1)
qqPlot(resid(lar_sc4_m1))
# no interaction is suggested
lar_sc4_m2 <- lmer(s_sc4_ddiff ~ co2 + s_logmoist+s_temp + s_logpar +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
Anova(lar_sc4_m2, test.statistic = "F")
plot(lar_sc4_m2)
qqPlot(resid(lar_sc4_m2))


# dominant c4 -------------------------------------------------------------
lar_dc4_m1 <- lmer(s_dc4_ddiff ~ co2 * (s_logmoist+s_temp + s_logpar) +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
lar_dc4_m1_full <- dredge(lar_dc4_m1, fixed = c("co2", "s_logmoist", "s_temp",  "s_logpar"))
lar_dc4_m1_full
plot(lar_dc4_m1)
qqPlot(resid(lar_dc4_m1))
# no interaction is suggested
lar_dc4_m2 <- lmer(s_dc4_ddiff ~ co2 + s_logmoist+s_temp + s_logpar +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
Anova(lar_dc4_m2, test.statistic = "F")




# subordinate c3 ----------------------------------------------------------

lar_sc3_m1 <- lmer(s_sc3_ddiff ~ co2 * (s_logmoist+s_temp + s_logpar) +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
lar_sc3_m1_full <- dredge(lar_sc3_m1, fixed = c("co2", "s_logmoist", "s_temp",  "s_logpar"))
lar_sc3_m1_full
plot(lar_sc3_m1)
qqPlot(resid(lar_sc3_m1))
# no interaction is suggested
lar_sc3_m2 <- lmer(s_sc3_ddiff ~ co2 + s_logmoist+s_temp + s_logpar +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
Anova(lar_sc3_m2, test.statistic = "F")
plot(lar_sc3_m2)
qqPlot(resid(lar_sc3_m2))


# dominant c3 ----------------------------------------------------------

lar_dc3_m1 <- lmer(s_dc3_ddiff ~ co2 * (s_logmoist+s_temp + s_logpar) +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
plot(lar_dc3_m1)
qqPlot(resid(lar_dc3_m1))
lar_dc3_m2 <- update(lar_dc3_m1, subset = -which.min(resid(lar_dc3_m1)))
plot(lar_dc3_m2)
qqPlot(resid(lar_dc3_m2))
lar_dc3_m1_full <- dredge(lar_dc3_m2, fixed = c("co2", "s_logmoist", "s_temp",  "s_logpar"))
lar_dc3_m1_full
# no interaction is suggested

lar_dc3_m3 <- lmer(s_dc3_ddiff ~ co2 + s_logmoist+s_temp + s_logpar +(1|ring)+(1|RY)+(1|id), data = grass_DS_acr)
plot(lar_dc3_m3)
qqPlot(resid(lar_dc3_m3))

# potential outliers are suggested
mcp.fnc(lar_dc3_m3)
oldf_dc3 <- romr.fnc(lar_dc3_m3, data = data.frame(grass_DS_acr))
dplyr::setdiff(oldf_dc3$data0, oldf_dc3$data)
  # three outliers are indicated

# remove the outlier
lar_dc3_m4 <- update(lar_dc3_m3, data = oldf_dc3$data)
plot(lar_dc3_m4)
qqPlot(resid(lar_dc3_m4))
Anova(lar_dc3_m4, test.statistic = "F")




# summary -----------------------------------------------------------------
lar_sd_c34_ml <- list(subordinate_c4 = lar_sc4_m2, 
                      subordinate_c3 = lar_sc3_m2, 
                      dominant_c4    = lar_dc4_m2, 
                      dominant_c3    = lar_dc3_m3)

lar_sd_c34_aov <- ldply(lar_sd_c34_ml, function(x) tidy(Anova(x, test.statistic = "F"))) %>% 
  mutate(term = mapvalues(term, c("s_logmoist", "s_temp", "s_logpar"), c("Moist", "Temp", "PAR")),
         Df.res = round(Df.res, 0),
         p.value = round(p.value, 3),
         statistic = round(statistic, 2)) %>% 
  rename(Fval = statistic)
