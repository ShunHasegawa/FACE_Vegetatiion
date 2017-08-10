
# here, I analyse the relationship between soil nutrients (IEM-extracted N and
# P) and LAR of dominant/subordinate C3/C4 grass


# merge with IEM data
iem <- read.xlsx2(file = "Data/newIem_smmry_tbl.xlsx", sheetIndex = 1, stringsAsFactors = FALSE)
iem <- iem %>% 
  select(-Moist, -Temp_Mean) %>% 
  mutate_each(funs(as.Date), insertion, sampling, date) %>% 
  mutate_each(funs(as.numeric), nitrate, ammonium, phosphate) %>%
  mutate_each(funs(as.factor), ring, plot) %>% 
  filter(time %in% c(6, 7, 13, 14, 16, 19)) %>% 
  mutate(year = ifelse(time %in% c(6, 7), "Year0",
                       ifelse(time %in% c(13, 14), "Year1", 
                              ifelse(time == 16, "Year2", "Year3"))))


iem %>% 
  select(insertion, sampling, time, year) %>% 
  gather(variable, value, insertion, sampling) %>% 
  mutate(day2 = yday(value)) %>% 
  distinct() %>% 
  arrange(as.numeric(as.character(time))) %>%
  mutate(day2 = ifelse(day2 < 100, 365 + day2, day2)) %>% 
  ggplot(., aes(x = day2, y = 1, group = time, col = time))+
  geom_point()+
  geom_path()+
  facet_grid(year ~ .)


iem_raw <- iem %>% 
  rename(no = nitrate, nh = ammonium, p = phosphate) %>% 
  group_by(year, time, ring) %>% 
  summarise_each(funs(mean), no, nh, p) %>% 
  group_by(year, ring) %>% 
  summarise_each(funs(mean), no, nh, p) %>% 
  ungroup() %>% 
  mutate(nitr = no + nh,
         np   = nitr / p)

plot(nitr ~ as.numeric(factor(year)), data = iem_raw, col = ring, pch = 19)
plot(p ~ as.numeric(factor(year)), data = iem_raw, col = ring, pch = 19)


# merge annual change rats of grass species and iem data
str(grass_DS_acr)
vi_dd <- grass_DS_acr %>% 
  group_by(year, ring, co2) %>% 
  summarise_each(funs(mean(.)), ends_with("_ddiff")) %>% 
  left_join(iem_raw) %>%
  ungroup() %>% 
  mutate(nitr      = no + nh,
         np        = nitr / p,
         s_nitr    = scale(nitr)[, 1],
         s_logp    = scale(log(p))[, 1],
         s_lognitr = scale(log(nitr))[, 1],
         s_p    = scale(p)[, 1]
         ) %>% 
  .[complete.cases(.), ]
names(vi_dd)
names(grass_DS_acr)




# analysis ----------------------------------------------------------------

# subordinate C4
plot(S_c4_ddiff ~ s_nitr, data = vi_dd, col = factor(year), pch = 19)
plot(S_c4_ddiff ~ s_lognitr, data = vi_dd, col = factor(year), pch = 19)
plot(S_c4_ddiff ~ s_p, data = vi_dd, col = factor(year), pch = 19)
plot(S_c4_ddiff ~ s_logp, data = vi_dd, col = factor(year), pch = 19)
plot(S_c4_ddiff ~ log(np), data = vi_dd, col = factor(year), pch = 19)


sc4_m1 <- lmer(s_sc4_ddiff ~ s_nitr + s_p + (1|ring) + (1|year), data = vi_dd)
sc4_m1 <- lmer(s_sc4_ddiff ~ s_nitr + s_p + (1|ring) +(1|year), data = vi_dd)
summary(sc4_m1)
Anova(sc4_m1, test.statistic = "F")

plot(sc4_m1)
qqPlot(resid(sc4_m1))
par(mfrow = c(2, 2), mar = c(4, 3, 1, 1))
plot(sc4_m1)


# subordinate C3
plot(S_c3_ddiff ~ log(p), data = vi_dd, col = factor(year), pch = 19)
plot(S_c3_ddiff ~ log(nitr), data = vi_dd[-4, ], col = factor(year), pch = 19)
plot(S_c3_ddiff ~ log(np), data = vi_dd[-4, ], col = factor(year), pch = 19)
sc3_m1 <- lmer(s_sc3_ddiff ~ log(p) + log(nitr) + (1|ring)+ (1|year), data = vi_dd[-4, ])
Anova(sc3_m1, test.statistic = "F")
summary(sc3_m1)
plot(sc3_m1)
qqPlot(resid(sc3_m1))
visreg(sc3_m1, xvar = "nitr", xtrans = log)


# dominant C4
plot(s_dc4_ddiff ~ log(p), data = vi_dd, col = co2, pch = 19)
plot(s_dc4_ddiff ~ log(nitr), data = vi_dd, col = co2, pch = 19)
dc4_m1 <- lmer(s_dc4_ddiff ~ log(p) + log(nitr) + (1|year) + (1|ring), data = vi_dd[-4, ])
Anova(dc4_m1, test.statistic = "F")
plot(dc4_m1)
qqPlot(resid(dc4_m1))
which.min(resid(dc4_m1))
qqPlot(resid(dc4_m1))


# dominant C3
plot(s_dc3_ddiff ~ log(p), data = vi_dd, col = co2, pch = 19)
plot(s_dc3_ddiff ~ p, data = vi_dd, col = co2, pch = 19)
plot(s_dc3_ddiff ~ log(nitr), data = vi_dd, col = co2, pch = 19)
plot(s_dc3_ddiff ~ nitr, data = vi_dd[-3, ], col = co2, pch = 19)
dc3_m1 <- lmer(s_dc3_ddiff ~ p + nitr + (1|ring) + (1|year), data = vi_dd[-3, ])
Anova(dc3_m1, test.statistic = "F")
plot(dc3_m1)
qqPlot(resid(dc3_m1))
which.min(resid(dc3_m1))



# summary -----------------------------------------------------------------

lar_sd_c34_soil_ml <- list(subordinate_c4 = sc4_m1, 
                           dominant_c4    = dc4_m1,
                           subordinate_c3 = sc3_m1,
                           dominant_c3    = dc3_m1)

lar_sd_c34_soil_aov <- ldply(lar_sd_c34_soil_ml, function(x) tidy(Anova(x, test.statistic = "F"))) %>% 
  mutate(term = mapvalues(term, c("log(p)", "log(nitr)"), c("IEM-P", "IEM-N")),
         p.value = round(p.value, 3),
         statistic = round(statistic, 2)) %>% 
  rename(Fval = statistic)