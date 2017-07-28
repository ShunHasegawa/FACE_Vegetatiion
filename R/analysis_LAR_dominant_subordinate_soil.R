
# here, I analyse the relationship between soil nutrients (IEM-extracted N and
# P) and LAR of dominant/subordinate C3/C4 grass


# merge with IEM data
load("Data/SoilVariables/FACE_IEM.RData")
iem_raw <- iem %>% 
  filter(time %in% c(5, 6, 7, 12, 13, 14)) %>%                     # use the data only from summer months when eCO2-induced P enhancement was observed
  mutate(Year = ifelse(year(date) == 2012, "Year0", "Year1")) %>% 
  group_by(Year, time, ring) %>% 
  summarise_each(funs(mean), no, nh, p) %>% 
  group_by(Year, ring) %>% 
  summarise_each(funs(mean), no, nh, p) %>% 
  ungroup() %>% 
  gather(key = variable, value, no, nh, p) %>% 
  mutate(yv = paste(Year, variable, sep = "_")) %>%
  select(-Year, -variable) %>% 
  spread(key = yv, value = value)

# merge annual change rats of grass species and iem data by probe coordinates
names(grass_DS_acr)
vi_dd <- grass_DS_acr %>% 
  filter(year == "Year1") %>% 
  group_by(year, ring, co2) %>% 
  summarise_each(funs(mean(.)), ends_with("_ddiff")) %>% 
  left_join(iem_raw) %>%
  ungroup() %>% 
  mutate(Year1_n  = Year1_no + Year1_nh,
         Year1_np = Year1_n / Year1_p)
names(vi_dd)
names(grass_DS_acr)




# analysis ----------------------------------------------------------------

# subordinate C4
plot(S_c4_ddiff ~ log(Year1_p), data = vi_dd, col = co2, pch = 19)
plot(S_c4_ddiff ~ log(Year1_n), data = vi_dd, col = co2, pch = 19)
sc4_m1 <- lm(s_sc4_ddiff ~ log(Year1_p) + log(Year1_n), data = vi_dd)
Anova(sc4_m1, test.statistic = "F")
par(mfrow = c(2, 2), mar = c(4, 3, 1, 1))
plot(sc4_m1)


# subordinate C3
plot(S_c3_ddiff ~ log(Year1_p), data = vi_dd, col = co2, pch = 19)
plot(S_c3_ddiff ~ log(Year1_n), data = vi_dd, col = co2, pch = 19)
sc3_m1 <- lm(s_sc3_ddiff ~ log(Year1_p) + log(Year1_n), data = vi_dd)
Anova(sc3_m1, test.statistic = "F")
plot(sc3_m1)


# dominant C4
plot(s_dc4_ddiff ~ log(Year1_p), data = vi_dd, col = co2, pch = 19)
plot(s_dc4_ddiff ~ log(Year1_n), data = vi_dd, col = co2, pch = 19)
dc4_m1 <- lmer(s_dc4_ddiff ~ log(Year1_p) + log(Year1_n) + (1|ring), data = vi_dd)
Anova(dc4_m1, test.statistic = "F")
plot(dc4_m1)
qqPlot(resid(dc4_m1))


# dominant C3
plot(s_dc3_ddiff ~ log(Year1_p), data = vi_dd, col = co2, pch = 19)
plot(s_dc3_ddiff ~ log(Year1_n), data = vi_dd, col = co2, pch = 19)
dc3_m1 <- lm(s_dc3_ddiff ~ log(Year1_p) + log(Year1_n), data = vi_dd[-9, ])
Anova(dc3_m1, test.statistic = "F")
plot(dc3_m1)


