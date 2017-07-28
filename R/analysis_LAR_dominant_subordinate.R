
# here, I analyse the relationship between soil nutrients (IEM-extracted N and
# P) and LAR of dominant/subordinate C3/C4 grass


# assign closest IEMs to vegetation plots ---------------------------------

# probe coordinate
plotcor <- read.csv("Data/ProbeCoordinate.csv")
plotcor_vi <- plotcor %>% 
  filter(Sample %in% c("Ion exchange resin", "vegetation")) %>% 
  mutate(Sample = ifelse(Sample == "vegetation", "vegetation", "IEM"))

iem_pos <- plotcor_vi %>% 
  filter(Sample == "IEM") %>% 
  mutate(iemplot = paste("IEM", Plot, sep = "_")) %>% 
  select(-Sample, -Plot) %>% 
  rename(iem_N = Northing,
         iem_E = Easting)
veg_pos <- filter(plotcor_vi, Sample == "vegetation") 


vi_pos <- right_join(veg_pos, iem_pos) %>% 
  select(Ring, Plot, Sample, iemplot, everything()) %>% 
  arrange(Ring, Plot) %>% 
  mutate(d = sqrt((Northing - iem_N)^2 + (Easting - iem_E)^2)) %>%  # didstance bewteen the vetation plots and each probe
  group_by(Ring, Plot) %>% 
  mutate(mind = d == min(d)) %>%                                    # choose the closest pro e
  ungroup() %>% 
  filter(mind) %>% 
  rename(ring = Ring, plot = Plot) %>% 
  select(ring, plot, iemplot, d) %>% 
  mutate(ring = factor(ring), plot = factor(plot))


# merge with IEM data
load("Data/SoilVariables/FACE_IEM.RData")
iem_raw <- iem %>% 
  filter(time %in% c(5, 6, 7, 12, 13, 14)) %>%                     # use the data only from summer months when eCO2-induced P enhancement was observed
  mutate(Year = ifelse(year(date) == 2012, "Year0", "Year1")) %>% 
  group_by(Year, time, ring, plot) %>% 
  summarise_each(funs(mean), no, nh, p) %>% 
  group_by(Year, ring, plot) %>% 
  summarise_each(funs(mean), no, nh, p) %>% 
  ungroup() %>% 
  gather(key = variable, value, no, nh, p) %>% 
  mutate(yv = paste(Year, variable, sep = "_")) %>%
  select(-Year, -variable) %>% 
  spread(key = yv, value = value) %>% 
  mutate(iemplot = paste("IEM", plot, sep = "_")) %>% 
  select(-plot)


# merge annual change rats of grass species and iem data by probe coordinates
vi_dd <- grass_DS_acr %>% 
  filter(year == "Year1") %>% 
  left_join(vi_pos) %>% 
  left_join(iem_raw) %>%
  ungroup() %>% 
  mutate(Year1_n  = Year1_no + Year1_nh,
         Year1_np = Year1_n / Year1_p)
names(vi_dd)
names(grass_DS_acr)



# analysis ----------------------------------------------------------------

# subordinate C4
plot(s_sc4_ddiff ~ log(Year1_p), data = vi_dd, col = co2, pch = 19)
plot(s_sc4_ddiff ~ log(Year1_n), data = vi_dd, col = co2, pch = 19)
sc4_m1 <- lmer(s_sc4_ddiff ~ log(Year1_p) + log(Year1_n) + (1|ring), data = vi_dd)
Anova(sc4_m1, test.statistic = "F")
plot(sc4_m1)
qqPlot(resid(sc4_m1))


# subordinate C3
plot(s_sc3_ddiff ~ log(Year1_p), data = vi_dd, col = co2, pch = 19)
plot(s_sc3_ddiff ~ log(Year1_n), data = vi_dd, col = co2, pch = 19)
sc3_m1 <- lmer(s_sc3_ddiff ~ log(Year1_p) + log(Year1_n) + (1|ring), data = vi_dd)
Anova(sc3_m1, test.statistic = "F")
plot(sc3_m1)
qqPlot(resid(sc3_m1))


# dominant C4
plot(s_dc4_ddiff ~ log(Year1_p), data = vi_dd, col = co2, pch = 19)
plot(s_dc4_ddiff ~ log(Year1_n), data = vi_dd, col = co2, pch = 19)
dc4_m1 <- lmer(s_dc4_ddiff ~ log(Year1_p) + log(Year1_n) + (1|ring), data = vi_dd)
Anova(dc4_m1, test.statistic = "F")
plot(dc4_m1)
qqPlot(resid(dc4_m1))


# dominant C3
plot(s_dc3_ddiff ~ log(Year1_p), data = vi_dd[-9, ], col = co2, pch = 19)
plot(s_dc3_ddiff ~ log(Year1_n), data = vi_dd[-9, ], col = co2, pch = 19)
dc3_m1 <- lmer(s_dc3_ddiff ~ log(Year1_p) + log(Year1_n) + (1|ring), data = vi_dd[-9, ])
Anova(dc3_m1, test.statistic = "F")
plot(dc3_m1)
qqPlot(resid(dc3_m1))


