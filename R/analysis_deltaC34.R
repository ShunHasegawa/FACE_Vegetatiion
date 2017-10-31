# Here we will:
# 1) download moisture and understorey air temperature data from HIEv
# 2) identify growing season for each of C3 and C4 species according to Murphy 2007
# 3) Fit moisture and temperature to delta C3 and C4 (annual change in abundance)



# download from HIEv ------------------------------------------------------


# # setToken(tokenfile = "Data/token.txt")
# # 
# # airvar_raw <- downloadTOA5(filename = "FACE_.*_AirVars_.*dat",
# #                           topath    = "Data/hievdata/raw_data/",
# #                           maxnfiles = 999)
# # 
# # save(airvar_raw, file = "output/Data/FACE_airvar_raw.RData")
# load("output/Data/FACE_airvar_raw.RData")
# 
# 
# 
# 
# # process HIEv data -------------------------------------------------------
# 
# 
# # get daily mean, min and max temparater at two layers (at a height of 2m and 15m)
# boxplot(airvar_raw$AirTC_1_Avg)
# airvar_day <- airvar_raw %>%
#   filter(Date        > as.Date("2012-09-18"),         # when co2 treatment was commenced
#          Date        < as.Date("2016-2-15"),
#          AirTC_1_Avg > -10) %>%                       # remove obviously weird value
#   distinct() %>%                                      # remove duplicates
#   mutate(ring = factor(substring(Source, 7, 7)),      # add ring number
#          year = factor(year(Date))) %>%
#   group_by(year, Date, ring) %>% 
#   summarise_each(funs(Mean = mean(., na.rm = TRUE), 
#                       Min  = min(., na.rm = TRUE),
#                       Max  = max(., na.rm = TRUE),
#                       N    = sum(!is.na(.))), 
#                  airtemp2m = AirTC_1_Avg, airtemp15m = AirTC_2_Avg) %>% 
#   group_by(year, ring) %>% 
#   mutate(annual_temp2m = mean(airtemp2m_Mean),
#          T = max(27, 1.745 * annual_temp2m  + 11.143)[1]) %>% 
#   ungroup() %>% 
#   mutate(c3growth = airtemp2m_Min >= -1 & airtemp2m_Max >= 10 & airtemp2m_Max < 24,  # c3 growing dates
#          c4growth = airtemp2m_Max >= 21 & airtemp2m_Max < T & month(Date))           # c4 growtin dates
# save(airvar_day, file = "output/Data/FACE_air_variables.RData")

load("output/Data/FACE_air_variables.RData")  # airvar_day; run the above to obtain up-to-date data from HIEv

c34growthdate <- airvar_day %>% 
  select(year, Date, ring, c3growth, c4growth, annual_temp2m)


## understorey temperature
annualtemp <- c34growthdate %>% 
  select(year, ring, annual_temp2m) %>% 
  distinct() %>% 
  filter(!year %in% c("2012", "2016")) %>% 
  mutate(year = mapvalues(year, c(2013, 2014, 2015), paste0("Year", 1:3)))


## understorey PAR
load("output/Data/FACE_FloorPAR.RData")

Lightdf$PAR <- rowMeans(Lightdf[, c("PAR_Den_1_Avg", "PAR_Den_2_Avg", "PAR_Den_3_Avg")], na.rm = TRUE)
underpar <- Lightdf %>% 
  mutate(year = year(DateTime)) %>% 
  filter(year %in% c(2013:2015)) %>% 
  group_by(year, ring) %>% 
  summarise(PAR = mean(PAR, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(year = factor(year, labels = paste0("Year", 1:3)))

## Murphy 2007. For C 3 grasses, growth was considered possible in any month 
## where the daily minimum temperature was ≥ −1 oC, and the daily maximum 
## temperature was ≥ 10 oC and < 24 oC. For C 4 grasses, growth was considered 
## possible in months where the daily maximum temperature was ≥ 21 °C and less 
## than the tem- perature define (Murphy and Bowman 2007).




# assign moisture data to each vegetation plot based on their corrdinates --------


# soil moisture for each vegetation plot assined based probe coordinates
load("Data/FACE_TDR_ProbeDF.RData")
veg_moist <- FACE_TDR_ProbeDF %>%
  filter(Date        > as.Date("2012-09-18"),         # when co2 treatment was commenced
         Date        < as.Date("2016-2-15"),
         Sample == "vegetation") %>% 
  mutate(year = factor(year(Date)), plot = factor(plot)) %>% 
  select(year, Date, ring, plot, Moist) 


# merge
summary(c34growthdate)
summary(veg_moist)
summary(Lightdf)


## survey dates: 2012-12-15 (Year0), 2014-1-15 (Year1), 2015-1-30 (Year2), 2016 (Year3)

c34growth_moist <- veg_moist %>% 
  filter(year %in% 2013:2015) %>% 
  mutate(year = factor(year, labels = paste0("Year", 1:3))) %>% 
  group_by(year, ring) %>%
  summarise(totalmoist = mean(Moist)) %>%
  ungroup() %>%
  left_join(annualtemp) %>% 
  left_join(underpar)




# abundance change ------------------------------------------------------------


c34sum <- C3grassC4 %>%
  group_by(year, block, ring, plot, co2, id, PFG, RY, variable) %>% 
  summarise(value = sum(value)) %>%
  mutate(value = log(value + 1)) %>% 
  group_by(year, block, ring, plot, co2, id, PFG, RY) %>% 
  summarise(value = sum(value)) %>% 
  spread(key = PFG, value = value) %>% 
  ungroup() %>% 
  arrange(id, year) %>%
  group_by(id) %>%
  mutate_each(funs(ddiff = . - lag(., 1)),             # Year1-Year0 and etc.
                   c3, c4) %>%
  filter(year != "Year0") %>%
  group_by(year, ring, co2) %>% 
  summarise_each(funs(mean), ends_with("ddiff")) %>% 
  left_join(c34growth_moist) %>% 
  ungroup() %>% 
  mutate(s_c4_ddiff = scale(c4_ddiff)[, 1],
         s_c3_ddiff = scale(c3_ddiff)[, 1],
         s_logmoist = scale(log(totalmoist))[, 1],
         s_temp     = scale(annual_temp2m)[, 1],
         s_logpar   = scale(log(PAR))[, 1])

plot(s_c4_ddiff ~ log(totalmoist), data = c34sum, pch = 19, col = co2)
plot(s_c4_ddiff ~ annual_temp2m, data = c34sum, pch = 19, col = co2)
plot(s_c4_ddiff ~ log(PAR), data = c34sum, pch = 19, col = co2)




# c4 abundance ------------------------------------------------------------

c4d_m1      <- lmer(s_c4_ddiff ~ co2 * (s_logmoist+s_temp+s_logpar) + (1|ring), data = c34sum)
c4d_m1_full <- dredge(c4d_m1, REML = F, fixed = c("co2", "s_logmoist", "s_temp",  "s_logpar"))
plot(c4d_m1)
qqPlot(resid(c4d_m1))
c4d_m2      <- lmer(s_c4_ddiff ~ co2 + s_logmoist + s_temp + s_logpar + (1|ring), data = c34sum)
summary(c4d_m2)
Anova(c4d_m2, test.statistic = "F")
VarCorr(c4d_m2)
c4d_m3      <- lm(s_c4_ddiff ~ co2 + s_logmoist + s_temp + s_logpar, data = c34sum)
c4_coef    <- confint(c4d_m3)
c4_coef_90 <- confint(c4d_m3, level = .9)
c4d_m2_full <- dredge(c4d_m2, REML = F)
c4_coef_imp <- importance(c4d_m2_full)
c4_coef_imp

# partial residual plot

pdf(file = "output/figs/lar_C4_partial_regression_plot.pdf", width = 5, height = 5)
create_resplot(c4d_m2,
               ylab = expression(Adj.~annual~rates~of~change~"in"~C[4]),
               ylim = c(-1.5, 2.5))
dev.off()


png("output/figs/lar_C4_partial_regression_plot.png", width = 5, height = 5, res = 600, units = "in")
create_resplot(c4d_m2,
               ylab = expression(Adj.~annual~rates~of~change~"in"~C[4]),
               ylim = c(-1.5, 2.5))
dev.off()




# C3 abundance ------------------------------------------------------------


c3d_m1      <- lmer(s_c3_ddiff ~ co2 * (s_logmoist+s_temp+s_logpar) + (1|ring), data = c34sum)
c3d_m1_full <- dredge(c3d_m1, REML = F, fixed = c("co2", "s_logmoist", "s_temp",  "s_logpar"))
plot(c3d_m1)
qqPlot(resid(c3d_m1))
c3d_m2      <- lmer(s_c3_ddiff ~ co2 + s_logmoist + s_temp + s_logpar + (1|ring), data = c34sum)
summary(c3d_m2)
Anova(c3d_m2, test.statistic = "F")
VarCorr(c3d_m2)
plot(c3d_m2)
qqPlot(resid(c3d_m2))
mcp.fnc(c3d_m2)
# one outlier was suggested

# remove outlier
olrm_df <- romr.fnc(c3d_m2, data.frame(c34sum))
dplyr::setdiff(olrm_df$data0, olrm_df$data)
c3d_m3 <- update(c3d_m2, data = olrm_df$data)
plot(c3d_m3)
qqPlot(resid(c3d_m3))

VarCorr(c3d_m2)
VarCorr(c3d_m3)

confint(c3d_m2, method = "boot", nsim = 99, level = .9)
confint(c3d_m3, method = "boot", nsim = 99, level = .9)
# no major difference, so use the original one


c3_coef    <- confint(c3d_m2, method = "boot", nsim = 999)
c3_coef_90 <- confint(c3d_m2, method = "boot", nsim = 999, level = .9)
c3d_m2_full <- dredge(c3d_m2, REML = F)
c3_coef_imp <- importance(c3d_m2_full)
c3_coef_imp



# > partial residual plot ---------------------------------------------------

pdf(file = "output/figs/lar_C3_partial_regression_plot.pdf", width = 5, height = 5)
create_resplot(c3d_m2,
               ylab = expression(Adj.~annual~rates~of~change~"in"~C[3]),
               ylim = c(-2.5, 3))
dev.off()


png("output/figs/lar_C3_partial_regression_plot.pdf.png", width = 5, height = 5, res = 600, units = "in")
create_resplot(c3d_m2,
               ylab = expression(Adj.~annual~rates~of~change~"in"~C[3]),
               ylim = c(-2.5, 3))
dev.off()

gmai