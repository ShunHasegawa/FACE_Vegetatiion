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

c34growth_moist <- left_join(veg_moist, c34growthdate) %>% 
  filter(year %in% 2013:2015) %>% 
  mutate(year     = factor(year, labels = paste0("Year", 1:3)),
         plot     = factor(plot)) %>% 
  group_by(year, ring, plot) %>%
  summarise_each(funs(c3moist = sum(.[c3growth]),
                      c4moist = sum(.[c4growth]),
                      totalmoist = mean), 
                 Moist) %>% 
  ungroup() %>% 
  mutate(swa_c3 = c3moist / (c3moist + c4moist),  # seasonal water availability
         swa_c4 = c4moist / (c3moist + c4moist),
         plot = factor(plot)) %>% 
  left_join(annualtemp) %>% 
  left_join(underpar)




# abundance change ------------------------------------------------------------


c34sum <- C3grassC4 %>%
  group_by(year, block, ring, plot, co2, id, PFG, RY, variable) %>% 
  summarise(value = sum(value)) %>%
  # filter(variable != "Cynodon.dactylon") %>%
  mutate(value = log(value + 1)) %>% 
  group_by(year, block, ring, plot, co2, id, PFG, RY) %>% 
  summarise(value = sum(value)) %>% 
  spread(key = PFG, value = value) %>% 
  ungroup() %>% 
  arrange(id, year) %>%
  group_by(id) %>%
  mutate_each(funs(ddiff = . - lag(., 1),             # Year1-Year0 and etc.
                   dprop = (. + 1) / lag(. + 1, 1)),  # Year1/Year0 and etc.
                   c3, c4) %>%
  filter(year != "Year0") %>%
  left_join(c34growth_moist) %>% 
  ungroup() %>% 
  mutate(s_c4_ddiff = scale(c4_ddiff)[, 1],
         s_c3_ddiff = scale(c3_ddiff)[, 1],
         s_logmoist = scale(log(totalmoist))[, 1],
         s_c3       = scale(log(c3))[, 1],
         s_temp     = scale(annual_temp2m)[, 1],
         s_logpar   = scale(log(PAR))[, 1])

plot(s_c4_ddiff ~ log(totalmoist), data = c34sum, pch = 19, col = co2)
plot(s_c4_ddiff ~ annual_temp2m, data = c34sum, pch = 19, col = co2)
plot(s_c4_ddiff ~ log(PAR), data = c34sum, pch = 19, col = co2)

which.max(c34sum$c4_dprop)
plot(log(c4_dprop) ~ log(totalmoist), data = c34sum, pch = 19, col = co2, subset = -5)
plot(log(c4_dprop) ~ annual_temp2m, data = c34sum, pch = 19, col = co2, subset = -5)
plot(log(c4_dprop) ~ log(PAR), data = c34sum, pch = 19, col = co2, subset = -5)



# c4 abundance ------------------------------------------------------------



# > diff ------------------------------------------------------------------
names(c34sum)

c4d_m0 <- lmer(s_c4_ddiff ~ co2*(s_logmoist+s_temp+s_logpar)+(1|ring)+(1|RY)+(1|id), data = c34sum)
summary(c4d_m0)

tt <- getME(c4d_m0,"theta")
ll <- getME(c4d_m0,"lower")
min(tt[ll==0])
## it has a singularity problem

plot(c4d_m0)
qqnorm(resid(c4d_m0))
qqline(resid(c4d_m0))

c4d_m0_full <- dredge(c4d_m0, REML = F, extra = "r.squaredGLMM")
c4nest <- subset(c4d_m0_full, !nested(.))
c4d_m0_avg  <- model.avg(get.models(c4nest, subset = delta <= 2))
confint(c4d_m0_avg, full = TRUE)
c4d_m0_bs   <- get.models(c4d_m0_full, subset = 1)[[1]]
summary(c4d_m0_avg)

# remove potential outlier
plot(c4d_m0)
qqnorm(resid(c4d_m0))
qqline(resid(c4d_m0))
which.min(resid(c4d_m0))
boxplot(c34sum$s_c4_ddiff)
points(c34sum$s_c4_ddiff[40], col = "red", pch = 19)
c34sum[40, ]

c4d_m1      <- lmer(s_c4_ddiff ~ co2*(s_logmoist+s_temp+s_logpar)+(1|ring)+(1|RY)+(1|id), data = c34sum[-40, ])
c4d_m1_full <- dredge(c4d_m1, REML = F, extra = "r.squaredGLMM")
c4_m1_nest  <- subset(c4d_m1_full, !nested(.))
    # interaction was driven by the outlier
    # now there is no no indication of interaction, so refit the model only with main effects
c4d_m2      <- lmer(s_c4_ddiff ~ co2 + s_logmoist+s_temp+s_logpar+(1|ring)+(1|RY)+(1|id), data = c34sum[-40, ])
c4d_m3      <- lmer(s_c4_ddiff ~ co2 + s_logmoist+s_temp+s_logpar+(1|RY), data = c34sum[-40, ])
summary(c4d_m2)
c4d_m2_full <- dredge(c4d_m2, REML = F, extra = "r.squaredGLMM")
c4_coef <- confint(c4d_m2, method = "boot", nsim = 999)
c4_coef_90 <- confint(c4d_m2, method = "boot", level = .9, nsim = 999)
c4_coef_imp <- importance(c4d_m2_full)

# miosture required for delta C4 to be positive in eCO2 relative to ambient (ignoring temp and par as their coeeficients are close to 0)
c4_estimate <- summary(c4d_m2)$coeff[, "Estimate"]
exp(-c4_estimate[2] * sd(log(c34sum$totalmoist))/c4_estimate[3]) # See Docs/LAR_C4_and_water.html for calculation


write.csv(c4d_m0_full, file = "output/table/delta_c4_modelsel.csv", na = "-")





# . level plot ------------------------------------------------------------


# predict values
c4_envdf <- with(c34sum, expand.grid(s_logmoist = seq(-2.5, 2.5, length.out = 200),
                                     s_logpar   = seq(-2.5, 2.5, length.out = 200),
                                     s_temp     = median(s_temp),
                                     co2        = c("amb", "elev"), 
                                     ring = 1, id = "1:1", RY = "1:Year1"))


c4_pred <- predict(c4d_m2, c4_envdf, re.form = NA)


# reverse transform
moist_sd <- sd(log(c34sum$totalmoist))
moist_m  <- mean(log(c34sum$totalmoist))
temp_sd  <- sd(c34sum$annual_temp2m)
temp_m   <- mean(c34sum$annual_temp2m)
par_sd   <- sd(log(c34sum$PAR))
par_m    <- mean(log(c34sum$PAR))
c4d_sd   <- sd(c34sum$c4_ddiff)
c4d_m    <- mean(c34sum$c4_ddiff)


c4_pred_df <- cbind(fit = c4_pred, c4_envdf) %>% 
  mutate(r_moist = rev_ztrans(s_logmoist, xsd = moist_sd, xmean = moist_m),
         r_temp  = rev_ztrans(s_temp, xsd = temp_sd, xmean = temp_m),
         r_par   = rev_ztrans(s_logpar, xsd = par_sd, xmean = par_m),
         r_fit   = rev_ztrans(fit, c4d_sd, c4d_m),
         co2     = mapvalues(co2, c("amb", "elev"), c("Ambient", "eCO[2]")))

c34sum_temp <- c34sum %>% 
  mutate(co2     = mapvalues(co2, c("amb", "elev"), c("Ambient", "eCO[2]")))

c4_levelplot <- ggplot(c4_pred_df, aes(x = r_moist, y = r_par)) + 
  geom_tile(aes(fill = r_fit)) +
  scale_fill_gradient2("Log annual\nchange \nrates of C4", low = "blue",high = "red", mid = "white")+
  stat_ellipse(data = c34sum_temp, aes(x = log(totalmoist), y = log(PAR), linetype = year), 
               type = "norm", level = .7)+
  geom_point(data = c34sum_temp, aes(x = log(totalmoist), y = log(PAR), shape = year), size = 2)+
  scale_shape_manual("Year", values = c(0:2), label = 2013:2015)+
  scale_linetype_manual("Year", values = c(1:3), label = 2013:2015)+
  facet_grid(. ~ co2, labeller = label_parsed)+
  labs(x = expression(Log[e](Moist)), y = expression(Log[e](PAR,~mu*mol~s^'-1'~m^"-2")))+
  theme(legend.title = element_text(size = 8))
ggsavePP(filename = "output/figs/LARC4_levelplot_byMoistPAR", plot = c4_levelplot, 
         width = 6.5, height = 3)




# >partial residual plot ------------------------------------------------

deltac4_regplt <- function(){
  
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, .5, .5))
  visreg(c4d_m2, xvar = "s_logpar", 
         ylab = expression(Adj.~annual~change~rates~of~C[4]), 
         xlab = expression(Adj.~Log[e](PAR,~mu*mol~s^'-1'~m^"-2")))
  
  visreg(c4d_m2, xvar = "s_logmoist", 
         expression(Adj.~annual~change~rates~of~C[4]), 
         xlab = expression(Adj.~Log[e](Moist)))
  
  visreg(c4d_m2, xvar = "s_temp", 
         expression(Adj.~annual~change~rates~of~C[4]), 
         xlab = expression(Adj.~Temp~(degree*C)))
  
  visreg(c4d_m2, xvar = "co2", ylab = expression(Adj.~LAR[C4]), 
         xlab = expression(CO[2]))
  
  
}

pdf(file = "output/figs/deltaC4_partial_regression_plot.pdf", width = 5, height = 5)
deltac4_regplt()
dev.off()


png("output/figs/deltaC4_partial_regression_plot.png", width = 5, height = 5, res = 600, units = "in")
deltac4_regplt()
dev.off()




# C3 abundance ------------------------------------------------------------



# > diff ------------------------------------------------------------------

# moisture
xyplot(s_c3_ddiff ~ s_logmoist | ring, group = id, data = c34sum, type=c("p", "r"))
xyplot(s_c3_ddiff ~ s_logmoist, group = ring, data = c34sum, type=c("p", "r"))
xyplot(s_c3_ddiff ~ s_logmoist, group = RY, data = c34sum, type=c("p", "r"))
xyplot(s_c3_ddiff ~ s_logmoist, group = year, data = c34sum, type=c("p", "r"))

# temperature
xyplot(s_c3_ddiff ~ s_temp | ring, group = id, data = c34sum, type=c("p", "r"))
xyplot(s_c3_ddiff ~ s_temp, group = ring, data = c34sum, type=c("p", "r"))


# remove outlier
boxplot(c34sum$s_c3_ddiff)
c34sum2 <- c34sum %>% 
  filter(s_c3_ddiff != max(s_c3_ddiff))

c3d_m0     <- lmer(s_c3_ddiff ~ co2 * (s_logmoist+s_temp+s_logpar) + (1|ring) + (1|RY) +(1|id), data = c34sum2)
plot(c3d_m0)
qqnorm(resid(c3d_m0))
qqline(resid(c3d_m0))

c3d_m0full <- dredge(c3d_m0, REML = F, extra = "r.squaredGLMM")
confint(model.avg(get.models(c3d_m0full, subset = delta <= 4)))
c3nest   <- subset(c3d_m0full, !nested(.))
    # no interaction is indicated


# remove potential outlier
which.min(resid(c3d_m0))
c3d_m1     <- lmer(s_c3_ddiff ~ co2 * (s_logmoist+s_temp+s_logpar) + (1|ring) + (1|RY) +(1|id), data = c34sum2[-48, ])
plot(c3d_m1)
qqnorm(resid(c3d_m1))
qqline(resid(c3d_m1))
c3d_m1full <- dredge(c3d_m1, REML = F, extra = "r.squaredGLMM")
c3nest1   <- subset(c3d_m1full, !nested(.))
c3nest1
  # similar to the above


# coefficient
c3d_m2     <- lmer(s_c3_ddiff ~ co2 + s_logmoist+s_temp+s_logpar + (1|ring) + (1|RY) +(1|id), data = c34sum2)
c3d_m2_full <- dredge(c3d_m2, REML = F)

summary(c3d_m2)
c3_coef <- confint(c3d_m2, method = "boot", nsim = 999)
c3_coef_90 <- confint(c3d_m2, method = "boot", level = .9, nsim = 999)
c3_coef_impo <- importance(c3d_m2_full)


# > partial residual plot ---------------------------------------------------

deltac3_regplt <- function(){
  
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, .5, .5))
  visreg(c3d_m2, xvar = "s_logpar", 
         expression(Adj.~annual~change~rates~of~C[3]),
         xlab = expression(Adj.~Log[e](PAR,~mu*mol~s^'-1'~m^"-2")))
  
  visreg(c3d_m2, xvar = "s_logmoist", 
         expression(Adj.~annual~change~rates~of~C[3]), 
         xlab = expression(Adj.~Log[e](Moist)))
  
  visreg(c3d_m2, xvar = "s_temp", 
         expression(Adj.~annual~change~rates~of~C[3]), 
         xlab = expression(Adj.~Temp~(degree*C)))
  
  visreg(c3d_m2, xvar = "co2", ylab = expression(Adj.~LAR[C3]), 
         xlab = expression(CO[2]))
  
  
}

pdf(file = "output/figs/deltaC3_partial_regression_plot.pdf", width = 5, height = 5)
deltac3_regplt()
dev.off()


png("output/figs/deltaC3_partial_regression_plot.png", width = 5, height = 5, res = 600, units = "in")
deltac3_regplt()
dev.off()


write.csv(c3d_m0full, file = "output/table/delta_c3_modelsel.csv", na = "-")
