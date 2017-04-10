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




# C4/3 proprtion change ---------------------------------------------

c34_prop <- PfgRDF$C3vsC4 %>% 
  arrange(id, year) %>% 
  group_by(id) %>% 
  mutate_each(funs(rat_diff = . - lag(., 1),             # Year1-Year0 and etc.
                   rat_prop = (. + 1) / lag(. + 1, 1)),  # Year1/Year0 and etc.
              ratios) 
    
# move year0 to a new column so that it can be used as a covariate
c34_prop_year0 <- c34_prop %>% 
  select(year, id, ratios) %>% 
  filter(year == "Year0") %>% 
  rename(ratios0 = ratios) %>%
  select(-year) %>% 
  right_join(c34_prop) %>% 
  left_join(c34growth_moist) %>% 
  filter(year != "Year0")



# > rat_diff  --------------------------------------------------------------
names(c34_prop_year0)
summary(c34_prop_year0)
ratd_m1 <- lmer(rat_diff ~ co2 * (swa_c4     + annual_temp2m) + (1|ring) + (1|RY) + (1|id), data = c34_prop_year0)
ratd_m2 <- lmer(rat_diff ~ co2 * (c4moist    + annual_temp2m) + (1|ring) + (1|RY) + (1|id), data = c34_prop_year0)
ratd_m3 <- lmer(rat_diff ~ co2 * (totalmoist + annual_temp2m) + (1|ring) + (1|RY) + (1|id), data = c34_prop_year0)
ratd_init <- list(ratd_m1, ratd_m2, ratd_m3)
model.sel(llply(ratd_init, function(x) update(x, REML = F)))
plot(ratd_m2)
qqnorm(resid(ratd_m2))
qqline(resid(ratd_m2))

Anova(ratd_m2, test.statistic = "F")
ratd_m2_full <- dredge(ratd_m2, REML = F)
ratd_m2_avg <- model.avg(get.models(ratd_m2_full, subset = cumsum(weight) <= .95))
summary(ratd_m2_avg)
confint(ratd_m2_avg)
# delta AICc for null model is only 2.39




# > rat_prop  --------------------------------------------------------------
names(c34_prop_year0)
summary(c34_prop_year0)
ratp_m1 <- lmer(rat_prop ~ co2 * (swa_c4     + annual_temp2m) + (1|ring) + (1|RY) + (1|id), data = c34_prop_year0)
ratp_m2 <- lmer(rat_prop ~ co2 * (c4moist    + annual_temp2m) + (1|ring) + (1|RY) + (1|id), data = c34_prop_year0)
ratp_m3 <- lmer(rat_prop ~ co2 * (totalmoist + annual_temp2m) + (1|ring) + (1|RY) + (1|id), data = c34_prop_year0)
ratp_init <- list(ratp_m1, ratp_m2, ratp_m3)
model.sel(llply(ratp_init, function(x) update(x, REML = F)))
plot(ratp_m2)
qqnorm(resid(ratp_m2))
qqline(resid(ratp_m2))

Anova(ratp_m2, test.statistic = "F")
ratp_m2_full <- dredge(ratp_m2, REML = F)
# delta AICc for null model is only 1.69




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




# . random slopes ---------------------------------------------------------

# moisture
xyplot(s_c4_ddiff ~ s_logmoist | ring, group = id, data = c34sum, type=c("p", "r"))
  # this should be included
xyplot(s_c4_ddiff ~ s_logmoist, group = ring, data = c34sum, type=c("p", "r"))
  # this is probably required
xyplot(c4_ddiff ~ s_logmoist | ring, group = RY, data = c34sum, type=c("p", "r"))
  # this should be included 
xyplot(c4_ddiff ~ s_logmoist, group = year, data = c34sum, type=c("p", "r"))


# temperature
xyplot(s_c4_ddiff ~ s_temp | ring, group = id, data = c34sum, type=c("p", "r"))
xyplot(s_c4_ddiff ~ s_temp, group = ring, data = c34sum, type=c("p", "r"))
xyplot(s_c4_ddiff ~ s_temp | year, group = RY, data = c34sum, type=c("p", "r")) # not run
xyplot(s_c4_ddiff ~ s_temp, group = year, data = c34sum, type=c("p", "r"))      # not run

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
confint(c4d_m0_avg, level = .9, full = TRUE)
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
confint(c4d_m2, method = "boot", level = .9)


exp(-m4coef[2]/m4coef[3]) # miosture required for delta C4 to be positive in eCO2 relative to ambient (ignoring temp and par as their coeeficients are close to 0)

write.csv(c4d_m0_full, file = "output/table/delta_c4_modelsel.csv", na = "-")




# . predicted values ------------------------------------------------------


# > Use each year climate conditions --------------------------------------
names(c34sum)

# yearly climate
yearly_climate <- c34sum[-40, ] %>% 
  group_by(year) %>% 
  summarise_each(funs(mean), s_logmoist, s_temp, s_logpar, totalmoist, annual_temp2m, PAR)
yearly_climate_df <- rbind(cbind(yearly_climate, co2 = "amb"), cbind(yearly_climate, co2 = "elev"))


# number of C4 spp per plot
c34_spn <- C3grassC4 %>% 
  group_by(year, co2, PFG, id, variable) %>% 
  summarise(value = sum(value)) %>% 
  group_by(year, co2, PFG, id) %>% 
  summarise(sp_n = sum(value > 0)) %>% 
  group_by(year, co2, PFG) %>% 
  summarise(spn_plt = round(mean(sp_n), 0), sample_n = sum(!is.na(sp_n))) %>% 
  ungroup() %>% 
  select(-sample_n) %>% 
  spread(PFG, spn_plt) %>% 
  filter(year != "Year0") %>% 
  rename(c3_spn = c3, c4_spn = c4)




yearly_pred <- predict(c4d_m0_bs_lm, yearly_climate_df, interval = "confidence", level = .9)
yearly_pred_df <- data.frame(yearly_climate_df, yearly_pred) %>% 
  mutate(r_moist = exp(rev_ztrans(s_logmoist, xsd = moist_sd, xmean = moist_m)),
         r_temp  = rev_ztrans(s_temp, xsd = temp_sd, xmean = temp_m),
         r_par   = exp(rev_ztrans(s_logpar, xsd = par_sd, xmean = par_m)),
         r_lwr   = exp(rev_ztrans(lwr, c4d_sd, c4d_m)), 
         r_fit   = exp(rev_ztrans(fit, c4d_sd, c4d_m)),
         r_upr   = exp(rev_ztrans(upr, c4d_sd, c4d_m))) %>% 
  left_join(c34_spn) %>% 
  mutate(lwr_persp = r_lwr^(1/4), 
         upr_persp = r_upr^(1/4), 
         fit_persp = r_fit^(1/4)) %>% 
  arrange(year, co2)

ggplot(yearly_pred_df, aes(x = year, y = log(r_fit)))+
  geom_point(aes(col = co2), position = position_dodge(.5))+
  geom_point(data = c34sum, aes(y = c4_ddiff, shape = co2), col = "gray", 
             position = position_dodge(.5))+
  geom_errorbar(aes(ymin = log(r_lwr), ymax = log(r_upr), col = co2), width = .1, 
                position = position_dodge(.5))+
  geom_hline(yintercept = 0)

# > regression lines ------------------------------------------------------


par(mfrow = c(2, 2))
boxplot(totalmoist ~ year, c34sum, ylab = "Moist")
boxplot(annual_temp2m ~ year, c34sum, ylab = "Temp")
boxplot(PAR ~ year, c34sum, ylab = "PAR")





# .figure  -----------------------------------------------------------


newdf_moist <- data.frame(co2 = c("amb", "elev"), 
                    s_temp     = median(c34sum$s_temp),
                    s_logpar   = median(c34sum$s_logpar),
                    s_logmoist = seq(min(c34sum$s_logmoist), max(c34sum$s_logmoist), length.out = 1000))

newdf_temp <- data.frame(co2 = c("amb", "elev"), 
                    s_logmoist     = median(c34sum$s_logmoist),
                    s_logpar   = median(c34sum$s_logpar),
                    s_temp = seq(min(c34sum$s_temp), max(c34sum$s_temp), length.out = 1000))

newdf_par<- data.frame(co2 = c("amb", "elev"), 
                    s_temp     = median(c34sum$s_temp),
                    s_logmoist = median(c34sum$s_logmoist),
                    s_logpar   = seq(min(c34sum$s_logpar), max(c34sum$s_logpar), length.out = 1000))

newdf_l <- list(moist = newdf_moist, temp = newdf_temp, par = newdf_par)


# remove uninformative random factors to avoid convergence problem
c4d_m3 <- lmer(s_c4_ddiff ~ co2 + s_logmoist+s_temp+s_logpar+(1|RY), data = c34sum[-40, ])


# perform bootstrap and 90% CI

pred_df_l <- llply(names(newdf_l), function(x){
  print(x)
  rm(envd)
  envd <<- newdf_l[[x]]
  pred_bb <- bootMer(c4d_m3,
                     FUN = function(y) predict(y, envd, re.form = NA),
                     nsim = 999)
  pred_df <- cbind(get_ci(pred_bb, a = .1), envd)
  return(pred_df)
})
names(pred_df_l) <- names(newdf_l)


# reverse transform
moist_sd <- sd(log(c34sum$totalmoist))
moist_m  <- mean(log(c34sum$totalmoist))
temp_sd  <- sd(c34sum$annual_temp2m)
temp_m   <- mean(c34sum$annual_temp2m)
par_sd   <- sd(log(c34sum$PAR))
par_m    <- mean(log(c34sum$PAR))
c4d_sd   <- sd(c34sum$c4_ddiff)
c4d_m    <- mean(c34sum$c4_ddiff)



c4d_m2_preddf_rv_l <- llply(pred_df_l, function(x){
  x %>% 
    mutate(r_moist = rev_ztrans(s_logmoist, xsd = moist_sd, xmean = moist_m),
           r_temp  = rev_ztrans(s_temp    , xsd = temp_sd , xmean = temp_m),
           r_par   = rev_ztrans(s_logpar  , xsd = par_sd  , xmean = par_m),
           r_fit   = rev_ztrans(fit, c4d_sd, c4d_m),
           r_lwr   = rev_ztrans(lwr, c4d_sd, c4d_m), 
           r_upr   = rev_ztrans(upr, c4d_sd, c4d_m))
})
names(c4d_m2_preddf_rv_l)
c4d_m2_preddf_rv_l[[1]]$xval <- c4d_m2_preddf_rv_l[[1]]$r_moist
c4d_m2_preddf_rv_l[[2]]$xval <- c4d_m2_preddf_rv_l[[2]]$r_temp
c4d_m2_preddf_rv_l[[3]]$xval <- c4d_m2_preddf_rv_l[[3]]$r_par



llply(c4d_m2_preddf_rv_l, summary)
c4d_p_l <- llply(c4d_m2_preddf_rv_l, function(x){
  ggplot(x, aes(x = xval, y = r_fit, col = co2)) +
    geom_hline(yintercept = 0, col = "gray")+
    geom_line(size = 1)+
    geom_line(aes(y = r_lwr), linetype = "dashed") +
    geom_line(aes(y = r_upr), linetype = "dashed") +
    science_theme+
    theme(legend.position = "none")+
    scale_color_manual(values = c("blue", "red"),labels = c("Ambient", expression(eCO[2])))+
    labs(y = expression(Delta*C[4]))
})



c4d_p_l[[1]] <- c4d_p_l[[1]] + 
  geom_point(data = c34sum, aes(x = log(totalmoist), y = c4_ddiff), size = 2, alpha = .7)+
  labs(x = "ln(Soil moisture)") +
  theme(legend.position = c(.25, .85))
c4d_p_l[[2]] <- c4d_p_l[[2]] + 
  geom_point(data = c34sum, aes(x = annual_temp2m, y = c4_ddiff), size = 2, alpha = .7)+
  labs(x = expression(Temperature~(degree*C)))
c4d_p_l[[3]] <- c4d_p_l[[3]] + 
  geom_point(data = c34sum, aes(x = log(PAR), y = c4_ddiff), size = 2, alpha = .7)+
  labs(x = expression(ln(Understorey~PAR~(mu*mol~s^'-1'~m^"-2"))))


sublabels <- c("(a)", "(b)", "(c)")

for(i in 1:3){
  c4d_p_l[[i]] <- c4d_p_l[[i]] +
    annotate("text", label = sublabels[i], x = -Inf, y = Inf, hjust = -.5, vjust = 1.5, 
              fontface = "bold")
  }

c4d_p_l[[1]]
c4d_p_l[[3]]


# >partial residual plot ------------------------------------------------

deltac4_regplt <- function(){

  par(mfrow = c(2, 2), mar = c(4.5, 4.5, .5, .5))
  visreg(c4d_m2, xvar = "s_logpar", ylab = expression(Adj.~LAR[C4]), 
         xlab = expression(Adj.~ln(PAR,~mu*mol~s^'-1'~m^"-2")),
         alpha = .1)
  
  visreg(c4d_m2, xvar = "s_logmoist", ylab = expression(Adj.~LAR[C4]), 
         xlab = "Adj. ln(Soil moisture)",
         alpha = .1)
  
  visreg(c4d_m2, xvar = "s_temp", ylab = expression(Adj.~LAR[C4]), 
         xlab = expression(Adj.~Temperature~(degree*C)),
         alpha = .1)
  
  visreg(c4d_m2, xvar = "co2", ylab = expression(Adj.~LAR[C4]), 
         xlab = expression(CO[2]),
         alpha = .1)
  
  
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
summary(c3d_m2)
confint(c3d_m2, method = "boot", level = .9)



# > partial residual plot ---------------------------------------------------

deltac3_regplt <- function(){
  
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, .5, .5))
  visreg(c3d_m2, xvar = "s_logpar", ylab = expression(Adj.~LAR[C3]), 
         xlab = expression(Adj.~ln(PAR,~mu*mol~s^'-1'~m^"-2")),
         alpha = .1)
  
  visreg(c3d_m2, xvar = "s_logmoist", ylab = expression(Adj.~LAR[C3]), 
         xlab = "Adj. ln(Soil moisture)",
         alpha = .1)
  
  visreg(c3d_m2, xvar = "s_temp", ylab = expression(Adj.~LAR[C3]), 
         xlab = expression(Adj.~Temperature~(degree*C)),
         alpha = .1)
  
  visreg(c3d_m2, xvar = "co2", ylab = expression(Adj.~LAR[C3]), 
         xlab = expression(CO[2]),
         alpha = .1)
  
  
}

pdf(file = "output/figs/deltaC3_partial_regression_plot.pdf", width = 5, height = 5)
deltac3_regplt()
dev.off()


png("output/figs/deltaC3_partial_regression_plot.png", width = 5, height = 5, res = 600, units = "in")
deltac3_regplt()
dev.off()




# predicted values --------------------------------------------------------


# get 95% CI by parametric bootstrap
par_df <- data.frame(s_logpar = with(c3d_m2, seq(min(s_logpar), 
                                                  max(s_logpar), length.out = 1000)))
bb <- bootMer(c3d_m0bs,
              FUN=function(x) predict(x, par_df, re.form = NA),
              nsim=999)
lwr <- apply(bb$t, 2, quantile, 0.025)
upr <- apply(bb$t, 2, quantile, 0.975)
fit <- bb$t0
c3_pred_df <- cbind(lwr, upr, fit, par_df)


# reverse transform
par_sd   <- sd(log(c34sum$PAR))
par_m    <- mean(log(c34sum$PAR))
c3d_sd   <- sd(c34sum$c3_ddiff)
c3d_m    <- mean(c34sum$c3_ddiff)

c3d_preddf_rv <- c3_pred_df %>% 
  mutate(r_par = rev_ztrans(s_logpar, xsd = par_sd, xmean = par_m), 
         r_fit   = rev_ztrans(fit, c3d_sd, c3d_m),
         r_lwr   = rev_ztrans(lwr, c3d_sd, c3d_m), 
         r_upr   = rev_ztrans(upr, c3d_sd, c3d_m))


# create a plot
c3d_p <- ggplot(c3d_preddf_rv, aes(x = r_par, y = r_fit))+
  geom_hline(yintercept = 0, col = "gray")+
  geom_point(data = c34sum2, aes(x = log(PAR), y = s_c3_ddiff, col = co2), size = 2, alpha = .7)+
  geom_line(size = 1)+
  geom_line(aes(y = r_lwr), linetype = "dashed")+
  geom_line(aes(y = r_upr), linetype = "dashed")+
  science_theme+
  theme(legend.position = "none")+
  scale_color_manual(values = c("blue", "red"),labels = c("Ambient", expression(eCO[2])))+
  labs(x = expression(ln(Understorey~PAR*','~mu*mol~s^'-1'~m^"-2")), 
       y = expression(Delta*C[3]))+
  annotate("text", label = "(b)", x = -Inf, y = Inf, hjust = -.5, vjust = 1.5, 
           fontface = "bold")




# merge figures -----------------------------------------------------------

deltaC34_fig <- cbind(ggplotGrob(c4d_p_l[[1]]), ggplotGrob(c3d_p))
grid.newpage()
grid.draw(deltaC34_fig)

ggsavePP(filename = "output/figs/delta_C34_envvar", plot = deltaC34_fig,
         width = 6.5, height = 3)

