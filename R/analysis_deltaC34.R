# Here we will:
# 1) download moisture and understorey air temperature data from HIEv
# 2) identify growing season for each of C3 and C4 species according to Murphy 2007
# 3) Fit moisture and temperature to delta C3 and C4 (annual change in abundance)



# download from HIEv ------------------------------------------------------


# setToken(tokenfile = "Data/token.txt")
# 
# airvar_raw <- downloadTOA5(filename = "FACE_.*_AirVars_.*dat",
#                           topath    = "Data/hievdata/raw_data/",
#                           maxnfiles = 999)
# 
# save(airvar_raw, file = "output/Data/FACE_airvar_raw.RData")
load("output/Data/FACE_airvar_raw.RData")




# process HIEv data -------------------------------------------------------


# get daily mean, min and max temparater at two layers (at a height of 2m and 15m)
boxplot(airvar_raw$AirTC_1_Avg)
airvar_day <- airvar_raw %>%
  filter(Date        > as.Date("2012-09-18"),         # when co2 treatment was commenced
         Date        < as.Date("2016-2-15"),
         AirTC_1_Avg > -10) %>%                       # remove obviously weird value
  distinct() %>%                                      # remove duplicates
  mutate(ring = factor(substring(Source, 7, 7)),      # add ring number
         year = factor(year(Date))) %>%
  group_by(year, Date, ring) %>% 
  summarise_each(funs(Mean = mean(., na.rm = TRUE), 
                      Min  = min(., na.rm = TRUE),
                      Max  = max(., na.rm = TRUE),
                      N    = sum(!is.na(.))), 
                 airtemp2m = AirTC_1_Avg, airtemp15m = AirTC_2_Avg) %>% 
  group_by(year) %>% 
  mutate(annual_temp2m = mean(airtemp2m_Mean),
         T = max(27, 1.745 * annual_temp2m  + 11.143)[1]) %>% 
  ungroup() %>% 
  mutate(c3growth = airtemp2m_Min >= -1 & airtemp2m_Max >= 10 & airtemp2m_Max < 24, # c3 growing dates
         c4growth = airtemp2m_Max >= 21 & airtemp2m_Max < T & month(Date))                        # c4 growtin dates

head(airvar_day$c4growth)

airvar_day$Date[airvar_day$c4growth]


c34growthdate <- airvar_day %>% 
  select(year, Date, ring, c3growth, c4growth, annual_temp2m)
annualtemp <- c34growthdate %>% 
  select(year, annual_temp2m) %>% 
  distinct() %>% 
  filter(!year %in% c("2012", "2016")) %>% 
  mutate(year = mapvalues(year, c(2013, 2014, 2015), paste0("Year", 1:3)))

## Murphy 2007. For C 3 grasses, growth was considered possible in any month 
## where the daily minimum temperature was ≥ −1 oC, and the daily maximum 
## temperature was ≥ 10 oC and < 24 oC. For C 4 grasses, growth was considered 
## possible in months where the daily maximum temperature was ≥ 21 °C and less 
## than the tem- perature define (Murphy and Bowman 2007).




# assign moisture data to each vegetation plot based on their corr --------


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


## survey dates: 2012-12-15 (Year0), 2014-1-15 (Year1), 2015-1-30 (Year2), 2016 (Year3)

c34growth_moist <- left_join(veg_moist, c34growthdate) %>% 
  filter(Date >= as.Date("2012-12-15"), Date <= as.Date("2016-2-15")) %>% 
  mutate(year     = ifelse(Date < as.Date("2014-1-15"), "Year1", 
                           ifelse(Date < as.Date("2015-1-30"), "Year2", "Year3")),
         plot     = factor(plot)
         # ,
         # c4growth = ifelse(c4growth & year == "Year1" & Date <= as.Date("2013-12-31"), TRUE,
         #                   ifelse(c4growth & year == "Year2" & Date <= as.Date("2014-12-31"), TRUE, 
         #                          ifelse(c4growth & year == "Year3" & Date <= as.Date("2015-12-31"), TRUE, FALSE)))
         ) %>% 
  group_by(year, ring, plot) %>%
  summarise_each(funs(c3moist = sum(.[c3growth]),
                      c4moist = sum(.[c4growth]),
                      totalmoist = sum), 
                 Moist) %>% 
  ungroup() %>% 
  mutate(swa_c3 = c3moist / (c3moist + c4moist),  # seasonal water availability
         swa_c4 = c4moist / (c3moist + c4moist),
         plot = factor(plot)) %>% 
  left_join(annualtemp)

plot(c4moist ~ totalmoist, data = c34growth_moist)


plot(swa_c4 ~ totalmoist, data = c34growth_moist)
plot(swa_c3 ~ totalmoist, data = c34growth_moist)
c34growth_moist %>% 
  gather(key = category, value = moisture, c3moist, c4moist, totalmoist) %>%
  group_by(category, year, ring) %>% 
  summarise(moisture = mean(moisture)) %>% 
ggplot(., aes(x = year, y = moisture, col = category))+
  geom_path(aes(group =  category)) +
  facet_wrap(~ ring)


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
  left_join(c34growth_moist)





# c4 abundance ------------------------------------------------------------



# > diff ------------------------------------------------------------------
names(c34sum)

c34sum <- c34sum %>% 
  ungroup() %>% 
  mutate(s_c4_ddiff = scale(c4_ddiff)[, 1],
         s_c3_ddiff = scale(c3_ddiff)[, 1],
         s_logmiost = scale(log(totalmoist))[, 1],
         s_c3       = scale(log(c3))[, 1],
         s_temp     = scale(annual_temp2m)[, 1])

plot(c4_ddiff ~ sqrt(c4moist), c34sum, col = co2, pch = 19)
plot(c4_ddiff ~ log(c4moist), c34sum, col = co2, pch = 19)
plot(c4_ddiff ~ log(totalmoist), c34sum, col = co2, pch = 19)
plot(c4_ddiff ~ c3, c34sum, col = co2, pch = 19)



# . random slopes ---------------------------------------------------------

# moisture
xyplot(s_c4_ddiff ~ s_logmiost | ring, group = id, data = c34sum, type=c("p", "r"))
  # this should be included
xyplot(s_c4_ddiff ~ s_logmiost, group = ring, data = c34sum, type=c("p", "r"))
  # this is probably required
xyplot(c4_ddiff ~ s_logmiost, group = RY, data = c34sum, type=c("p", "r"))
  # this should be included 
xyplot(c4_ddiff ~ s_logmiost, group = year, data = c34sum, type=c("p", "r"))


# temperature
xyplot(s_c4_ddiff ~ s_temp | ring, group = id, data = c34sum, type=c("p", "r"))
xyplot(s_c4_ddiff ~ s_temp, group = ring, data = c34sum, type=c("p", "r"))
xyplot(s_c4_ddiff ~ s_temp | year, group = RY, data = c34sum, type=c("p", "r")) # not run
xyplot(s_c4_ddiff ~ s_temp, group = year, data = c34sum, type=c("p", "r"))      # not run


names(c34sum)
c4d_m0 <- lmer(s_c4_ddiff ~ co2*(s_logmiost+s_temp)+(1|ring)+(1|RY)+(1|id), data = c34sum)
c4d_m1 <- lmer(s_c4_ddiff ~ co2*(s_logmiost+s_temp)+(1|RY), data = c34sum)
c4d_m1 <- lmer(s_c4_ddiff ~ co2+s_logmiost+(1|RY), data = c34sum)
# c4d_m1 <- lmer(s_c4_ddiff ~ co2*(s_logmiost+s_temp)+(1+s_logmiost|ring)+(1+s_logmiost|year)+(1|RY)+(1|id), data = c34sum)
anova(c4d_m0, c4d_m1)
confint(c4d_m1, method = "boot")

c4d_m0 <- lmer(s_c4_ddiff ~ co2*(s_logmiost+s_temp)+(1|ring)+(1|RY)+(1|id), data = c34sum)
c4d_m0 <- lmer(s_c4_ddiff ~ co2+s_logmiost+(1|ring)+(1|RY)+(1|id), data = c34sum)

c4d_m0 <- lmer(s_c4_ddiff ~ 1+(1+s_logmiost|ring)+(1+s_logmiost|year), data = c34sum)
# c4d_m1 <- lmer(s_c4_ddiff ~ 1+(1|ring)+(1+s_logmiost|year), data = c34sum)
c4d_m2 <- lmer(s_c4_ddiff ~ 1+(1+s_logmiost|ring)+(1|year), data = c34sum)
c4d_m3 <- lmer(s_c4_ddiff ~ 1+(1|ring)+(1|year), data = c34sum)

c4d_m2 <- lmer(s_c4_ddiff ~ co2*s_logmiost+(1|ring/id), data = c34sum)
c4d_m3 <- lmer(s_c4_ddiff ~ co2*s_logmiost+(1+s_logmiost|ring/id), data = c34sum)
anova(c4d_m2, c4d_m3)


c4d_m4 <- lmer(s_c4_ddiff ~ 1+(1+s_logmiost|year/ring), data = c34sum)
summary(c4d_m3)
summary(c4d_m4)

ggplot(c34sum, aes(x = s_logmiost, y = s_c4_ddiff))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, aes(group = RY))+
  geom_smooth(method = "lm", se = FALSE, col = "red", size = 2)

ggplot(c34sum, aes(x = s_logmiost, y = s_c4_ddiff))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, aes(group = id), col = "blue")+
  geom_smooth(method = "lm", se = FALSE, aes(group = RY), col = "brown")+
  geom_smooth(method = "lm", se = FALSE, col = "red", size = 2)

rand_is_RY <- ddply(c34sum, .(year), function(x){
  slope <- coef(lm(s_c4_ddiff ~ s_logmiost, data = x))[2]
  inter <- mean(x$s_c4_ddiff)
  return(data.frame(slope, inter))
  })
plot(slope~inter, data = rand_is_RY)
boxplot(rand_is_RY$slope)
var(rand_is_RY$slope)
cor(rand_is_RY$inter, rand_is_RY$slope)

summary(c4d_m0)
r.squaredGLMM(c4d_m0)
r.squaredGLMM(c4d_m3)

Anova(c4d_m0, test.statistic = "F")
summary(c4d_m0)
dredge(c4d_m0, REML = F)
# c4d_m1 <- lmer(s_c4_ddiff ~ co2*(s_logmiost+s_temp)+(1+s_logmiost|ring)+(1+s_logmiost|year)+(1|id)+(1|RY), data = c34sum)

# c4d_m1 <- lmer(s_c4_ddiff ~ co2+s_logmiost+(1|ring)+(1|year)+(1|id)+(1|RY), data = c34sum)
# c4d_m1 <- lmer(s_c4_ddiff ~ co2+s_logmiost+(1|ring)+(1|year)+(1|RY), data = c34sum)
# c4d_m1 <- lmer(s_c4_ddiff ~ co2+s_logmiost+(1|ring)+(1|id)+(1|year), data = c34sum)
# c4d_m1 <- lmer(s_c4_ddiff ~ co2+s_logmiost+(1|ring)+(1|year)+(1|RY)+(1|id), data = c34sum)
# c4d_m1 <- lmer(s_c4_ddiff ~ co2+s_logmiost+(1|year)+(1|RY), data = c34sum)
Anova(c4d_m1, test.statistic = "F")
plot(c4d_m1)
qqnorm(resid(c4d_m1))
qqline(resid(c4d_m1))
confint(c4d_m1, method = "boot", level = .9)
summary(c4d_m1)


summary(c4d_m1)
a <- dredge(c4d_m1, REML = F, extra = "r.squaredGLMM")
aavg <- model.avg(get.models(a, subset = delta <= 2))
confint(aavg, full = TRUE)


Anova(c4d_m1, test.statistic = "F")
c4d_m1 <- lmer(s_c4_ddiff ~ co2*(s_logmiost+s_temp)  +(1+s_logmiost|ring)+(1|RY)+(1|id), data = c34sum, REML = F)


c4d_m1 <- lmer(s_c4_ddiff ~ co2*(s_logmiost+s_temp)  +(1+s_logmiost|ring)+(1+s_logmiost|year), data = c34sum, REML = F)

c4d_m1 <- lmer(s_c4_ddiff ~ co2*(s_logmiost+s_temp)  +(1+s_logmiost|ring)+(1+s_logmiost|year)+(1+s_logmiost|id), data = c34sum, REML = F)
c4d_m2 <- lmer(s_c4_ddiff ~ co2*s_logmiost  +(1+s_logmiost|ring)+(1+s_logmiost|year)+(1|id), data = c34sum, REML = F)
c4d_m2 <- lmer(s_c4_ddiff ~ co2+s_logmiost  +(1+s_logmiost|id)+(1+s_logmiost|year), data = c34sum, REML = F)
summary(c4d_m2)
dredge(c4d_m2, REML = F)
Anova(update(c4d_m2, REML = T), test.statistic = "F")

ggplot(c34sum, aes(x = s_logmiost, y = s_c4_ddiff))+
  geom_point()+
  geom_smooth(aes(group = ring),method = "lm", se = FALSE)+
  geom_smooth(method = "lm", se = FALSE, col = "red", size = 2)+
  facet_grid(. ~ co2)

ggplot(c34sum, aes(x = s_logmiost, y = s_c4_ddiff))+
  geom_point()+
  geom_smooth(aes(group = year),method = "lm", se = FALSE)+
  geom_smooth(method = "lm", se = FALSE, col = "red", size = 2)+
  facet_grid(. ~ co2)

ggplot(c34sum, aes(x = s_logmiost, y = s_c4_ddiff))+
  geom_point()+
  geom_smooth(aes(group = id),method = "lm", se = FALSE)+
  geom_smooth(method = "lm", se = FALSE, col = "red", size = 2)+
  facet_grid(. ~ co2)

c4d_m1 <- lmer(s_c4_ddiff ~ co2*(s_logmiost+s_temp)  +(1|ring)+(1|RY)+(1|id), data = c34sum, REML = F)
c4d_m2 <- lmer(s_c4_ddiff ~ co2+s_logmiost+s_temp+co2:s_logmiost+(1|ring)+(1|RY)+(1|id), data = c34sum, REML = F)
c4d_m3 <- lmer(s_c4_ddiff ~ co2+s_logmiost+s_temp+    co2:s_temp+(1+s_logmiost|ring)+(1|RY)+(1|id), data = c34sum, REML = F)
c4d_m4 <- lmer(s_c4_ddiff ~ co2+s_logmiost+s_temp+(1+s_logmiost|ring)+(1|RY)+(1|id), data = c34sum, REML = F)
?dredge

options(na.action = "na.fail")
af <- dredge(c4d_m1, fixed = c("s_logmiost", "s_temp", "co2"), REML = F)
ab <- get.models(af, subset = 1)[[1]]
summary(ab)
confint(ab, method = "boot", level = .9)

model.sel(c4d_m1, c4d_m2, c4d_m3)


summary(c4d_m0)
c4d_m0full <- dredge(c4d_m0, REML = F)
c4d_m0bsm <- get.models(c4d_m0full, subset = 1)[[1]]
c4d_m0abg <- model.avg(get.models(c4d_m0full, subset = delta <=2))
summary(c4d_m0abg)
Anova(c4d_m0bsm, test.statistic = "F")
coef(c4d_m0abg, full = TRUE)
confint(c4d_m0abg, full = TRUE, level = .9)




# . predicted values ------------------------------------------------------

sitedf <- c34sum %>% 
  select(ring, id, RY, co2) %>% 
  ungroup() %>% 
  distinct()
moistval <- seq(min(c34sum$s_logmiost), max(c34sum$s_logmiost), length.out = 1000)
tempval  <- median(c34sum$s_temp)


# .. prepare df -----------------------------------------------------------


# temp is median
c4d_m0_preddf <- ldply(1, function(x){
  
  newdf       <- ldply(1:10, function(y){ 
    cbind(sitedf, 
          s_logmiost  = moistval[sample(1000, nrow(sitedf), replace = TRUE)])
  })
  
  newdf <- newdf %>% 
    mutate(s_temp      = tempval[x])
  c4d_m0_pred    <- predict(c4d_m0abg, newdf, se.fit = TRUE, re.form = NA, full = TRUE)
  c4d_m0_pred_df <- cbind(c4d_m0_pred, newdf) %>% 
    mutate(lwr = fit - se.fit * 1.96,
           upr = fit + se.fit * 1.96) 
  return(c4d_m0_pred_df)
})

# reverse transform
rev_ztrans <- function(x, xsd, xmean) x * xsd + xmean

moist_sd <- sd(log(c34sum$totalmoist))
moist_m <- mean(log(c34sum$totalmoist))
c4d_sd <- sd(c34sum$c4_ddiff)
c4d_m <- mean(c34sum$c4_ddiff)

c4d_m0_preddf_rv <- c4d_m0_preddf %>% 
  mutate(r_moist = rev_ztrans(s_logmiost, xsd = moist_sd, xmean = moist_m), 
         r_fit   = rev_ztrans(fit, c4d_sd, c4d_m),
         r_lwr   = rev_ztrans(lwr, c4d_sd, c4d_m), 
         r_upr   = rev_ztrans(upr, c4d_sd, c4d_m))

c4d_p <- ggplot(c4d_m0_preddf_rv, aes(x = r_moist, y = r_fit, col = co2)) +
  geom_line()+
  geom_line(aes(y = r_lwr), linetype = "dashed") +
  geom_line(aes(y = r_upr), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  science_theme
c4d_p




# C3 abundance ------------------------------------------------------------



# > diff ------------------------------------------------------------------

# moisture
xyplot(s_c3_ddiff ~ s_logmiost | ring, group = id, data = c34sum, type=c("p", "r"))
xyplot(s_c3_ddiff ~ s_logmiost, group = ring, data = c34sum, type=c("p", "r"))
xyplot(s_c3_ddiff ~ s_logmiost, group = RY, data = c34sum, type=c("p", "r"))
xyplot(s_c3_ddiff ~ s_logmiost, group = year, data = c34sum, type=c("p", "r"))

# temperature
xyplot(s_c3_ddiff ~ s_temp | ring, group = id, data = c34sum, type=c("p", "r"))
xyplot(s_c3_ddiff ~ s_temp, group = ring, data = c34sum, type=c("p", "r"))

c3d_m0 <- lmer(s_c3_ddiff ~ co2 * (s_logmiost+s_temp) + (1|ring) + (1|RY) + (1|id), data = c34sum)
summary(c3d_m0)
c3d_m1 <- lmer(s_c3_ddiff ~ co2 * (s_logmiost+s_temp) + (1+s_logmiost|ring)+(1|RY)+(1|id), data = c34sum)
r.squaredGLMM(c3d_m0)
r.squaredGLMM(c3d_m1)
anova(c3d_m0, c3d_m1)
c3d_m1full <- dredge(c3d_m1, REML = F)
plot(c3d_m1)
qqnorm(resid(c3d_m1))
qqline(resid(c3d_m1))
summary(c3d_m1)
c3d_m1avg <- model.avg(get.models(c3d_m1full, subset = delta <= 2))
summary(c3d_m1avg)
confint(c3d_m1avg, full = TRUE)
coef(c3d_m1avg, full = TRUE)

ggplot(c34sum, aes(x = s_logmiost, y = s_c3_ddiff, col = ring))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)

