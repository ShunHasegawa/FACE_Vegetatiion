
# prepare data frame ------------------------------------------------------


# Compute total log annual change ratios (LAR) for C3 and C4 graminoids, and
# merge with environmental variables

lar_data <- graminoid_pfg_df %>% 
  mutate(ring = as.character(ring), plot = as.character(plot)) %>% 
  arrange(variable, ring, plot, year) %>%
  mutate(value = value + 1) %>%                 # some speceis were absent, so add 1
  group_by(co2, ring, plot, pfg, variable) %>%  
  mutate(lar = log(value / lag(value, 1))) %>%  # comput LAR for each species
  filter(year != "Year0") %>%
  group_by(year, co2, ring, plot, pfg) %>% 
  summarise(total_lar = sum(lar)) %>%           # compute total LAR for each PFG
  ungroup() %>% 
  spread(key = pfg, value = total_lar) %>% 
  rename(lar_c3 = c3, lar_c4 = c4) %>% 
  left_join(env_data) %>% 
  mutate(s_lar_c3    = scale(lar_c3)[, 1],         # Z-standardize numeric variables
         s_lar_c4    = scale(lar_c4)[, 1],
         s_logmoist  = scale(log(Moist))[, 1],
         s_temp      = scale(Temp)[, 1],
         s_logpar    = scale(log(PAR))[, 1])


# analysis  ------------------------------------------------------------

options(na.action = "na.fail")


# > LAR_C4 ----------------------------------------------------------------


# test interactions
lar_c4_m1 <- lmer(s_lar_c4 ~ co2 * (s_logmoist+s_temp+s_logpar) + 
                    (1|ring)+(1|ring:plot)+(1|ring:year), data = lar_data)

dredge(lar_c4_m1, REML = F, fixed = c("co2", "s_logmoist", "s_temp", "s_logpar"))
dredge(lar_c4_m1, REML = F)
  # no interaction is suggested


# test main effects
lar_c4_m2 <- lmer(s_lar_c4 ~ co2 + s_logmoist + s_temp +  s_logpar + 
                    (1|ring)+(1|ring:plot)+(1|ring:year), data = lar_data)


# model diagnosis
plot(lar_c4_m2)
qqnorm(resid(lar_c4_m2))
qqline(resid(lar_c4_m2))


# non-normality of the data was suggested
lar_c4_m3 <- update(lar_c4_m2, subset = -which.min(resid(lar_c4_m2)))
plot(lar_c4_m3)
qqnorm(resid(lar_c4_m3))
qqline(resid(lar_c4_m3))


# get 95% and 90% CI for coefficients
summary(lar_c4_m3)
r.squaredGLMM(lar_c4_m3)
confint(lar_c4_m3, method = "boot", nsim = 999)
confint(lar_c4_m3, method = "boot", level = .9, nsim = 999)




# . level plot ------------------------------------------------------------

# plot predicted values from the model above on a level plot to show the
# relations of LAR_C4 to soil moisture and PAR at a given temperature

# predict values
c4_envdf <- with(lar_data, expand.grid(s_logmoist = seq(-2.5, 2.5, length.out = 100),
                                     s_logpar   = seq(-2.5, 2.5, length.out = 100),
                                     s_temp     = median(s_temp),
                                     co2        = c("amb", "elev"), 
                                     ring = 1, id = "1:1", RY = "1:Year1"))


c4_pred <- predict(lar_c4_m3, c4_envdf, re.form = NA)


# reverse transform

rev_ztrans <- function(x, xsd, xmean){
  x * xsd + xmean
}

moist_sd <- sd(log(lar_data$Moist))
moist_m  <- mean(log(lar_data$Moist))
temp_sd  <- sd(lar_data$Temp)
temp_m   <- mean(lar_data$Temp)
par_sd   <- sd(log(lar_data$PAR))
par_m    <- mean(log(lar_data$PAR))
c4d_sd   <- sd(lar_data$lar_c4)
c4d_m    <- mean(lar_data$lar_c4)

c4_pred_df <- cbind(fit = c4_pred, c4_envdf) %>% 
  mutate(r_moist = rev_ztrans(s_logmoist, xsd = moist_sd, xmean = moist_m),
         r_temp  = rev_ztrans(s_temp, xsd = temp_sd, xmean = temp_m),
         r_par   = rev_ztrans(s_logpar, xsd = par_sd, xmean = par_m),
         r_fit   = rev_ztrans(fit, c4d_sd, c4d_m),
         co2     = mapvalues(co2, c("amb", "elev"), c("Ambient", "eCO[2]")))

lar_data_temp <- lar_data %>% 
  mutate(co2     = mapvalues(co2, c("amb", "elev"), c("Ambient", "eCO[2]")))

c4_levelplot <- ggplot(c4_pred_df, aes(x = r_moist, y = r_par)) + 
  geom_tile(aes(fill = r_fit)) +
  scale_fill_gradient2(expression(LAR[C4]), low = "blue",high = "red", mid = "white")+
  stat_ellipse(data = lar_data_temp, aes(x = log(Moist), y = log(PAR), linetype = year), 
               type = "norm", level = .7)+
  geom_point(data = lar_data_temp, aes(x = log(Moist), y = log(PAR), shape = year), size = 2)+
  scale_shape_manual("Year", values = c(0:2), label = 2013:2015)+
  scale_linetype_manual("Year", values = c(1:3), label = 2013:2015)+
  facet_grid(. ~ co2, labeller = label_parsed)+
  labs(x = "ln(Moist)", y = expression(ln(PAR,~mu*mol~s^'-1'~m^"-2")))


# > LAR_C3 ------------------------------------------------------------------


# test interaction
lar_c3_m1 <- lmer(s_lar_c3 ~ co2 * (s_logmoist+s_temp+s_logpar) 
                  + (1|ring) + (1|ring:plot) +(1|ring:year), data = lar_data)
plot(lar_c3_m1)
qqnorm(resid(lar_c3_m1))
qqline(resid(lar_c3_m1))
  # one large outlier
which.max(resid(lar_c3_m1))


lar_c3_m2 <- lmer(s_lar_c3 ~ co2 * (s_logmoist+s_temp+s_logpar) 
                  + (1|ring) + (1|ring:plot) +(1|ring:year),  
                  data = lar_data[-56, ])
plot(lar_c3_m2)
qqnorm(resid(lar_c3_m2))
qqline(resid(lar_c3_m2))


dredge(lar_c3_m1, REML = F, fixed = c("co2", "s_logmoist", "s_temp", "s_logpar"))
dredge(lar_c3_m2, REML = F, fixed = c("co2", "s_logmoist", "s_temp", "s_logpar"))
  # no interaction is suggsted


# test main effects
lar_c3_m3 <- lmer(s_lar_c3 ~ co2 + s_logmoist + s_temp + s_logpar 
                  + (1|ring) + (1|ring:plot) +(1|ring:year), 
                  data = lar_data[-56, ])


# get 95% and 90% CI for coefficients
summary(lar_c3_m3)
r.squaredGLMM(lar_c3_m3)
confint(lar_c3_m3, method = "boot", nsim = 999)
confint(lar_c3_m3, method = "boot", level = .9, nsim = 999)

