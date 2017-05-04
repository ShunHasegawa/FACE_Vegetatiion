
# prepare data frame ------------------------------------------------------

# reshape graminoid_data to a long-format and combine with plant functional groups

graminoid_pfg_df <- graminoid_data %>% 
  gather(key = variable, value = value, one_of(SppName_gram)) %>% 
  left_join(sp_pfg)


# data frame for microlaena and cynodon
mic_cyn <- graminoid_pfg_df %>% 
  filter(variable %in% c("Microlaena.stipoides", "Cynodon.dactylon")) %>% 
  select(-pfg)



# data frame for total abundance of C3 and C4
C34_abund <- graminoid_pfg_df %>%
  group_by(pfg, year, co2, ring, plot) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  mutate(pfg = factor(paste(pfg, "total", sep = "_"))) %>% 
  rename(variable = pfg)


# merge the data frames above
abund_data <- bind_rows(mic_cyn, C34_abund)


# move value in Year0 to a new column in order to use it as a covariate
abund_data_year0 <- abund_data %>%
  filter(year == "Year0") %>%
  select(ring, plot, value, variable) %>%
  rename(value0 = value) %>% 
  left_join(filter(abund_data, year != "Year0"), by = c("ring", "plot", "variable")) 



# analysis ----------------------------------------------------------------


abund_mcr <- filter(abund_data_year0, variable == "Microlaena.stipoides")
abund_cyn <- filter(abund_data_year0, variable == "Cynodon.dactylon")
abund_c3  <- filter(abund_data_year0, variable == "c3_total")
abund_c4  <- filter(abund_data_year0, variable == "c4_total")


# > Microlaena --------------------------------------------------------------


# model
plot(value ~ value0, data = abund_mcr)
plot(logit(value) ~ sqrt(value0), data = abund_mcr)
  # Microlaena is bounded to the maximam value (100) so transformed with logit.
  # value0 is transformed with sqrt for a better linearity

mcr_m1 <- lmer(I(car::logit(value / 100)) ~ co2 * year + I(sqrt(value0)) + (1|ring) + (1|ring:plot) + (1|ring:year), data = abund_mcr)


# model diagnosis
plot(mcr_m1)
qqPlot(residuals(mcr_m1))


# non-normality of the data is suggested
mcr_m2 <- update(mcr_m1, subset = -c(order(resid(mcr_m1))[1:3]))
plot(mcr_m2)
qqPlot(residuals(mcr_m2))


# F test
Anova(mcr_m1, test.statistic = "F")
Anova(mcr_m2, test.statistic = "F")
  # no difference so present the first one without removing the outliers


# 95% confidence intervals for covariate-adjusted means
mcr_lsmean <- lsmeans::lsmeans(mcr_m1, ~ co2 | year)
mr_95CI <- data.frame(summary(mcr_lsmean)) %>% 
  mutate_each(funs(r = boot::inv.logit(.) * 100 / 4), lsmean, lower.CL, upper.CL)  # response variable was devided by 100 and logit-transformed, so estimates and assocaited 95% CI were reverse-transformed. Then, the reverse-transformed values were devided by 4 to express on the basis of 1 x 1 m (plot = 2 x 2 m)
mr_95CI


# CO2 response ratios (RR) on covariate-adjusted means
mr_95CI %>% 
  group_by(co2) %>% 
  summarise(value = mean(lsmean_r)) %>%
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1)




# > Cynodon --------------------------------------------------------------


# model
plot(value ~ value0, data = abund_cyn)
plot(sqrt(value) ~ value0, data = abund_cyn)
  # transform with sqrt for a better linearity
cyn_m1 <- lmer(sqrt(value) ~ co2 * year + value0 + (1|ring) + (1|ring:plot) + (1|ring:year), data = abund_cyn)


# model diagnosis
plot(cyn_m1)
qqPlot(residuals(cyn_m1))


# F test
Anova(cyn_m1, test.statistic = "F")


# pairwise comparisons and 95% confidence intervals for covariate-adjusted means
cyn_lsmean <- lsmeans::lsmeans(cyn_m1, ~ co2 | year)
summary(pairs(cyn_lsmean)[1:3], adjust = "none")
cyn_95CI <- data.frame(summary(cyn_lsmean)) %>% 
  mutate_each(funs(.^2 / 4), lsmean, lower.CL, upper.CL)
cyn_95CI


# CO2 response ratios (RR) on covariate-adjusted means
cyn_95CI %>% 
  group_by(co2) %>% 
  summarise(value = mean(lsmean)) %>%
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1)




# > C3 abundance --------------------------------------------------------------


# model
plot(value ~ value0, abund_c3)
plot(value ~ log(value0), abund_c3)
  # value0 is log-transformed for a better linearity

c3_m1 <- lmer(value ~ co2 * year + I(log(value0)) + (1|ring) + (1|ring:plot) + (1|ring:year), data = abund_c3)


# model diagnosis
plot(c3_m1)
qqPlot(residuals(c3_m1))


# non-normality of the data is sugested
c3_m2 <- update(c3_m1, subset = -which(abs(resid(c3_m1)) > 20))
plot(c3_m2)
qqPlot(residuals(c3_m2))


# F test
Anova(c3_m1, test.statistic = "F")
Anova(c3_m2, test.statistic = "F")
  # no major difference, so present the first model without removing the outliers


# 95% confidence intervals for covariate-adjusted means
c3_lsmean <- lsmeans::lsmeans(c3_m1, ~ co2 | year)
c3_95CI <- data.frame(summary(c3_lsmean)) %>% 
  mutate_each(funs(. / 4), lsmean, lower.CL, upper.CL)
c3_95CI


# CO2 response ratios (RR) on covariate-adjusted means
c3_95CI %>% 
  group_by(co2) %>% 
  summarise(value = mean(lsmean)) %>%
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1)




# > C4 abundance --------------------------------------------------------------


# model
c4_m1 <- lmer(value ~ co2 * year + value0 + (1|ring) + (1|ring:plot) + (1|ring:year), data = abund_c4)


# model diagnosis
plot(c4_m1)
qqPlot(residuals(c4_m1))


# F test
Anova(c4_m1, test.statistic = "F")


# pairwise comparisons and 95% confidence intervals for covariate-adjusted means
c4_lsmean <- lsmeans::lsmeans(c4_m1, ~ co2 | year)
summary(pairs(c4_lsmean)[1:3], adjust = "none")
c4_95CI <- data.frame(summary(c4_lsmean)) %>% 
  mutate_each(funs(. / 4), lsmean, lower.CL, upper.CL)
c4_95CI


# CO2 response ratios (RR) on covariate-adjusted means
c4_95CI %>% 
  group_by(co2) %>% 
  summarise(value = mean(lsmean)) %>%
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1)

