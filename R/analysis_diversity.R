
# prepare data frames -----------------------------------------------------


# compute diversity indices
graminoid_diversity <- mutate(site_data, 
                              H = diversity(graminoid_data[, SppName_gram]),  # Shannon's index
                              S = specnumber(graminoid_data[, SppName_gram]), # number of spp
                              J = H/log(S))                                   # Pielou's evenness

forb_diversity <- mutate(site_data, 
                         H = diversity(forb_data[ , SppName_forb]),
                         S = specnumber(forb_data[ , SppName_forb]),
                         J = H/log(S))



# Move Year0 value to a new column to be used as covariate for the analysis
diversity_list <- list(graminoid = graminoid_diversity, forb = forb_diversity)


diversity_year0_list <- llply(diversity_list, function(x){
  
  DivDF_mlt <- gather(x, variable, value, H, S, J)
  
  DivDF_year0 <- DivDF_mlt %>% # Year0
    filter(year == "Year0") %>%
    select(ring, plot, value, variable) %>%
    rename(value0 = value) %>% 
    left_join(filter(DivDF_mlt, year != "Year0"), by = c("ring", "plot", "variable")) %>% 
    filter(!is.na(value))
  
  return(DivDF_year0)
})

graminoid_diversity <- diversity_year0_list$graminoid
forb_diversity      <- diversity_year0_list$forb




# analysis ----------------------------------------------------------------


# > Graminoid -------------------------------------------------------------

graminoid_h <- filter(graminoid_diversity, variable == "H")
graminoid_j <- filter(graminoid_diversity, variable == "J")
graminoid_s <- filter(graminoid_diversity, variable == "S")


# . H ---------------------------------------------------------------------


# model
gr_h_m1 <- lmer(value ~ co2 * year + value0 + (1|ring) + (1|ring:plot) + (1|ring:year), data = graminoid_h)


# model diagnosis
plot(gr_h_m1)
qqnorm(resid(gr_h_m1))
qqline(resid(gr_h_m1))


# F test
Anova(gr_h_m1, test.statistic = "F")


# Pairwise comparisons and 95% confidence intervals for covariate-adjusted means
gr_h_lsmean <- lsmeans::lsmeans(gr_h_m1, ~ co2 | year)
summary(pairs(gr_h_lsmean)[1:3], adjust = "none")
summary(gr_h_lsmean)


# CO2 response ratios (RR) on covariate-adjusted means
data.frame(summary(gr_h_lsmean)) %>% 
  group_by(co2) %>% 
  summarise(value = mean(lsmean)) %>%
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1)





# . J -----------------------------------------------------------------------

# model
gr_j_m1 <- lmer(value ~ co2 * year + value0 + (1|ring) + (1|ring:plot) + (1|ring:year), data = graminoid_j)


# model diagnosis
plot(gr_j_m1)
qqnorm(resid(gr_j_m1))
qqline(resid(gr_j_m1))


# non-normality of the data is suggested
gr_j_m2 <- update(gr_j_m1, subset = -which.min(resid(gr_j_m1)))
plot(gr_j_m2)
qqnorm(resid(gr_j_m2))
qqline(resid(gr_j_m2))


# F test
Anova(gr_j_m1, test.statistic = "F")
Anova(gr_j_m2, test.statistic = "F")
  # not majoer difference so present the first one withought removing the outlier


# 95% confidence intervals for covariate-adjusted means
gr_j_lsmean <- lsmeans::lsmeans(gr_j_m1, ~ co2 | year)
summary(gr_j_lsmean)


# CO2 response ratios (RR) on covariate-adjusted means
data.frame(summary(gr_j_lsmean)) %>% 
  group_by(co2) %>% 
  summarise(value = mean(lsmean)) %>%
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1)




# . S ---------------------------------------------------------------------


# model
gr_s_m1 <- lmer(value ~ co2 * year + value0 + (1|ring) + (1|ring:plot) + (1|ring:year), data = graminoid_s)


# model diagnosis
plot(gr_s_m1)
qqnorm(resid(gr_s_m1))
qqline(resid(gr_s_m1))


# F test
Anova(gr_s_m1, test.statistic = "F")


# 95% confidence intervals for covariate-adjusted means
gr_s_lsmean <- lsmeans::lsmeans(gr_s_m1, ~ co2 | year)
summary(gr_s_lsmean)


# CO2 response ratios (RR) on covariate-adjusted means
data.frame(summary(gr_s_lsmean)) %>% 
  group_by(co2) %>% 
  summarise(value = mean(lsmean)) %>%
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1)




# > Forb --------------------------------------------------------------------


forb_h <- filter(forb_diversity, variable == "H")
forb_j <- filter(forb_diversity, variable == "J")
forb_s <- filter(forb_diversity, variable == "S")


# . H ---------------------------------------------------------------------


# model
fo_h_m1 <- lmer(value ~ co2 * year + value0 + (1|ring) + (1|ring:plot) + (1|ring:year), data = forb_h)


# model diagnosis
plot(fo_h_m1)
qqnorm(resid(fo_h_m1))
qqline(resid(fo_h_m1))


# non-normality of the data is suggested
fo_h_m2 <- update(fo_h_m1, subset = -which.max(resid(fo_h_m1)))
plot(fo_h_m2)
qqnorm(resid(fo_h_m2))
qqline(resid(fo_h_m2))


# F test
Anova(fo_h_m1, test.statistic = "F")
Anova(fo_h_m2, test.statistic = "F")
# no major difference, so present the first one without removing the outlier


# 95% confidence intervals for covariate-adjusted means
fo_h_lsmean <- lsmeans::lsmeans(fo_h_m1, ~ co2 | year)
summary(fo_h_lsmean)


# CO2 response ratios (RR) on covariate-adjusted means
data.frame(summary(fo_h_lsmean)) %>% 
  group_by(co2) %>% 
  summarise(value = mean(lsmean)) %>%
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1)




# . J -----------------------------------------------------------------------


# model
fo_j_m1 <- lmer(value ~ co2 * year + value0 + (1|ring) + (1|ring:plot) + (1|ring:year), data = forb_j)


# model diagnosis
plot(fo_j_m1)
qqnorm(resid(fo_j_m1))
qqline(resid(fo_j_m1))


# F test
Anova(fo_j_m1, test.statistic = "F")


# 95% confidence intervals for covariate-adjusted means
fo_j_lsmean <- lsmeans::lsmeans(fo_j_m1, ~ co2 | year)
summary(fo_j_lsmean)


# CO2 response ratios (RR) on covariate-adjusted means
data.frame(summary(fo_j_lsmean)) %>% 
  group_by(co2) %>% 
  summarise(value = mean(lsmean)) %>%
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1)




# . S ---------------------------------------------------------------------


# model
fo_s_m1 <- lmer(value ~ co2 * year + value0 + (1|ring) + (1|ring:plot) + (1|ring:year), data = forb_s)


# model diagnosis
plot(fo_s_m1)
qqnorm(resid(fo_s_m1))
qqline(resid(fo_s_m1))


# F test
Anova(fo_s_m1, test.statistic = "F")


# 95% confidence intervals for covariate-adjusted means
fo_s_lsmean <- lsmeans::lsmeans(fo_s_m1, ~ co2 | year)
summary(fo_s_lsmean)


# CO2 response ratios (RR) on covariate-adjusted means
data.frame(summary(fo_s_lsmean)) %>% 
  group_by(co2) %>% 
  summarise(value = mean(lsmean)) %>%
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1)

