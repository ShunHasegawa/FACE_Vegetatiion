# Prepare data frames -----------------------------------------------------

summary(veg_FullVdf)

# proportion of each form for each plot
PropDF <- veg_FullVdf %>%
  filter(form %in% c("Forb", "Grass")) %>% 
  group_by(year, co2, block, ring, plot, id, form) %>% 
  summarise(value = sum(value)) %>% 
  spread(form, value) %>% 
  mutate(Total      = Forb + Grass, 
         grass_prop = Grass / Total) %>% 
  ungroup()
  
# Move Year0 value to a new column to be used as covariate in the analyssis
# below

# subsequent years
subyear_dd <- filter(PropDF, year != "Year0") 

# merge with Year0
PropDF_year0 <- PropDF %>% 
  filter(year == "Year0") %>%
  select(id, grass_prop) %>%
  rename(value0 = grass_prop) %>%
  left_join(subyear_dd, by = "id") %>% 
  mutate(logitv0 = logit(value0),
         sqrtv0  = sqrt(value0),
         logv0   = log(value0 + 1)) # logit transformation

# plot
par(mfrow = c(1, 2))
plot(grass_prop ~ value0, pch = 19, col = year, data = PropDF_year0, main = "raw")
plot(logit(grass_prop) ~ logitv0, pch = 19, col = year, data = PropDF_year0, main = "logit")




# Analysis ----------------------------------------------------------------

m1 <- lmer(logit(grass_prop) ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), data = PropDF_year0)
m2 <- lmer(logit(grass_prop) ~ co2 * year + logitv0 + (1|block) + (1|ring) + (1|id), data = PropDF_year0)
m3 <- lmer(logit(grass_prop) ~ co2 * year + sqrtv0 + (1|block) + (1|ring) + (1|id), data = PropDF_year0)
m4 <- lmer(logit(grass_prop) ~ co2 * year + logv0 + (1|block) + (1|ring) + (1|id), data = PropDF_year0)
AICc(m1, m2, m3, m4)
# m4 is slightly better

m5 <- update(m4, ~ . - (1|block))
AICc(m4, m5)




# model diagnosis ---------------------------------------------------------

plot(m5)
qqnorm(resid(m5))
qqline(resid(m5))

which.max(resid(m5))
m6 <- update(m2, subset = -43)
plot(m6)
llply(list(m5, m6), function(x) Anova(x, test.statistic = "F"))
# CO2xYear interaction was driven by outlier so remove
grassprop_m <- m6

# CI and post-hoc ---------------------------------------------------------


# compute 95 CI and post-hoc test
lsmeans_dd <- lsmeans::lsmeans(m6, ~ co2 | year)

# 95% CI
CI_dd <- data.frame(summary(lsmeans_dd))

# post-hoc test
contrast_dd <- data.frame(summary(pairs(lsmeans_dd)[1:3], adjust = "fdr")) %>% 
  mutate(co2  = factor("amb", levels = c("amb", "elev")),
         star = get_star(p.value)) %>% 
  select(year, co2, p.value, star)

# merge
grassprop_ci_dd <- left_join(CI_dd, contrast_dd, by = c("year", "co2")) %>% 
  mutate(year = factor(year, levels = paste0("Year", 0:3)),
         rlsmean = boot::inv.logit(lsmean), # reverse transform
         rlowerCL = boot::inv.logit(lower.CL),
         rupperCL = boot::inv.logit(upper.CL),
         year = factor(year, levels = paste0("Year", 0:3)),
         form = "Grass")
grassprop_ci_dd$star[is.na(grassprop_ci_dd$star)] <- ""

# df for observed values
grassprop_d <- PropDF %>%
  group_by(year, co2, ring) %>% 
  summarise_each(funs(sum), Forb, Grass) %>% 
  ungroup() %>% 
  mutate(grass_prop = Grass / (Grass + Forb),
         year = factor(year, levels = paste0("Year", 0:3)))
