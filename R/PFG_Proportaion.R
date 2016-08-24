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
  mutate(logitv0 = logit(value0)) # logit transformation

# plot
par(mfrow = c(1, 2))
plot(grass_prop ~ value0, pch = 19, col = year, data = PropDF_year0, main = "raw")
plot(logit(grass_prop) ~ logitv0, pch = 19, col = year, data = PropDF_year0, main = "logit")

# Analysis ----------------------------------------------------------------

m1 <- lmer(logit(grass_prop) ~ co2 * year + logitv0 + (1|block) + (1|ring) + (1|id), 
           data = PropDF_year0)
m2 <- update(m1, ~ . - (1|block))
AICc(m1, m2)
grassprop_m <- m2

# compute 95 CI and post-hoc test
lsmeans_dd <- summary(lsmeans::lsmeans(m2, pairwise ~ co2 | year))

# 95% CI
CI_dd <- data.frame(lsmeans_dd$lsmeans)

# post-hoc test
contrast_dd <- data.frame(lsmeans_dd$contrast) %>% 
  mutate(co2 = factor("amb", levels = c("amb", "elev")),
         star = cut(p.value, right = FALSE,
                    breaks = c(0, .1, .05, .01, .001, 1),  
                    labels = c("***", "**", "*", "\u2020", ""))) %>% 
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

# df for Year0
grassprop_d <- PropDF %>%
  filter(year == "Year0") %>% 
  group_by(year, co2, ring) %>% 
  summarise_each(funs(sum), Forb, Grass) %>% 
  ungroup() %>% 
  mutate(grass_prop = Grass / (Grass + Forb),
         year = factor(year, levels = paste0("Year", 0:3)))
