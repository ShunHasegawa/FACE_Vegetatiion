head(C3grassC4)

c43_ratio <- C3grassC4 %>% 
  group_by(year, co2, ring, plot, PFG) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  spread(key = PFG, value) %>% 
  mutate(c43_r = c4 / c3)


c43_ratio_year0 <- c43_ratio %>%
  filter(year == "Year0") %>%
  select(ring, plot, c43_r) %>%
  rename(ratios0 = c43_r) %>% 
  left_join(filter(c43_ratio, year != "Year0"), by = c("ring", "plot")) %>% 
  droplevels()

plot(log(c43_r + .01) ~ log(ratios0 + .01), data = c43_ratio_year0)
summary(c43_ratio_year0)
m1 <- lmer(I(log(c43_r + .01)) ~ co2 * year + I(log(ratios0 + .01)) + (1|ring) + (1|ring:plot) + (1|ring:year), data = c43_ratio_year0)
plot(m1)
qqPlot(residuals(m1))
m2 <- update(m1, subset = -order(resid(m1))[1:2])
plot(m2)
qqPlot(residuals(m2))

Anova(m2, test.statistic = "F")
visreg(m2, xvar = "co2")
visreg(m1, xvar = "year", by = "co2", overlay = TRUE)
