
# response ratios ---------------------------------------------------------

head(grass_DS)

rr_grass <- grass_DS %>% 
  group_by(year, variable, co2, PFG, type) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  mutate(value = value + 1) %>% 
  spread(co2, value) %>% 
  mutate(rr = log(elev / amb), 
         pfgtype = paste(PFG, type, sep = "_"))

m1 <- lmer(rr ~ (1|year) + (1|type) + (1|PFG/variable), data = rr_grass)
summary(m1)
boxplot(rr ~ year * pfgtype, data = rr_grass)
Anova(m1, test.statistic = "F")


rr_sc4 <- filter(rr_grass, pfgtype == "c4_S")
m1 <- lmer(rr ~ year + (1|variable), data = rr_sc4)
Anova(m1, test.statistic = "F")
summary(m1)
lsmeans::lsmeans(m1, year)


rr_dc4 <- filter(rr_grass, pfgtype == "c4_D")
m1 <- lmer(rr ~ year + (1|variable), data = rr_dc4)
anova(m1)
Anova(m1, test.statistic = "F")
summary(m1)
plot(m1)
qqPlot(resid(m1))



rr_sc3 <- filter(rr_grass, pfgtype == "c3_S")
m1 <- lmer(rr ~ year + (1|variable), data = rr_sc3)
Anova(m1, test.statistic = "F")
summary(m1)
lsmeans::lsmeans(m1, year)




ggplot(data = rr_grass, aes(x = year, y = rr))+
  geom_boxplot()+
  facet_grid(. ~ pfgtype)
boxplot(rr ~ pfgtype * year, rr_grass)
abline(h = 0, lty = 2)
