# species richness

c34_S <- C3grassC4 %>% 
  group_by(year, block, ring, co2, plot, id, variable, PFG) %>% 
  summarise_each(funs(sum), value) %>%
  group_by(year, block, ring, co2, plot, id, PFG) %>% 
  summarise_each(funs(sum(. > 0)), S = value) %>% 
  ungroup() %>% 
  spread(key = PFG, value = S) %>% 
  mutate(c3r = c3 / (c3 + c4)) %>% 
  arrange(id, year) %>%
  group_by(id) %>%
  mutate_each(funs(delta = c(NA, diff(.))), c3, c4, c3r) %>% 
  filter(year != "Year0") %>% 
  left_join(c34growth_moist)

  
plot(c3_delta ~ c3moist, data = c34_S, col = co2, pch = 19)
plot(c4_delta ~ c4moist, data = c34_S, col = co2, pch = 19)
plot(c4_delta ~ log(c3), data = c34_S, col = co2, pch = 19)
plot(c3r ~ swa_c3, data = c34_S, col = co2, pch = 19)
boxplot(S ~ PFG, c34_S)

mm1 <- lmer(c3_delta ~ co2 * c3moist + (1|ring) + (1|id), data = c34_S, REML = F)
dredge(mm1)

c34_prop <- PfgRDF$C3vsC4 %>% 
  arrange(id, desc(year)) %>% 
  group_by(id) %>% 
  mutate(dratios = c(NA, diff(ratios))) %>% 
  ungroup() %>% 
  filter(!is.na(dratios)) %>%
  left_join(c34growth_moist)


?specnumber
