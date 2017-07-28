
# identify dominant-subordinate spp ---------------------------------------

# dominant or subordinate spp were deifnied as those that occur across all the
# study years (but no need to be in all rings)

ds_spp <- PlotSumVeg[, c(SppName_grass, 
                         "year", "ring", "plot", "block", "co2", "id", "RY")] %>% 
  gather(key = species, abund, one_of(SppName_grass)) %>% 
  group_by(species, year) %>% 
  summarise(abund = sum(abund)) %>% 
  group_by(species) %>% 
  summarise(DS = !any(abund == 0)) %>% 
  ungroup() %>% 
  filter(DS)


# set an arbitrary threshold between dominant and subordinate species at
# relative abundance of 0.1 (see the fig below)

grass_df <- PlotSumVeg[, c(SppName_grass, 
                           "year", "ring", "plot", "block", "co2", "id", "RY")] %>% 
  gather(key = species, abund, one_of(ds_spp$species)) %>% 
  group_by(species) %>% 
  summarise(abund = mean(abund)) %>% 
  ungroup() %>% 
  arrange(abund) %>% 
  mutate(r_abund = abund / max(abund),
         cumcov = cumsum(abund),
         cumcov = cumcov * 100 / max(cumcov),
         type = ifelse(r_abund > .1, "D", "S"),
         species = factor(species, 
                          levels = species[order(abund, decreasing = TRUE)]))

ggplot(grass_df, aes(x = as.numeric(species), y = r_abund, label = type))+
  geom_path()+
  geom_text(aes(col = type)) +
  geom_hline(yintercept = .1, col = "blue", linetype = "dashed") +
  scale_x_continuous(breaks = seq(1, nrow(grass_df), 1),
                     labels =  as.character(grass_df$species[order(grass_df$abund, decreasing = TRUE)]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = NULL, y = "Relative abundance")




# analysis ----------------------------------------------------------------

# df for dominant and subordinate spp
grass_DS <- veg_FullVdf %>% 
  filter(variable %in% ds_spp$species) %>% 
  left_join(grass_df[, c("species", "type")], by = c("variable" = "species")) %>% 
  droplevels(.)


# . subordinate: dominant ratio -------------------------------------------

# prepare df (use Year0 as a covariate)
DS_ratio <- grass_DS %>% 
  group_by(year, ring, co2, plot, id, RY, type) %>% 
  summarise(abund = sum(value)) %>% 
  ungroup() %>% 
  spread(key = type, value = abund) %>% 
  mutate(total = D + S,
         sd_ratio = S / D)

DS_ratio_y0 <- DS_ratio %>% 
  filter(year == "Year0") %>% 
  rename(sd_ratio0 = sd_ratio,
         d0        = D,
         s0        = S) %>% 
  select(id, sd_ratio0, d0, s0)

DS_ratio_ed <- DS_ratio %>% 
  filter(year != "Year0") %>% 
  left_join(DS_ratio_y0) %>% 
  mutate(logit_ratio = logit(sd_ratio),
         logit_ratio0 = logit(sd_ratio0))

plot(sqrt(sd_ratio) ~ sqrt(sd_ratio0), DS_ratio_ed, pch = 19, col = co2)
m1 <- lmer(sqrt(sd_ratio) ~ co2 * year + sqrt(sd_ratio0) + (1|ring) + (1|id) + (1|RY), data = DS_ratio_ed)
Anova(m1, test.statistic = "F")
summary(m1)
plot(m1)
qqPlot(resid(m1))
r.squared(m1)
visreg(m1, xvar = "sd_ratio0", by = "co2", overlay = TRUE, xtrans = sqrt)

# Dominant
plot(D ~ d0, DS_ratio_ed, pch = 19, col = co2)
d_m1 <- lmer(D ~ co2 * year + d0 + (1|ring) + (1|id) + (1|RY), data = DS_ratio_ed)
Anova(d_m1, test.statistic = "F")
plot(d_m1)
qqPlot(resid(d_m1))


# Subordinate 
plot(S ~ s0, DS_ratio_ed, pch = 19, col = co2)
plot(sqrt(S) ~ sqrt(s0), DS_ratio_ed, pch = 19, col = co2)
s_m1 <- lmer(sqrt(S) ~ co2 * year + sqrt(s0) + (1|ring) + (1|id) + (1|RY), data = DS_ratio_ed)
Anova(s_m1, test.statistic = "F")
plot(s_m1)
qqPlot(resid(s_m1))



# Paspalidium.distans -----------------------------------------------------


grass_DS %>% 
  filter(type == "S") %>% 
  group_by(year, co2, plot, id, RY, PFG, variable) %>% 
  summarise(value = sum(value))

ddd <- filter(grass_DS, variable == "Paspalidium.distans") %>% 
  group_by(year, co2, ring, plot, id, RY) %>% 
  summarise(value = sum(value)) %>% 
  ungroup()
dd0 <- ddd %>% 
  filter(year == "Year0") %>% 
  rename(value0 = value) %>% 
  select(id, value0)
ddd2 <- ddd %>% 
  filter(year != "Year0") %>% 
  left_join(dd0) %>% 
  mutate(YC = factor(paste(year, co2, sep = "")))

plot(value ~ value0, data = ddd2, col = co2, pch = 19)
plot(sqrt(value) ~ value0, data = ddd2, col = co2, pch = 19)
plot(log(value) ~ log(value0), data = ddd2, col = co2, pch = 19)
m1 <- lmer(sqrt(value) ~ co2 * year + value0 + (1|ring) + (1|id) + (1|RY), data = ddd2)  
summary(m1)
Anova(m1, test.statistic = "F")
plot(m1)
qqPlot(resid(m1))
which.max(resid(m1))


grass_DS_sum <- grass_DS %>% 
  group_by(year, ring, co2, plot, id, RY, PFG, type) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  mutate(pfg_type = paste(type, PFG, sep = "_"))



 

# analysis on Dminant and subordinate C3 and C4 species -------------------


# > prepare data frame ----------------------------------------------------

# year0
grass_DS_sum_y0 <- grass_DS_sum %>% 
  filter(year == "Year0") %>% 
  rename(value0 = value) %>% 
  select(id, value0, pfg_type, PFG, type)

# merge
grass_DS_dd <- grass_DS_sum %>% 
  filter(year != "Year0") %>% 
  left_join(grass_DS_sum_y0)

ggplot(grass_DS_dd, aes(x = sqrt(value0 + 1), y = sqrt(value + 1)))+
  geom_point() +
  facet_wrap(PFG ~ type, scale = "free")

# dominant C3
dc3_dd <- filter(grass_DS_dd, pfg_type == "D_c3")
dc3_dd_ed <- dc3_dd %>% 
  rename(dc3 = value) %>% 
  select(year, id, dc3)


# dominant C4
dc4_dd <- filter(grass_DS_dd, pfg_type == "D_c4")
dc4_dd_ed <- dc4_dd %>% 
  rename(dc4 = value) %>% 
  select(year, id, dc4)


# subordinate C3
sc3_dd <- grass_DS_dd %>% 
  filter(pfg_type == "S_c3") %>%
  left_join(dc3_dd_ed) %>% 
  left_join(dc4_dd_ed)

  
# subordinate C4  
sc4_dd <- grass_DS_dd %>% 
  filter(pfg_type == "S_c4") %>% 
  left_join(dc3_dd_ed) %>% 
  left_join(dc4_dd_ed)




# > analysis --------------------------------------------------------------

dc3_m1 <- lmer(value ~ co2 * year + value0 + (1|ring) + (1|id) + (1|RY), dc3_dd)
Anova(dc3_m1, test.statistic = "F")

dc4_m1 <- lmer(sqrt(value + 1) ~ co2 * year + sqrt(value0 + 1) + (1|ring) + (1|id) + (1|RY), dc4_dd)
Anova(dc4_m1, test.statistic = "F")
summary(dc4_m1)

sc3_m1 <- glmer(value/100 ~ co2 * year + I(value0/100) + (1|ring) + (1|id) + (1|RY), 
                family = "binomial", sc3_dd)

sc3_m1 <- lmer(value ~ value0 + dc3 + dc4 + (1|ring) + (1|id) + (1|RY), sc3_dd)
Anova(sc3_m1, test.statistic = "F")

sc4_m1 <- lmer(sqrt(value + 1) ~ co2 * year + sqrt(value0 + 1) + (1|ring) + (1|id) + (1|RY), sc4_dd)
qqPlot(resid(sc4_m1))
plot(sc4_m1)

Anova(sc4_m1, test.statistic = "F")
Anova(get.models(sc4_f, subset = 1)[[1]], test.statistic = "F")


summary(sc4_m1)