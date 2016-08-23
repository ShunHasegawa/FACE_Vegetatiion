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

Anova(m2, test.statistic = "F")

plot(m2)
qqnorm(resid(m2))  
qqline(resid(m2))  

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
ci_dd <- left_join(CI_dd, contrast_dd, by = c("year", "co2")) %>% 
  mutate(year = factor(year, levels = paste0("Year", 0:3)),
         rlsmean = boot::inv.logit(lsmean), # reverse transform
         rlowerCL = boot::inv.logit(lower.CL),
         rupperCL = boot::inv.logit(upper.CL),
         year = factor(year, levels = paste0("Year", 0:3)),
         form = "Grass")
ci_dd$star[is.na(ci_dd$star)] <- ""

# Create fig --------------------------------------------------------------

# df for Year0
d <- PropDF %>%
  filter(year == "Year0") %>% 
  mutate(year = factor(year, levels = paste0("Year", 0:3)))

# df for median of Year0
d_med <- d %>%  
  summarise(M = median(grass_prop)) %>% 
  mutate(Med = "Md[Year0]")

# fig
dodgeval <- .4
p <- ggplot(ci_dd, aes(x = year, y = rlsmean, fill = co2, group = co2)) +
  
  geom_errorbar(aes(ymin = rlowerCL, ymax = rupperCL), width = 0, 
                position = position_dodge(width = dodgeval)) +
  geom_line(aes(linetype = co2), position = position_dodge(width = dodgeval)) +
  geom_boxplot(data = d, aes(x = year, y = grass_prop),  alpha = .6, 
               position = position_dodge(.7), outlier.shape = 21, width = .7, 
               show.legend = FALSE) +
  geom_point(shape = 21, size = 3, position = position_dodge(width = dodgeval)) +
  geom_hline(data = d_med, aes(yintercept = M, col = Med), alpha = .8) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_text(aes(label = star, y = rupperCL), fontface = "bold", vjust = 0) +
  
  scale_fill_manual(values = c("black", "white"), 
                    labels = c("Ambient", expression(eCO[2]))) +
  scale_linetype_manual(values = c("solid", "dashed"), 
                        labels = c("Ambient", expression(eCO[2]))) +
  scale_color_manual(values = "grey50", labels = expression(Md[Year0])) +
  scale_x_discrete("Year", labels = 0:4, drop = FALSE) +
  scale_y_continuous(breaks = seq(.3, .9, by = .1), 
                     labels = paste(seq(.3, .9, by = .1), 
                                    1 - seq(.3, .9, by = .1),
                                    sep = " | ")) +
  
  science_theme +
  theme(legend.position = c(.8, .85),
        legend.margin = unit(-.5, "line"),
        strip.text.x = element_text(face = "italic")) +
  guides(linetype = guide_legend(order = 1),
         fill     = guide_legend(order = 1),
         col      = guide_legend(order = 2)) +
  
  labs(y = "Adjusted proportion (Grass | Forb)")

ggsavePP(filename = "output/figs/adjusted_Grass-Forb_proportion",
         plot = p, width = 4, height = 3)
