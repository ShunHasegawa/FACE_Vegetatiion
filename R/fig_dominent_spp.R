
# prepare data frame ------------------------------------------------------

# create models to be tests
dominentsp <- veg_FullVdf %>% 
  filter(variable %in% DmSpp) %>% 
  group_by(variable, year, block, co2, ring, plot, id) %>% 
  summarise(value = sum(value)) %>% 
  ungroup()

dominentsp_year0 <- dominentsp %>%
  filter(year == "Year0") %>%
  select(id, value, variable) %>%
  rename(value0 = value) %>% 
  left_join(filter(dominentsp, year != "Year0"), by = c("id", "variable")) %>% 
  mutate(logitv = logit(value), 
         logitv0 = logit(value0))

par(mfrow = c(2, 4), mar = c(2, 2, 1, 1))
d_ply(dominentsp_year0, .(variable), 
      function(x) plot(logit(value) ~ logitv0, data = x, pch = 19, col = factor(year), main = "logit"))  
d_ply(dominentsp_year0, .(variable), 
      function(x) plot(value ~ value0, data = x, pch = 19, col = factor(year), main = "raw"))  

m_list <- dlply(dominentsp_year0, .(variable), function(x){
  m1 <- lmer(logit(value) ~ co2 * year + logitv0 + (1|block) + (1|ring) + (1|id), data = x)
  m2 <- update(m1, ~ . - (1|block))
  if (AICc(m1) >= AICc(m2)) return(m2) else return(m1)
})

# compute 95 CI and post-hoc test
lsmeans_list <- llply(m_list, function(x) {
  summary(lsmeans::lsmeans(x, pairwise ~ co2 | year))
})

# 95% CI
CI_dd <- ldply(lsmeans_list, function(x) data.frame(x$lsmeans)) 

# post-hoc test
contrast_dd <- ldply(lsmeans_list, function(x) data.frame(x$contrast)) %>% 
  mutate(co2 = factor("amb", levels = c("amb", "elev")),
         star = cut(p.value, right = FALSE,
                    breaks = c(0, .1, .05, .01, .001, 1),  
                    labels = c("***", "**", "*", "\u2020", ""))) %>% 
  select(.id, year, co2, p.value, star)

# merge
ci_dd <- left_join(CI_dd, contrast_dd, by = c(".id", "year", "co2")) %>% 
  mutate(year = factor(year, levels = paste0("Year", 0:3)))
ci_dd$star[is.na(ci_dd$star)] <- ""


# create fig --------------------------------------------------------------


# df for Year0
d <- dominentsp %>%
  filter(year == "Year0") %>% 
  mutate(year = factor(year, levels = paste0("Year", 0:3)),
         value = value/4, # standardise for 1mx1m plot
         .id = gsub("[.]", " ", variable))

# df for median of Year0
d_med <- d %>%  
  group_by(.id) %>% 
  summarise(M = median(value)) %>% 
  mutate(Med = "Md[Year0]")


ci_dd <- ci_dd %>% 
  mutate(rlsmean = boot::inv.logit(lsmean) * 25, # reverse transform and standardise for 1mx1m plot
         rlowerCL = boot::inv.logit(lower.CL) * 25,
         rupperCL = boot::inv.logit(upper.CL) * 25,
         .id = gsub("[.]", " ", .id),
         year = factor(year, levels = paste0("Year", 0:3)))

# fig
dodgeval <- .4
p <- ggplot(ci_dd, aes(x = year, y = rlsmean, fill = co2, group = co2))
p2 <- p +
  
  geom_errorbar(aes(ymin = rlowerCL, ymax = rupperCL), width = 0, 
                position = position_dodge(width = dodgeval)) +
  geom_line(aes(linetype = co2), position = position_dodge(width = dodgeval)) +
  geom_boxplot(data = d, aes(x = year, y = value),  alpha = .6, 
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
  scale_y_continuous(limits = c(0, 25)) +
  scale_x_discrete("Year", labels = 0:4, drop = FALSE) +
  
  science_theme +
  theme(legend.position = c(.87, .91),
        legend.margin = unit(-.5, "line"),
        strip.text.x = element_text(face = "italic")) +
  guides(linetype = guide_legend(order = 1),
         fill     = guide_legend(order = 1),
         col      = guide_legend(order = 2)) +
  
  
  facet_wrap( ~ .id) +
  labs(y = expression(Abundance~(Counts~m^'-1')))
p2

ggsavePP(filename = "output/figs/adjusted_abundance_dominenetSPP", 
         plot = p2, width = 5, height = 4)
