# merge data frames for grass and pfg proportions

# 95% CI
all_pfg_prop_ci_dd <- grassprop_ci_dd %>% 
  rename(.id = form) %>% 
  bind_rows(pfgprop_ci_dd) %>% 
  mutate(value_type = "adjusted")

# observed values
all_pfg_d <- grassprop_d %>%
  rename(ratios = grass_prop) %>% 
  mutate(.id = "Grass") %>% 
  bind_rows(pfgprop_d) %>% 
  select(.id, year, co2, ratios) %>% 
  mutate(value_type = "observed")

# create a fig

# modify labels for facet_wrap subplots
facet_labels <- c(Grass       = "Grass~(Grass~vs.~Forb)",
                  C3vsC4      = "C[3*'_'*grass]~(C[3*'_'*grass]~vs.~C[4*'_'*grass])",
                  LegvsNonleg = "Legume~(Legume~vs.~Non*-legume)",
                  NatvsIntr   = "Native~plant~(Native~vs.~Introduced)")

all_pfg_prop_ci_dd$.id <- factor(all_pfg_prop_ci_dd$.id, 
                                 levels = c("Grass", "C3vsC4", "LegvsNonleg", "NatvsIntr"),
                                 labels = facet_labels)
all_pfg_d$.id <- factor(all_pfg_d$.id, 
                        levels = c("Grass", "C3vsC4", "LegvsNonleg", "NatvsIntr"),
                        labels = facet_labels)

# create fig
dodgeval <- .3
fig_pfgprop <- ggplot(all_pfg_prop_ci_dd, 
                      aes(x = year, y = rlsmean, shape = co2, group = co2, 
                          col = value_type)) +
  facet_wrap(~ .id, scales = "free_y", ncol = 2, labeller = label_parsed) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +

  
  # observed values
  geom_point(data = all_pfg_d, aes(x = year, y = ratios),  fill = "grey80", size = 2,
             position = position_dodge(dodgeval)) +
  
  
  # adjusted values
  geom_errorbar(aes(ymin = rlowerCL, ymax = rupperCL), width = 0,
                position = position_dodge(width = dodgeval)) +
  geom_line(aes(linetype = co2), position = position_dodge(width = dodgeval)) +
  geom_point(size = 2.5, position = position_dodge(width = dodgeval)) +
  geom_text(aes(label = star, y = rupperCL), fontface = "bold", vjust = .4) +

  
  # scaling
  scale_shape_manual(values = c(16, 17),
                    labels = c("Ambient", expression(eCO[2]))) +
  scale_linetype_manual(values = c("solid", "dashed"),
                        labels = c("Ambient", expression(eCO[2]))) +
  scale_color_manual(values = c("black", "grey80"),
                     guide  = guide_legend(override.aes = list(linetype = "blank",
                                                               size     = 2))) +
  scale_x_discrete("Year", labels = 0:4, drop = FALSE) +

  
  # legend and theme
  science_theme +
  theme(legend.position   = "top",
        legend.box        = "horizontal", 
        legend.direction  = "vertical", 
        legend.text.align = 0,
        strip.text.x      = element_text(size = 8)) +

  labs(y = expression("Proportion"))
fig_pfgprop

ggsavePP(filename = "output/figs/adjusted_PFG_proportion", width = 5, height = 5.5,
         plot = fig_pfgprop)
