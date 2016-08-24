# merge data frames for grass and pfg proportions

# 95% CI
all_pfg_prop_ci_dd <- grassprop_ci_dd %>% 
  rename(.id = form) %>% 
  bind_rows(pfgprop_ci_dd)

# Year0
all_pfg_d <- grassprop_d %>%
  rename(ratios = grass_prop) %>% 
  mutate(.id = "Grass") %>% 
  bind_rows(pfgprop_d) %>% 
  select(.id, year, co2, ratios)

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
dodgeval <- .4
fig_pfgprop <- ggplot(all_pfg_prop_ci_dd, aes(x = year, y = rlsmean, fill = co2, group = co2)) +
  facet_wrap(~ .id, scales = "free_y", ncol = 2, labeller = label_parsed) +

  geom_errorbar(aes(ymin = rlowerCL, ymax = rupperCL), width = 0,
                position = position_dodge(width = dodgeval)) +
  geom_line(aes(linetype = co2), position = position_dodge(width = dodgeval)) +
  geom_point(data = all_pfg_d, aes(x = year, y = ratios),  alpha = .7, shape = 21, size = 3,
             position = position_dodge(dodgeval), show.legend = FALSE) +
  geom_point(shape = 21, size = 3, position = position_dodge(width = dodgeval)) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_text(aes(label = star, y = rupperCL), fontface = "bold", vjust = .4) +

  scale_fill_manual(values = c("black", "white"),
                    labels = c("Ambient", expression(eCO[2]))) +
  scale_linetype_manual(values = c("solid", "dashed"),
                        labels = c("Ambient", expression(eCO[2]))) +
  scale_color_manual(values = "grey50", labels = expression(Md[Year0])) +
  scale_x_discrete("Year", labels = 0:4, drop = FALSE) +

  science_theme +
  theme(legend.position = c(.3, .95),
        strip.text.x    = element_text(size = 8)) +

  labs(y = expression("Adjusted proportion"))
fig_pfgprop

ggsavePP(filename = "output/figs/adjusted_PFG_proportion", width = 6, height = 5,
         plot = fig_pfgprop)
