# figure ------------------------------------------------------------------

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


# modify labels for facet_wrap subplots
facet_labels <- c(Grass       = "Graminoid",
                  C3vsC4      = "C[4]~graminoid",
                  LegvsNonleg = "Legume",
                  NatvsIntr   = "Native~species")

all_pfg_prop_ci_dd <- mutate(all_pfg_prop_ci_dd,
                             .id        = factor(all_pfg_prop_ci_dd$.id, 
                                                 levels = c("Grass", "C3vsC4","LegvsNonleg", "NatvsIntr"),
                                                 labels = facet_labels),
                             plot_lab   = as.character(factor(.id, labels = c("(a)", "", "(b)", "(c)"))))  # sub-plot label


all_pfg_d$.id          <- factor(all_pfg_d$.id, 
                                 levels = c("Grass", "C3vsC4", "LegvsNonleg", "NatvsIntr"),
                                 labels = facet_labels)


# df for plot labels and response ratios
all_pfg_plab_d <- all_pfg_prop_ci_dd %>% 
  group_by(.id, plot_lab, co2, co2star) %>% 
  summarise(value = mean(rlsmean)) %>% 
  group_by(.id, plot_lab, co2star) %>% 
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1) %>% 
  mutate(rr = ifelse(rr >= 0,
                     paste0("RR= +", format(rr, digits = 0, nsmall = 2), co2star), 
                     paste0("RR= ", format(rr, digits = 0, nsmall = 2), co2star)))



# create fig
dodgeval <- .3




# fig for C4 proportion -------------------------------------------------


fig_c4prop <- ggplot(filter(all_pfg_prop_ci_dd, .id == "C[4]~graminoid"), 
                     aes(x = year, y = rlsmean)) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  
  
  # observed values
  geom_point(data = filter(all_pfg_d, .id == "C[4]~graminoid"), 
             aes(x = year, y = ratios, shape = co2, group = co2, col = value_type),  
             fill = "grey80", size = 2,
             position = position_dodge(dodgeval)) +
  
  
  # adjusted values
  geom_errorbar(aes(ymin = rlowerCL, ymax = rupperCL, 
                    shape = co2, group = co2, col = value_type), 
                width = 0,
                position = position_dodge(width = dodgeval)) +
  geom_line(aes(linetype = co2, shape = co2, group = co2, col = value_type), 
            position = position_dodge(width = dodgeval)) +
  geom_point(aes(shape = co2, group = co2, col = value_type), 
             size = 2.5, position = position_dodge(width = dodgeval)) +
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
  
  labs(y = expression(C[4]~proportion)) +
  geom_text(data = filter(all_pfg_plab_d, .id == "C[4]~graminoid"), 
            aes(label = rr), x =  Inf, y = Inf, 
            hjust = 1.1, vjust = 1.5, size = 3) +

  geom_blank(data = data.frame(rlsmean = 0.8, year = "Year0"),
             aes(x = year, y = rlsmean))
  
fig_c4prop

ggsavePP(filename = "output/figs/adjusted_c4_proportion", width = 3, height = 3,
         plot = fig_c4prop)




# fig for other PFGs ------------------------------------------------------


fig_PFGprop <- ggplot(filter(all_pfg_prop_ci_dd, .id != "C[4]~graminoid"), 
                      aes(x = year, y = rlsmean)) +
  facet_wrap(~ .id, scales = "free_y", ncol = 2, labeller = label_parsed) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  
  
  # observed values
  geom_point(data = filter(all_pfg_d, .id != "C[4]~graminoid"), 
             aes(x = year, y = ratios, shape = co2, group = co2, col = value_type),  
             fill = "grey80", size = 2,
             position = position_dodge(dodgeval)) +
  
  
  # adjusted values
  geom_errorbar(aes(ymin = rlowerCL, ymax = rupperCL, 
                    shape = co2, group = co2, col = value_type), 
                width = 0,
                position = position_dodge(width = dodgeval)) +
  geom_line(aes(linetype = co2, shape = co2, group = co2, col = value_type), 
            position = position_dodge(width = dodgeval)) +
  geom_point(aes(shape = co2, group = co2, col = value_type), 
             size = 2.5, position = position_dodge(width = dodgeval)) +
  geom_text(aes(label = star, y = rupperCL), fontface = "bold", vjust = .4) +
  geom_blank(data = all_pfg_d %>%
               filter(.id != "C[4]~graminoid") %>% 
               select(.id) %>%
               distinct() %>%
               mutate(rlsmean = c(NA, NA, 1.03), year = "Year0"),
             aes(x = year, y = rlsmean)) +
  
  
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
  
  labs(y = expression("Proportion")) +
  geom_text(data = filter(all_pfg_plab_d, .id != "C[4]~graminoid"), 
            aes(label = plot_lab), x = -Inf, y = Inf, 
            hjust = -.1, vjust = 1.5, size = 3, fontface = "bold") +
  geom_text(data = filter(all_pfg_plab_d, .id != "C[4]~graminoid"), 
            aes(label = rr), x =  Inf, y = Inf, 
            hjust = 1.1, vjust = 1.5, size = 3)

ggsavePP(filename = "output/figs/adjusted_PFG_proportion", width = 4, height = 4.5,
         plot = fig_PFGprop)





# summary table -----------------------------------------------------------


# . table for observed values ---------------------------------------------

# treatment mean

# grass prop
grassprop_co2_d <- grassprop_d %>% 
  group_by(year, co2) %>% 
  summarise(M  = sum(Grass) / (sum(Grass) + sum(Forb))) %>% 
  mutate(.id = "Graminoid")

# other pfg prop
pfgprop_co2_d <- pfgprop_d %>% 
  mutate(pfg_value = Total * ratios) %>%          # get pfg abundance before getting proportion
  group_by(.id, year, co2) %>% 
  summarise(M = sum(pfg_value) / sum(Total))

# table
obs_tbl <- bind_rows(grassprop_co2_d, pfgprop_co2_d) %>%
  mutate(value_type = "observed",
         .id = mapvalues(.id, c("C3vsC4", "LegvsNonleg", "NatvsIntr"),
                         c("C4graminoid", "Legume", "Native")))




# . bind with adjusted values ---------------------------------------------

pfgprop_adjMean_tble <- all_pfg_prop_ci_dd %>% 
  mutate(.id = dplyr::recode(.id, "C[4]~graminoid" = "C4graminoid",
         "Native~species" = "Native")) %>%           # change .id labels 
  select(.id, co2, year, rlsmean, value_type) %>% 
  rename(M = rlsmean) %>% 
  bind_rows(obs_tbl) %>% 
  mutate(variable = paste(value_type, co2, sep = "_")) %>% 
  select(-value_type, -co2) %>% 
  spread(key  = variable, value = M) %>% 
  mutate(resp = adjusted_elev / adjusted_amb - 1) %>%                             # response ratios
  group_by(.id, year) %>% 
  summarise_each(funs(round(., 2)), everything(), -.id) %>% 
  select(year, starts_with("observed"), starts_with("adjusted"), resp) %>% 
  ungroup()


# split df by .id
pfgprop_tbl_l <- dlply(pfgprop_adjMean_tble, .(.id), function(x) select(x, -.id))


# save as excel
writeWorksheetToFile(file        = "output/table/summary_tbl_pfg_prop.xlsx",  # define file name to be saved
                     data        = pfgprop_tbl_l,                             # writeWorksheetToFile doesn't take dplyr object so turn them into data frames using as.data.frame
                     sheet       = names(pfgprop_tbl_l),                      # sheet names in excel are defined by object names a list
                     clearSheets = TRUE)
