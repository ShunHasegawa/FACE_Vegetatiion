
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
sd_m1 <- lmer(sqrt(sd_ratio) ~ co2 * year + sqrt(sd_ratio0) + (1|ring) + (1|id) + (1|RY), data = DS_ratio_ed)
Anova(sd_m1, test.statistic = "F")
summary(sd_m1)
plot(sd_m1)
qqPlot(resid(sd_m1))
r.squared(sd_m1)
visreg(sd_m1, xvar = "sd_ratio0", by = "co2", overlay = TRUE, xtrans = sqrt)

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



# . fig -------------------------------------------------------------------

# CI and postdoc test -----------------------------------------------------


# compute 95 CI
sd_lsmeans <- lsmeans::lsmeans(sd_m1, ~ co2 | year)


# CO2 effect
sd_co2_pval <- tidy(Anova(sd_m1, test.statistic = "F")) %>%
  filter(term == "co2") %>% 
  mutate(co2star = get_star(p.value)) %>% 
  select(co2star)


sd_CI_dd <- data.frame(summary(sd_lsmeans, type = "response")) %>% 
  mutate(year       = factor(year, levels = paste0("Year", 0:3)),
         co2star    = sd_co2_pval[,1],
         value_type = "adjusted")


# df for observed values
sd_obs_d <- DS_ratio %>%
  group_by(year, co2, ring) %>%
  summarise(value = mean(sd_ratio)) %>%
  ungroup() %>%
  mutate(year = factor(year, levels = paste0("Year", 0:3)),
         value_type = "observed")


# co2 response ratio
sd_rr_d <- sd_CI_dd %>% 
  group_by(co2, co2star) %>% 
  summarise(value = mean(response)) %>% 
  group_by(co2star) %>% 
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1) %>% 
  mutate(rr = ifelse(rr >= 0,
                     paste0("RR= +", format(rr, digits = 0, nsmall = 2), co2star), 
                     paste0("RR= ", format(rr, digits = 0, nsmall = 2), co2star)))




# create a figure ---------------------------------------------------------

# create fig
dodgeval <- .3


fig_sd <- ggplot(sd_CI_dd, aes(x = year, y = response)) +
  
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  
  # observed values
  geom_point(data = sd_obs_d, 
             aes(x = year, y = value, shape = co2, group = co2, col = value_type),  
             fill = "grey80", size = 2,
             position = position_dodge(dodgeval)) +
  
  # adjusted values
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, 
                    shape = co2, group = co2, col = value_type), 
                width = 0,
                position = position_dodge(width = dodgeval)) +
  geom_line(aes(linetype = co2, shape = co2, group = co2, col = value_type), 
            position = position_dodge(width = dodgeval)) +
  geom_point(aes(shape = co2, group = co2, col = value_type), 
             size = 2.5, position = position_dodge(width = dodgeval)) +
  
  # scaling
  scale_shape_manual(values = c(16, 15),
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
  
  labs(y = "Subordinate:Dominant ratios") +
  geom_text(data = sd_rr_d, 
            aes(label = rr), x =  Inf, y = Inf, 
            hjust = 1.1, vjust = 1.5, size = 3)


ggsavePP(filename = "output/figs/adjusted_sd_ratio", width = 3, height = 3,
         plot = fig_sd)




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

# dominant c3
dc3_m1 <- lmer(value ~ co2 * year + value0 + (1|ring) + (1|id) + (1|RY), dc3_dd)
Anova(dc3_m1, test.statistic = "F")
plot(dc3_m1)
qqPlot(resid(dc3_m1))
dc3_m2 <- lmer(value ~ co2 * year + value0 + (1|ring) + (1|id) + (1|RY), dc3_dd[-which.min(resid(dc3_m1)), ])
qqPlot(resid(dc3_m2))
Anova(dc3_m2, test.statistic = "F")


# dominant c4
dc4_m1 <- lmer(sqrt(value + 1) ~ co2 * year + sqrt(value0 + 1) + (1|ring) + (1|id) + (1|RY), dc4_dd)
Anova(dc4_m1, test.statistic = "F")
summary(dc4_m1)
plot(dc4_m1)
qqPlot(resid(dc4_m1))


# subordinate c3
boxplot(log(value + 1) ~ co2 * year, data = sc3_dd)
sc3_m1 <- glmer(value ~ co2 * year + I(value0/100) +  (1|id) + (1|RY), 
                family = "poisson", sc3_dd)
summary(sc3_m1)
plot(sc3_m1)
qqPlot(resid(sc3_m1))

plot(log(value + 1)  ~ log(value0 + 1), data = sc3_dd, col = co2, pch = 19)
sc3_m1 <- lmer(log(value + 1) ~ co2 * year + log(value0 + 1) + (1|ring) + (1|id) + (1|RY), sc3_dd)
Anova(sc3_m1, test.statistic = "F")
plot(sc3_m1)
qqPlot(resid(sc3_m1))


# subordinate c4
sc4_m1 <- lmer(sqrt(value + 1) ~ co2 * year + sqrt(value0 + 1) + (1|ring) + (1|id) + (1|RY), sc4_dd)
qqPlot(resid(sc4_m1))
plot(sc4_m1)

Anova(sc4_m1, test.statistic = "F")
summary(sc4_m1)


# . fig -------------------------------------------------------------------


# CI and post-hoc test ----------------------------------------------------


# compute 95 CI and post-hoc test (only need grass and forb species)

sd_m_list <- list(Dominant_C3    = dc3_m1, 
                  Dominant_C4    = dc4_m1, 
                  Subordinate_C3 = sc3_m1, 
                  Subordinate_C4 = sc4_m1)
nl <- names(sd_m_list)

# defint transformation
sapply(lsmeans_list, str)
sd_m_list_tr <- sd_m_list
sd_m_list_tr[[2]] <- update(ref.grid(sd_m_list_tr[[2]]), tran = make.tran("power", c(.5, -1))) # (y + 1)^.5
sd_m_list_tr[[3]] <- update(ref.grid(sd_m_list_tr[[3]]), tran = make.tran("genlog", 1))        # log(y + 1)
sd_m_list_tr[[4]] <- update(ref.grid(sd_m_list_tr[[4]]), tran = make.tran("power", c(.5, -1))) # (y + 1)^.5



lsmeans_list <- llply(sd_m_list_tr, function(x) {
  lsmeans::lsmeans(x, ~ co2 | year)
})

# 95% CI
CI_dd <- ldply(lsmeans_list, function(x) data.frame(summary(x, type = "response"))) %>% 
  mutate(lsmean = ifelse(is.na(lsmean), response, lsmean)) %>%
  select(-response)


# post-hoc test
contrast_dd <- ldply(lsmeans_list, function(x) {
  data.frame(summary(pairs(x)[1:3], adjust = "none"))
}) %>% 
  mutate(co2 = factor("amb", levels = c("amb", "elev")),
         star = get_star(p.value)) %>% 
  select(.id, year, co2, p.value, star)


# CO2 effect
sd_aov_df <- ldply(sd_m_list, function(x) tidy(Anova(x, test.statistic = "F")),  # Anova result of models
                   .progress = "text")
sd_co2_pval <- sd_aov_df %>%                                                     # get p-values for CO2 term
  filter(term == "co2") %>% 
  mutate(co2star = get_star(p.value)) %>% 
  select(.id, co2star)


# merge
ci_dd <- left_join(CI_dd, contrast_dd, by = c(".id", "year", "co2")) %>% 
  mutate(type       = tstrsplit(.id, "_")[[1]], 
         PFG        = tstrsplit(.id, "_")[[2]],
         year       = factor(year, levels = paste0("Year", 0:3)),
         value_type = "adjusted",
         plot_lab   = factor(.id, labels = paste0("(", letters[1:4], ")"))) %>%  # sub-plot label
  left_join(sd_co2_pval, by = ".id")                                             # merge with pvalues for CO2 term

ci_dd$star[is.na(ci_dd$star)]           <- ""  # turn NA to ""
ci_dd$star[ci_dd$.id == "Dominant_C4"]  <- ""  # no CO2xTime interaction, so don't show post-hoc results


# Observed vlaues for each variable
sd_obs_dd <- grass_DS_sum %>%
  group_by(pfg_type, PFG, type, year, co2, ring) %>% 
  summarise(value = mean(value)) %>% 
  ungroup() %>% 
  mutate(value_type = "observed",
         type = mapvalues(type, c("D", "S"), c("Dominant", "Subordinate")),
         PFG  = mapvalues(PFG, c("c3", "c4"), c("C3", "C4")))


# plot label
plab_d <- ci_dd %>% 
  group_by(type, PFG, plot_lab, co2, co2star) %>% 
  summarise(value = mean(lsmean)) %>% 
  group_by(type, PFG, plot_lab, co2star) %>% 
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1) %>% 
  mutate(rr = ifelse(rr >= 0,
                     paste0("RR= +", format(rr, digits = 0, nsmall = 2), co2star), 
                     paste0("RR= ", format(rr, digits = 0, nsmall = 2), co2star)))


# generate a fugure
ci_dd$PFG     <- factor(ci_dd$PFG, labels = c("C[3]", "C[4]"))
sd_obs_dd$PFG <- factor(sd_obs_dd$PFG, labels = c("C[3]", "C[4]"))
plab_d$PFG    <- factor(plab_d$PFG, labels = c("C[3]", "C[4]"))
dodgeval <- .4
sd_abund_fig <- ggplot(ci_dd, aes(x = year, y = lsmean)) +
  
  facet_grid(type ~ PFG, labeller = label_parsed) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  
  
  # observed
  geom_point(data = sd_obs_dd, 
             aes(x = year, y = value, shape = co2, col = value_type), 
             size = 2, fill = "grey80", 
             position = position_dodge(dodgeval)) +
  
  
  # adjusted
  geom_line(aes(linetype = co2, group = co2), 
            position = position_dodge(width = dodgeval)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, group = co2), width = 0, 
                position = position_dodge(width = dodgeval)) +
  geom_point(aes(shape = co2, col = value_type), size = 2.5, position = position_dodge(width = dodgeval)) +
  geom_text(aes(y = upper.CL, label = star), fontface = "bold", vjust = -.1) +
  
  
  # scaaling
  scale_shape_manual(values = c(16, 15), 
                     labels = c("Ambient", expression(eCO[2]))) +
  scale_linetype_manual(values = c("solid", "dashed"), 
                        labels = c("Ambient", expression(eCO[2]))) +
  scale_color_manual(values = c("black", "grey80"),
                     guide = guide_legend(override.aes = list(linetype = "blank",size = 2))) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20)) +
  scale_x_discrete("Year", labels = 0:4, drop = FALSE) +
  
  
  # legend and theme
  science_theme +
  theme(legend.position   = "top",
        legend.box        = "horizontal", 
        legend.direction  = "vertical", 
        legend.text.align = 0) +
  
  labs(y = expression(Abundance~(Counts~plot^'-1'))) +
  
  geom_text(data = plab_d, aes(label = plot_lab), x = -Inf, y = Inf, 
            hjust = -.1, vjust = 1.5, size = 3, fontface = "bold") +
  geom_text(data = plab_d, aes(label = rr), x =  Inf, y = Inf, hjust = 1.1, 
            vjust = 1.5, size = 3)

sd_abund_fig

ggsavePP(filename = "output/figs/adjusted_SD_C34_abund", 
         plot = sd_abund_fig, width = 4.5, height  = 4.5)



# summary table -----------------------------------------------------------

# table for observed values
obs_tbl <- sd_obs_dd %>% 
  select(PFG, type, year, co2, value_type, ring, value) %>% 
  group_by(PFG, type, year, co2, value_type) %>% 
  summarise(M = mean(value, na.rm = TRUE))


# bind with adjusted values
sd_adjMean_tble <- ci_dd %>%
  select(PFG, type, co2, year, lsmean, value_type) %>% 
  rename(M = lsmean) %>% 
  bind_rows(obs_tbl) %>% 
  mutate(variable = paste(value_type, co2, sep = "_")) %>% 
  select(-value_type, -co2) %>% 
  spread(key = variable, value = M) %>% 
  mutate(resp     = adjusted_elev / adjusted_amb - 1) %>% 
  group_by(type, PFG, year) %>% 
  summarise_each(funs(round(., 2)), everything(), -PFG, -year, -type) %>% 
  select(PFG, type, year, starts_with("observed"), starts_with("adjusted"), 
         resp) %>% 
  ungroup() %>% 
  arrange(PFG, type)


# save as csv
write.table(sd_adjMean_tble, file = "output/table/summary_tbl_dom_sub_c34_abund.csv", 
            row.names = FALSE)
