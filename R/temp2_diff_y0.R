diff_fh_m1  <- lmer(I(value - value0)       ~ co2 * year + (1|ring), data = fh_d)
diff_fj_m1  <- lmer(I(value - value0)       ~ co2 * year + (1|ring), data = fj_d)
diff_fs_m1  <- lmer(I(value - value0)       ~ co2 * year + (1|ring), data = fs_d)
diff_gh_m1  <- lmer(I(value - value0)       ~ co2 * year + (1|ring), data = gh_d)
diff_gj_m1  <- lmer(I(value - value0)       ~ co2 * year + (1|ring), data = gj_d)
diff_dc3_m1 <- lmer(I(value - value0)       ~ co2 * year + (1|ring), dc3_dd)
diff_gs_m1  <- lmer(I(value - value0)       ~ co2 * year + (1|ring), data = gs_d)
diff_dc4_m1 <- lmer(I(value - value0)       ~ co2 * year + (1|ring), dc4_dd)
diff_sc3_m1 <- lmer(I(value - value0)       ~ co2 * year + (1|ring), sc3_dd)
diff_sc4_m1 <- lmer(I(value - value0)       ~ co2 * year + (1|ring), sc4_dd)
diff_c43_m1 <- lmer(I(c43_r - ratios0)      ~ co2 * year + (1|ring), data = c43_ratio_year0)
diff_sd_m1  <- lmer(I(sd_ratio - sd_ratio0) ~ co2 * year + (1|ring), data = DS_ratio_ed)
summary(diff_sd_m1)

diffy0_m_list <- list('forb_H'         = diff_fh_m1,
                      'forb_J'         = diff_fj_m1,
                      'forb_S'         = diff_fs_m1,
                      'grass_H'        = diff_gh_m1,
                      'grass_J'        = diff_gj_m1,
                      'grass_S'        = diff_gs_m1,
                      'dominant_c3'    = diff_dc3_m1,
                      'dominant_c4'    = diff_dc4_m1,
                      'subordinate_c3' = diff_sc3_m1,
                      'subordinate_c4' = diff_sc4_m1,
                      'c43_ratio'      = diff_c43_m1,
                      'sd_ratio'       = diff_sd_m1)
diffy0_aov <- ldply(diffy0_m_list, function(x) tidy(Anova(x, test.statistic = "F")))
diffy0_aov_tbl <- diffy0_aov %>% 
  mutate(p.value = round(p.value, 3),
         Fval = paste0("F(", df, ",", round(Df.res, 0), ")=",round(statistic, 2))) %>% 
  select(.id, term, Fval, p.value) %>% 
  gather(variable, value, Fval, p.value) %>% 
  mutate(CF = paste(term, variable, sep = "_")) %>% 
  select(.id, CF, value) %>% 
  spread(CF, value) %>% 
  select(.id, starts_with("co2_"), starts_with("year_"), starts_with("co2:year")) %>% 
  mutate(.id = factor(.id, levels = c("grass_H", "grass_J","grass_S", "forb_H", "forb_J", "forb_S",
                                      "dominant_c3", "subordinate_c3", "dominant_c4", "subordinate_c4",
                                      "c43_ratio", "sd_ratio"))) %>% 
  arrange(.id)
write.csv(diffy0_aov_tbl, "output/table/diff_y0_analysis_anova.csv", row.names = FALSE)












par(mfrow = c(3, 4))
boxplot(I(value - value0)       ~ co2 * year, data = fh_d, main = "forbH")
boxplot(I(value - value0)       ~ co2 * year, data = fj_d, main = "forbJ")
boxplot(I(value - value0)       ~ co2 * year, data = fs_d, main = "forbS")
boxplot(I(value - value0)       ~ co2 * year, data = gh_d, main = "grassH")
boxplot(I(value - value0)       ~ co2 * year, data = gj_d, main = "grassJ")
boxplot(I(value - value0)       ~ co2 * year, data = gs_d, main = "grassS")
boxplot(I(value - value0)       ~ co2 * year, dc3_dd, main = "dc3")
boxplot(I(value - value0)       ~ co2 * year, dc4_dd, main = "dc4")
boxplot(I(value - value0)       ~ co2 * year, sc3_dd, main = "sc3")
boxplot(I(value - value0)       ~ co2 * year, sc4_dd, main = "sc4")


boxplot(log(I(c43_r - ratios0) + .6) ~ co2 * year, data = c43_ratio_year0, main = "c43")
diff_c43_m1 <- lmer(log(I(c43_r - ratios0) + .6) ~ co2 * year + (1|ring), data = c43_ratio_year0)
diff_c43_m1 <- lmer(log(I(c43_r - ratios0) + .6) ~ co2 * year + (1|ring), data = c43_ratio_year0)
Anova(diff_c43_m1, test.statistic = "F")
plot(diff_c43_m1)
qqPlot(resid(diff_c43_m1))



boxplot(sqrt(I(c43_r - ratios0) + 0.5)      ~ co2 * year, data = c43_ratio_year0, main = "c43")
boxplot(I(c43_r - ratios0)      ~ co2 * year, data = c43_ratio_year0, main = "c43")
boxplot(log(I(sd_ratio - sd_ratio0) + .9) ~ co2 * year, data = DS_ratio_ed, main = "sd")
boxplot(I(sd_ratio - sd_ratio0) ~ co2 * year, data = DS_ratio_ed, main = "sd")
