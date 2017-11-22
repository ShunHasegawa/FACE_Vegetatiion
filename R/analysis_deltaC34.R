# Here we will:
# 1) download moisture and understorey air temperature data from HIEv
# 2) identify growing season for each of C3 and C4 species according to Murphy 2007
# 3) Fit moisture and temperature to delta C3 and C4 (annual change in abundance)

# source("R/process_env.R")


# abundance change ------------------------------------------------------------


c34sum <- C3grassC4 %>%
  group_by(year, block, ring, plot, co2, id, PFG, RY, variable) %>% 
  summarise(value = sum(value)) %>%
  mutate(value = log(value + 1)) %>% 
  group_by(year, block, ring, plot, co2, id, PFG, RY) %>% 
  summarise(value = sum(value)) %>% 
  spread(key = PFG, value = value) %>% 
  ungroup() %>% 
  arrange(id, year) %>%
  group_by(id) %>%
  mutate_each(funs(ddiff = . - lag(., 1)),             # Year1-Year0 and etc.
                   c3, c4) %>%
  filter(year != "Year0") %>%
  group_by(year, ring, co2) %>% 
  summarise_each(funs(mean), ends_with("ddiff")) %>% 
  left_join(c34growth_moist) %>% 
  ungroup() %>% 
  mutate(s_c4_ddiff = scale(c4_ddiff)[, 1],
         s_c3_ddiff = scale(c3_ddiff)[, 1],
         s_logmoist = scale(log(totalmoist))[, 1],
         s_temp     = scale(annual_temp2m)[, 1],
         s_logpar   = scale(log(PAR))[, 1])

plot(s_c4_ddiff ~ log(totalmoist), data = c34sum, pch = 19, col = co2)
plot(s_c4_ddiff ~ annual_temp2m, data = c34sum, pch = 19, col = co2)
plot(s_c4_ddiff ~ log(PAR), data = c34sum, pch = 19, col = co2)




# c4 abundance ------------------------------------------------------------

c4d_m1      <- lmer(s_c4_ddiff ~ co2 * (s_logmoist+s_temp+s_logpar) + (1|ring), data = c34sum)
c4d_m1_full <- dredge(c4d_m1, REML = F, fixed = c("co2", "s_logmoist", "s_temp",  "s_logpar"))
plot(c4d_m1)
qqPlot(resid(c4d_m1))
c4d_m2      <- lmer(s_c4_ddiff ~ co2 + s_logmoist + s_temp + s_logpar + (1|ring), data = c34sum)
summary(c4d_m2)
Anova(c4d_m2, test.statistic = "F")
VarCorr(c4d_m2)
c4d_m3      <- lm(s_c4_ddiff ~ co2 + s_logmoist + s_temp + s_logpar, data = c34sum)
c4_coef    <- confint(c4d_m3)
c4_coef_90 <- confint(c4d_m3, level = .9)
c4d_m2_full <- dredge(c4d_m2, REML = F)
c4_coef_imp <- importance(c4d_m2_full)
c4_coef_imp

# partial residual plot

pdf(file = "output/figs/lar_C4_partial_regression_plot.pdf", width = 5, height = 5)
create_resplot(c4d_m2,
               ylab = expression(Adj.~annual~rates~of~change~"in"~C[4]),
               ylim = c(-1.5, 2.5))
dev.off()


png("output/figs/lar_C4_partial_regression_plot.png", width = 5, height = 5, res = 600, units = "in")
create_resplot(c4d_m2,
               ylab = expression(Adj.~annual~rates~of~change~"in"~C[4]),
               ylim = c(-1.5, 2.5))
dev.off()




# C3 abundance ------------------------------------------------------------


c3d_m1      <- lmer(s_c3_ddiff ~ co2 * (s_logmoist+s_temp+s_logpar) + (1|ring), data = c34sum)
c3d_m1_full <- dredge(c3d_m1, REML = F, fixed = c("co2", "s_logmoist", "s_temp",  "s_logpar"))
plot(c3d_m1)
qqPlot(resid(c3d_m1))
c3d_m2      <- lmer(s_c3_ddiff ~ co2 + s_logmoist + s_temp + s_logpar + (1|ring), data = c34sum)
summary(c3d_m2)
Anova(c3d_m2, test.statistic = "F")
VarCorr(c3d_m2)
plot(c3d_m2)
qqPlot(resid(c3d_m2))
mcp.fnc(c3d_m2)
# one outlier was suggested

# remove outlier
olrm_df <- romr.fnc(c3d_m2, data.frame(c34sum))
dplyr::setdiff(olrm_df$data0, olrm_df$data)
c3d_m3 <- update(c3d_m2, data = olrm_df$data)
plot(c3d_m3)
qqPlot(resid(c3d_m3))

VarCorr(c3d_m2)
VarCorr(c3d_m3)

confint(c3d_m2, method = "boot", nsim = 99, level = .9)
confint(c3d_m3, method = "boot", nsim = 99, level = .9)
# no major difference, so use the original one


c3_coef    <- confint(c3d_m2, method = "boot", nsim = 999)
c3_coef_90 <- confint(c3d_m2, method = "boot", nsim = 999, level = .9)
c3d_m2_full <- dredge(c3d_m2, REML = F)
c3_coef_imp <- importance(c3d_m2_full)
c3_coef_imp



# > partial residual plot ---------------------------------------------------

pdf(file = "output/figs/lar_C3_partial_regression_plot.pdf", width = 5, height = 5)
create_resplot(c3d_m2,
               ylab = expression(Adj.~annual~rates~of~change~"in"~C[3]),
               ylim = c(-2.5, 3))
dev.off()


png("output/figs/lar_C3_partial_regression_plot.pdf.png", width = 5, height = 5, res = 600, units = "in")
create_resplot(c3d_m2,
               ylab = expression(Adj.~annual~rates~of~change~"in"~C[3]),
               ylim = c(-2.5, 3))
dev.off()
