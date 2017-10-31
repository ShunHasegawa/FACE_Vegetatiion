
# Here, I analyse LAR of dominant/subordinate C3/C4 species against
# environmental variables



# prepare data frame ------------------------------------------------------

summary(c34growth_moist)


# compute annual change rate (acr) for each group
grass_DS_acr <- grass_DS %>%
  mutate(type_pfg = paste(type, PFG, sep = "_")) %>% 
  group_by(year, block, ring, co2, plot, id, RY, type_pfg, variable) %>% 
  summarise(value = sum(value)) %>%
  mutate(value = log(value + 1)) %>% 
  group_by(year, block, ring, co2, plot, id, RY, type_pfg) %>% 
  summarise(value = sum(value)) %>% 
  spread(key = type_pfg, value = value) %>% 
  ungroup() %>% 
  arrange(id, year) %>%
  group_by(id) %>%
  mutate_each(funs(ddiff = . - lag(., 1)),
              D_c3, D_c4, S_c3, S_c4) %>%
  filter(year != "Year0") %>%
  group_by(year, ring, co2) %>%
  summarise_each(funs(mean), ends_with("_ddiff")) %>%
  left_join(c34growth_moist) %>% 
  ungroup() %>% 
  mutate(s_dc4_ddiff = scale(D_c4_ddiff)[, 1],
         s_dc3_ddiff = scale(D_c3_ddiff)[, 1],
         s_sc4_ddiff = scale(S_c4_ddiff)[, 1],
         s_sc3_ddiff = scale(S_c3_ddiff)[, 1],
         s_logmoist  = scale(log(totalmoist))[, 1],
         s_temp      = scale(annual_temp2m)[, 1],
         s_logpar    = scale(log(PAR))[, 1])


# subordinate c4  -------------------------------------------------------------

lar_sc4_m1 <- lmer(s_sc4_ddiff ~ co2 * (s_logmoist+s_temp + s_logpar) +(1|ring), data = grass_DS_acr)
lar_sc4_m1_full <- dredge(lar_sc4_m1, fixed = c("co2", "s_logmoist", "s_temp",  "s_logpar"), REML = F)
lar_sc4_m1_full
plot(lar_sc4_m1)
qqPlot(resid(lar_sc4_m1))
# no interaction is suggested
lar_sc4_m2 <- lmer(s_sc4_ddiff ~ co2 + s_logmoist + s_temp + s_logpar + (1|ring), data = grass_DS_acr)
Anova(lar_sc4_m2, test.statistic = "F")
plot(lar_sc4_m2)
qqPlot(resid(lar_sc4_m2))
summary(lar_sc4_m2)

# coefficients and RI

VarCorr(lar_sc4_m2) # no variaation was explained by ring, so remove to avaoid convergence issues

lar_sc4_m3 <- lm(s_sc4_ddiff ~ co2 + s_logmoist+s_temp + s_logpar, data = grass_DS_acr)
sc4_coef       <- confint(lar_sc4_m3)
sc4_coef_90    <- confint(lar_sc4_m3, level = .9)
save(sc4_coef, sc4_coef_90, file = "output/Data/coef_sc4_env.RData")
load("output/Data/coef_sc4_env.RData")
lar_sc4_m2_full <- dredge(lar_sc4_m2, REML = F)
sc4_coef_imp    <- importance(lar_sc4_m2_full)


# partial residual plot
pdf(file = "output/figs/lar_sC4_partial_regression_plot.pdf", width = 5, height = 5)
create_resplot(lar_sc4_m2, 
               ylab = expression(Adj.~annual~rates~of~change~"in"~subordinate~C[4]),
               ylim <- c(-2, 2.6))

dev.off()

png(file = "output/figs/lar_sC4_partial_regression_plot.png", width = 5, height = 5, res = 600, units = "in")
create_resplot(lar_sc4_m2, 
               ylab = expression(Adj.~annual~rates~of~change~"in"~subordinate~C[4]),
               ylim <- c(-2, 2.6))

dev.off()




# dominant c4 -------------------------------------------------------------
lar_dc4_m1      <- lmer(s_dc4_ddiff ~ co2 * (s_logmoist+s_temp + s_logpar) +(1|ring), data = grass_DS_acr)
lar_dc4_m1_full <- dredge(lar_dc4_m1, fixed = c("co2", "s_logmoist", "s_temp",  "s_logpar"), REML = F)
lar_dc4_m1_full
plot(lar_dc4_m1)
qqPlot(resid(lar_dc4_m1))
# no interaction is suggested
lar_dc4_m2 <- lmer(s_dc4_ddiff ~ co2 + s_logmoist+s_temp + s_logpar +(1|ring), data = grass_DS_acr)
Anova(lar_dc4_m2, test.statistic = "F")
summary(lar_dc4_m2)
VarCorr(lar_dc4_m2)
plot(lar_dc4_m2)
qqPlot(resid(lar_dc4_m2))

# no variaation was explained by ring, so remove to avaoid convergence issues
lar_dc4_m3 <- lm(s_dc4_ddiff ~ co2 + s_logmoist+s_temp + s_logpar, data = grass_DS_acr)
# coefficients and RI
dc4_coef        <- confint(lar_dc4_m3)
dc4_coef_90     <- confint(lar_dc4_m3, level = .9)
save(dc4_coef, dc4_coef_90, file  = "output/Data/coef_dc4_env.RData")
load("output/Data/coef_dc4_env.RData")
lar_dc4_m2_full <- dredge(lar_dc4_m2, REML = F)
dc4_coef_imp    <- importance(lar_dc4_m2_full)
dc4_coef_imp


# partial residual plot
range(visreg(lar_dc4_m2, xvar = "s_temp")$res$visregRes)

pdf(file = "output/figs/lar_dC4_partial_regression_plot.pdf", width = 5, height = 5)
create_resplot(lar_dc4_m2, 
               ylab = expression(Adj.~annual~rates~of~change~"in"~dominant~C[4]),
               ylim = c(-2, 2.5))

dev.off()

png(file = "output/figs/lar_dC4_partial_regression_plot.png", width = 5, height = 5, res = 600, units = "in")
create_resplot(lar_dc4_m2, 
               ylab = expression(Adj.~annual~rates~of~change~"in"~dominant~C[4]),
               ylim = c(-2, 2.5))

dev.off()




# subordinate c3 ----------------------------------------------------------

lar_sc3_m1 <- lmer(s_sc3_ddiff ~ co2 * (s_logmoist+s_temp + s_logpar) +(1|ring), data = grass_DS_acr)
lar_sc3_m1_full <- dredge(lar_sc3_m1, fixed = c("co2", "s_logmoist", "s_temp",  "s_logpar"), REML = F)
lar_sc3_m1_full
plot(lar_sc3_m1)
qqPlot(resid(lar_sc3_m1))
# no interaction is suggested
lar_sc3_m2 <- lmer(s_sc3_ddiff ~ co2 + s_logmoist+s_temp + s_logpar +(1|ring), data = grass_DS_acr)
Anova(lar_sc3_m2, test.statistic = "F")
plot(lar_sc3_m2)
qqPlot(resid(lar_sc3_m2))
summary(lar_sc3_m2)
# mcp.fnc(lar_sc3_m2)
# lar_sc3_m3 <- update(lar_sc3_m2, subset = -which.min(resid(lar_sc3_m2)))
# plot(lar_sc3_m3)
# qqPlot(resid(lar_sc3_m3))
# confint(lar_sc3_m2, method = "boot", nsim = 99, level = .9)
# confint(lar_sc3_m3, method = "boot", nsim = 99, level = .9)


# coefficients and RI
# sc3_coef        <- confint(lar_sc3_m2, method = "boot", nsim = 999)
# sc3_coef_90     <- confint(lar_sc3_m2, method = "boot", level = .9, nsim = 999)
# save(sc3_coef, sc3_coef_90, file = "output/Data/coef_sc3_env.RData")
load("output/Data/coef_sc3_env.RData")
lar_sc3_m2_full <- dredge(lar_sc3_m2, REML = F)
sc3_coef_imp    <- importance(lar_sc3_m2_full)
sc3_coef_imp


# partial residual plot
range(visreg(lar_sc3_m2, xvar = "s_temp")$res$visregRes)

pdf(file = "output/figs/lar_sC3_partial_regression_plot.pdf", width = 5, height = 5)
create_resplot(lar_sc3_m2, 
               ylab = expression(Adj.~annual~rates~of~change~"in"~subordinate~C[3]),
               ylim = c(-2.1, 3.0))

dev.off()

png(file = "output/figs/lar_sC3_partial_regression_plot.png", width = 5, height = 5, res = 600, units = "in")
create_resplot(lar_sc3_m2, 
               ylab = expression(Adj.~annual~rates~of~change~"in"~subordinate~C[3]),
               ylim = c(-2.1, 2.8))

dev.off()




# dominant c3 ----------------------------------------------------------

lar_dc3_m1 <- lmer(s_dc3_ddiff ~ co2 * (s_logmoist+s_temp + s_logpar) +(1|ring), data = grass_DS_acr)
plot(lar_dc3_m1)
qqPlot(resid(lar_dc3_m1))
lar_dc3_m1_full <- dredge(lar_dc3_m1, fixed = c("co2", "s_logmoist", "s_temp",  "s_logpar"), REML = F)
lar_dc3_m1_full
  # no interaction is suggested

lar_dc3_m2 <- lmer(s_dc3_ddiff ~ co2 + s_logmoist+s_temp + s_logpar +(1|ring), data = grass_DS_acr)
Anova(lar_dc3_m2, test.statistic = "F")
plot(lar_dc3_m2)
qqPlot(resid(lar_dc3_m2))

VarCorr(lar_dc3_m2) # no variation was explained by random factor

lar_dc3_m3 <- lm(s_dc3_ddiff ~ co2 + s_logmoist+s_temp + s_logpar, data = grass_DS_acr)
# coefficients and RI
dc3_coef    <- confint(lar_dc3_m3)
dc3_coef_90 <- confint(lar_dc3_m3, level = .9)
save(dc3_coef, dc3_coef_90, file = "output/Data/coef_dc3_env.RData")
load("output/Data/coef_dc3_env.RData")
lar_dc3_m2_full <- dredge(lar_dc3_m2, REML = F)
dc3_coef_imp    <- importance(lar_dc3_m2_full)
dc3_coef_imp


# partial residual plot
range(visreg(lar_dc3_m2, xvar = "s_temp")$res$visregRes)

pdf(file = "output/figs/lar_dC3_partial_regression_plot.pdf", width = 5, height = 5)
create_resplot(lar_dc3_m2, 
               ylab = expression(Adj.~annual~rates~of~change~"in"~dominant~C[3]),
               ylim = c(-2.4, 2))

dev.off()

png(file = "output/figs/lar_dC3_partial_regression_plot.png", width = 5, height = 5, res = 600, units = "in")
create_resplot(lar_dc3_m2, 
               ylab = expression(Adj.~annual~rates~of~change~"in"~dominant~C[3]),
               ylim = c(-2.4, 2))
dev.off()




# summary -----------------------------------------------------------------
lar_sd_c34_ml <- list(subordinate_c4 = lar_sc4_m2, 
                      subordinate_c3 = lar_sc3_m2, 
                      dominant_c4    = lar_dc4_m2, 
                      dominant_c3    = lar_dc3_m2)

lar_sd_c34_aov <- ldply(lar_sd_c34_ml, function(x) tidy(Anova(x, test.statistic = "F"))) %>% 
  mutate(term = mapvalues(term, c("s_logmoist", "s_temp", "s_logpar"), c("Moist", "Temp", "PAR")),
         Df.res = round(Df.res, 0),
         p.value = round(p.value, 3),
         statistic = round(statistic, 2)) %>% 
  rename(Fval = statistic)



# . summary table of multiple regression analysis on LAR ------------------

coef_imp_df <-ldply(list(c4.subordinate = sc4_coef_imp,
                         c3.subordinate = sc3_coef_imp,
                         c4.dominant    = dc4_coef_imp,
                         c3.dominant    = dc3_coef_imp,
                         c4.all         = c4_coef_imp,
                         c3.all         = c3_coef_imp), 
                    function(x) data.frame(RI = x) %>% 
                      mutate(predictor = dplyr::recode(row.names(.), co2 = 'co2elev'),
                             RI = round(RI, 3))) %>% 
  mutate(type = tstrsplit(.id, split = "[.]")[[2]],
         PFG  = tstrsplit(.id, split = "[.]")[[1]]) %>% 
  select(-.id)

create_smmry_coeftbl <- function(coefaval, model, coefpos1, coefpos2){
  d <- cbind(summary(model)$coef[, "Estimate"], coefaval[coefpos1:coefpos2, ])
  colnames(d) <- c("estimate", "lwr", "upr")
  dd <- data.frame(d) %>% 
    mutate(predictor = row.names(.)) %>% 
    mutate_each(funs(round(., 3)), estimate, lwr, upr) %>% 
    transmute(predictor, coefs = paste0(estimate, " [", lwr, ", ", upr, "]"))
  return(dd)
}


summary(c4d_m2)$coef[, "Estimate"]
coeftbl_c4  <- ldply(list('coef95CI' = c4_coef,  'coef90CI' = c4_coef_90),  create_smmry_coeftbl, coefpos1 = 1, coefpos2 = 5, model = c4d_m2)
coeftbl_sc4 <- ldply(list('coef95CI' = sc4_coef, 'coef90CI' = sc4_coef_90), create_smmry_coeftbl, coefpos1 = 1, coefpos2 = 5, model  = lar_sd_c34_ml$subordinate_c4)
coeftbl_dc4 <- ldply(list('coef95CI' = dc4_coef, 'coef90CI' = dc4_coef_90), create_smmry_coeftbl, coefpos1 = 1, coefpos2 = 5, model  = lar_sd_c34_ml$dominant_c4)
coeftbl_c3  <- ldply(list('coef95CI' = c3_coef,  'coef90CI' = c3_coef_90),  create_smmry_coeftbl, coefpos1 = 3, coefpos2 = 7, model = c3d_m2)
coeftbl_sc3 <- ldply(list('coef95CI' = sc3_coef, 'coef90CI' = sc3_coef_90), create_smmry_coeftbl, coefpos1 = 3, coefpos2 = 7, model  = lar_sd_c34_ml$subordinate_c3)
coeftbl_dc3 <- ldply(list('coef95CI' = dc3_coef, 'coef90CI' = dc3_coef_90), create_smmry_coeftbl, coefpos1 = 1, coefpos2 = 5, model  = lar_sd_c34_ml$dominant_c3)



coeftbl_all <- ldply(list('c4.all'         = coeftbl_c4, 
                          'c4.subordinate' = coeftbl_sc4, 
                          'c4.dominant'    = coeftbl_dc4,
                          'c3.all'         = coeftbl_c3, 
                          'c3.subordinate' = coeftbl_sc3,
                          'c3.dominant'    = coeftbl_dc3),
                     .id = "type") %>% 
  mutate(PFG  = tstrsplit(type, split = "[.]")[[1]],
         type = tstrsplit(type, split = "[.]")[[2]]) %>%
  spread(.id, coefs) %>% 
  left_join(coef_imp_df) %>% 
  arrange(PFG, type, predictor) %>% 
  select(PFG, type, predictor, everything())

coeftbl_all