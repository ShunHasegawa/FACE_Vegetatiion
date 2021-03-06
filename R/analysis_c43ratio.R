
# prepare df --------------------------------------------------------------

head(C3grassC4)

c43_ratio <- C3grassC4 %>% 
  group_by(year, co2, ring, plot, PFG) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  spread(key = PFG, value) %>% 
  mutate(c43_r = c4 / c3) %>% 
  group_by(year, co2, ring) %>% 
  summarise(c43_r = mean(c43_r)) %>% 
  ungroup()


c43_ratio_year0 <- c43_ratio %>%
  filter(year == "Year0") %>%
  select(ring, c43_r) %>%
  rename(ratios0 = c43_r) %>% 
  left_join(filter(c43_ratio, year != "Year0"), by = "ring") %>% 
  droplevels() %>% 
  group_by(year, ring, co2) %>% 
  summarise_each(funs(mean), c43_r, ratios0)




# analysis ----------------------------------------------------------------


plot(log(c43_r) ~ log(ratios0), data = c43_ratio_year0, pch = 19, col = co2)

c43_m1 <- lmer(I(log(c43_r)) ~ co2 * year + I(log(ratios0)) + (1|ring), data = c43_ratio_year0)
plot(c43_m1)
qqPlot(residuals(c43_m1))
Anova(c43_m1, test.statistic = "F")
VarCorr(c43_m1)


# CI and postdoc test -----------------------------------------------------


# compute 95 CI and post-hoc test
c43_lsmeans <- lsmeans::lsmeans(c43_m1, ~ co2 | year)


# CO2 effect
c43_co2_pval <- tidy(Anova(c43_m1, test.statistic = "F")) %>%
  filter(term == "co2") %>% 
  mutate(co2star = get_star(p.value)) %>% 
  select(co2star)

c43_CI_dd <- data.frame(summary(c43_lsmeans)) %>% 
  mutate(year       = factor(year, levels = paste0("Year", 0:3)),
         rlsmean    = exp(lsmean),
         rlowerCL   = exp(lower.CL),
         rupperCL   = exp(upper.CL),
         co2star    = c43_co2_pval[,1],
         value_type = "adjusted")


# df for observed values
c43_obs_d <- c43_ratio %>%
  group_by(year, co2, ring) %>%
  summarise(value = mean(c43_r)) %>%
  ungroup() %>% # grouping informaiton is not required later
  mutate(year = factor(year, levels = paste0("Year", 0:3)),
         value_type = "observed")


# co2 response ratio
c43_rr_d <- c43_CI_dd %>% 
  group_by(co2, co2star) %>% 
  summarise(value = mean(rlsmean)) %>% 
  group_by(co2star) %>% 
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1) %>% 
  mutate(rr = ifelse(rr >= 0,
                     paste0("RR= +", format(rr, digits = 0, nsmall = 2), co2star), 
                     paste0("RR= ", format(rr, digits = 0, nsmall = 2), co2star)))



# create a figure ---------------------------------------------------------

# create fig
dodgeval <- .3


fig_c43 <- ggplot(c43_CI_dd, aes(x = year, y = rlsmean)) +
  
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  
  # observed values
  geom_point(data = c43_obs_d, 
             aes(x = year, y = value, shape = co2, group = co2, col = value_type),  
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
  
  labs(y = expression(C[4]:C[3]~ratios)) +
  geom_text(data = c43_rr_d, 
            aes(label = rr), x =  Inf, y = Inf, 
            hjust = 1.1, vjust = 1.5, size = 3)


ggsavePP(filename = "output/figs/adjusted_c43_ratio", width = 3, height = 3,
         plot = fig_c43)




# fit soil N and P (IEM) and other env vars------------------------------------------

# fit soil nutrient data from Hasegawa et al. 2016 and Juan & Raul's data



# . merge with env  data-------------------------------------------------------------


source("R/process_env.R") # env vars
head(c34growth_moist)
head(iem_raw)

# merge with C34 ratiodata
c43_ratio_iem <- c43_ratio %>% 
  left_join(iem_raw) %>% 
  filter(year != "Year0") %>% 
  left_join(c34growth_moist) %>% 
  mutate(s_c43_r     = scale(log(c43_r))[, 1],    # Z-standardize for multiple regression
         s_n         = scale(log(nitr))[, 1],
         s_p         = scale(log(p))[, 1],
         probe       = ifelse(year %in% c("Year0", "Year1"), "IEM", "PRS"),
         s_moist     = scale(totalmoist)[, 1],
         s_logmoist  = scale(log(totalmoist))[, 1],
         s_temp      = scale(annual_temp2m)[, 1],
         s_par       = scale(PAR)[, 1],
         s_logpar    = scale(log(PAR))[, 1]) %>% 
  arrange(year, ring)
summary(c43_ratio_iem)
plot(s_c43_r ~ log(nitr), c43_ratio_iem, col = factor(year), pch = 19)
plot(s_c43_r ~ log(nitr), c43_ratio_iem, col = co2, pch = 19)
plot(s_c43_r ~ totalmoist, c43_ratio_iem, col = co2, pch = 19)
plot(log(c43_r + .1) ~ log(p), c43_ratio_iem, col = co2, pch = 19)


# anlaysis

car::vif(lm(s_c43_r ~ s_n + s_p + s_moist + s_temp + s_par, data = c43_ratio_iem))

c43_soil_m1 <- lmer(s_c43_r ~ s_n + s_p + s_moist + s_temp + s_par + (1|ring) + (1|probe) + (1|year), data = c43_ratio_iem)
summary(c43_soil_m1)
# probe type doeesn't explain any variation (probe is redundant with year anyway so remove)
c43_soil_m2 <- lmer(s_c43_r ~ s_n + s_p + s_moist + s_temp + s_par + (1|ring) + (1|year), data = c43_ratio_iem)
c43_soil_m3 <- lmer(s_c43_r ~ s_n + s_p + s_logmoist + s_temp + s_logpar + (1|ring) + (1|year), data = c43_ratio_iem)
r.squared(c43_soil_m2)
r.squared(c43_soil_m3)
summary(c43_soil_m2)
Anova(c43_soil_m2, test.statistic = "F")
plot(c43_soil_m2)
qqPlot(resid(c43_soil_m2))

# coefficient and RI
# c43_soil_coef      <- confint(c43_soil_m2, method = "boot", nsim = 999)
# c43_soil_coef_90   <- confint(c43_soil_m2, method = "boot", nsim = 999, level = .9)
# save(c43_soil_coef, c43_soil_coef_90, file = "output/Data/coef_c43ratio.RData")
load("output/Data/coef_c43ratio.RData")
c43_soil_full      <- dredge(c43_soil_m2, REML = F)
c43_soil_coef_impo <- importance(c43_soil_full)
c43_soil_coef_impo




# . fig -------------------------------------------------------------------


# Nitrogen ----------------------------------------------------------------
vireg_obs_n <- visreg(c43_soil_m2,  xvar = "s_n")$res  # adjusted observed values
vireg_fit_n <- visreg(c43_soil_m2,  xvar = "s_n")$fit  # fitted values
vireg_obs_n$co2 <- c43_ratio_iem$co2

nlim <- c(-2.5, 1.5)
resplot_n <- function(){
  plot(visregRes ~ s_n, data = vireg_obs_n, type = "n", 
       ylim = c(-1.7, 2.6),
       xlim = nlim,
       axes = F)
  axis(side = 1)
  axis(side = 2, las = 2)
  box()
  polygon(c(rev(vireg_fit_n$s_n), vireg_fit_n$s_n), 
          c(rev(vireg_fit_n$visregLwr), vireg_fit_n$visregUpr), 
          col = "gray80", border = NA)
  lines(visregFit ~ s_n, vireg_fit_n, lwd = 3)
  points(visregRes ~ s_n, data = vireg_obs_n, subset = co2 == "amb", pch = 19)
  points(visregRes ~ s_n, data = vireg_obs_n, subset = co2 == "elev", pch = 0)
}
resplot_n()

# N boxplot
n_boxplot <- function(){
  boxplot(s_n ~ co2 * year, c43_ratio_iem, horizontal = TRUE, 
          ylim = nlim,
          axes = F, col = c("gray80", "white"))
  axis(side = 1, labels = FALSE, tck = .03)
  axis(side = 2, at = c(1.5, 3.5, 5.5), labels = paste0("Year", 1:3), 
       las = 2)
  box()
}



# Phosphorus ----------------------------------------------------------------
vireg_obs_p <- visreg(c43_soil_m2,  xvar = "s_p")$res  # adjusted observed values
vireg_fit_p <- visreg(c43_soil_m2,  xvar = "s_p")$fit  # fitted values
vireg_obs_p$co2 <- c43_ratio_iem$co2

plim <- c(-1.8, 1.4)
resplot_p <- function(){
  plot(visregRes ~ s_p, data = vireg_obs_p, type = "n",
       ylim = c(-1.7, 2.6), axes = F, ann = F,
       xlim = plim)
  axis(side = 1, at = seq(-1.5, 1.5, 1), labels = seq(-1.5, 1.5, 1))
  axis(side = 2, las = 2, labels = F)
  box()
  # polygon(c(rev(vireg_fit_p$s_p), vireg_fit_p$s_p), 
  #         c(rev(vireg_fit_p$visregLwr), vireg_fit_p$visregUpr), 
  #         col = "gray80", border = NA)
  # lines(visregFit  ~ s_p, vireg_fit_p, lwd = 3)
  points(visregRes ~ s_p, data = vireg_obs_p, subset = co2 == "amb", pch = 19)
  points(visregRes ~ s_p, data = vireg_obs_p, subset = co2 == "elev", pch = 0)
}
resplot_p()


# P boxplot
p_boxplot <- function(){
  boxplot(s_p ~ co2 * year, c43_ratio_iem, horizontal = TRUE, axes = F, 
          col = c("gray80", "white"),
          ylim = plim)
  axis(side = 1, at = seq(-1.5, 1.5, 1), labels = F, tck = .03)
  axis(side = 2, at = c(1.5, 3.5, 5.5), labels = F, las = 2)
  box()
}


# merge figs --------------------------------------------------------------

source("http://www.math.mcmaster.ca/bolker/R/misc/legendx.R") # allow to change legend box sizes when defined using fill

plot_c43r_np <- function(){
  par(mfrow = c(2, 2), tck = .03, oma = c(4, 5, .5, .5), mgp = c(3, .1, 0))
  par(mar = c(0, 0, .5, .5))
  n_boxplot()
  text(par("usr")[1], par("usr")[4], expression(bold((a))), adj = c(0, 1))
  legend("bottomleft", legend = c(expression(eCO[2]), "Amb"), 
         fill = c("white", "gray80"), bty = "n", box.cex = c(1.1, .8))
  p_boxplot()
  text(par("usr")[1], par("usr")[4], expression(bold((b))), adj = c(0, 1))
  
  # partial regression plots
  resplot_n()
  text(par("usr")[1], par("usr")[4], expression(bold((c))), adj = c(0, 1))
  mtext(side = 1, text = expression(Adj.~Log[e](soil~N,~ng~cm^"-2"~d^"-1")), 
        line = 2, cex = .9)
  mtext(side = 2, text = expression(Adj.~Log[e](C[4]:C[3])), line = 2)
  legend("bottomleft", legend = c(expression(eCO[2]), "Amb"), pch = c(0, 19),
         pt.cex = 2, bty = "n")
  resplot_p()
  text(par("usr")[1], par("usr")[4], expression(bold((d))), adj = c(0, 1))
  mtext(side = 1, text = expression(Adj.~Log[e](soil~P,~ng~cm^"-2"~d^"-1")), 
        line = 2, cex = .9)
}


pdf(file = "output/figs/C43ratio_NP_partial_regression_plot.pdf", width = 5.5, height = 5)
plot_c43r_np()
dev.off()

png(file = "output/figs/C43ratio_NP_partial_regression_plot.png", width = 5.5, height = 5, res = 600, units = "in")
plot_c43r_np()
dev.off()

