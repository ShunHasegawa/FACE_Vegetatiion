
# prepare df --------------------------------------------------------------

head(C3grassC4)

c43_ratio <- C3grassC4 %>% 
  group_by(year, co2, ring, plot, PFG) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  spread(key = PFG, value) %>% 
  mutate(c43_r = c4 / c3)


c43_ratio_year0 <- c43_ratio %>%
  filter(year == "Year0") %>%
  select(ring, plot, c43_r) %>%
  rename(ratios0 = c43_r) %>% 
  left_join(filter(c43_ratio, year != "Year0"), by = c("ring", "plot")) %>% 
  droplevels()




# analysis ----------------------------------------------------------------


plot(log(c43_r + .01) ~ log(ratios0 + .01), data = c43_ratio_year0)
summary(c43_ratio_year0)


c43_m1 <- lmer(I(log(c43_r + .01)) ~ co2 * year + I(log(ratios0 + .01)) + (1|ring) + (1|ring:plot) + (1|ring:year), data = c43_ratio_year0)
plot(c43_m1)
qqPlot(residuals(c43_m1))
c43_m2 <- update(m1, subset = -order(resid(m1))[1:2])
plot(c43_m2)
qqPlot(residuals(c43_m2))

Anova(c43_m1, test.statistic = "F")
Anova(c43_m2, test.statistic = "F")



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
         rlsmean    = exp(lsmean) - 0.01,
         rlowerCL   = exp(lower.CL) - 0.01,
         rupperCL   = exp(upper.CL) - 0.01,
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
  
  labs(y = expression(C[4]:C[3]~ratios)) +
  geom_text(data = c43_rr_d, 
            aes(label = rr), x =  Inf, y = Inf, 
            hjust = 1.1, vjust = 1.5, size = 3)


ggsavePP(filename = "output/figs/adjusted_c43_ratio", width = 3, height = 3,
         plot = fig_c43)