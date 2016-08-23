# Prepare data frames -----------------------------------------------------

summary(veg_FullVdf)

# proportion of each form for each plot
PropDF <- veg_FullVdf %>%
  group_by(year, co2, block, ring, plot, id) %>% 
  mutate(Total = sum(value)) %>% # total sum for each plot
  group_by(form, Total, add = TRUE) %>% 
  summarise(value = sum(value)) %>% # sum for each form
  mutate(value = value / Total) %>% 
  filter(form %in% c("Forb", "Grass")) %>% 
  ungroup()


# Move Year0 value to a new column to be used as covariate in the analyssis
# below

# subsequent years
subyear_dd <- filter(PropDF, year != "Year0") 

# merge with Year0
PropDF_year0 <- PropDF %>% 
  filter(year == "Year0") %>%
  select(form, id, value) %>%
  rename(value0 = value) %>%
  left_join(subyear_dd, by = c("id", "form")) %>% 
  mutate(logitv0 = logit(value0)) # logit transformation

# plot
unique(PropDF_year0$form)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 0.5), cex = .5)
d_ply(PropDF_year0, .(form), function(x) {
  plot(value ~ value0, pch = 19, col = year, data = x, main = unique(x$form))
  plot(logit(value) ~ logitv0, pch = 19, col = year, data = x, 
       main = unique(x$form))
})

# Analysis ----------------------------------------------------------------

# models to be tested
m_list <- dlply(PropDF_year0, .(form), function(x){
  m1 <- lmer(logit(value) ~ co2 * year + logitv0 + (1|block) + (1|ring) + (1|id), data = x)
  m2 <- update(m1, ~ . - (1|block))
  if (AICc(m1) >= AICc(m2)) return(m2) else return(m1)
})
llply(m_list, function(x) Anova(x, test.statistic = "F"))
  
PropDF_year0 %>%
  select(id, year, form, value) %>% 
  spread(form, value) %>% 
  qplot(x = Forb, y = Grass, geom = "point", data = .)

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



# Create fig --------------------------------------------------------------

# df for Year0
d <- PropDF %>%
  filter(year == "Year0") %>% 
  mutate(year = factor(year, levels = paste0("Year", 0:3)))

# df for median of Year0
d_med <- d %>%  
  group_by(form) %>% 
  summarise(M = median(value)) %>% 
  mutate(Med = "Md[Year0]")


ci_dd <- ci_dd %>% 
  mutate(rlsmean = boot::inv.logit(lsmean), # reverse transform
         rlowerCL = boot::inv.logit(lower.CL),
         rupperCL = boot::inv.logit(upper.CL),
         year = factor(year, levels = paste0("Year", 0:3)))

# fig
dodgeval <- .4


p <- ggplot(ci_dd, aes(x = year, y = rlsmean, fill = co2)) +
  geom_bar(stat = "identity")+
  scale_y_continuous()
  facet_grid(. ~ .id)
p
  
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








  # > Grass proportion -------------------------------------------------------

gProp_dd <- filter(PropDF_year0, variable == "Grassprop")

par(mfrow = c(1, 3))
plot(value ~ value0, pch = 19, col = year,  data = gProp_dd)
plot(asin(value) ~ value0, pch = 19, col = year, data = gProp_dd)
plot(logit(value) ~ value0, pch = 19, col = year, data = gProp_dd)

gm1 <- lmer(value ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
            data = gProp_dd)
gm2 <- lmer(logit(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
            data = gProp_dd)
gm3 <- lmer(asin(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
            data = gProp_dd)
ldply(list(gm1, gm2, gm3), r.squared)

plot(gm2)
qqnorm(resid(gm2))
qqline(resid(gm2))
AnvF_Grass <- Anova(gm2, test.statistic = "F")
AnvF_Grass

# . post-hoc test ---------------------------------------------------------

plot(lmerTest::lsmeans(gm2))
lsmeans::lsmeans(gm2, pairwise ~ co2 | year)

# Forb --------------------------------------------------------------------

########
# Forb #
########

forbProp_dd <- filter(PropDF_year0, variable == "Forbprop")

par(mfrow = c(1, 3))
plot(value ~ value0, pch = 19, col = year, data = forbProp_dd)
plot(asin(value) ~ value0,  pch = 19, col = year, data = forbProp_dd)
plot(logit(value) ~ value0, pch = 19, col = year, data = forbProp_dd)

fm1 <- lmer(value ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
            data = forbProp_dd)
fm2 <- lmer(logit(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
            data = forbProp_dd)
fm3 <- lmer(asin(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
            data = forbProp_dd)
ldply(list(fm1, fm2, fm3), r.squared)

Anova(fm2)
AnvF_forb <- Anova(fm2, test.statistic = "F")
AnvF_forb

# model diagnosis
plot(fm2)
qqnorm(resid(fm2))
qqline(resid(fm2))


# . Legume proportion -----------------------------------------------------
legProp_dd <- filter(PropDF_year0, variable == "Legume")

par(mfrow = c(1, 3), cex = .8)
plot(value ~ value0, pch = 19, col = year, data = legProp_dd)
plot(asin(value) ~ value0,  pch = 19, col = year, data = legProp_dd)
plot(logit(value) ~ value0, pch = 19, col = year, data = legProp_dd)

lpm1 <-  lmer(value ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
              data = legProp_dd)
lpm2 <-  lmer(logit(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
              data = legProp_dd)
lpm3 <-  lmer(asin(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
              data = legProp_dd)
ldply(list(lpm1, lpm2, lpm3), r.squared)

plot(lpm2)
qqnorm(resid(lpm2))
qqline(resid(lpm2))

AnvF_legume <- Anova(lpm2, test.statistic = "F")
AnvF_legume

# . Non_legume proportion -----------------------------------------------------
nonlegProp_dd <- filter(PropDF_year0, variable == "Non_leg")

par(mfrow = c(1, 3), cex = .8, mar = c(4, 4, 1, .5))
plot(value ~ value0, pch = 19, col = year, data = nonlegProp_dd)
plot(asin(value) ~ value0,  pch = 19, col = year, data = nonlegProp_dd)
plot(logit(value) ~ value0, pch = 19, col = year, data = nonlegProp_dd)

nlpm1 <-  lmer(value ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
               data = nonlegProp_dd)
nlpm2 <-  lmer(logit(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
               data = nonlegProp_dd)
nlpm3 <-  lmer(asin(value) ~ year * co2 + value0 + (1|block) +  (1|ring) + (1|id), 
               data = nonlegProp_dd)
ldply(list(nlpm1, nlpm2, nlpm3), r.squared)

plot(nlpm1)
qqnorm(resid(nlpm1))
qqline(resid(nlpm1))

AnvF_nonlegume <- Anova(nlpm1, test.statistic = "F")
AnvF_nonlegume


# Summary -----------------------------------------------------------------

###########
# Summary #
###########

# F from LMM
SummaryAnvF_PFG <- ldply(list(Forb       = AnvF_forb,
                              Legume     = AnvF_legume,
                              Non_legume = AnvF_nonlegume,
                              Grass      = AnvF_Grass, 
                              C3total    = AnvF_PropC3_total, 
                              C3grass    = AnvF_C3grass, 
                              C4grass    = AnvF_C4grass), 
                         function(x) data.frame(x, terms = row.names(x)), 
                         .id = "variable")

# summary of anova
PFGResAnvF <- SummaryAnvF_PFG
names(PFGResAnvF)[5] <- "Pr"
PFGResAnvF <- within(PFGResAnvF, {
              F      <- round(F, 2)
              Df.res <- round(Df.res, 0)
              Pr     <- round(Pr, 3)
            })
PFGResAnvF$terms <- factor(PFGResAnvF$terms, 
                           levels = c("value0", "co2", "year", "year:co2"))
PFGResAnvF <- PFGResAnvF[order(PFGResAnvF$variable, PFGResAnvF$terms), ]
write.csv(PFGResAnvF, "output/table/FACE_PFG_AnvF.csv", row.names = FALSE)

# summary fig with 95% CI
models <- list(Forb = fm2, 
               Legume = lpm2,
               Non_legume = nlpm1,
               Grass = gm2,
               C3_grass = c3gm3,
               C4_grass = m3Lmer)

sapply(models, function(x) x@call$data)

CI_dd <- ldply(models, function(x) 
  summary(lsmeans::lsmeans(x, pairwise ~ co2 | year, type = "response")$lsmeans),
  .id = "variable")

CI_dd <- CI_dd %>% 
  mutate(value = ifelse(is.na(response), lsmean, response)) %>% 
  select(-response, -lsmean)

CI_grass <- filter(CI_dd, variable %in% c("Grass", "C3_grass", "C4_grass"))
CI_forb <- filter(CI_dd, !variable %in% c("Grass", "C3_grass", "C4_grass"))

p <- ggplot(CI_grass, aes(x = as.numeric(year), y = value, fill = co2, 
                          shape = variable, group = co2:variable, 
                          linetype = co2))

p2 <- p +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0, 
                position = position_dodge(.1), linetype = "solid") +
  geom_line(position = position_dodge(.1)) + 
  geom_point(position = position_dodge(.1)) +
  scale_shape_manual(values = c(21:23)) +
  scale_linetype_manual(values = c("solid", "dashed"),
                        labels = c("Ambient", expression(eCO[2]))) +
  scale_fill_manual(values = c("black", "white"), 
                    labels = c("Ambient", expression(eCO[2])),
                    guide = guide_legend(override.aes = list(shape = 21))) +
  science_theme +
  theme(legend.position = c(.5, .9), 
        legend.key.height = unit(.7, "line"),
        legend.box = "horizontal", 
        legend.box.just = "top") +
  scale_x_continuous(breaks = 1:3, labels = 1:3) +
  labs(x = "Year", y = "Proportion (adjusted by Year0 value)") +
  scale_y_continuous(limits = c(0, .65))
p2

p3 <- ggplot(CI_forb, aes(x = as.numeric(year), y = value, fill = co2, 
                          shape = variable, group = co2:variable, 
                          linetype = co2))

p4 <- p3 +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0, 
                position = position_dodge(.1), linetype = "solid") +
  geom_line(position = position_dodge(.1)) + 
  geom_point(position = position_dodge(.1)) +
  scale_shape_manual(values = c(21:23)) +
  scale_linetype_manual(values = c("solid", "dashed"),
                        labels = c("Ambient", expression(eCO[2]))) +
  scale_fill_manual(values = c("black", "white"), 
                    labels = c("Ambient", expression(eCO[2])),
                    guide = guide_legend(override.aes = list(shape = 21))) +
  science_theme +
  theme(legend.position = c(.4, .9), 
        legend.key.height = unit(.7, "line"),
        legend.box = "horizontal", 
        legend.box.just = "top") +
  scale_x_continuous(breaks = 1:3, labels = 1:3) +
  labs(x = "Year", y = NULL) +
  scale_y_continuous(limits = c(0, .65))
p4

# combine two plots
grid.arrange(p2, p4, ncol = 2)
p <- recordPlot() # save

pdf(file = "output/figs/PFG_proportion_95CI.pdf", width = 6, height = 3)
replayPlot(p)
dev.off()

save_png600(filename = "output/figs/PFG_proportion_95CI.png", width = 6, height = 3)
replayPlot(p)
dev.off()

ggsavePP(filename = "output/figs/PFG_proportion_95CI_grass", 
         width = 3, height = 3, plot = p2)
