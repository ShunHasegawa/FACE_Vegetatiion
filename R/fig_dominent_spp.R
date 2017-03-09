
# prepare data frame ------------------------------------------------------

# create models to be tested


# microlaena and cynodon
mic_cyn <- veg_FullVdf %>% 
  filter(variable %in% c("Microlaena.stipoides", "Cynodon.dactylon")) %>% 
  group_by(variable, year, block, co2, ring, plot, id, RY) %>% 
  summarise(value = sum(value)) %>% 
  ungroup()


# C3 and C4 grass abundance
C3grassC4_abund <- veg_FullVdf %>%
  filter(form == "Grass") %>% 
  group_by(PFG, year, block, ring, co2, plot, id, RY) %>% 
  summarise_each(funs(sum), value) %>% 
  ungroup() %>% 
  mutate(PFG = factor(paste(PFG, "total", sep = "_"))) %>% 
  rename(variable = PFG)


# merge the above
dominentsp <- bind_rows(mic_cyn, C3grassC4_abund)


# Year0
dominentsp_year0 <- dominentsp %>%
  filter(year == "Year0") %>%
  select(id, value, variable) %>%
  rename(value0 = value) %>% 
  left_join(filter(dominentsp, year != "Year0"), by = c("id", "variable")) 


# transformations for x value
xtrans <- c("value0", "logitv0", "logv0", "sqrtv0")




# analysis ----------------------------------------------------------------


# > Microlaena --------------------------------------------------------------

mic_df <- dominentsp_year0 %>%
  filter(variable == "Microlaena.stipoides") %>% 
  mutate(logitv  = logit(value/100),
         logitv0 = logit(value0/100),
         logv0   = log(value0 + 1), 
         sqrtv0  = sqrt(value0))

mic_m1 <- lmer(logitv ~ co2 * year + logitv0 + (1|block) + (1|ring) + (1|id) + (1|RY), data = mic_df, REML = F)
mic_m2 <- lmer(logitv ~ co2 * year + logv0   + (1|block) + (1|ring) + (1|id) + (1|RY), data = mic_df, REML = F)
mic_m3 <- lmer(logitv ~ co2 * year + value0  + (1|block) + (1|ring) + (1|id) + (1|RY), data = mic_df, REML = F)
mic_m4 <- lmer(logitv ~ co2 * year + sqrtv0  + (1|block) + (1|ring) + (1|id) + (1|RY), data = mic_df, REML = F)
model.sel(mic_m1, mic_m2, mic_m3, mic_m4, extra = "r.squaredGLMM")

mic_m5 <- lmer(logitv ~ co2 * year + sqrtv0  + (1|block) + (1|ring) + (1|id) + (1|RY), data = mic_df)
mic_m6 <- lmer(logitv ~ co2 * year + sqrtv0  + (1|ring) + (1|id) + (1|RY), data = mic_df)
model.sel(mic_m5, mic_m6)
plot(mic_m6)
qqnorm(resid(mic_m6))
qqline(resid(mic_m6))
mic_m_fin <- mic_m6




# > Cynodon --------------------------------------------------------------


cyn_df <- dominentsp_year0 %>%
  filter(variable == "Cynodon.dactylon") %>% 
  mutate(logitv  = logit(value/100),
         logitv0 = logit(value0/100),
         logv0   = log(value0 + 1), 
         sqrtv0  = sqrt(value0))


# . possible transformation -----------------------------------------------


cyn_p <- ldply(xtrans, function(x){
  f <- formula(paste("I(value + 1) ~ co2 * year +", x))
  a <- boxcox(f, data = cyn_df, plotit = FALSE)
  p <- a$x[which.max(a$y)]
  return(data.frame(xt = x, p))
})


# inspect the above-suggested transformation
cyn_ms <- mlply(cyn_p, function(xt, p){
  f <- formula(paste("value^(", p, ") ~ co2 * year +", xt, 
                     "+ (1|block) + (1|ring) + (1|id) + (1|RY)"))
  m <- lmer(f, data = cyn_df)
  return(m)
})
llply(cyn_ms, function(x) Anova(x, test.statistic = "F"))

plot(cyn_ms[[1]])
plot(cyn_ms[[2]])
plot(cyn_ms[[3]])
plot(cyn_ms[[4]])
par(mfrow = c(2, 2))
l_ply(cyn_ms, function(x){
  qqnorm(resid(x))
  qqline(resid(x))
})



# . model selection -------------------------------------------------------


# check model4; there was an outlier
cyn_m1 <- cyn_ms[[4]]
which.min(qqnorm(resid(cyn_m1))$y)
cyn_m2 <- update(cyn_ms[[4]], subset = -72)
Anova(cyn_m1, test.statistic = "F")
Anova(cyn_m2, test.statistic = "F")
 # maybe interaction was driven by the outlier


# check model1
cyn_m3 <- cyn_ms[[1]]
cyn_m4 <- update(cyn_m3, ~ . - (1|block))
AICc(cyn_m3, cyn_m4)
Anova(cyn_m4, test.statistic = "F")
plot(cyn_m4)
qqnorm(resid(cyn_m4))
qqline(resid(cyn_m4))
### use this model
cyn_m_fin <- cyn_m4




# > C4 abundance ------------------------------------------------------------



c4ab_df <- dominentsp_year0 %>%
  filter(variable == "c4_total") %>% 
  mutate(logitv  = logit(value/100),
         logitv0 = logit(value0/100),
         logv0   = log(value0 + 1),
         sqrtv0  = sqrt(value0))


# determin the "best" xtransformation for each of variaous y tranformtion to
# express linearity
c4ab_df_p <- ldply(xtrans, function(x){
  f <- formula(paste("I(value + 1) ~ co2 * year +", x))
  a <- boxcox(f, data = c4ab_df, plotit = FALSE)
  p <- a$x[which.max(a$y)]
  return(data.frame(xt = x, p))
})


# inspect the above-suggested transformation
c4ab_df_ms <- mlply(c4ab_df_p, function(xt, p){
  f <- formula(paste("value^(", p, ") ~ co2 * year +", xt, 
                     "+ (1|block) + (1|ring) + (1|id) + (1|RY)"))
  m <- lmer(f, data = c4ab_df)
  return(m)
})
llply(c4ab_df_ms, function(x) Anova(x, test.statistic = "F"))

plot(c4ab_df_ms[[1]])  # not bad
plot(c4ab_df_ms[[2]])  # bad
plot(c4ab_df_ms[[3]])  # bad
plot(c4ab_df_ms[[4]])  # not bad
par(mfrow = c(2, 2))
l_ply(c4ab_df_ms, function(x){
  qqnorm(resid(x))
  qqline(resid(x))
})
## model 1 or 4

# check model1
c4ab_m1 <- c4ab_df_ms[[1]]
c4ab_m2 <- update(c4ab_df_ms[[1]], ~ . - (1|block))
AICc(c4ab_m1, c4ab_m2)
Anova(c4ab_m2, test.statistic = "F")
plot(c4ab_m2)
qqnorm(resid(c4ab_m2))
qqline(resid(c4ab_m2))
## this looks good
c4ab_m_fin <- c4ab_m2



# > C3 abundance ------------------------------------------------------------


c3ab_df <- dominentsp_year0 %>%
  filter(variable == "c3_total") %>% 
  mutate(logitv  = logit(value/100),
         logitv0 = logit(value0/100),
         logv0   = log(value0 + 1),
         sqrtv0  = sqrt(value0))


# check transformation
c3ab_df_p <- ldply(c(xtrans), function(x){
  f <- formula(paste("I(value + 1) ~ co2 * year +", x))
  a <- boxcox(f, data = c3ab_df, plotit = FALSE, lambda = seq(-10, 10, .2))
  p <- a$x[which.max(a$y)]
  return(data.frame(xt = x, p))
})
c3ab_df_p
  # not improved


plot(value ~ log(value0), data = c3ab_df)
plot(value ~ value0, data = c3ab_df)
plot(log(value) ~ log(value0), data = c3ab_df)
plot(logit(value / 10) ~ value0, data = c3ab_df)

c3ab_m1 <- lmer(value ~ co2 * year + log(value0) + (1|block) + (1|ring) + (1|id) + (1|RY), data = c3ab_df)
c3ab_m2 <- lmer(value ~ co2 * year + log(value0) + (1|ring) + (1|id) + (1|RY), data = c3ab_df)
model.sel(c3ab_m1, c3ab_m2, extra = "r.squaredGLMM")
plot(c3ab_m2)
qqnorm(resid(c3ab_m2))
qqline(resid(c3ab_m2))
Anova(c3ab_m2, test.statistic = "F")
c3ab_m_fin <- c3ab_m2




# CI and post-hoc test ----------------------------------------------------

dom_m_list <- list('Microlaena.stipoides' = mic_m_fin, 
                   'Cynodon.dactylon'     = cyn_m_fin, 
                   'C3 abundance'         = c3ab_m_fin,
                   'C4 abundance'         = c4ab_m_fin)


# compute 95 CI and post-hoc test
lsmeans_list <- llply(dom_m_list, function(x) {
  lsmeans::lsmeans(x, ~ co2 | year)
})


# 95% CI
CI_ll <- llply(lsmeans_list, function(x) data.frame(summary(x)))
CI_dd

# reverse transformation and standerdise for 1m x 1m
sapply(dom_m_list, function(x) x@call)
  ## Microlaena: logit
  ## Cynodon: 0.5
  ## C3: no transforamtion
  ## C4: 0.5

CI_ll[["Microlaena.stipoides"]] <- CI_ll[["Microlaena.stipoides"]] %>% 
  mutate_each(funs(r = boot::inv.logit(.)*100/4), lsmean, lower.CL, upper.CL)

CI_ll[["Cynodon.dactylon"]]     <- CI_ll[["Cynodon.dactylon"]] %>% 
  mutate_each(funs(r = (.^2)/4), lsmean, lower.CL, upper.CL)

CI_ll[["C4 abundance"]]         <- CI_ll[["C4 abundance"]] %>% 
  mutate_each(funs(r = (.^2)/4), lsmean, lower.CL, upper.CL)

CI_ll[["C3 abundance"]]         <- CI_ll[["C3 abundance"]] %>% 
  mutate_each(funs(r = ./4), lsmean, lower.CL, upper.CL)

CI_dd <- ldply(CI_ll)

# post-hoc test
contrast_dd <- ldply(lsmeans_list, function(x) {
  data.frame(summary(pairs(x)[1:3], adjust = "fdr"))
}) %>% 
  mutate(co2 = factor("amb", levels = c("amb", "elev")),
         star = cut(p.value, right = FALSE,
                    breaks = c(0, .1, .05, .01, .001, 1),  
                    labels = c("***", "**", "*", "\u2020", ""))) %>% 
  select(.id, year, co2, p.value, star)


# CO2 effect
dom_aov_df <- ldply(dom_m_list, function(x) tidy(Anova(x, test.statistic = "F")),  # Anova result of models
                    .progress = "text")
dom_co2_pval <- dom_aov_df %>%                                                     # get p-values for CO2 term
  filter(term == "co2") %>% 
  mutate(co2star = get_star(p.value)) %>% 
  select(.id, co2star) 


# merge
ci_dd <- left_join(CI_dd, contrast_dd, by = c(".id", "year", "co2")) %>% 
  rename(rlsmean = lsmean_r, rlowerCL = lower.CL_r, rupperCL = upper.CL_r) %>% 
  mutate(year       = factor(year, levels = paste0("Year", 0:3)),
         value_type = "adjusted") %>% 
  left_join(dom_co2_pval, by = ".id") %>%                                         # merge with pvalues for CO2 term
  mutate(.id        = gsub("[.]", " ", .id),
         .id        = factor(.id, levels = c("C3 abundance", "Microlaena stipoides", 
                                             "C4 abundance", "Cynodon dactylon")),
         plot_lab   = as.character(factor(.id, labels = paste0("(", letters[1:4], ")"))))                 # sub-plot label
         
ci_dd$star[is.na(ci_dd$star)] <- ""




# create fig --------------------------------------------------------------


# df for observed vlaues
obs_d <- dominentsp %>%
  group_by(year, co2, ring, variable) %>% 
  summarise(N = sum(!is.na(value)), value = mean(value)) %>% 
  ungroup() %>% 
  mutate(year = factor(year, levels = paste0("Year", 0:3)),
         value = value/4,                                    # standardise for 1mx1m plot
         .id = mapvalues(variable, paste0(c("c3", "c4"), "_total"), 
                         paste(c("C3", "C4"), "abundance")),
         .id = gsub("[.]", " ", .id),
         .id = factor(.id, levels = c("C3 abundance", "Microlaena stipoides", 
                                      "C4 abundance", "Cynodon dactylon")),
         value_type = "observed")


# df for plot labels and response ratios
plab_d <- ci_dd %>% 
  group_by(.id, plot_lab, co2, co2star) %>% 
  summarise(value = mean(rlsmean)) %>% 
  group_by(.id, plot_lab, co2star) %>% 
  summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1) %>% 
  mutate(rr = ifelse(rr >= 0,
                     paste0("RR= +", format(rr, digits = 0, nsmall = 2), co2star), 
                     paste0("RR= ",  format(rr, digits = 0, nsmall = 2), co2star)))


# fig
dodgeval <- .3
fig_domspp <- ggplot(ci_dd, aes(x = year, y = rlsmean)) +
  
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  facet_wrap( ~ .id) +
  
  
  # observed
  geom_point(data = obs_d, 
             aes(x = year, y = value, shape = co2, group = co2, col = value_type),
             size = 2, fill = "grey80", position = position_dodge(dodgeval)) +
  
  
  # adjusted
  geom_errorbar(aes(ymin = rlowerCL, ymax = rupperCL, shape = co2, group = co2, col = value_type), 
                width = 0, position = position_dodge(width = dodgeval)) +
  geom_line(aes(linetype = co2, shape = co2, group = co2, col = value_type), 
            position = position_dodge(width = dodgeval)) +
  geom_point(aes(shape = co2, group = co2, col = value_type),size = 2.5, position = position_dodge(width = dodgeval)) +
  geom_text(aes(label = star, y = rupperCL), fontface = "bold", vjust = 0) +
  
  
  # scaling
  scale_shape_manual(values = c(16, 17), 
                    labels = c("Ambient", expression(eCO[2]))) +
  scale_linetype_manual(values = c("solid", "dashed"), 
                        labels = c("Ambient", expression(eCO[2]))) +
  scale_color_manual(values = c("black", "grey80"),
                     guide = guide_legend(override.aes = list(linetype = "blank",size = 2))) +
  scale_y_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30)) +
  scale_x_discrete("Year", labels = 0:4, drop = FALSE) +
  
  
  # theme and legend
  science_theme +
  theme(legend.position   = "top",
        legend.box        = "horizontal", 
        legend.direction  = "vertical", 
        legend.text.align = 0, 
        strip.text.x      = element_text(face = "italic")) +
  
  labs(y = expression(Abundance~(Counts~m^'-1'))) +
  
  geom_text(data = plab_d, aes(label = plot_lab), x = -Inf, y = Inf, hjust = -.1, 
            fontface = "bold", vjust = 1.5, size = 3) +
  geom_text(data = plab_d, aes(label = rr), x =  Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3)


fig_domspp

ggsavePP(filename = "output/figs/adjusted_abundance_dominenetSPP", 
         plot = fig_domspp, width = 4.5, height = 4.5)


# summary table -----------------------------------------------------------

# table for observed values
obs_tbl <- obs_d %>% 
  select(.id, year, co2, value_type, ring, value) %>% 
  group_by(.id, year, co2, value_type) %>% 
  summarise(M = mean(value))


# bind with adjusted values
domsp_adjMean_tble <- ci_dd %>% 
  select(.id, co2, year, rlsmean, value_type) %>% 
  rename(M = rlsmean) %>% 
  bind_rows(obs_tbl) %>% 
  mutate(variable = paste(value_type, co2, sep = "_")) %>% 
  select(-value_type, -co2) %>% 
  spread(key  = variable, value = M) %>% 
  mutate(resp = adjusted_elev / adjusted_amb - 1) %>% 
  group_by(.id, year) %>% 
  summarise_each(funs(round(., 2)), everything(), -.id) %>% 
  select(year, starts_with("observed"), starts_with("adjusted"), resp) %>% 
  ungroup()


# split df by .id
domsp_tbl_l <- dlply(domsp_adjMean_tble, .(.id), function(x) select(x, -.id))


# save as excel
writeWorksheetToFile(file  = "output/table/summary_tbl_dom_spp.xlsx",  # define file name to be saved
                     data  = domsp_tbl_l,                              # writeWorksheetToFile doesn't take dplyr object so turn them into data frames using as.data.frame
                     sheet = names(domsp_tbl_l))                       # sheet names in excel are defined by object names a list

