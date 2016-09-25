
# prepare data frame ------------------------------------------------------

# create models to be tests
dominentsp <- veg_FullVdf %>% 
  filter(variable %in% DmSpp) %>% 
  group_by(variable, year, block, co2, ring, plot, id) %>% 
  summarise(value = sum(value)) %>% 
  ungroup()

dominentsp_year0 <- dominentsp %>%
  filter(year == "Year0") %>%
  select(id, value, variable) %>%
  rename(value0 = value) %>% 
  left_join(filter(dominentsp, year != "Year0"), by = c("id", "variable")) %>% 
  mutate(logitv  = logit(value), 
         logitv0 = logit(value0),
         logv0   = log(value0 + 1),
         sqrtv0  = sqrt(value0))




# analysis ----------------------------------------------------------------


# independent variable transformation
par(mfrow = c(2, 4), mar = c(2, 2, 1, 1))

# raw
d_ply(dominentsp_year0, .(variable), 
      function(x) plot(logit(value) ~ value0, data = x, pch = 19, col = factor(year), 
                       main = "raw"))  

# sqrt
d_ply(dominentsp_year0, .(variable), 
      function(x) plot(logit(value) ~ sqrtv0, data = x, pch = 19, col = factor(year), 
                       main = "sqrt"))  

# log
d_ply(dominentsp_year0, .(variable), 
      function(x) plot(logit(value) ~ logv0, data = x, pch = 19, col = factor(year), 
                       main = "log"))  

# logit
d_ply(dominentsp_year0, .(variable), 
      function(x) plot(logit(value) ~ logitv0, data = x, pch = 19, col = factor(year), 
                       main = "logit"))  


dom_m_list <- dlply(dominentsp_year0, .(variable), function(x){
  m1 <- lmer(logit(value) ~ co2 * year + value0  + (1|block) + (1|ring) + (1|id), data = x)
  m2 <- lmer(logit(value) ~ co2 * year + sqrtv0  + (1|block) + (1|ring) + (1|id), data = x)
  m3 <- lmer(logit(value) ~ co2 * year + logv0   + (1|block) + (1|ring) + (1|id), data = x)
  m4 <- lmer(logit(value) ~ co2 * year + logitv0 + (1|block) + (1|ring) + (1|id), data = x)
  ml <- list(m1, m2, m3, m4)
  bm <- ml[[which.min(AICc(m1, m2, m3, m4)$AICc)]]
  m2 <- update(bm, ~ . - (1|block))
  if (AICc(bm) >= AICc(m2)) return(m2) else return(m1)
})




# model diagnosis ---------------------------------------------------------

pdf("output/figs/mod_diag_domspp.pdf", onefile = TRUE, width = 4, height = 4)
l_ply(names(dom_m_list), function(x){
  m <- dom_m_list[[x]]
  print(plot(m, main = x))
  qqnorm(resid(m, main = x))
  qqline(resid(m, main = x))
})
dev.off()





# CI and post-hoc test ----------------------------------------------------


# compute 95 CI and post-hoc test
lsmeans_list <- llply(dom_m_list, function(x) {
  lsmeans::lsmeans(x, ~ co2 | year)
})


# 95% CI
CI_dd <- ldply(lsmeans_list, function(x) data.frame(summary(x))) 


# post-hoc test
contrast_dd <- ldply(lsmeans_list, function(x) {
  data.frame(summary(pairs(x)[1:3], adjust = "fdr"))
}) %>% 
  mutate(co2 = factor("amb", levels = c("amb", "elev")),
         star = cut(p.value, right = FALSE,
                    breaks = c(0, .1, .05, .01, .001, 1),  
                    labels = c("***", "**", "*", "\u2020", ""))) %>% 
  select(.id, year, co2, p.value, star)


# merge
ci_dd <- left_join(CI_dd, contrast_dd, by = c(".id", "year", "co2")) %>% 
  mutate(year = factor(year, levels = paste0("Year", 0:3)),
         rlsmean = boot::inv.logit(lsmean) * 25, # reverse transform and standardise for 1mx1m plot
         rlowerCL = boot::inv.logit(lower.CL) * 25,
         rupperCL = boot::inv.logit(upper.CL) * 25,
         .id = gsub("[.]", " ", .id),
         value_type = "adjusted")
ci_dd$star[is.na(ci_dd$star)] <- ""


# create fig --------------------------------------------------------------


# df for observed vlaues
obs_d <- dominentsp %>%
  group_by(year, co2, ring, variable) %>% 
  summarise(N = sum(!is.na(value)), value = mean(value)) %>% 
  ungroup() %>% 
  mutate(year = factor(year, levels = paste0("Year", 0:3)),
         value = value/4,                                    # standardise for 1mx1m plot
         .id = gsub("[.]", " ", variable),
         value_type = "observed")

# fig
dodgeval <- .3
fig_domspp <- ggplot(ci_dd, aes(x = year, y = rlsmean, shape = co2, group = co2, 
                                col = value_type)) +
  
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  facet_wrap( ~ .id) +
  
  
  # observed
  geom_point(data = obs_d, aes(x = year, y = value),  size = 2, fill = "grey80",
               position = position_dodge(dodgeval)) +
  
  
  # adjusted
  geom_errorbar(aes(ymin = rlowerCL, ymax = rupperCL), width = 0, 
                position = position_dodge(width = dodgeval)) +
  geom_line(aes(linetype = co2), position = position_dodge(width = dodgeval)) +
  geom_point(size = 2.5, position = position_dodge(width = dodgeval)) +
  geom_text(aes(label = star, y = rupperCL), fontface = "bold", vjust = 0) +
  
  
  # scaling
  scale_shape_manual(values = c(16, 17), 
                    labels = c("Ambient", expression(eCO[2]))) +
  scale_linetype_manual(values = c("solid", "dashed"), 
                        labels = c("Ambient", expression(eCO[2]))) +
  scale_color_manual(values = c("black", "grey80"),
                     guide = guide_legend(override.aes = list(linetype = "blank",size = 2))) +
  scale_y_continuous(limits = c(0, 25)) +
  scale_x_discrete("Year", labels = 0:4, drop = FALSE) +
  
  
  # theme and legend
  science_theme +
  theme(legend.position   = "top",
        legend.box        = "horizontal", 
        legend.direction  = "vertical", 
        legend.text.align = 0, 
        strip.text.x      = element_text(face = "italic")) +
  
  labs(y = expression(Abundance~(Counts~m^'-1')))
fig_domspp

ggsavePP(filename = "output/figs/adjusted_abundance_dominenetSPP", 
         plot = fig_domspp, width = 5, height = 5.5)


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

