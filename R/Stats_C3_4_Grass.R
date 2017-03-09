
# prepare data frames -----------------------------------------------------

# Create C3grass:C4, legume:non-legume and Native:introduced
C3grassC4 <- veg_FullVdf %>% 
  filter(form == "Grass") %>% 
  mutate(yval = factor(ifelse(PFG == "c4", "p", "q")))

  
legumeR <- veg_FullVdf %>% 
  filter(form == "Forb") %>% 
  mutate(yval = factor(ifelse(PFG == "legume", "p", "q")))

NativeR <- veg_FullVdf %>% 
  filter(!is.na(origin)) %>% 
  mutate(yval = factor(ifelse(origin == "native", "p", "q")))
summary(NativeR)  

dfList <- list('C3vsC4'      = C3grassC4, 
               'LegvsNonleg' = legumeR, 
               'NatvsIntr'   = NativeR)

# compute ratios and total number
PfgRDF <- llply(dfList, function(x){
  x %>% 
    group_by(year, co2, block, ring, plot, id, RY) %>%
    summarise(Total = sum(value),
              ratios = sum(value[yval == "p"]/Total)) %>%
    ungroup() # grouping informaiton is not required later
  })

# Move Year0 value to a new column to be used as covariate in the analyssis
# below
PfgRDF_year0 <- llply(PfgRDF, function(x) {

  newyear_dd <- x %>% # year0
    filter(year == "Year0") %>%
    select(id, ratios) %>%
    rename(ratios0 = ratios) %>% 
    left_join(filter(x, year != "Year0"), by = "id")  # merge with subsequent years
  
  newyear_dd$obs <- 1:nrow(newyear_dd) # add id
  return(newyear_dd)
})

# assess lineariry
plot(logit(ratios) ~ ratios0, pch = 19, col = year, data = PfgRDF_year0[[1]])
plot(logit(ratios) ~ sqrt(ratios0), pch = 19, col = year, data = PfgRDF_year0[[2]])
plot(logit(ratios) ~ sqrt(ratios0), pch = 19, col = year, data = PfgRDF_year0[[3]])

PfgRDF_year0[[1]]$t_value0 <- PfgRDF_year0[[1]]$ratios0
PfgRDF_year0[[2]]$t_value0 <- sqrt(PfgRDF_year0[[2]]$ratios0)
PfgRDF_year0[[3]]$t_value0 <- sqrt(PfgRDF_year0[[3]]$ratios0)


# Aanlysis ----------------------------------------------------------------


# models to be tested
pfgprop_m_list <- llply(PfgRDF_year0, function(x) {
  m1 <- lmer(logit(ratios) ~ co2 * year + t_value0 + (1|block) + (1|ring) + (1|id) + (1|RY), data = x)
  m2 <- lmer(logit(ratios) ~ co2 * year + t_value0 + (1|ring) + (1|id) + (1|RY), data = x)
  if (AICc(m1) >= AICc(m2)) return(m2) else return(m1)
})


llply(pfgprop_m_list, function(x) Anova(x, test.statistic = "F"))

m1 <- pfgprop_m_list[[1]]
m2 <- update(m1, ~ . -co2:year)
m3 <- update(m2, ~ . -year)
Anova(m3, test.statistic = "F")

# model diagnosis ---------------------------------------------------------


pdf("output/figs/mod_diag_PFGprop.pdf", onefile = TRUE, width = 4, height = 4)
l_ply(names(pfgprop_m_list), function(x){
  m <- pfgprop_m_list[[x]]
  print(plot(m, main = x))
  qqnorm(resid(m, main = x))
  qqline(resid(m, main = x))
})
dev.off()


# insepct legume
d <- PfgRDF_year0[["LegvsNonleg"]]
m1 <- lmer(logit(ratios) ~ co2 * year + t_value0 + (1|block) + (1|ring) + (1|id), data = d)
which.min(resid(m1))
m2 <- update(m1, subset = -6)
llply(list(m1, m2), function(x) Anova(x, test.statistic = "F"))
# no difference so leave it as it is

# inspect Native
d <- PfgRDF_year0[["NatvsIntr"]]
m1 <- lmer(logit(ratios) ~ co2 * year + t_value0 + (1|block) + (1|ring) + (1|id), data = d)
which.max(resid(m1))
m2 <- update(m1, subset = -72)
plot(m2)
llply(list(m1, m2), function(x) Anova(x, test.statistic = "F"))
# not much difference, so leave it as it is




# CI and postdoc test -----------------------------------------------------


# compute 95 CI and post-hoc test
lsmeans_list <- llply(pfgprop_m_list, function(x) {
  lsmeans::lsmeans(x, ~ co2 | year)
})

# 95% CI
CI_dd <- ldply(lsmeans_list, function(x) data.frame(summary(x))) 

# post-hoc test
contrast_dd <- ldply(lsmeans_list, function(x) {
  data.frame(summary(pairs(x)[1:3], adjust = "none"))
}) %>% 
  mutate(co2  = factor("elev", levels = c("amb", "elev")),
         star = get_star(p.value)) %>% 
  select(.id, year, co2, p.value, star)
contrast_dd$star[contrast_dd$.id == "C3vsC4"] <- ""  # no co2:year interaction

# CO2 effect
pfgprop_aov_df <- ldply(pfgprop_m_list, function(x) tidy(Anova(x, test.statistic = "F")),  # Anova result of models
                    .progress = "text")
pfgprop_co2_pval <- pfgprop_aov_df %>%                                                     # get p-values for CO2 term
  filter(term == "co2") %>% 
  mutate(co2star = get_star(p.value)) %>% 
  select(.id, co2star)

# merge
pfgprop_ci_dd <- left_join(CI_dd, contrast_dd, by = c(".id", "year", "co2")) %>% 
  mutate(year     = factor(year, levels = paste0("Year", 0:3)),
         rlsmean  = boot::inv.logit(lsmean),                                      # reverse transform and standardise for 1mx1m plot
         rlowerCL = boot::inv.logit(lower.CL),
         rupperCL = boot::inv.logit(upper.CL)) %>% 
  left_join(pfgprop_co2_pval, by = ".id")                                         # merge with pvalues for CO2 term
pfgprop_ci_dd$star[is.na(pfgprop_ci_dd$star)] <- ""

# df for observed values
pfgprop_d <- ldply(dfList) %>%
  group_by(year, .id, co2, block, ring) %>%
  summarise(Total = sum(value),
            ratios = sum(value[yval == "p"]/Total)) %>%
  ungroup() %>% # grouping informaiton is not required later
  mutate(year = factor(year, levels = paste0("Year", 0:3)))

