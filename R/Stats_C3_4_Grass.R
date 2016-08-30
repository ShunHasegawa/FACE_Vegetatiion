
# prepare data frames -----------------------------------------------------

# Create C3grass:C4, legume:non-legume and Native:introduced
C3grassC4 <- veg_FullVdf %>% 
  filter(form == "Grass") %>% 
  mutate(yval = factor(ifelse(PFG == "c3", "p", "q")))
  
legumeR <- veg_FullVdf %>% 
  filter(form == "Forb") %>% 
  mutate(yval = factor(ifelse(PFG == "legume", "p", "q")))

NativeR <- veg_FullVdf %>% 
  filter(!is.na(origin)) %>% 
  mutate(yval = factor(ifelse(origin == "native", "p", "q")))
summary(NativeR)  

dfList <- list('C3vsC4' = C3grassC4, 
               'LegvsNonleg' = legumeR, 
               'NatvsIntr' = NativeR)

# compute ratios and total number
PfgRDF <- llply(dfList, function(x){
  x %>% 
    group_by(year, co2, block, ring, plot, id) %>%
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
  m1 <- lmer(logit(ratios) ~ co2 * year + t_value0 + (1|block) + (1|ring) + (1|id), data = x)
  m2 <- update(m1, ~ . - (1|block))
  if (AICc(m1) >= AICc(m2)) return(m2) else return(m1)
}
)

# compute 95 CI and post-hoc test
lsmeans_list <- llply(pfgprop_m_list, function(x) {
  summary(lsmeans::lsmeans(x, pairwise ~ co2 | year))
})

# 95% CI
CI_dd <- ldply(lsmeans_list, function(x) data.frame(x$lsmeans)) 

# post-hoc test
contrast_dd <- ldply(lsmeans_list, function(x) data.frame(x$contrast)) %>% 
  mutate(co2  = factor("elev", levels = c("amb", "elev")),
         star = get_star(p.value)) %>% 
  select(.id, year, co2, p.value, star)

# merge
pfgprop_ci_dd <- left_join(CI_dd, contrast_dd, by = c(".id", "year", "co2")) %>% 
  mutate(year     = factor(year, levels = paste0("Year", 0:3)),
         rlsmean  = boot::inv.logit(lsmean), # reverse transform and standardise for 1mx1m plot
         rlowerCL = boot::inv.logit(lower.CL),
         rupperCL = boot::inv.logit(upper.CL))
pfgprop_ci_dd$star[is.na(pfgprop_ci_dd$star)] <- ""

# df for observed values
pfgprop_d <- ldply(dfList) %>%
  group_by(year, .id, co2, block, ring) %>%
  summarise(Total = sum(value),
            ratios = sum(value[yval == "p"]/Total)) %>%
  ungroup() %>% # grouping informaiton is not required later
  mutate(year = factor(year, levels = paste0("Year", 0:3)))

