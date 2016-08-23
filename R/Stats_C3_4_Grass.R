
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

dfList <- list('C3 in grass: C3 vs. C4' = C3grassC4, 
               'Legume in forb: Legume vs. Non-legume' = legumeR, 
               'Native plants: native vs. introduced' = NativeR)

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
m_list <- llply(PfgRDF_year0, function(x) {
  m1 <- lmer(logit(ratios) ~ co2 * year + t_value0 + (1|block) + (1|ring) + (1|id), data = x)
  m2 <- update(m1, ~ . - (1|block))
  if (AICc(m1) >= AICc(m2)) return(m2) else return(m1)
}
)

# compute 95 CI and post-hoc test
lsmeans_list <- llply(m_list, function(x) {
  summary(lsmeans::lsmeans(x, pairwise ~ co2 | year))
})

# 95% CI
CI_dd <- ldply(lsmeans_list, function(x) data.frame(x$lsmeans)) 

# post-hoc test
contrast_dd <- ldply(lsmeans_list, function(x) data.frame(x$contrast)) %>% 
  mutate(co2 = factor("elev", levels = c("amb", "elev")),
         star = cut(p.value, right = FALSE,
                    breaks = c(0, .1, .05, .01, .001, 1),  
                    labels = c("***", "**", "*", "\u2020", ""))) %>% 
  select(.id, year, co2, p.value, star)

# merge
ci_dd <- left_join(CI_dd, contrast_dd, by = c(".id", "year", "co2")) %>% 
  mutate(year = factor(year, levels = paste0("Year", 0:3)),
         rlsmean = boot::inv.logit(lsmean), # reverse transform and standardise for 1mx1m plot
         rlowerCL = boot::inv.logit(lower.CL),
         rupperCL = boot::inv.logit(upper.CL))
ci_dd$star[is.na(ci_dd$star)] <- ""

# fig ---------------------------------------------------------------------

# df for Year0
d <- ldply(dfList) %>%
  filter(year == "Year0") %>% 
  group_by(year, .id, co2, block, ring) %>%
  summarise(Total = sum(value),
            ratios = sum(value[yval == "p"]/Total)) %>%
  ungroup() %>% # grouping informaiton is not required later
  mutate(year = factor(year, levels = paste0("Year", 0:3)))

dodgeval <- .4
p <- ggplot(ci_dd, aes(x = year, y = rlsmean, fill = co2, group = co2)) +
  facet_wrap(~ .id, scales = "free_y", ncol = 1) +
  
  geom_errorbar(aes(ymin = rlowerCL, ymax = rupperCL), width = 0, 
                position = position_dodge(width = dodgeval)) +
  geom_line(aes(linetype = co2), position = position_dodge(width = dodgeval)) +
  geom_point(data = d, aes(x = year, y = ratios),  alpha = .7, shape = 21, size = 3,
             position = position_dodge(dodgeval), show.legend = FALSE) +
  geom_point(shape = 21, size = 3, position = position_dodge(width = dodgeval)) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_text(aes(label = star, y = rupperCL), fontface = "bold", vjust = .4) +
  
  scale_fill_manual(values = c("black", "white"), 
                    labels = c("Ambient", expression(eCO[2]))) +
  scale_linetype_manual(values = c("solid", "dashed"), 
                        labels = c("Ambient", expression(eCO[2]))) +
  scale_color_manual(values = "grey50", labels = expression(Md[Year0])) +
  scale_x_discrete("Year", labels = 0:4, drop = FALSE) +
  
  science_theme +
  theme(legend.position = c(.7, .75),
        strip.text.x    = element_text(size = 8)) +
 
  labs(y = expression("Adjusted proportion"))
p

ggsavePP(filename = "output/figs/adjusted_PFG_proportion", width = 3, height = 6,
         plot = p)
