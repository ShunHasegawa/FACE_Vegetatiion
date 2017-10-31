DivDF_list2 <- llply(vegDF_list, function(x) {
  siteDF %>% 
    mutate( 
      H = diversity(x),  # Shannon's index
      S = specnumber(x), # number of spp
      J = H/log(S)       # Pielou's evenness
    ) %>% 
    group_by(year, ring, co2, plot) %>% 
    summarise_each(funs(mean(., na.rm = TRUE)), H, S, J) %>% 
    ungroup()
})



DivDF_year0_list2 <- llply(DivDF_list2, function(x){
  
  DivDF_mlt <- gather(x, variable, value, H, S, J)
  
  DivDF_year0 <- DivDF_mlt %>% # Year0
    filter(year == "Year0") %>%
    select(ring, plot, value, variable) %>%
    rename(value0 = value) %>% 
    left_join(filter(DivDF_mlt, year != "Year0"), by = c("ring", "plot","variable")) %>% 
    filter(!is.na(value)) %>% 
    mutate(value0_log  = log(value0 + 1),
           value0_sqrt = sqrt(value0),
           cy0 = value - value0) %>% 
    group_by(year, ring, co2, variable) %>% 
    summarise_each(funs(mean), cy0) %>% 
    ungroup() %>% 
    spread(variable, cy0)
  
  return(DivDF_year0)
})

boxplot(sqrt(H + 1) ~ co2 * year, data = DivDF_year0_list2$grass_spp)
boxplot(J ~ co2 * year, data = DivDF_year0_list2$grass_spp)
boxplot(sqrt(S + 10) ~ co2 * year, data = DivDF_year0_list2$grass_spp)

mm1 <- lmer(H ~ co2 * year + (1|ring), data = DivDF_year0_list2$grass_spp)
mm1 <- lmer(J ~ co2 * year + (1|ring), data = DivDF_year0_list2$grass_spp)
mm1 <- lmer(S ~ co2 * year + (1|ring), data = DivDF_year0_list2$grass_spp)
r.squared(mm1)
r.squared(gs_m1)
mm1 <- lmer(sqrt(S + 10) ~ co2 * year + (1|ring), data = DivDF_year0_list2$grass_spp)
r.squared(mm1)
VarCorr(mm1)
Anova(mm1, test.statistic = "F")
summary(pairs(lsmeans(mm1, ~ co2 | year))[1:3], adjust = "none")

plot(mm1)
qqPlot(resid(mm1))
visreg(mm1, xvar = "year", by = "co2", overlay = TRUE)



boxplot(I(value - value0) ~ co2 * year, data = gh_d)
mm2 <- lmer(I(value - value0) ~ co2 * year+ (1|ring), data = gh_d)






boxplot(I(log(value - value0 + 20)) ~ co2 * year, sc4_dd)
mmm2 <- lmer(I(log(value - value0 + 20)) ~ co2 * year + (1|ring), sc4_dd)
plot(mmm2)
qqPlot(resid(mmm2))
Anova(mmm2, test.statistic = "F")
summary(pairs(lsmeans(mmm2, ~ co2 | year))[1:3], adjust = "none")
VarCorr(mmm2)





library(nlme)
library(varComp)
data(Oxide)
lmef = lme(Thickness~Source, Oxide, ~1|Lot/Wafer)
vcf = varComp(Thickness~Source, Oxide, ~Lot/Wafer)
VarCorr(lmef)
coef(vcf, 'varComp') ## same values as above
varComp.test(vcf, test = "RLRT") ## test against linear model
varComp.test(vcf, null=1) ## test against model with only Lot random effect







