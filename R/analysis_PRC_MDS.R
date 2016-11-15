
#  Prepare dataframe ------------------------------------------------------

prc_sp   <- decostand(PlotSumVeg[, SppName], method = "log") # log-transformed sp df
prc_site <- PlotSumVeg %>%                                   # site df 
  select(year, ring, plot, co2) %>% 
  mutate(id = ring:plot)


# analysis ----------------------------------------------------------------


# > PRC ---------------------------------------------------------------------


## principal response curve anlsyis with bray-curtis dissimilarity

prc_all <- capscale(prc_sp ~  year * co2 + Condition(year), data = prc_site, 
                    distance = "bray")




# . permutation test ------------------------------------------------------


## define permutation for this experiment; year is nested within plot within
## ring.

## 1. plot is a random factor, should never be exchanged. 
## 2. year is a time series measurement. This can be exchanged, but its order 
## should be kept to take autocorrelation into account (i.e. [1, 2, 3], [2, 3, 
## 1], [3, 1, 2]). Also this order should be consistent among all plots.
## 3. ring is an experimental unit so fully exchangable


## 1) define permutation for year within each plot

cntrol_year <- how(within = Within(type = "series"),                     # year within plot is a time seiries
                   plot   = Plots(strata = prc_site$id, type = "none"),  # plot should not be exchanged
                   nperm = 4999)  
perm_year   <- shuffleSet(nrow(prc_site), control = cntrol_year)
nrow(perm_year)

## 2) the above permutation did not exchange between rings. do so as below

perm_year_or <- rbind(1:nrow(prc_site), perm_year)                               # add the original (i.e. 1:96) as this should be also permuted between rings
perm_ring_l  <- llply(levels(prc_site$ring), function(x) perm_year_or[, prc_site$ring == x])  # split each row of the permutation by ring into a list
perm6        <- rbind(1:6, allPerms(6))                                          # potential permutation for 6 rings (i.e. 6!)
perm_by_ring <- alply(perm6, 1, function(x) do.call(cbind, perm_ring_l[x]))      # permute by ring, exchanging perm_year_or between rings (e.g. origin[1, 2, 3, 4, 5, 6], 2nd[1, 2, 3, 4, 6, 5] etc.)
comp_perm    <- do.call(rbind, perm_by_ring)[-1, ]                               # row bind all results and remove the original
nrow(comp_perm)


## check permutation; ring is freely exchanged. Plot is not exchanged. Year is
## exchanged within plots, but the order (time series) is preserved

llply(split(comp_perm[sample(nrow(comp_perm), 1), ], prc_site$ring), function(x){
  m <- matrix(x, ncol = 4)
  apply(m, c(1, 2),  function(x) paste(prc_site$id, prc_site$year, sep = "-")[x])
})


# run permutation test
prc_res <- anova(prc_all, permutations = comp_perm[sample(4999), ], by = "axis")
prc_res <- anova(prc_all, permutations = alperms_bind[sample(4999), ], by = "axis")
prc_res <- anova(prc_all, permutations = alperms_bind[sample(999), ], by = "axis")
prc_res


class(summary_prc$sites[, "CAP1"])

?anova.cca
prc_res <- anova(prc_all, permutations = comp_perm[sample(999), ], by = "axis")
prc_res

EnvDF2 <- ldply(1:4, function(x) data.frame(EnvDF_3df, plot = x)) %>% 
  arrange(ring, plot, year) %>% 
  mutate(plot =factor(plot))
identical(data.frame(prc_site[, c("year","ring", "plot")]), 
          data.frame(EnvDF2[, c("year", "ring", "plot")]))

envf2 <- envfit(res_prc_site[, c("MDS1", "MDS2")], 
               EnvDF2[, c("Depth_HL", "TotalC", "Drysoil_ph", "moist", 
                             "gapfraction", "temp", "sand", "silt", "clay")],
               permutations = alperms_bind[sample(4999), ])
envf
envf2
?envfit
# . summarise -------------------------------------------------------------


## prepare df for a figure from the result of prc
summary_prc       <- summary(prc_all, scaling = 3)                            # use scaling = 3, it doesn't matter too much, but this is used in the original prc function
MDSprop           <- summary_prc$cont$importance[2, c("MDS1", "MDS2")]        # proportion explained by the first 2 MDS axies
MDSaxes           <- paste0("MDS", 1:2, " (", round(MDSprop * 100, 0), "%)")  # MDS axis with explained proportion
res_prc_site      <- cbind(summary_prc$sites, prc_site)                       # site score
res_prc_site_ring <- res_prc_site %>%                                         # site scores by ring
  group_by(year, ring, co2) %>% 
  summarise_each(funs(mean), CAP1, MDS1, MDS2) %>%                            # CAP1 is for PRC plot and MDSs for visualisation
  arrange(ring, year)



# . fit envronmental variable ---------------------------------------------


load("output/Data/EucFACE_understorey_env_vars_2012-2016.RData")                  # load environmental variables, generated in "R/FitEnvironmentalVars.R"
identical(EnvDF_3df[, c("year", "ring")],                                         # check EnvDf_3df and species df are in the same order for site
          data.frame(res_prc_site_ring[, c("year", "ring")]))  
cntr       <- how(within = Within(type = "series"),                               # define permutatoin. ring is exchangable; autocorrelation is taken into consideration by series
                  plots = Plots(strata = res_prc_site_ring$ring, type = "free"),
                  nperm = 4999)
mds_envfit <- envfit(res_prc_site_ring[, c("MDS1", "MDS2")],                      # fit environmental variables
                     EnvDF_3df[, c("Depth_HL", "TotalC", "Drysoil_ph", "moist", 
                                   "gapfraction", "temp", "sand", "silt", "clay")],
                     permutations = cntr)


mds_arrw_d <- data.frame(scores(mds_envfit, "vectors"), 
                         pval = mds_envfit$vector$pvals) %>%                      # arrows for environmental variables 
  mutate(env = recode(row.names(.), 
                      Depth_HL    = "HL", 
                      TotalC      = "Total C",
                      Drysoil_ph  = "pH",
                      moist       = "Moist",
                      gapfraction = "Light",
                      temp        = "Temp",
                      sand        = "Sand",
                      silt        = "Silt",
                      clay        = "Clay"),
         co2 = factor("amb", levels = c("amb", "elev")),
         year = factor("Year0", levels = paste0("Year", 0:3))) %>% 
  filter(pval <= 0.1)

res_prc_site_co2  <- res_prc_site_ring %>%                                    # canonical coefficients for CO2 treatment (i.e. treatment difference for each year)
  group_by(year, co2) %>% 
  summarise_each(funs(mean), CAP1) %>% 
  group_by(year) %>% 
  summarise_each(funs(-diff(.)), CAP1) %>% 
  mutate(co2 = "elev") %>% 
  bind_rows(data.frame(year = levels(.$year), CAP1 = 0, co2 = "amb")) %>% 
  mutate(fig = "site")



# figure ------------------------------------------------------------------


# ggplot theme for prc
science_theme_prc <- science_theme +
  theme(legend.title      = element_text(size = 8),
        legend.text       = element_text(size = 8),
        legend.key.width  = unit(1.7, "lines"),
        axis.title        = element_text(size = 9),
        axis.text         = element_text(size = 8))


# . canonical coefficient -------------------------------------------------


fig_prc_site <- ggplot(res_prc_site_co2, aes(x = as.numeric(year), y = CAP1, 
                                             fill = co2, linetype = co2)) +
  labs(x = "", y = expression(Canonical~coefficient~(C[dt]))) +
  geom_hline(yintercept = 0, linetype = "dotted", col = "gray") +
  
  geom_line() +
  geom_point(size = 2, shape = 21) +
  
  science_theme_prc +
  theme(legend.position = "none") +
  scale_x_continuous(labels = paste0("Year", 0:3)) +
  scale_fill_manual(values = c("black", "white"),
                     labels = c("Ambient", expression(CO[2]))) +
  scale_linetype_manual(values = c("solid", "dashed"), 
                        labels = c("Ambient", expression(CO[2])))
fig_prc_site




# . Sp weight plot -----------------------------------------------------------


## merge sp score with pfg and propotion in Year0


### sp score given by PRC
res_prc_sp <- data.frame(summary_prc$species,
                         variable = row.names(summary_prc$species))

### PFG, form table
pfg_spp_tbl <- veg_FullVdf %>%
  select(variable, form, PFG, origin) %>% 
  distinct() %>% 
  mutate(form    = as.character(form),
         PFG     = as.character(PFG),
         pfgform = ifelse(form == "Grass", paste0(PFG, form),
                          ifelse(form == "Forb", PFG, form))) %>% 
  select(variable, pfgform, origin)


## proportion for spp, PFG or origin in Year0
year0_df        <- PlotSumVeg %>%                                               # df for Year0
  filter(year == "Year0") %>% 
  gather(variable, value, one_of(SppName)) %>% 
  group_by(variable) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  left_join(pfg_spp_tbl, by = "variable")                                       # merge with PFG and forb table
year0_sp_prop   <-  transmute(year0_df, variable, spprop = value / sum(value))  # sp proportion in Year0
year0_pfg_prop  <- year0_df %>%                                                 # PFG proportion in Year0
  group_by(pfgform) %>%  
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  transmute(pfgform, pfgprop = value / sum(value))
year0_orgn_prop <- year0_df %>%                                                 # origin proportion in Year0
  group_by(origin) %>% 
  filter(!is.na(origin)) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  transmute(origin, orgnprop = value / sum(value))


### merge the above dfs
res_pric_sp_d <- res_prc_sp %>%
  left_join(pfg_spp_tbl    , by = "variable") %>% 
  left_join(year0_sp_prop  , by = "variable") %>% 
  left_join(year0_pfg_prop , by = "pfgform") %>% 
  left_join(year0_orgn_prop, by = "origin")
summary(res_pric_sp_d)


### re-order pfgform by median
pfgform_lev <- res_pric_sp_d %>% 
  group_by(pfgform) %>% 
  summarise(med = median(CAP1)) %>% 
  arrange(-med) %>% 
  .$pfgform
res_pric_sp_d$pfgform <- factor(res_pric_sp_d$pfgform, levels = pfgform_lev)



## create a plot
head(res_pric_sp_d)
res_pric_sp_d %>% 
  select(pfgform, pfgprop) %>% 
  distinct()
boxplot(CAP1 ~ pfgform, width = unique(res_pric_sp_d$pfgprop), data = res_pric_sp_d)
hist(log(res_pric_sp_d2$spprop, base = 10))
hist(res_pric_sp_d2$spprop)
nrow(res_pric_sp_d)
res_pric_sp_d2 <- res_pric_sp_d %>%
  gather(measure, value = type, pfgform, origin) %>% 
  mutate(measure = factor(measure, levels = c("pfgform", "origin"), labels = c("PFG", "Origin")),
         type    = factor(type, levels = c(pfgform_lev, "native", "naturalised")),
         type    = mapvalues(type, c(pfgform_lev, "native", "naturalised"),
                             c("C4 grass", "Fern", "C3 grass", "Forb",  "Moss", 
                               "Wood", "Legume", "Native", "Introduced")),
         year    = factor("Year0", levels = paste0("Year", 0:4)),
         co2     = factor("Ambient", levels = c("Ambient", "eCO[2]"))) %>% 
  filter(!is.na(type))




fig_prc_spp_byPfg <- ggplot(res_pric_sp_d2, aes(x = type, y = CAP1)) +
  facet_grid(. ~ measure, space = "free_x", scale = "free_x") +
  geom_hline(yintercept = 0, linetype = "dotted", size = .5) +
  geom_boxplot(outlier.size = 1, na.rm = TRUE) +
  geom_jitter(aes(size = spprop), alpha = .5, width = .6, col = "gray20") +
  scale_size_continuous("Proportion\nin Year0", range = c(1, 5), 
                        breaks = c(0, 0.05, 0.1, 0.3)) +
  science_theme_prc +
  theme(legend.position  = "right",
        axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)) +
  labs(y = "Species weight", x = "") +
  ylim(c(-.5, .5))
fig_prc_spp_byPfg





# . PCoA (MDS) plot -------------------------------------------------------


fig_prc_mds <- ggplot(res_prc_site_ring, aes(x = MDS1, y = MDS2, shape = year, 
                                             fill = co2)) +
  
  labs(x = MDSaxes[1], y = MDSaxes[2]) +
  geom_hline(yintercept = 0, linetype = "dotted", col = "gray") +
  geom_vline(xintercept = 0, linetype = "dotted", col = "gray") +
  
  geom_path(aes(group = ring, linetype = co2)) +
  geom_point(size = 2, alpha = .6) +
  
  geom_segment(data = mds_arrw_d,
               aes(x = 0, y = 0, xend = MDS1 * .7, yend = MDS2 * .7),
               arrow = arrow(length = unit(.2, "cm")),  alpha = .7) +
  geom_text(data = mds_arrw_d,  aes(x = MDS1 * .8, y = MDS2* .8, label = env),  
            size = 2, fontface = "italic") +
  

  science_theme_prc +
  theme(legend.position = "right") +
  scale_fill_manual(name   = expression(CO[2]),
                    values = c("black", "white"), 
                    labels = c("Ambient", expression(eCO[2])),
                    guide  = guide_legend(override.aes = list(shape = 21, alpha = 1))) +
  scale_shape_manual(name = "Year", 
                     values = c(21, 22, 23, 24)) +
  scale_linetype_manual(name = expression(CO[2]),
                        values = c("solid", "dashed"), 
                        labels = c("Ambient", expression(eCO[2])))
fig_prc_mds




# . merge figs ------------------------------------------------------------


prc_figs <- list(fig_prc_site, fig_prc_spp_byPfg, fig_prc_mds)


# ## set margins
# prc_margins <- llply(list(c(1, .5, 1, .5),  c(1, .5, 6, .5), c(0, .5, 1, .5)),
#                      function(x) unit(x, "line"))
# 
# for (i in 1:3){
#   prc_figs[[i]] <-  prc_figs[[i]]  + theme(plot.margin = prc_margins[[i]],
#                                            axis.text = element_text(size = 9))
# }


## define layout
lo <- matrix(c(1, 1, 3, 3, 3,
               1, 1, 3, 3, 3,
               2, 2, 2, 2, 2,
               2, 2, 2, 2, 2,
               2, 2, 2, 2, 2), 
             ncol = 5, byrow = TRUE)


## save as pdf
pdf(file = "output/figs/prc_merge_fig.pdf", width = 6.5, height = 6.5)
multiplot(plotlist = prc_figs, layout = lo)
dev.off()
