
#  Prepare dataframe ------------------------------------------------------

prc_sp   <- decostand(PlotSumVeg[, SppName], method = "log") # log-transformed sp df
prc_site <- select(PlotSumVeg, year, ring, plot, co2)        # site df



# analysis ----------------------------------------------------------------

# ## environmental variable
# load("output/Data/EucFACE_understorey_env_vars_2012-2016.RData")            # load environmental variables, generated in "R/FitEnvironmentalVars.R"
# identical(EnvDF_3df[, c("year", "ring")], RingSumVeg[, c("year", "ring")])  # check EnvDf_3df and species df are in the same order for site

## Principal response curve anlsyis with bray-curtis dissimilarity
prc_all <- capscale(prc_sp ~  year * co2 + Condition(year), data = prc_site, 
                    distance = "bray")


## prepare df for a figure from the result of prc
summary_prc       <- summary(prc_all, scaling = 3)       # use scaling = 3, it doesn't matter too much, but this is used in the original prc function
res_prc_site      <- cbind(summary_prc$sites, prc_site)  # site score
res_prc_site_ring <- res_prc_site %>%                    # site scores by ring
  group_by(year, ring, co2) %>% 
  summarise_each(funs(mean), CAP1, MDS1, MDS2)           # CAP1 is for PRC plot and MDSs for visualisation
res_prc_site_co2  <- res_prc_site_ring %>%               # canonical coefficients for CO2 treatment (i.e. treatment difference for each year)
  group_by(year, co2) %>% 
  summarise_each(funs(mean), CAP1) %>% 
  group_by(year) %>% 
  summarise_each(funs(-diff(.)), CAP1) %>% 
  mutate(co2 = "elev") %>% 
  bind_rows(data.frame(year = levels(.$year), CAP1 = 0, co2 = "amb")) %>% 
  mutate(fig = "site")



# figure ------------------------------------------------------------------


# . canonical coefficient -------------------------------------------------


fig_prc_site <- ggplot(res_prc_site_co2, aes(x = as.numeric(year), y = CAP1, 
                                             fill = co2, linetype = co2)) +
  labs(x = "", y = expression(Canonical~coefficient~(C[dt]))) +
  geom_hline(yintercept = 0, linetype = "dotted", col = "gray") +
  
  geom_line() +
  geom_point(size = 2, shape = 21) +
  
  science_theme +
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
  scale_size_continuous("Proportion\nin Year0", range = c(1, 5), breaks = c(0, 0.05, 0.1, 0.3)) +
  science_theme +
  theme(legend.title     = element_text(size = 10),
        legend.position  = "right",
        axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.key.width = unit(.1, "inches")) +
  labs(y = "Species weight", x = "") +
  ylim(c(-.5, .5))
fig_prc_spp_byPfg





# . PCoA (MDS) plot -------------------------------------------------------


fig_prc_mds <- ggplot(res_prc_site_ring, aes(x = MDS1, y = MDS2, shape = year, 
                                             fill = co2)) +
  geom_hline(yintercept = 0, linetype = "dotted", col = "gray") +
  geom_vline(xintercept = 0, linetype = "dotted", col = "gray") +
  
  geom_point(size = 2, alpha = .8) +
  geom_path(aes(group = ring, linetype = co2)) +
  
  science_theme +
  theme(legend.title = element_text(size = 10)) +
  scale_fill_manual(name   = expression(CO[2]),
                    values = c("black", "white"), 
                    labels = c("Ambient", expression(eCO[2])),
                    guide  = guide_legend(override.aes = list(shape = 21))) +
  scale_shape_manual(name = "Year", 
                     values = c(21, 22, 23, 24)) +
  scale_linetype_manual(name = expression(CO[2]),
                        values = c("solid", "dashed"), 
                        labels = c("Ambient", expression(eCO[2]))) +
  theme(legend.position = "right")
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
