
#  Prepare dataframe ------------------------------------------------------

prc_sp   <- decostand(PlotSumVeg[, SppName], method = "log") # log-transformed sp df
prc_site <- select(PlotSumVeg, year, ring, plot, co2)                      # site df



# analysis ----------------------------------------------------------------


## Principal response curve anlsyis with bray-curtis dissimilarity
prc_all <- capscale(prc_sp ~  year * co2 + Condition(year), data = prc_site, 
                    distance = "bray")
summary(prc_all)
a$species
plot(prc_all, scaling = 3)

# ## environmental variable
# load("output/Data/EucFACE_understorey_env_vars_2012-2016.RData")            # load environmental variables, generated in "R/FitEnvironmentalVars.R"
# identical(EnvDF_3df[, c("year", "ring")], RingSumVeg[, c("year", "ring")])  # check EnvDf_3df and species df are in the same order for site




#####

summary_prc <- summary(prc_all, scaling = 3)
res_prc_site <- cbind(summary_prc$sites, prc_site)
res_prc_site_ring <- res_prc_site %>% 
  group_by(year, ring, co2) %>% 
  summarise_each(funs(mean), CAP1, MDS1)
res_prc_site_co2 <- res_prc_site_ring %>% 
  group_by(year, co2) %>% 
  summarise_each(funs(mean), CAP1, MDS1)

fig_prc_site <- ggplot(res_prc_site_ring, aes(x = year, y = CAP1, col = co2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(position = position_dodge(.2)) +
  geom_line(data = res_prc_site_co2, aes(x = as.numeric(year), y = CAP1, linetype = co2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  science_theme
# +
  # theme(legend.position = "top") +
  # ylim(c(-1.4, 1.4))
fig_prc_site

fig_prc_mds <- ggplot(res_prc_site_ring, aes(x = MDS1, y = CAP1, col = co2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point() +
  geom_path(aes(group = ring)) +
  science_theme +
  # theme(legend.position = "top") +
  # ylim(c(-1.4, 1.4)) +
  labs(y = "")
fig_prc_mds

# merge sp score with pfg and propotion in Year0
res_prc_sp <- data.frame(summary_prc$species,                          # sp score given by PRC
                         variable = row.names(summary_prc$species))
pfg_spp_tbl <- veg_FullVdf %>%                                         # PFG, form table
  select(variable, form, PFG) %>% 
  distinct() %>% 
  mutate(form    = as.character(form),
         PFG     = as.character(PFG),
         pfgform = ifelse(form == "Grass", paste0(PFG, form),
                          ifelse(form == "Forb", PFG, form))) %>% 
  select(variable, pfgform)
  
year0_prop <- PlotSumVeg %>%                                            # proportion of each sp in Year0
  filter(year == "Year0") %>% 
  gather(variable, value, one_of(SppName)) %>% 
  group_by(variable) %>% 
  summarise(value = sum(value)) %>% 
  mutate(prop = value / sum(value))
res_pric_sp_d <- Reduce(function(...) merge(..., by = "variable"),
                        list(res_prc_sp, pfg_spp_tbl,year0_prop))
summary(res_pric_sp_d)

# re-order pfgform by median
pfgform_lev <- res_pric_sp_d %>% 
  group_by(pfgform) %>% 
  summarise(med = median(CAP1)) %>% 
  arrange(-med) %>% 
  .$pfgform
res_pric_sp_d$pfgform <- factor(res_pric_sp_d$pfgform, levels = pfgform_lev)


fig_prc_spp_byPfg <- ggplot(res_pric_sp_d, aes(x = pfgform, y = CAP1)) +
  geom_hline(yintercept = 0, linetype = "dotted", size = .5) +
  geom_boxplot(size = .3, outlier.size = 1, outlier.color = "grey") +
  geom_jitter(aes(size = prop), shape = 4, alpha = .6, width = .4) +
  scale_size_continuous(range = c(0, 1)) +
  science_theme +
  theme(legend.position = c(.1, .2),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.key.width = unit(.1, "inches")) +
  labs(y = "") +
  ylim(c(-.5, .5))
fig_prc_spp_byPfg

prc_figs <- list(fig_prc_site, fig_prc_mds, fig_prc_spp_byPfg)

# set margins

prc_margins <- llply(list(c(1, -1, 0, .5), c(0, -1, 0, 0), c(0, .5, 0, 0)),
                     function(x) unit(x, "line"))

for (i in 1:3){
  prc_figs[[i]] <-  prc_figs[[i]] + theme(plot.margin = prc_margins[[i]],
                                          axis.text   = element_text(size = 6),
                                          legend.text = element_text(size = 6))
}


fig_prc_merged <- cbind(ggplotGrob(prc_figs[[1]]),
                        ggplotGrob(prc_figs[[2]]), 
                        ggplotGrob(prc_figs[[3]]))
grid.newpage()
grid.draw(fig_prc_merged)
ggsavePP("output/figs/prc_spp_byPfg", fig_prc_merged, width = 6.5, height = 3)
