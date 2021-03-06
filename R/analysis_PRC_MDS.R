
#  Prepare dataframe ------------------------------------------------------

prc_sp       <- log(PlotSumVeg[, SppName_grass] + 1)                  # log-transformed grass sp df
prc_sp_all   <- log(PlotSumVeg[, SppName] + 1)                  # log-transformed all sp df
prc_site     <- PlotSumVeg[, c("year", "ring", "plot", "co2", "id")]  # site df 




# analysis ----------------------------------------------------------------


# > PRC ---------------------------------------------------------------------


prc_all <- capscale(prc_sp ~  year * co2 + Condition(year), data = prc_site,  # principal response curve anlsyis
                    distance = "euclidean")
mds_all <- cmdscale(d = vegdist(prc_sp_all, method = "bray"), eig = TRUE, k = 2)  # MDS with bray-curtis dissimilarity


# . define permutation ------------------------------------------------------


## define permutation for this experiment; year is nested within plot within
## ring.

## 1. plot is a random factor, should never be exchanged. 
## 2. year is a time series measurement. This can be exchanged, but its order 
## should be kept to take autocorrelation into account. Also this order should 
## be consistent among four plots within each ring, whereas it doesn't need to
## be consistent between rings. (e.g. Ring1 [0, 1, 2, 3],  Ring2 [2, 3, 0, 1], 
## Ring3 [3, 0, 1, 2]). 
## 3. ring is an experimental unit so fully exchangable


## 1) define permutation for year within each plot for each ring
ctrl_year <- how(within = Within(type = "series", constant = TRUE),    # year within plot is a time seiries and the order should be consistent aamong plots
                 plot   = Plots(strata = prc_site$id, type = "none"))  # plot should not be exchanged
perm_year <- shuffleSet(nrow(prc_site), control = ctrl_year)           # define permutation for year
perm_year <- rbind(1:nrow(prc_site), perm_year)                        # add the original (i.e. 1:96) as this should be also permuted between rings
nrow(perm_year)


## 2) get permutation for each ring
perm_year_byRing   <- alply(perm_year, 1, function(x) split(x, prc_site$ring))  # split the above defined permutaiton by ring and get permutation for each year. there are four permutations for each ring
perm_year_byRing_l <- llply(1:6, function(x){                                   # store permutations for each ring in a list
  laply(perm_year_byRing, function(y) y[[x]])
})
names(perm_year_byRing_l) <- paste0("ring", 1:6)


## 3) define the order of rings to combine
allperm6   <- rbind(1:6, allPerms(6))          # get all possible permutaiton for 6 rings (6! = 720); the above permutaitons for each ring will be exchanged between rings according to these combination; this defines the order of rings to combine (i.e., [1, 2, 3, 4, 5, 6], [1, 3, 2, 4, 5, 6])
perm6_ring <- alply(allperm6, 1, function(x){  # reorder the above ring-permutation for each of combination defined in allperm6 (720 combination) 
  perm_year_byRing_l[x]
}, .progress = "text")
llply(perm6_ring, summary)


## 4) combine rings; there are 4 permutaitons for each ring x 6 rings x 720 combination = 4^6 * 720 = 2949120
allcomb4    <- as.matrix(expand.grid(1:4, 1:4, 1:4, 1:4, 1:4, 1:4))  # define all combinations of 4 pemutaions to be combined among the 6 rings (4^6 = 4096); e.g. [1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 2], [1, 1, 1, 1, 1, 3]

get_allperm <- function(ring_perm_list, ncomb = 99){
  
  # ring_perm_list: a list of permutations for each ring  (e.g. perm6_ring[[10]], this contains 6 object (i.e. rings) which stores permutaitons of plots (4 permutations))
  # ncomb         : number of combinations to be used out of 46096 in allcmob4. One can use all combinations, although it take a while
  
  r1_perm <- alply(allcomb4[sample(nrow(allcomb4), ncomb), ], 1, function(x){  # randomly select ncomb rows from allcomb4 and apply the function below for each row
    do.call(c, llply(1:6, function(y){
      ring_perm_list[[y]][x[y], ]  
      ## 1) ring_perm_list[[y]] is yth object in ring_perm_list (e.g. ring1)
      ## 2) x is xth row in  allcomb4[sample(nrow(allcomb4), ncomb), ] (e.g. [1, 3, 2, 4, 2, 3])
      ## 3) x[y] is yth objct in x (e.g. 4)
      ## 4) e.g. x = [1, 3, 2, 1, 2, 3], y = 4; ring_perm_list[[4]][x[4], ] extracts the 1st row (x[4] = 1) of the 4th object in ring_perm_list
      ## 5) ring_perm_list contains 6 object, so repeate 4) for 6 times. do.call(c, ...) cmobines the resulted 6 objects
    })) 
    
  })
  r1_perm <- do.call(rbind, r1_perm)  # rbind the all results from ncomb rows of allcomb4
  return(r1_perm)
}

## get all permutation; this requires computatino power, so use multicores and carry out parallel processing to save time
detectCores()                # number of cores in the current machine
registerDoParallel(cores = 3)  # register parallel background

alperms <- llply(perm6_ring, get_allperm, ncomb = 45, .parallel = TRUE,  
                 .paropts = list(.export = "allcomb4"))  # allcomb4 is difined outside of this function, so export it
# llply(alperms, dim)
alperms_bind <- do.call(rbind, alperms) # ncomb x 720
any(apply(alperms_bind, 1, function(x) identical(as.vector(x), c(1:96))))  # check if this contains the original order (ie.e 1:96). if so remove
dim(alperms_bind)


## check permutation; ring is freely exchanged. Plot is not exchanged. Year is
## exchanged within plots, but the order (time series) is preserved

llply(split(alperms_bind[sample(nrow(alperms_bind), 1), ], prc_site$ring), function(x){
  m <- matrix(x, ncol = 4)
  apply(m, c(1, 2),  function(x) paste(prc_site$id, prc_site$year, sep = "-")[x])
})

## also one can check how permutation function generates permutation here
# source("R/check_permut.R")

# . run permutation test --------------------------------------------------

allperm_n <- nrow(alperms_bind)  # number of all permutations
prc_res <- anova(prc_all, permutations = alperms_bind[sample(allperm_n, 9999), ], # randomly select 9999 rows
                 by = "axis", parallel = 3)
prc_res




# . summarise -------------------------------------------------------------


## prepare df for a figure from the result of prc
summary_prc       <- summary(prc_all, scaling = 3)                              # use scaling = 3, it doesn't matter too much, but this is used in the original prc function
MDSprop           <- mds_all$eig[1:2] * mds_all$GOF[2]/(sum(mds_all$eig[1:2]))  # proportion explained by the first 2 MDS axies
MDSaxes           <- paste0("MDS", 1:2, " (", round(MDSprop * 100, 0), "%)")    # MDS axis with explained proportion
MDSsite           <- scores(mds_all)
colnames(MDSsite) <- paste0("MDS", 1:ncol(MDSsite))                             # modify column names
res_prc_site      <- cbind(CAP1 = summary_prc$sites[, 1], MDSsite, prc_site)    # site score



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
                  nperm = 9999)
mds_envfit <- envfit(res_prc_site_ring[, c("MDS1", "MDS2")],                      # fit environmental variables
                     EnvDF_3df[, c("co2", "Depth_HL", "TotalC", "Drysoil_ph", "moist", 
                                   "gapfraction", "temp", "sand", "silt", "clay")],
                     permutations = cntr)

mds_vect_d <- data.frame(scores(mds_envfit, "vectors"),   # scores for vecorts
                         pval = mds_envfit$vector$pvals,
                         r2  = mds_envfit$vector$r) %>% 
  mutate(env = recode(row.names(.),
                      Depth_HL    = "HL", 
                      TotalC      = "Total C",
                      Drysoil_ph  = "pH",
                      moist       = "Moist",
                      gapfraction = "Light",
                      temp        = "Temp",
                      sand        = "Sand",
                      silt        = "Silt",
                      clay        = "Clay"))

mds_fact_d <- data.frame(scores(mds_envfit, "factors"),   # scores for facctors (CO2)
                         pval = mds_envfit$factors$pvals,
                         r2   = mds_envfit$factors$r,
                         env  = "CO2")[1, ] 

mds_arrw_d <-  mds_vect_d %>% 
  mutate(co2 = factor("amb", levels = c("amb", "elev")),
         year = factor("Year0", levels = paste0("Year", 0:3))) %>% 
  filter(pval <= 0.1)

mds_envfit_d <- bind_rows(mds_vect_d, mds_fact_d) %>% 
  mutate(pval = round(pval, 3)) %>% 
  mutate_each(funs(format(., digit = 0, nsmall = 2)), r2, MDS1, MDS2) %>% 
  mutate(variable = recode_factor(env, 
                                  CO2 = "CO2", sand = "Sand", silt = "Silt", 
                                  clay = "Clay", Drysoil_ph = "pH", 
                                  'Total C' = "Total C", moist = "Moist", 
                                  temp = "Temp", Depth_HL = "HL", 
                                  gapfraction = "Light")) %>%
  select(variable, MDS1, MDS2, r2, pval) %>% 
  arrange(variable)
write.csv(mds_envfit_d, file = "output/table/result_MDS_envfit.csv", row.names = FALSE)


res_prc_site_co2  <- res_prc_site_ring %>%                                    # canonical coefficients for CO2 treatment (i.e. treatment difference for each year)
  group_by(year, co2) %>% 
  summarise_each(funs(mean), CAP1) %>% 
  group_by(year) %>% 
  summarise_each(funs(diff(.)), CAP1) %>% 
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
                                             linetype = co2)) +
  labs(x = "Year", y = expression(Canonical~coefficient~(C[dt]))) +
  geom_hline(yintercept = 0, linetype = "dotted", col = "gray") +
  
  geom_line() +
  geom_point(size = 2.5, aes(shape = co2)) +
  
  science_theme_prc +
  theme(legend.position = c(.15, .75),
        legend.title    = element_blank()) +
  scale_x_continuous(labels = 0:3) +
  scale_linetype_manual(values = c("solid", "dashed"), 
                        labels = c("Ambient", expression(eCO[2]))) +
scale_shape_manual(values = c(16, 15), 
                        labels = c("Ambient", expression(eCO[2]))) +
  annotate("text", x = -Inf, y = Inf, label = "(a)", hjust = -.5, vjust = 1.5, fontface = "bold")
fig_prc_site




# . Sp weight plot -----------------------------------------------------------


## merge sp score with pfg and CO2 response ratios in Year3


### sp score given by PRC
res_prc_sp <- data.frame(summary_prc$species,
                         variable = row.names(summary_prc$species))

### PFG, form table
pfg_spp_tbl <- veg_FullVdf %>%
  filter(form == "Grass") %>% 
  select(variable, PFG) %>% 
  distinct() %>% 
  mutate(PFG = ifelse(PFG == "c3", "C3", "C4")) %>% 
  select(variable, PFG)


### merge species score and pfg and origin table
res_pric_sp_d <- res_prc_sp %>% 
  left_join(pfg_spp_tbl, by = "variable") %>% 
  mutate(year    = factor("Year0", levels = paste0("Year", 0:4)),
         co2     = factor("Ambient", levels = c("Ambient", "eCO[2]")))

### save species scores
unique(res_pric_sp_d$pfgform)
res_pric_sp_tbl <- res_pric_sp_d %>% 
  select(variable, PFG, CAP1) %>%
  transmute(Species   = gsub("[.]", " ", as.character(variable)),
            PRC.score = round(CAP1, 2),
            PFG       = PFG) %>%
  arrange(PRC.score)
write.csv(res_pric_sp_tbl, "output/table/summary_PRCscore.csv", row.names = FALSE)




## create a plot
fig_prc_spp_byPfg <- ggplot(res_pric_sp_d, aes(x = PFG, y = CAP1)) +
  geom_hline(yintercept = 0, linetype = "dotted", size = .5) +
  geom_boxplot(outlier.size = 2, na.rm = TRUE) +
  science_theme_prc +
  labs(y = "Species weight", x = "PFG") +
  geom_text(aes(label = "(b)"), x = -Inf, y = Inf, 
            hjust = -.5, vjust = 1.5, fontface = "bold")
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
               aes(x = 0, y = 0, xend = MDS1 * .23, yend = MDS2 * .23),
               arrow = arrow(length = unit(.2, "cm")),  alpha = .7) +
  geom_text(data = mds_arrw_d,  aes(x = MDS1 * .25, y = MDS2* .25, label = env),  
            size = 2.5, fontface = "italic") +
  

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
lo <- matrix(c(1, 1, 2), ncol = 3)


## save as pdf and png
pdf(file = "output/figs/prc_merge_fig.pdf", width = 5.5, height = 3)
multiplot(plotlist = list(fig_prc_site, fig_prc_spp_byPfg), layout = lo)
dev.off()

save_png600(file = "output/figs/prc_merge_fig.png", width = 5.5, height = 3)
multiplot(plotlist = list(fig_prc_site, fig_prc_spp_byPfg), layout = lo)
dev.off()

setEPS()
postscript(file = "output/figs/prc_merge_fig.eps", width = 5.5, height = 3)
multiplot(plotlist = list(fig_prc_site, fig_prc_spp_byPfg), layout = lo)
dev.off()
