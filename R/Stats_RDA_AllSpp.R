
# prepare df --------------------------------------------------------------


# . transformation --------------------------------------------------------
# try transfomation as species abundance data are highly skewed


# raw
pdf(file = "output/figs/histgram_all_spp.pdf", onefile = TRUE, width = 6, height = 6)
par(mfrow = c(3, 3))
l_ply(SppName, function(x) hist(RingSumVeg[, x], main = x))
dev.off()


# log
pdf(file = "output/figs/histgram_all_spp_log.pdf", onefile = TRUE, width = 6, height = 6)
par(mfrow = c(3, 3))
l_ply(SppName, function(x) hist(log(RingSumVeg[, x] + 1), main = x))
dev.off()


# log but return 0 for 0
log2_d <- decostand(RingSumVeg[, SppName], method = "log", logbase = 10)
pdf(file = "output/figs/histgram_all_spp_logAnderson.pdf", onefile = TRUE, width = 6, height = 6)
par(mfrow = c(3, 3))
l_ply(SppName, function(x) hist(log2_d[, x], main = x))
dev.off()


# hellinger
d_hel <- decostand(RingSumVeg[, SppName], method = "hellinger")
pdf(file = "output/figs/histgram_all_spp_hellinger.pdf", onefile = TRUE, width = 6, height = 6)
par(mfrow = c(3, 3))
l_ply(SppName, function(x) hist(d_hel[, x], main = x))
dev.off()


# after visual inspection of histgram hellinger seems bettter, especially for
# abundant spp
tRingSumVeg <- RingSumVeg
tRingSumVeg[, SppName] <- decostand(RingSumVeg[, SppName], method = "hellinger")


# . combine with environmental vars ---------------------------------------


# combine environment and spp df, then split df for each year

SiteName_rda <- c("year", "ring", "block", "co2")

seDFs <- llply(list('all' = SppName, 'grass' = SppName_grass, 'forb' = SppName_forb),
               function(x) {
                 
                 # df for each form
                 d <- RingSumVeg[, c(x, SiteName_rda)]
                 
                 # split by year
                 dy <- split(d, d$year)
                 
                 # remove un-observed spp and transform data
                 new_d_list <- llply(dy, function(y){
                   
                   # sp sum
                   spsum <- colSums(y[, x])                             
                   
                   # spp to be used
                   sp_to_use <- names(spsum)[spsum != 0]   
                   
                   # subset spp to be used
                   new_d <- y[, c(sp_to_use, SiteName_rda)]
                   
                   # transform using hellinger 
                   new_d[, sp_to_use] <- decostand(new_d[, sp_to_use], method = "hellinger")
                   
                   # merge with environmental variables
                   new_d_merge <- left_join(new_d, EnvDF_3df, by = SiteName_rda)
                   
                   
                   return(new_d_merge)
                 })
                 
                 return(new_d_list)
                 
               })

llply(seDFs, summary)
seDFs <- unlist(seDFs, recursive = FALSE)
seDFs <- new_d_list
summary(seDFs)
llply(seDFs, names)



  
# single term -----------------------------------------------------------


# R2adj for each single term

single_adj_r2 <- ldply(seDFs, function(x){
  get_adjR_singl(x, ignoredd = TRUE, SiteName_rda = SiteName_rda, expl = expl)
})

# terms with positive adjR2
pos_adjr <- filter(single_adj_r2, adjR > 0)
  



# full models -----------------------------------------------------------


# create combinations of 1-4 terms and make formulas to be tested
full_formulas <- llply(1:4, function(x){
  dlply(pos_adjr, .(.id), function(y) {
    get_full_formula(y$variable, n = x)
  })
})  
names(full_formulas) <- paste0("term_n", 1:4)  # number of terms used in models
summary(full_formulas)


# combind formulas, df, environmental vars (expl) and site (SiteName_rda)

termn <- names(full_formulas)               # number of terms
fy    <- names(full_formulas[['term_n1']])  # form.year

f_df_list <- list()
for(i in termn){  # for each number of terms (i.e. 1-4)
  for(j in fy){   # for each form by year
    l <- paste(i, j, sep = ".")
    f_df_list[[l]] <- list(formula_list = full_formulas[[i]][[j]],
                           df           = seDFs[[j]],
                           expl         = expl,
                           SiteName_rda = SiteName_rda)
  }
}


# rda summary
rda_summary <- ldply(f_df_list, function(x) do.call("get_rda_summary", args = x))


# acceptable models
rda_accept <- rda_summary %>% 
  filter(vif_less10, p_value < .1) %>% 
  group_by(.id) %>%
  filter(adjr == max(adjr)) %>%                                # choose the one with the largest adjr
  ungroup() %>% 
  mutate(term_n  = as.numeric(tstrsplit(.id, "_n|[.]")[[2]]),  # number of terms
         dataset = gsub(".*_n.[.]", "", .id)) %>%              # form by year
  group_by(dataset) %>%
  filter(term_n == max(term_n)) %>%                            # choose the one with largest number of terms
  arrange(dataset)
    



# model simplification --------------------------------------------------


# get model formulas and associated dataframe to be simplified
models_to_simplify <- mlply(rda_accept[, c(".id", "f_id")],  # form.year and formula  
                            function(.id, f_id) {
                              f  <- f_df_list[[.id]]$formula_list[f_id]
                              df <- f_df_list[[.id]]$df
                              return(list(f = f, df = df, expl = expl, 
                                          SiteName_rda = SiteName_rda))
                            })

names(models_to_simplify) <- rda_accept$.id


simple_rdas <- llply(models_to_simplify, function(x) do.call("get_simple_rda", x))   
llply(simple_rdas, summary)




# Summary ---------------------------------------------------------------



# . Summary of each term in RDA -------------------------------------------


# get P and F for each term in the final model
rda_margin_term <- ldply(simple_rdas, function(x) tidy(x$anova_final_mod)) %>% 
  mutate(.id = gsub("term_n.[.]", "", .id)) %>% 
  group_by(.id) %>%  # degree of freedom for denominator for each id 
  mutate(DFden     = df[term == "Residual"],
         df        = df[term != "Residual"][1],
         statistic = round(statistic, 3),
         p.value   = get_star(p.value, dagger = FALSE)) %>%   
  ungroup() %>%
  filter(!term %in% c("Residual", "Model")) %>% 
  mutate(Fdf       = paste0("F(", df, ", ", DFden, ")=", statistic, p.value)) %>% 
  select(-df, -DFden, -Variance, -statistic, -p.value) %>% 
  rename(variable = term)


# get full and final models from the list
rda_mods <- llply(simple_rdas, function(x) {
  x$anova_final_mod <- NULL
  x
})


# get terms removed by model simplificaiton
rm_terms_dd <- ldply(rda_mods, function(x){
  t1 <- attr(terms(x$full_mod), "term.labels")   # terms in the full model
  t2 <- attr(terms(x$final_mod), "term.labels")  # terms in the final model
  
  rmval <- setdiff(t1, t2)
  
  if(length(rmval) > 0){ # where no term was removed, return NA
    removed <- TRUE  
  } else {
    rmval   <- NA
    removed <- NA
  }  
  
  data.frame(variable = rmval, removed = removed)
}) %>% 
  mutate(.id = gsub("term_n.[.]", "", .id)) %>% 
  filter(!is.na(variable))


# merge single term R2, terms in tha final model, and removed terms by simplificaiton
all_model_results <- Reduce(function(...) merge(..., by = c(".id", "variable"), all = TRUE),
                            list(single_adj_r2, rda_margin_term, rm_terms_dd))
all_model_results_tbl <- all_model_results %>% 
  mutate(Form      = tstrsplit(.id, split = "[.]")[[1]],
         Year      = tstrsplit(.id, split = "[.]")[[2]],
         adjR      = ifelse(adjR < 0, "<0", round(adjR, 3)),
         removed   = as.logical(removed),
         Fdf       = ifelse((is.na(Fdf) & is.na(removed)), "-", 
                            ifelse(is.na(removed), Fdf, "rm")),
         Term      = factor(variable, 
                            levels = c("co2", "TotalC", "moist", "Drysoil_ph", 
                                       "Depth_HL", "gapfraction", "temp"),
                            labels = c("CO2", "Total C", "Moist", "pH", "HL depth", 
                                       "Canopy transmittance", "Temp"))) %>% 
  select(-.id, -removed, -variable) %>%
  gather(variable, value, adjR, Fdf) %>% 
  mutate(Year_var = paste(Year, variable, sep = "_")) %>% 
  select(-Year, -variable) %>% 
  spread(key = Year_var, value) %>% 
  select(Form, Term, everything()) %>% 
  arrange(Form, Term)


write.csv(file = "output/table/RDA_result_table.csv", 
          all_model_results_tbl, row.names = FALSE)




# . Results of anova for full/final models -------------------------


# get model summary
rda_model_summary <- ldply(unlist(rda_mods, recursive = FALSE), get_rda_model_summary) %>% 
  mutate(.id = gsub("term_n.[.]", "", .id))


# format the summary
rda_model_tbl <- rda_model_summary %>%
  mutate(Form     = tstrsplit(.id, split = "[.]")[[1]],
         Year     = tstrsplit(.id, split = "[.]")[[2]],
         Mod      = tstrsplit(.id, split = "[.]")[[3]],
         Year_Mod = paste(Year, Mod, sep = "_"),
         mod_p    = get_star(mod_p, dagger = FALSE),
         mod_adjr = paste0(round(mod_adjr, 3), mod_p)) %>% 
  select(-.id, -Year, -Mod, -mod_p) %>% 
  spread(key = Year_Mod, value = mod_adjr) %>% 
  select(Form,  # re-ordering columns 
         Year0_full_mod, Year0_final_mod, 
         Year1_full_mod, Year1_final_mod,
         Year3_full_mod, Year3_final_mod)

write.csv(file = "output/table/RDA_model_summary.csv", rda_model_tbl, row.names = FALSE)
  


# 4-year data set ---------------------------------------------------------

# From the above analysis, moist and temp are determied to be imporatnt
# driver

# list of dfs for 4-year data
seDF_4y_list <- llply(list('all' = SppName, 'grass' = SppName_grass, 'forb' = SppName_forb),
                      function(x) {
                        tRingSumVeg %>% 
                          select(one_of(x, SiteName_rda)) %>% 
                          left_join(EnvDF_3df, by = SiteName_rda)
                      })

rda_4y_list <- llply(seDF_4y_list, function(x) rda(log(seDF[, SppName] + 1) ~  moist + temp + year, x))


  
rda_all <- 



rda_all <- rda(log(seDF[, SppName] + 1) ~ TotalC + moist + as.numeric(year), seDF)
  
  # can't run anova as it is. cause different perumutation units need to be
  # defined for year and co2. so anyway create a triplot and see the pattern.
  rda_all
  
  # plot
  p <- TriPlot(MultValRes = rda_all, env = seDF, yaxis = "RDA axis", axispos = c(1, 2, 3), centcons = 2)
  ggsavePP(filename = "output//figs/FACE_RDA_EnvVar_Year0_3", plot = p, width = 6, height = 6)


# Fig for thesis ----------------------------------------------------------

  RdaAllRes <- summary(rda_all)
  seDF$year <- factor(seDF$year, labels = paste0("Year", 0:3))
  sitedd <- data.frame(RdaAllRes$site, seDF)
  
  bipldd <- data.frame(RdaAllRes$biplot, co2 = "amb", year = "Year0", 
                       variable = row.names(RdaAllRes$biplot))
  bipldd <- subsetD(bipldd, variable %in% c("TotalC", "moist", "as.numeric(year)"))
  bipldd$variable <- factor(bipldd$variable,
                            levels = c("TotalC", "moist", "as.numeric(year)"),
                            labels = c("Moist", "Total C", "Year"))
  
  VarProp <- RdaAllRes$cont$importance["Eigenvalue",] / RdaAllRes$tot.chi
  axislabs <- paste0(c("RDA1", "RDA2"), "(", round(VarProp[c(1, 2)] * 100, 2), "%)")
  
  # make a plot
  p <- ggplot(data = sitedd, aes(x = RDA1, y = RDA2, shape = year))
  p2 <- p + 
    geom_path(aes(group = ring), col = "black") +
    geom_point(aes(fill = co2), size = 4) + 
    scale_fill_manual(values = c("black", "white"), 
                      labels = c("Ambient", expression(eCO[2])),
                      guide = guide_legend(override.aes = list(shape = 21))) +
    # need to add overrisde here to make white circle in the legend with shape of
    # 21
    scale_shape_manual(values = c(21, 22, 23, 24)) + 
    geom_segment(data = bipldd,
                 aes(x = 0, y = 0, xend = RDA1 * 2, yend = RDA2 * 2), 
                 arrow = arrow(length = unit(.2, "cm")), 
                 color = "red") +
    geom_text(data = bipldd, 
              aes(x = RDA1 * 2.3 , y = RDA2 * 2.3, label = variable), 
              lineheight = .7, 
              color = "red", size = 4, 
              fontface = "bold") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    science_theme +
    theme(legend.position = c(.17, .15), 
          legend.box = "horizontal", 
          legend.box.just = "top") +
    labs(x = axislabs[1], y = axislabs[2])
  RDA_Plot_AllSpp <- p2
  RDA_Plot_AllSpp
  ggsavePP(plot     = RDA_Plot_AllSpp, 
           filename = "output/figs/Fig_Thesis/RDA_3yr_AllSpp", 
           width    = 6, 
           height   = 4)
