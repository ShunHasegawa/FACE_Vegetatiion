
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




# . combine with environmental vars ---------------------------------------


# combine environment and spp df, then split df for each year for each form
SiteName_rda <- c("year", "ring", "block", "co2")
seDF_b_year <- llply(list('all' = SppName, 'grass' = SppName_grass, 'forb' = SppName_forb),
                     function(x) {
                       d <- RingSumVeg %>%                        # df for each form 
                         select(one_of(x, SiteName_rda)) %>% 
                         left_join(EnvDF_3df, by = SiteName_rda)
                       d <- split(d, d$year)                      # split by year
                       return(d)
                     })

seDF_b_year <- unlist(seDF_b_year, recursive = FALSE)


# remove un-observed spp and transform data
seDFs <- llply(seDF_b_year, function(x){
  spsum     <- summarise_each(x, funs(sum), -one_of(SiteName_rda, expl))     # sp sum
  sp_to_use <- names(spsum)[spsum != 0]                                      # spp to be used
  new_d     <- select(x, one_of(sp_to_use, SiteName_rda, expl))              # subset spp to be used
  new_d[, sp_to_use] <- decostand(new_d[, sp_to_use], method = "hellinger")  # transform using hellinger
  return(new_d)
})
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


write.csv(file = "output/table/RDA_result_table.csv", all_model_results_tbl, 
          row.names = FALSE)




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
         Year2_full_mod, Year2_final_mod,
         Year3_full_mod, Year3_final_mod)

write.csv(file = "output/table/RDA_model_summary.csv", rda_model_tbl, row.names = FALSE)
  


# 4-year data set ---------------------------------------------------------

# From the above analysis, the following variables will be fitted for each form
# All: moist, Drysoil_ph
# Forb: TotalC, moist
# Grass: Drysoil_pH

# renemae env variable for plotting purposes
EnvDF_3df2 <- rename(EnvDF_3df, Moist = moist, pH = Drysoil_ph, Total_C = TotalC)

# list of df for each form
seDF_4y_list <- llply(list('all' = SppName, 'grass' = SppName_grass, 'forb' = SppName_forb),
                      function(x) {
                        d <- RingSumVeg %>% 
                          select(one_of(x, SiteName_rda)) %>%              # select spp for each form
                          left_join(EnvDF_3df2, by = SiteName_rda)         # merge with environmental variables
                        d[, x] <- decostand(d[, x], method = "hellinger")  # transform
                        return(d)
                      })

all_4y_d   <- seDF_4y_list[["all"]]
forb_4y_d  <- seDF_4y_list[["forb"]]
grass_4y_d <- seDF_4y_list[["grass"]]


# rda for each form
r_all   <- rda(all_4y_d[, SppName] ~ as.numeric(year) + Moist + pH, data = all_4y_d)
r_forb  <- rda(forb_4y_d[, SppName_forb] ~ as.numeric(year) + Total_C + Moist, data = forb_4y_d)
r_grass <- rda(all_4y_d[, SppName_grass] ~ as.numeric(year) + Moist + pH, data = all_4y_d)

rda_4y <- list('All' = r_all, 'Forb' = r_forb, 'Grass' = r_grass)
par(mfrow = c(2, 2))
l_ply(rda_4y, plot)

summary(r_all)$biplot  
  
  

# figure ------------------------------------------------------------------


# list of rda scores to create plots
rda_plot_par <- llply(rda_4y, get_rda_scores)


# add plant forms and constant for rescaling predictors on triplot
bc_cons <- data.frame(form   = c("All", "Forb", "Grass"),
                      b_cons = c(.7, .6, .5))              # constant for biplot

rda_plot_par <- mlply(bc_cons, function(form, b_cons){
  l <- rda_plot_par[[form]]
  l$sitedd <- mutate(l$sitedd, Form = form)
  l$bipldd <- mutate(l$bipldd, Form = form)
  l$b_cons <- b_cons
  return(l)
})
names(rda_plot_par) <- bc_cons$form
llply(rda_plot_par, summary)


# create plots
rda_plots <- llply(rda_plot_par, function(x) do.call("create_rda_plots", x))


rda_plots[[3]] <- rda_plots[[3]] +
  theme(legend.title      = element_text(size = 7),
        legend.position   = c(1.5, .6),
        legend.text.align = 0,
        legend.box.just   = "left",
        legend.margin     = unit(-.1, "line"),
        legend.text       = element_text(size = 7),
        legend.key.height = unit(.13, "in"),
        legend.key.width  = unit(.3, "in"))


rda_plots2 <- llply(rda_plots, function(x) x + 
                      theme(axis.title = element_text(size = 8),
                            axis.text  = element_text(size = 6)))

# blank plot
bp <- ggplot(rda_plot_par$All$sitedd, aes(x = RDA1, y = RDA2)) +
  geom_blank(inherit.aes = FALSE) +
  labs(x = " ", y = " ") +
  theme(panel.border     = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_blank()) +
  facet_wrap(~Form)


# merge plots
rda_plot_merged <- rbind(cbind(ggplotGrob(rda_plots2[[1]]), ggplotGrob(rda_plots2[[2]])), 
                         cbind(ggplotGrob(rda_plots2[[3]]), ggplotGrob(bp)), size = "last")

grid.newpage()
grid.draw(rda_plot_merged)

ggsavePP(filename = "output/figs/RDA_3forms", plot = rda_plot_merged,
         width = 6.5, height = 6)
