
# combine all models ------------------------------------------------------


all_m_list <- list('diversity'   = div_m_list,
                   'sd_abund'    = sd_m_list,                                     # abundance of dominant/subordiante C3/C4 ratios
                   'pfg_prop'    = list('sd_radio' = sd_m1, 'c43_ratio' = c43_m1))
all_m_list <- unlist(all_m_list)



# combine and organise anova table ----------------------------------------


anova_df <- ldply(all_m_list, function(x) tidy(Anova(x, test.statistic = "F")))

anova_df_ed <- anova_df %>% 
  mutate(term  = factor(ifelse(term %in% c("co2", "year", "co2:year"), term, "Baseline"),
                        levels = c("co2", "year", "co2:year", "Baseline"),
                        labels = c("CO2", "Year", "CO2xYear", "Baseline")),
         statistic = round(statistic, 2),
         Df.res    = round(Df.res, 0), 
         Fv        = paste0("F(", df, ",", Df.res, ")=", statistic),
         P         = round(p.value, 3),
         P         = ifelse(P < 0.001, "<0.001", P),
         Type      = tstrsplit(.id, "[.]")[[1]]) %>% 
  select(Type, .id, term, Fv, P) %>% 
  gather(variable, value, Fv, P)
anova_df_ed




# create tables -----------------------------------------------------------



# diversity indices
div_tbl <- anova_df_ed %>% 
  filter(Type == "diversity") %>% 
  mutate(Form       = tstrsplit(.id, "[.]")[[2]],
         Ind        = tstrsplit(.id, "[.]")[[3]],
         response   = paste(Ind, Form, sep =  "_"),
         report     = ifelse(Form == "grass_spp", "main", "appendix")) %>% 
  select(Type, term, variable, value, response, report) 
  

# abundance of each functional groups
abund_tbl <- anova_df_ed %>% 
  filter(Type == "sd_abund") %>%
  mutate(response = tstrsplit(.id, "[.]| ")[[2]],
         report   = "main") %>% 
  select(Type, term, variable, value, response, report)


# PFG proportion
pfg_tbl <- anova_df_ed %>% 
  filter(Type == "pfg_prop") %>% 
  mutate(response = tstrsplit(.id, "[.]| ")[[2]],
         report   = "main") %>% 
  select(Type, term, variable, value, response, report) 


# mege the above tables
all_tbl <- rbind.fill(div_tbl, abund_tbl, pfg_tbl) %>% 
  mutate(tv       = paste(term, variable, sep = "_")) %>% 
  select(-term, -variable) %>% 
  spread(tv, value) %>% 
  select(response, starts_with("Baseline"), starts_with("CO2_"), starts_with("Year"),
         starts_with("CO2x"), report, -Type)
all_tbl


# split by report type (main or appendix document)
all_tbl_l <- dlply(all_tbl, .(report), function(x) select(x, -report))


# save tables
writeWorksheetToFile(file        = "output/table/summary_rptdAov_tbl.xlsx",    # define file name to be saved
                     data        = all_tbl_l,                                  # writeWorksheetToFile doesn't take dplyr object so turn them into data frames using as.data.frame
                     sheet       = names(all_tbl_l),                           # sheet names in excel are defined by object names a list
                     clearSheets = TRUE)



# . summary table with VC -------------------------------------------------

# variance components
vc_df <- ldply(all_m_list, function(x){
  vc <- data.frame(VarCorr(x, comp="Variance")) %>% 
    select(grp, vcov) %>% 
    mutate(vcov = round(vcov, 4)) %>% 
    spread(grp, vcov) %>% 
    dplyr::rename(vc_residual = Residual,
                  vc_ring = ring) %>% 
    select(vc_ring, vc_residual)
  return(vc)
})

anova_vc_df <- anova_df %>% 
  mutate(term  = factor(ifelse(term %in% c("co2", "year", "co2:year"), term, "Baseline"),
                        levels = c("co2", "year", "co2:year", "Baseline"),
                        labels = c("CO2", "Year", "CO2xYear", "Baseline")),
         statistic = round(statistic, 2),
         Df.res    = round(Df.res, 0), 
         dfs        = paste0("F(", df, ",", Df.res, ")"),
         P         = round(p.value, 3),
         P         = ifelse(P < 0.001, "<0.001", P), 
         Fval = paste(term, dfs, sep = "_")) %>%
  select(-df, -Df.res, -p.value, term, -dfs) %>% 
  gather(variable, value, statistic, P) %>% 
  mutate(Fval = ifelse(variable == "statistic", Fval, paste(term, "p", sep = "_"))) %>% 
  select(-term, -variable) %>% 
  spread(Fval, value) %>% 
  left_join(vc_df) %>% 
  mutate(Type = tstrsplit(.id, "[.]")[[1]],
         Type = factor(Type, levels = c("diversity", "pfg_prop", "sd_abund"))) %>%
  arrange(Type) %>% 
  select(.id, starts_with("Baseline"), starts_with("CO2_"), starts_with("Year_"), 
         starts_with("CO2xYear_"), starts_with("vc_"))
anova_vc_df
write.csv(anova_vc_df, "output/table/summary_tbl_ringmean_analysis.csv", row.names = FALSE)


save.image(file = "output/Data/summary_analysis.RData")
