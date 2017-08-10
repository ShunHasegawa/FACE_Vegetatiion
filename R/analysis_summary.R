
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




save.image(file = "output/Data/summary_analysis.RData")
