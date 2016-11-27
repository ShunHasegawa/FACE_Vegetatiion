# combine all models
all_m_list <- list('diversity'   = div_m_list,
                   'dominentSpp' = dom_m_list,
                   'grass_prop'  = grassprop_m,
                   'pfg_prop'    = pfgprop_m_list)
all_m_list <- unlist(all_m_list)

# combine and organise anova table
anova_df <- ldply(all_m_list, function(x) tidy(Anova(x, test.statistic = "F")))

anova_df_ed <- anova_df %>% 
  mutate(term  = factor(ifelse(term %in% c("co2", "year", "co2:year"), term, "BL"),
                        levels = c("co2", "year", "co2:year", "BL"),
                        labels = c("CO2", "Year", "CO2xYear", "BL")),
         statistic = round(statistic, 2),
         Df.res    = round(Df.res, 0), 
         Fv        = paste0("F(", df, ",", Df.res, ")=", statistic),
         P         = round(p.value, 3),
         P         = ifelse(P < 0.001, "<0.001", P),
         Type      = tstrsplit(.id, "[.]")[[1]]) %>% 
  select(Type, .id, term, Fv, P) %>% 
  gather(variable, value, Fv, P)
anova_df_ed

# create tables

# diversity indices
div_tbl <- anova_df_ed %>% 
  filter(Type == "diversity") %>% 
  mutate(Form = tstrsplit(.id, "[.]")[[2]],
         Ind  = tstrsplit(.id, "[.]")[[3]],
         response = paste(Ind, Form, sep =  "_")) %>% 
  select(Type, term, variable, value, response)

# dominent species
dom_tbl <- anova_df_ed %>% 
  filter(Type == "dominentSpp") %>%
  mutate(response = paste(tstrsplit(.id, "[.]")[[2]], tstrsplit(.id, "[.]")[[3]])) %>% 
  select(Type, term, variable, value, response)

# PFG proportion
pfg_tbl <- anova_df_ed %>% 
  filter(Type %in% c("grass_prop", "pfg_prop")) %>% 
  mutate(response = ifelse(.id == "grass_prop", "grass_prop", tstrsplit(.id, "[.]")[[2]])) %>% 
  select(Type, term, variable, value, response) 

# mege the above tables
all_tbl <- rbind.fill(div_tbl, dom_tbl, pfg_tbl) %>% 
  mutate(tv = paste(term, variable, sep = "_")) %>% 
  select(-term, -variable) %>% 
  spread(tv, value) %>% 
  select(-Type)
all_tbl


# save tables
write.csv(all_tbl, "output/table/summary_rptdAov_tbl.csv", row.names = FALSE)


save.image(file = "output/Data/summary_analysis.RData")
