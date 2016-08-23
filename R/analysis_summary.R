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
         Fv        = paste0("F(", df, ", ", Df.res, ")=", statistic),
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
         Ind  = tstrsplit(.id, "[.]")[[3]]) %>% 
  select(term, variable, value, Form, Ind) %>% 
  mutate(Form_var = paste(Form, variable, sep = "_")) %>% 
  select(-Form, -variable) %>% 
  spread(key = Form_var, value) %>% 
  select(Ind, everything()) %>% 
  arrange(Ind, term)

# dominent species
dom_tbl <- anova_df_ed %>% 
  filter(Type == "dominentSpp") %>%
  mutate(Spp = paste(tstrsplit(.id, "[.]")[[2]], tstrsplit(.id, "[.]")[[3]])) %>% 
  select(term, variable, value, Spp) %>% 
  mutate(Spp_var = paste(Spp, variable, sep = "_")) %>% 
  select(-Spp, -variable) %>% 
  spread(key = Spp_var, value) %>% 
  arrange(term)

# PFG proportion
pfg_tbl <- anova_df_ed %>% 
  filter(Type %in% c("grass_prop", "pfg_prop")) %>% 
  mutate(PFG = ifelse(.id == "grass_prop", "grass_prop", tstrsplit(.id, "[.]")[[2]])) %>% 
  select(term, variable, value, PFG) %>% 
  mutate(PFG_var = paste(PFG, variable, sep = "_")) %>% 
  select(-PFG, -variable) %>% 
  spread(key = PFG_var, value) %>% 
  arrange(term)

# save tables
write.csv(div_tbl, file = "output/table/diversity_indices_anova.csv", row.names = FALSE)
write.csv(dom_tbl, file = "output/table/dominentSpp_anova.csv", row.names = FALSE)
write.csv(pfg_tbl, file = "output/table/PFG_anova.csv", row.names = FALSE)

save.image(file = "output/Data/summary_analysis.RData")
