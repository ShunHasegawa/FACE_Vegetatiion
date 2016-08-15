# load data ---------------------------------------------------------------

# Raw data for multi variate analysis (matrix)
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S1.RData")
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S2.RData")
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S3.RData")
veg_matrix_list <- list(s1 = FullVdf, s2 = uniqueYear0_Vdf,
                        s3 = uniqueYear0_plot_Vdf)
llply(veg_matrix_list, summary)

# dfs with plant functional groups (df)
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S1_PFG.RData")
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S2_PFG.RData")
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S3_PFG.RData")
veg_df_list <- list(s1 = veg_FullVdf, s2 = veg_uniqueYear0_Vdf,
                    s3 = veg_uniqueYear0_plot_Vdf)
llply(veg_df_list, summary)

# spp
SppName_list <- llply(veg_df_list, function(x) unique(as.character((x$variable))))

# number of spp in Year0 for each solution
veg_df <- ldply(veg_df_list)
S_Year0 <-  veg_df %>% 
  group_by(.id, variable, year, form) %>% 
  summarise_each(funs(sum), value) %>% 
  filter(year == "Year0" & value != 0) %>% # rmeove spp not observed in Year0
  group_by(.id, form) %>% 
  summarise(S = n()) %>% # number of rows (i.e. spp)
  spread(.id, S) # reshape for table

# organise for further analysis -------------------------------------------

# > grass and forb spp ----------------------------------------------------

gfspp_list <- llply(veg_df_list, function(x){
  x %>% 
    filter(form %in% c("Grass", "Forb")) %>% 
    select(variable, form) %>%
    mutate(variable = as.character(variable)) %>% 
    distinct()
  })
llply(gfspp_list, some)

SppName_grass_list <- llply(gfspp_list, function(x) x$variable[x$form == "Grass"])
SppName_forb_list <- llply(gfspp_list, function(x) x$variable[x$form == "Forb"])

# > plot sum --------------------------------------------------------------

PlotSumVeg_list <- llply(1:3, function(x) {
  s <- SppName_list[[x]]
  d <- veg_matrix_list[[x]]
  d %>% 
    group_by(year, block, ring, plot, id, co2) %>% 
    summarise_each_(funs(sum), s) %>% 
    ungroup()
})
  
# > diversity indices -------------------------------------------------------

# create list of dfs containing only species (not site values)
vegDF_list <- llply(1:3, function(x){
  d <- PlotSumVeg_list[[x]]
  d_list <- list(all_spp   = d[, SppName_list[[x]]],       # all spp
                 grass_spp = d[, SppName_grass_list[[x]]], # grass spp
                 forb_spp  = d[, SppName_forb_list[[x]]])  # forb spp
  return(d_list)
})
names(vegDF_list) <- paste0("s", 1:3)
llply(vegDF_list, summary)

# compute diversity indices
siteDF <- select(PlotSumVeg_list[[1]], -one_of(SppName_list[[1]]))

DivDF_list <- llply(vegDF_list, function(x) {
  llply(x, function(y) mutate(siteDF, 
                              H = diversity(y),  # Shannon's index
                              S = specnumber(y), # number of spp
                              J = H/log(S)       # Pielou's evenness
                              ))})
llply(DivDF_list, summary)


# . summary plots ---------------------------------------------------------

# ring summary
summary_ring_list <- llply(DivDF_list, function(x){
  d <- ldply(x)
  d2 <- d %>% 
    gather(variable, value, H, S, J) %>% 
    group_by(.id, year, co2, ring, variable) %>% 
    summarise_each(funs(mean(., na.rm = TRUE)), value)
  return(d2)
})

# co2 summary
summary_co2_list <- llply(summary_ring_list, function(x){
  x %>% 
    group_by(.id, year, co2, variable) %>% 
    summarise_each(funs(M = mean(., na.rm = TRUE), 
                        SE = ci(., na.rm = TRUE)[4],
                        N = sum(!is.na(.))), 
                   value)
  })
llply(summary_co2_list, summary)

# plot
summary_plot_list <- llply(summary_co2_list, function(x){
  d <- x
  p <- ggplot(d, aes(x = as.numeric(year), y = M, group = co2, col = co2))
  p2 <- p +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 0.2) +
    facet_grid(variable ~  .id, scale = "free") +
    labs(y = "Value", x = "Year")
  return(p2)
  })

# . df for ANCOVA ---------------------------------------------------------

# Move Year0 value to a new column to be used as covariate for the analysis
DivDF_year0_list <- llply(DivDF_list, function(x){
  llply(x, function(y){
    
    DivDF_mlt <- gather(y, variable, value, H, S, J)
    
    # year0
    year0_dd <- DivDF_mlt %>%
      filter(year == "Year0") %>%
      select(id, value, variable) %>%
      rename(value0 = value)
    
    # subseqent years
    subyear_dd <- filter(DivDF_mlt, year != "Year0")
    
    # merge
    DivDF_year0 <- merge(subyear_dd, year0_dd, by = c("id", "variable")) 
    
    return(DivDF_year0)
  })
})

# plot against Year0
  
plot_vsYear0 <- llply(DivDF_year0_list, function(x) {
  d <- ldply(x)
  p <- ggplot(d, aes(x = value0, y = value, col = year))
  p2 <- p + 
    geom_point(alpha = .9) +
    facet_wrap(variable ~ .id , scales = "free")
  return(p2)
})

plot_vsYear0[[1]]
plot_vsYear0[[2]]
plot_vsYear0[[3]]

# analysis ----------------------------------------------------------------

# > two-way ANOVA ---------------------------------------------------------

# formulas for 2-way anova and ancova 
f_2anv <- llply(paste(c("H", "log(S)", "J"), 
                      "~ co2 * year + (1|block) + (1|ring) + (1|id)"), 
                formula)
names(f_2anv) <- c("H", "S", "J")

# list of 2-way anova models for each index
two_anv_m_list <- llply(DivDF_list, function(x){ # s1, s2, s3
  llply(x, function(y){ # all_spp, grass_spp, forb_spp
    llply(f_2anv, function(z) lmer(z, data = y)) # H, S, J
  })
})  
llply(two_anv_m_list, function(x) llply(x, summary))

# they are list of list of list, so unclass
two_anv_m_list <- unlist(two_anv_m_list, recursive = TRUE)

# F test
two_anv_ftest <- ldply(two_anv_m_list, function(x) tidy(Anova(x, test.statistic = "F")))


# model diagnosis 

pdf(file = "output/figs/diagnostic_plot_two.way.anova.pdf", onefile = TRUE, 
    width = 3, height = 3)
llply(names(two_anv_m_list), function(x) {
  m <- two_anv_m_list[[x]]
  print(plot(m, main = x))
  qqnorm(resid(m), main = x)
  qqline(resid(m))
  })
dev.off()

# > ancova with Year0 value -------------------------------------------------

f_ancv <- formula("value ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id)")

# list of ancova models for each index
ancov_m_list <- llply(DivDF_year0_list, function(x){ # s1, s2, s3
  llply(x, function(y){ # all_spp, grass_spp, forb_spp
    dlply(y, .(variable), function(z) lmer(f_ancv, data = z)) # H, S, J
  })
})  
# they are list of list of list, so unclass
ancov_m_list <- unlist(ancov_m_list, recursive = TRUE)

# F test
ancov_ftest <- ldply(ancov_m_list, function(x) tidy(Anova(x, test.statistic = "F")))

# model diagnosis 

pdf(file = "output/figs/diagnostic_plot_ancova.pdf", onefile = TRUE, 
    width = 3, height = 3)
llply(names(ancov_m_list), function(x) {
  m <- ancov_m_list[[x]]
  print(plot(m, main = x))
  qqnorm(resid(m), main = x)
  qqline(resid(m))
})
dev.off()

# > merge results --------------------------------------------------------------

# merge two test results and create a table
anova_df <- bind_rows(Tanova = two_anv_ftest, ancova = ancov_ftest, .id = "test") %>% 
  rename(F = statistic, DFnum = df, DFden = Df.res, P = p.value) %>% 
  filter(!grepl("^s2.grass|^s3.grass", .id)) %>% # grass is identical for all
  mutate(term  = factor(term, 
                        levels = c("co2", "year", "co2:year", "value0"),
                        labels = c("CO2", "Year", "CO2xYear", "BL")),
         F     = round(F, 2),
         DFden = round(DFden, 0), 
         P     = round(P, 3),
         S     = ldply(strsplit(.id, "[.]"))[, 1],
         Form  = ldply(strsplit(.id, "[.]"))[, 2],
         Ind   = ldply(strsplit(.id, "[.]"))[, 3]) %>% 
  gather(variable, value, F, DFnum, DFden, P) %>% # reshape to make a tabale
  mutate(test_var = paste(test, variable, sep = "_")) %>% 
  select(S, Form, Ind, everything(), -test, -variable, -.id) %>% 
  spread(test_var, value, fill = "-") %>% 
  arrange(S, Form, Ind, term) 
names(anova_df)[5:12] <- gsub("_", "\n", names(anova_df)[5:12])


# save --------------------------------------------------------------------

save.image(file = "output/Data/solve_Year0_issues.RData")
