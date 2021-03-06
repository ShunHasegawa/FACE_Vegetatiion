
# ANCOVA ------------------------------------------------------------------


# . prepare data frames -----------------------------------------------------

# Move Year0 value to a new column to be used as covariate for the analysis
DivDF_year0_list <- llply(DivDF_list, function(x){
  
  DivDF_mlt <- gather(x, variable, value, H, S, J)
  
  DivDF_year0 <- DivDF_mlt %>% # Year0
    filter(year == "Year0") %>%
    select(ring, value, variable) %>%
    rename(value0 = value) %>% 
    left_join(filter(DivDF_mlt, year != "Year0"), by = c("ring", "variable")) %>% 
    filter(!is.na(value)) %>% 
    mutate(value0_log  = log(value0 + 1),
           value0_sqrt = sqrt(value0))
  
  return(DivDF_year0)
})

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

# scatter plot agianst Year0
l_ply(names(DivDF_year0_list), function(x){
  d_ply(DivDF_year0_list[[x]], .(variable), function(y){
    figtitle <- paste(x, unique(y$variable), sep = "-")
    plot(value ~ value0, pch = 19, col = year, data = y, main = figtitle)
  })
})




# . analysis ----------------------------------------------------------------


# prepare data frames
f_d  <- DivDF_year0_list$forb_spp
fh_d <- filter(f_d, variable == "H")
fj_d <- filter(f_d, variable == "J")
fs_d <- filter(f_d, variable == "S")

g_d  <- DivDF_year0_list$grass_spp
gh_d <- filter(g_d, variable == "H")
gj_d <- filter(g_d, variable == "J")
gs_d <- filter(g_d, variable == "S")


# forb
plot(value ~ value0, data = fh_d, col = year, pch = 19)
plot(value ~ value0, data = fj_d, col = year, pch = 19)
plot(value ~ value0, data = fs_d, col = year, pch = 19)

fh_m1 <- lmer(value ~ co2 * year + value0 + (1|ring), data = fh_d)
fj_m1 <- lmer(value ~ co2 * year + value0 + (1|ring), data = fj_d)
fs_m1 <- lmer(value ~ co2 * year + value0 + (1|ring), data = fs_d)

# grass
plot(value ~ value0, data = gh_d, col = year, pch = 19)
plot(value ~ value0, data = gj_d, col = year, pch = 19)
plot(sqrt(value) ~ sqrt(value0), data = gs_d, col = year, pch = 19)

gh_m1 <- lmer(value ~ co2 * year + value0 + (1|ring), data = gh_d)
gj_m1 <- lmer(value ~ co2 * year + value0 + (1|ring), data = gj_d)
gs_m1 <- lmer(sqrt(value) ~ co2 * year + value0_sqrt + (1|ring), data = gs_d)


div_m_list <- list(grass_spp.H = gh_m1,
                   grass_spp.J = gj_m1,
                   grass_spp.S = gs_m1,
                   forb_spp.H  = fh_m1,
                   forb_spp.J  = fj_m1,
                   forb_spp.S  = fs_m1)

names(div_m_list)
summary(div_m_list) 
llply(div_m_list, VarCorr)


# . model diagnosis -------------------------------------------------------


# diagnosing plot
pdf("output/figs/mod_diag_divind.pdf", onefile = TRUE, width = 4, height = 4)
l_ply(names(div_m_list), function(x){
  m <- div_m_list[[x]]
  print(plot(m, main = x))
  qqPlot(resid(m, main = x))
})
dev.off()


# inspect forb S
div_m_list[["forb_spp.S"]]
d_gsj <- filter(DivDF_year0_list[["forb_spp"]], variable == "S")
m1 <- lmer(value ~ co2 * year + value0 + (1|ring), data = d_gsj)
plot(m1)
qqPlot(resid(m1))
which.min(resid(m1))
m2 <- update(m1, subset = -15)
plot(m2)
qqPlot(resid(m2))
# model is improved
llply(list(m1, m2), function(x) Anova(x, test.statistic = "F"))
  # no difference so keep the original




# . CI and post-hoc test ----------------------------------------------------


# compute 95 CI and post-hoc test (only need grass and forb species)
nl <- names(div_m_list)
div_m_list <- div_m_list[grepl("grass|forb", nl)]  # extract grass and forb species

lsmeans_list <- llply(div_m_list, function(x) {
  lsmeans::lsmeans(x, ~ co2 | year, type = "response")
  })

# 95% CI
CI_dd <- ldply(lsmeans_list, function(x){
  d <- data.frame(summary(x))
  if(any(names(d) == "response")) d <- rename(d, lsmean = response)
  return(d)
  }) 


# post-hoc test
contrast_dd <- ldply(lsmeans_list, function(x) {
  data.frame(summary(pairs(x)[1:3], adjust = "none"))
}) %>% 
  mutate(co2 = factor("amb", levels = c("amb", "elev")),
         star = get_star(p.value)) %>% 
  select(.id, year, co2, p.value, star)

# CO2 effect
div_aov_df <- ldply(div_m_list, function(x) tidy(Anova(x, test.statistic = "F")),  # Anova result of models
                    .progress = "text")
div_co2_pval <- div_aov_df %>%                                                     # get p-values for CO2 term
  filter(term == "co2") %>% 
  mutate(co2star = get_star(p.value)) %>% 
  select(.id, co2star)


# merge
ci_dd <- left_join(CI_dd, contrast_dd, by = c(".id", "year", "co2")) %>% 
  mutate(Type       = tstrsplit(.id, "[.]")[[1]], 
         variable   = tstrsplit(.id, "[.]")[[2]],
         Type       = dplyr::recode_factor(Type, grass_spp = "Graminoid", forb_spp = "Forb"),
         year       = factor(year, levels = paste0("Year", 0:3)),
         value_type = "adjusted",
         vt         = paste(Type, variable, sep = "_"),
         plot_lab   = mapvalues(vt, 
                                c("Graminoid_H", "Forb_H", 
                                  "Graminoid_J", "Forb_J",
                                  "Graminoid_S", "Forb_S"), 
                                paste0("(", letters[1:6], ")"))) %>%             # sub-plot label
  select(-vt) %>% 
  left_join(div_co2_pval, by = ".id")                                            # merge with pvalues for CO2 term
ci_dd$star[is.na(ci_dd$star)]            <- ""                                   # turn NA to ""
ci_dd$star[ci_dd$.id == "grass_spp.S"]   <- ""                                   # no CO2xTime interaction, so don't show post-hoc results
ci_dd$star[ci_dd$.id == "grass_spp.J"]   <- ""                                   # no CO2xTime interaction, so don't show post-hoc results
# ci_dd$star[ci_dd$.id == "grass_spp.S"] <- ""                                       # no CO2xTime interaction, so don't show post-hoc results

         
# Observed vlaues for each variable
div_obs_dd <- ldply(DivDF_list) %>%
  filter(.id != "all_spp") %>% 
  group_by(.id, year, co2, ring) %>% 
  summarise_each(funs(mean), H, S, J) %>% 
  ungroup() %>% 
  mutate(Type       = dplyr::recode_factor(.id, grass_spp = "Graminoid", forb_spp = "Forb"),
         value_type = "observed") %>% 
  gather(variable, value, H, S, J)


# . create fig --------------------------------------------------------------


div_plots <- dlply(ci_dd, .(variable), function(x){
  
  # df for observed values
  d <- filter(div_obs_dd, Type %in% unique(x$Type) & variable == unique(x$variable))    
  
  # df for plot labels and response ratios
  plab_d <- x %>% 
    group_by(Type, plot_lab, co2, co2star) %>% 
    summarise(value = mean(lsmean)) %>% 
    group_by(Type, plot_lab, co2star) %>% 
    summarise(rr = value[co2 == "elev"] / value[co2 == "amb"] - 1) %>% 
    mutate(rr = ifelse(rr >= 0,
                       paste0("RR= +", format(rr, digits = 0, nsmall = 2), co2star), 
                       paste0("RR= ", format(rr, digits = 0, nsmall = 2), co2star)))
  
  # fig
  dodgeval <- .4
  yscale   <- 1.2
  p <- ggplot(x, aes(x = year, y = lsmean)) +
    
    facet_grid(. ~ Type) +
    geom_vline(xintercept = 1.5, linetype = "dashed") +
    
  
  # observed
    geom_point(data = d, 
               aes(x = year, y = value, shape = co2, col = value_type), 
               size = 2, fill = "grey80", 
               position = position_dodge(dodgeval)) +
    
    
    # adjusted
    geom_line(aes(linetype = co2, group = co2), 
              position = position_dodge(width = dodgeval)) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, group = co2), width = 0, 
                  position = position_dodge(width = dodgeval)) +
    geom_point(aes(shape = co2, col = value_type), size = 2.5, position = position_dodge(width = dodgeval)) +
    geom_text(aes(y = upper.CL, label = star), fontface = "bold", vjust = -.1) +
    
    
    # scaaling
    scale_shape_manual(values = c(16, 15), 
                       labels = c("Ambient", expression(eCO[2]))) +
    scale_linetype_manual(values = c("solid", "dashed"), 
                          labels = c("Ambient", expression(eCO[2]))) +
    scale_color_manual(values = c("black", "grey80"),
                       guide = guide_legend(override.aes = list(linetype = "blank",size = 2))) +
    scale_x_discrete("Year", labels = 0:4, drop = FALSE) +
    ylim(min(c(d$value, x$lower.CL), na.rm = TRUE), 
         max(c(d$value, x$upper.CL), na.rm = TRUE) * yscale) +
    
    
    # legend and theme
    science_theme +
    theme(legend.position = "none") +
    
    geom_text(data = plab_d, aes(label = plot_lab), x = -Inf, y = Inf, 
              hjust = -.1, vjust = 1.5, size = 3, fontface = "bold") +
    geom_text(data = plab_d, aes(label = rr), x =  Inf, y = Inf, hjust = 1.1, 
              vjust = 1.5, size = 3)
    
  
  return(p)
  })

div_plots[[3]]
names(div_plots)


# . fine tuning of figure ---------------------------------------------------


# add ylabels
div_ylabs <- c(expression(italic("H'")~(plot^'-1')), 
               expression(italic("J'")~(plot^'-1')),
               expression(italic(S)~(plot^'-1')))
div_ylabs <- rep(div_ylabs, 2)

for (i in 1:3){
  div_plots[[i]] <- div_plots[[i]] + 
    labs(x = NULL, y = div_ylabs[i]) +
    theme(axis.title.y = element_text(size = 9))
  }


# add legend 
div_plots[[1]] <- div_plots[[1]] +
  theme(legend.position   = "top",
        legend.direction  = "vertical", 
        legend.box        = "horizontal",
        legend.text.align = 0,
        legend.margin     = unit(0, "line"),
        legend.text       = element_text(size = 8))


# remove facet_grid labels
for (i in 2:3){
  div_plots[[i]] <- div_plots[[i]] +
    theme(strip.text.x = element_blank())
}

# remove x axis label
for (i in 1:2){
  div_plots[[i]] <- div_plots[[i]] +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank())
}

# # set margins
# div_margins <- llply(list(c(1, 1, -.5, .5), c(0, 1, -.5, .5), c(0, 1, 0, .5)),
#                      function(x) unit(x, "line"))
# 
# for (i in 1:3){
#   div_plots[[i]] <-  div_plots[[i]] + theme(plot.margin = div_margins[[i]])
# }




# . merge figures ---------------------------------------------------------


div_mplot <- do.call(rbind, llply(div_plots, ggplotGrob))
grid.newpage()
grid.draw(div_mplot)

ggsavePP(filename = "output/figs/adjusted_diversity_indices", 
         plot = div_mplot, width = 4.5, height  = 6)




# . summary table -----------------------------------------------------------

# table for observed values
obs_tbl <- div_obs_dd %>% 
  mutate(.id = paste0(tolower(Type), "_spp.", variable)) %>% 
  select(.id, year, co2, value_type, ring, value) %>% 
  group_by(.id, year, co2, value_type) %>% 
  summarise(M = mean(value, na.rm = TRUE))


# bind with adjusted values

div_pval <- ci_dd %>% 
  mutate(.id = paste0(tolower(Type), "_spp.", variable),
         p.value = round(p.value, 4)) %>% 
  filter(!is.na(p.value)) %>% 
  select(.id, year, p.value)


div_adjMean_tble <- ci_dd %>%
  mutate(.id = paste0(tolower(Type), "_spp.", variable)) %>% 
  select(.id, co2, year, lsmean, value_type) %>% 
  rename(M = lsmean) %>% 
  bind_rows(obs_tbl) %>% 
  mutate(variable = paste(value_type, co2, sep = "_")) %>% 
  select(-value_type, -co2) %>% 
  spread(key = variable, value = M) %>% 
  left_join(div_pval) %>% 
  mutate(resp     = adjusted_elev / adjusted_amb - 1,
         Type     = gsub("_.*", "", .id),
         variable = str_sub(.id, -1, -1)) %>% 
  group_by(Type, variable, year) %>% 
  mutate_each(funs(round(., 2)), everything(), -.id, -p.value) %>% 
  select(variable, Type, year, starts_with("observed"), starts_with("adjusted"), 
         resp, p.value) %>% 
  ungroup() %>% 
  arrange(variable, Type)



# split df by variable
div_tbl_l <- dlply(div_adjMean_tble, .(variable), function(x) select(x, -variable))

# save as excel
writeWorksheetToFile(file        = "output/table/summary_tbl_diversity.xlsx",  # define file name to be saved
                     data        = div_tbl_l,                                  # writeWorksheetToFile doesn't take dplyr object so turn them into data frames using as.data.frame
                     sheet       = names(div_tbl_l),                           # sheet names in excel are defined by object names a list
                     clearSheets = TRUE)




# Forb pre & post CO2 -----------------------------------------------------

forb_pre_d <- filter(DivDF_list$forb_spp, year == "Year0")
forb_post_d <- filter(DivDF_list$forb_spp, year != "Year0")

# . Pre CO2 ---------------------------------------------------------------

create_trans_boxplot(H ~ co2, data = forb_pre_d)
forb_pre_H_m1 <- lm(H ~ co2, data = forb_pre_d)
summary(forb_pre_H_m1)
par(mfrow = c(2, 2))
plot(forb_pre_H_m1)

create_trans_boxplot(S ~ co2, data = forb_pre_d)
forb_pre_S_m1 <- lm(S ~ co2, data = forb_pre_d)
par(mfrow = c(2, 2))
plot(forb_pre_S_m1)
summary(forb_pre_S_m1)

create_trans_boxplot(J ~ co2, data = forb_pre_d)
forb_pre_J_m1 <- lm(J ~ co2, data = forb_pre_d)
summary(forb_pre_J_m1)


# . Post CO2 ----------------------------------------------------------------

create_trans_boxplot(H ~ co2 * year, data = forb_post_d)
forb_post_H_m1 <- lmer(H ~ co2 * year + (1|ring), data = forb_post_d)
Anova(forb_post_H_m1, test.statistic = "F")
plot(forb_post_H_m1)
qqPlot(resid(forb_post_H_m1))


create_trans_boxplot(S ~ co2 * year, data = forb_post_d)
forb_post_S_m1 <- lmer(S ~ co2 * year + (1|ring), data = forb_post_d)
plot(forb_post_S_m1)
qqPlot(resid(forb_post_S_m1))
Anova(forb_post_S_m1, test.statistic = "F")


create_trans_boxplot(J ~ co2 * year, data = forb_post_d)
forb_post_J_m1 <- lmer(J ~ co2 * year + (1|ring), data = forb_post_d)
plot(forb_post_J_m1)
qqPlot(resid(forb_post_J_m1))
Anova(forb_post_J_m1, test.statistic = "F")




# . summary table -Forb ---------------------------------------------------

forb_pre_m_l <- list(pre_h  = forb_pre_H_m1, 
                     pre_j  = forb_pre_J_m1,
                     pre_S  = forb_pre_S_m1)
forb_pre_anv <- ldply(forb_pre_m_l, function(x) tidy(Anova(x))) %>% 
  select(-sumsq) %>% 
  filter(term == "co2") %>% 
  mutate(Df.res = 4)



forb_post_m_l <- list(post_h = forb_post_H_m1, 
                      post_j = forb_post_J_m1,
                      post_s = forb_post_S_m1)
forb_post_anv <- ldply(forb_post_m_l, function(x) tidy(Anova(x, test.statistic = "F")))



forb_anv <- bind_rows(forb_pre_anv, forb_post_anv) %>%
  rename(Fval = statistic) %>% 
  mutate(treat   = tstrsplit(.id, "_")[[1]],
         ind     = tstrsplit(.id, "_")[[2]],
         dfs     = paste0("F", "(", df, ", ", round(Df.res, 0), ")"),
         Fval    = round(Fval, 3),
         p.value = round(p.value, 3)) %>% 
  select(-.id, -df, -Df.res) %>%
  gather(variable, value, Fval, p.value, dfs) %>% 
  mutate(term_var = paste(term, variable, sep = "_")) %>% 
  select(-term, -variable) %>% 
  spread(term_var, value) %>% 
  select(ind, treat, starts_with("co2_"), starts_with("year_"), starts_with("co2:year")) %>% 
  arrange(ind, desc(treat))
forb_anv[is.na(forb_anv)] <- "-" 
write.csv(forb_anv, file = "output/table/analysis_PrePost_forb_diveristy.csv", 
          row.names = F)
