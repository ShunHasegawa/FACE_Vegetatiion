
# prepare data frames -----------------------------------------------------

# Move Year0 value to a new column to be used as covariate for the analysis
DivDF_year0_list <- llply(DivDF_list, function(x){
  
  DivDF_mlt <- gather(x, variable, value, H, S, J)
  
  DivDF_year0 <- DivDF_mlt %>% # Year0
    filter(year == "Year0") %>%
    select(id, value, variable) %>%
    rename(value0 = value) %>% 
    left_join(filter(DivDF_mlt, year != "Year0"), by = c("id", "variable")) %>% 
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




# analysis ----------------------------------------------------------------


# create models to be tested
div_m_list <- llply(DivDF_year0_list, function(x){
  dlply(x, .(variable), function(y){
    m1 <- lmer(value ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), data = y)
    m2 <- update(m1, ~ . - (1|block))
    if (AICc(m1) >= AICc(m2)) return(m2) else return(m1)
    })
  })

div_m_list <- unlist(div_m_list, recursive = FALSE)
summary(div_m_list)  



# . model diagnosis -------------------------------------------------------


# diagnosing plot
pdf("output/figs/mod_diag_divind.pdf", onefile = TRUE, width = 4, height = 4)
l_ply(names(div_m_list), function(x){
  m <- div_m_list[[x]]
  print(plot(m, main = x))
  qqnorm(resid(m, main = x))
  qqline(resid(m, main = x))
})
dev.off()


# inspect grass_spp.J
div_m_list[["grass_spp.J"]]
d_gsj <- filter(DivDF_year0_list[["grass_spp"]], variable == "J")
m1 <- lmer(value ~ co2 * year + value0 + (1|ring) + (1|id), data = d_gsj)
which.min(resid(m1))
m2 <- update(m1, subset = -8)
plot(m2)
# model is improved
llply(list(m1, m2), function(x) Anova(x, test.statistic = "F"))
# co2 x time interaction seems to be driven by outlier so remove

div_m_list[["grass_spp.J"]] <- m2




# CI and post-hoc test ----------------------------------------------------


# compute 95 CI and post-hoc test
lsmeans_list <- llply(div_m_list, function(x) {
  lsmeans::lsmeans(x, ~ co2 | year)
  })

# 95% CI
CI_dd <- ldply(lsmeans_list, function(x) data.frame(summary(x))) 

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
         Type       = factor(Type, labels = c("All", "Forb", "Grass")),
         year       = factor(year, levels = paste0("Year", 0:3)),
         value_type = "adjusted",
         plot_lab   = as.character(factor(.id, 
                                          labels = paste0("(", letters[c(1, 4, 7, 2, 5, 8,3, 6, 9)], ")")))  # sub-plot label
         ) %>% 
  left_join(div_co2_pval, by = ".id")                                             # merge with pvalues for CO2 term
ci_dd$star[is.na(ci_dd$star)]        <- ""                                        # turn NA to ""
ci_dd$star[ci_dd$.id == "all_spp.H"] <- ""                                        # no CO2xTime interaction, so don't show post-hoc results

         
# Observed vlaues for each variable
div_obs_dd <- ldply(DivDF_list) %>% 
  group_by(.id, year, co2, ring) %>% 
  summarise_each(funs(mean), H, S, J) %>% 
  ungroup() %>% 
  mutate(.id        = factor(.id, labels = c("All", "Forb", "Grass")),
         value_type = "observed") %>% 
  rename(Type = .id) %>% 
  gather(variable, value, H, S, J)


# create fig --------------------------------------------------------------


div_plots <- dlply(ci_dd, .(variable), function(x){
  
  # df for observed values
  d <- div_obs_dd %>%
    filter(variable == unique(x$variable))
    
  
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
  p <- ggplot(x, aes(x = year, y = lsmean)) +
    
    
    geom_vline(xintercept = 1.5, linetype = "dashed") +
    facet_grid(. ~ Type) +
   
    
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
    scale_shape_manual(values = c(16, 17), 
                       labels = c("Ambient", expression(eCO[2]))) +
    scale_linetype_manual(values = c("solid", "dashed"), 
                          labels = c("Ambient", expression(eCO[2]))) +
    scale_color_manual(values = c("black", "grey80"),
                       guide = guide_legend(override.aes = list(linetype = "blank",size = 2))) +
    scale_x_discrete("", labels = NULL, drop = FALSE) +
    
    
    # legend and theme
    science_theme +
    theme(legend.position = "none") +
    
    geom_text(data = plab_d, aes(label = plot_lab), x = -Inf, y = Inf, hjust = -.1, vjust = 1.5, size = 3) +
    geom_text(data = plab_d, aes(label = rr), x =  Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3)
    
  
  return(p)
  })
div_plots[[1]]


# . fine tuning of figure ---------------------------------------------------


# add ylabels
div_ylabs <- c(expression(Diversity~(italic("H'"))), 
               expression(Evenness~(italic("J'"))),
               expression(Species~richness~(italic(S))))


for (i in 1:3){
  div_plots[[i]] <- div_plots[[i]] + 
    labs(x = NULL, y = div_ylabs[i]) +
    theme(axis.title.y = element_text(size = 9))
  }

# remove facet_gird labels
for (i in 2:3) {
  div_plots[[i]] <-  div_plots[[i]] + theme(strip.background = element_blank(), 
                                            strip.text.x = element_blank())
  }

# x lab for the bottom plot
div_plots[[3]] <- div_plots[[3]] + scale_x_discrete("Year", labels = c(0:3))

# add legend in the top plot
div_plots[[1]] <- div_plots[[1]] + theme(legend.position = "top",
                                         legend.box        = "horizontal", 
                                         legend.direction  = "vertical", 
                                         legend.text.align = 0,
                                         legend.margin     = unit(0, "line"),
                                         legend.text = element_text(size = 8))


# set margins
div_margins <- llply(list(c(1, 1, 0, 0), c(0, 1, 0, 0), c(0, 1, 0, 0)),
                     function(x) unit(x, "line"))

for (i in 1:3){
  div_plots[[i]] <-  div_plots[[i]] + theme(plot.margin = div_margins[[i]])
}

# merge plots
div_plot_merged <- rbind(ggplotGrob(div_plots[[1]]), 
                         ggplotGrob(div_plots[[2]]), 
                         ggplotGrob(div_plots[[3]]), 
                         size = "last")

grid.newpage()
grid.draw(div_plot_merged)

ggsavePP(filename = "output/figs/adjusted_diversity_indices", 
         plot = div_plot_merged, width = 6, height  = 6)




# summary table -----------------------------------------------------------

# table for observed values
obs_tbl <- div_obs_dd %>% 
  mutate(.id = paste0(tolower(Type), "_spp.", variable)) %>% 
  select(.id, year, co2, value_type, ring, value) %>% 
  group_by(.id, year, co2, value_type) %>% 
  summarise(M = mean(value))


# bind with adjusted values
div_adjMean_tble <- ci_dd %>% 
  select(.id, co2, year, lsmean, value_type) %>% 
  rename(M = lsmean) %>% 
  bind_rows(obs_tbl) %>% 
  mutate(variable = paste(value_type, co2, sep = "_")) %>% 
  select(-value_type, -co2) %>% 
  spread(key = variable, value = M) %>% 
  mutate(resp     = adjusted_elev / adjusted_amb - 1,
         Type     = gsub("_.*", "", .id),
         variable = str_sub(.id, -1, -1)) %>% 
  group_by(Type, variable, year) %>% 
  summarise_each(funs(round(., 2)), everything(), -.id) %>% 
  select(variable, Type, year, starts_with("observed"), starts_with("adjusted"), 
         resp) %>% 
  ungroup() %>% 
  arrange(variable, Type)


# split df by variable
div_tbl_l <- dlply(div_adjMean_tble, .(variable), function(x) select(x, -variable))


# save as excel
writeWorksheetToFile(file  = "output/table/summary_tbl_diversity.xlsx",  # define file name to be saved
                     data  = div_tbl_l,                                  # writeWorksheetToFile doesn't take dplyr object so turn them into data frames using as.data.frame
                     sheet = names(div_tbl_l))                           # sheet names in excel are defined by object names a list

