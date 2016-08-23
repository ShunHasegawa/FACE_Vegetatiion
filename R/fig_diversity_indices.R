
# prepare data frames -----------------------------------------------------


# create models to be tests
m_list <- llply(DivDF_year0_list, function(x){
  dlply(x, .(variable), function(y){
    m1 <- lmer(value ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), data = y)
    m2 <- update(m1, ~ . - (1|block))
    if (AICc(m1) >= AICc(m2)) return(m2) else return(m1)
    })
  })

m_list <- unlist(m_list, recursive = FALSE)
summary(m_list)  

# compute 95 CI and post-hoc test
lsmeans_list <- llply(m_list, function(x) {
  summary(lsmeans::lsmeans(x, pairwise ~ co2 | year))
  })

# 95% CI
CI_dd <- ldply(lsmeans_list, function(x) data.frame(x$lsmeans)) 

# post-hoc test
contrast_dd <- ldply(lsmeans_list, function(x) data.frame(x$contrast)) %>% 
  mutate(co2 = factor("amb", levels = c("amb", "elev")),
         star = cut(p.value, right = FALSE,
                    breaks = c(0, .1, .05, .01, .001, 1),  
                    labels = c("***", "**", "*", "\u2020", ""))) %>% 
  select(.id, year, co2, p.value, star)

# merge
ci_dd <- left_join(CI_dd, contrast_dd, by = c(".id", "year", "co2")) %>% 
  mutate(Type = tstrsplit(.id, "[.]")[[1]], 
         variable = tstrsplit(.id, "[.]")[[2]],
         Type = factor(Type, labels = c("All", "Forb", "Grass")),
         year = factor(year, levels = paste0("Year", 0:3)))
ci_dd$star[is.na(ci_dd$star)] <- ""
         
# Year0 vlaue to create a box-whisker plot
div_Year0_dd <- ldply(DivDF_list) %>% 
  filter(year == "Year0") %>% 
  mutate(year = factor(year, levels = paste0("Year", 0:3)), 
         .id = factor(.id, labels = c("All", "Forb", "Grass"))) %>% 
  rename(Type = .id) %>% 
  gather(variable, value, H, S, J)


# create fig --------------------------------------------------------------


div_plots <- dlply(ci_dd, .(variable), function(x){
  
  # df for Year0
  d <- div_Year0_dd %>%
    filter(variable == unique(x$variable))
    
  # df for median of Year0
  d_med <- d %>%  
    group_by(Type) %>% 
    summarise(M = median(value)) %>% 
    mutate(Med = "Md[Year0]")
  
  # fig
  dodgeval <- .4
  p <- ggplot(x, aes(x = year, y = lsmean, fill = co2, group = co2))
  p2 <- p +
   
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0, 
                  position = position_dodge(width = dodgeval)) +
    geom_line(aes(linetype = co2), position = position_dodge(width = dodgeval)) +
    geom_boxplot(data = d, aes(x = year, y = value), alpha = .6, 
                 position = position_dodge(.7),  outlier.shape = 21, width = .7, 
                 show.legend = FALSE) +
    geom_point(shape = 21, size = 3, position = position_dodge(width = dodgeval)) +
    geom_hline(data = d_med, aes(yintercept = M, col = Med), alpha = .8) +
    geom_vline(xintercept = 1.5, linetype = "dashed") +
    
    scale_fill_manual(values = c("black", "white"), 
                      labels = c("Ambient", expression(eCO[2]))) +
    scale_linetype_manual(values = c("solid", "dashed"), 
                          labels = c("Ambient", expression(eCO[2]))) +
    scale_color_manual(values = "grey50", labels = expression(Md[Year0])) +
    scale_x_discrete("", labels = NULL, drop = FALSE) +
    geom_text(aes(label = star), fontface = "bold", vjust = -1) +
    
    science_theme +
    theme(legend.position = c(.88, .82),
          legend.margin = unit(-.5, "line")) +
    guides(linetype = guide_legend(order = 1),
           fill     = guide_legend(order = 1),
           col      = guide_legend(order = 2)) +
    
    facet_grid(. ~ Type)
  
  return(p2)
  })
div_plots[[1]]


# fine tuning of figure ---------------------------------------------------


# add ylabels
div_ylabs <- c(expression(Adjusted~diversity~(italic("H'"))), 
               expression(Adjusted~evenness~(italic("J'"))),
               expression(Adjusted~species~richness~(italic(S))))


for (i in 1:3){
  div_plots[[i]] <- div_plots[[i]] + 
    labs(x = NULL, y = div_ylabs[i]) +
    theme(axis.title.y = element_text(size = 9))
  }

# remove facet_gird labels and legends from two plots
for (i in 2:3) {
  div_plots[[i]] <-  div_plots[[i]] + theme(strip.background = element_blank(), 
                                            strip.text.x = element_blank(), 
                                            legend.position="none")
  }

# x lab for the bottom plot
div_plots[[3]] <- div_plots[[3]] + scale_x_discrete("Year", labels = c(0:3))

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
