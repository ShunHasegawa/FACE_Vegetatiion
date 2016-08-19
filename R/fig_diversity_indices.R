
m_list <- llply(DivDF_year0_list, function(x){
  dlply(x, .(variable), function(y){
    lmer(value ~ co2 * year + value0 + (1|block) + (1|ring) + (1|id), data = y)
    })
  })

m_list <- unlist(m_list, recursive = FALSE)
summary(m_list)  

lsmeans_list <- llply(m_list, function(x) {
  summary(lsmeans::lsmeans(x, pairwise ~ co2 | year))
  })


CI_dd <- ldply(lsmeans_list, function(x) data.frame(x$lsmeans)) %>% 
  mutate(Type     = tstrsplit(.id, "[.]")[[1]], 
         variable = tstrsplit(.id, "[.]")[[2]],
         Type = factor(Type, labels = c("All", "Forb", "Grass")))

contrast_dd <- ldply(lsmeans_list, function(x) data.frame(x$contrast)) %>% 
  mutate(Type     = tstrsplit(.id, "[.]")[[1]], 
         variable = tstrsplit(.id, "[.]")[[2]])

div_plots <- dlply(CI_dd, .(variable), function(x){
  p <- ggplot(x, aes(x = year, y = lsmean, fill = co2, group = co2))
  p2 <- p +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0, 
                  position = position_dodge(width = dodgeval)) +
    geom_line(aes(linetype = co2),
              position = position_dodge(width = dodgeval)) +
    geom_point(shape = 21, size = 3, 
               position = position_dodge(width = dodgeval)) +
    scale_fill_manual(values = c("black", "white"), 
                      labels = c("Ambient", expression(eCO[2]))) +
    scale_linetype_manual(values = c("solid", "dashed"), 
                          labels = c("Ambient", expression(eCO[2]))) +
    science_theme +
    theme(legend.position = c(.85, .85)) + 
    scale_x_discrete("", labels = c("", "", "")) +
    facet_grid(. ~ Type)
  return(p2)
  })

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
div_plots[[3]] <- div_plots[[3]] + scale_x_discrete("Year", labels = 1:3)

# set margins
div_margins <- llply(list(c(1, 1, -.5, 1), c(-.5, 1, -.5, 1), c(-.5, 1, 1, 1)),
                     function(x) unit(x, "line"))

for (i in 1:3){
  div_plots[[i]] <-  div_plots[[i]] + theme(plot.margin = div_margins[[i]])
}

# combine plots
grid.newpage()
div_plot_merged <- grid.draw(rbind(ggplotGrob(div_plots[[1]]), 
                                   ggplotGrob(div_plots[[2]]), 
                                   ggplotGrob(div_plots[[3]]), 
                                   size = "last"))

div_plot_merged <- rbind(ggplotGrob(div_plots[[1]]), 
                                   ggplotGrob(div_plots[[2]]), 
                                   ggplotGrob(div_plots[[3]]), 
                                   size = "last")
ggsavePP(filename = "output/figs/adjusted_diversity_indices", 
         plot = div_plot_merged, width = 6, height  = 6)
