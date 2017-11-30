
# num per 1mx1m -----------------------------------------------------------


# . prepare df ------------------------------------------------------------


# number of count per m2 for each plot
Veg_Plot <- veg_FullVdf %>% 
  group_by(variable, year, ring, plot, co2, form, PFG, origin) %>%
  summarise(value_m2 = sum(value)/4) %>% 
  ungroup() %>% 
  mutate(form       = car::recode(form, "c('Fern', 'Moss') = 'Moss/Fern'"),
         OrginalVar = variable,
         variable   = gsub("[.]", " ", as.character(variable)),                # e.g.  Araujia.sericifera ->  Araujia sericifera
         co2        = factor(co2, labels = c("Ambient", expression(eCO[2]))),
         PFG        = ifelse(form == "Grass", paste(PFG, form, sep = "_"),     # relabel PFG: (C3grass, C4grass, legume, non-legeume, wood, Moss/Fern)
                             ifelse(form == "Forb", as.character(PFG),
                                    ifelse(form == "Moss/Fern", "Moss/Fern", "Woody"))),
         PFG        = factor(PFG, 
                             levels = c("c3_Grass", "c4_Grass", "legume", 
                                        "Non_legume", "Woody", "Moss/Fern"),
                             labels = c("C[3]~graminoid","C[4]~graminoid",
                                        "Legume", "Non*-legume", "Woody", "Moss/Fern")))
summary(Veg_Plot)
  
  
# Ring mean
Veg_Ring <- Veg_Plot %>% 
  group_by(variable, OrginalVar, year, ring, co2, form, PFG, origin) %>% 
  summarise_each(funs(value = mean(.), N = sum(!is.na(.))), value_m2) %>% 
  ungroup()
summary(Veg_Ring)    


# Treatment Mean
veg_co2 <- Veg_Ring %>% 
  group_by(variable, OrginalVar, year, co2, form, PFG, origin) %>% 
  summarise_each(funs(Mean = mean(.), SE = ci(.)[4], N = sum(!is.na(.))), 
                 value) %>% 
  ungroup()
                 
                 


# . create a fig ----------------------------------------------------------


posdos <- 1
veg_co2 <- mutate(veg_co2, variable2 = factor(variable, levels = rev(unique(variable))))
                  
                  
p <- ggplot(veg_co2, aes(y = variable2, x = log10(Mean + 1),  col = year)) + 
  
  
  facet_grid(PFG ~ co2, scale = "free_y", space = "free",
             labeller = label_parsed) +


  geom_errorbarh(aes(xmin = log10(Mean-SE + 1), 
                     xmax = log10(Mean+SE + 1), 
                     col  = year), 
                height = .6, size = .5, alpha = .6) +
  geom_point(alpha = .8) +
  
  
  theme(axis.text.y        = element_text(face = "italic", size = 6.6), 
        strip.text.y       = element_text(size = 7, angle = 0),
        legend.title       = element_blank(),
        legend.position    = c(1.05, .57), 
        legend.key         = element_blank(),
        legend.background  = element_rect(colour = "black"),
        panel.grid.major   = element_line(colour = "grey90", size = .2),
        panel.grid.major.x = element_blank(),
        panel.margin       = unit(0, "lines")) +
  labs(y = NULL, x = expression(log[10](Cover+1)~(Count~m^'-2'))) 
p
ggsavePP(filename = "output/figs/FACE_vegetation_CO2_Scatter", plot = p, 
         width= 6, height = 8)




# PFG composition ---------------------------------------------------------

# total abundance

# ring sum
Veg_RingSum <- veg_FullVdf %>% 
  mutate(PFG = ifelse(form == "Grass", paste(PFG, form, sep = "_"),     # relabel PFG: (C3grass, C4grass, legume, non-legeume, Woody, Moss/Fern)
                      ifelse(form == "Forb", as.character(PFG), 
                             as.character(form))),
         PFG = factor(PFG, 
                      levels = c("c3_Grass", "c4_Grass", "legume", "Non_legume",
                                 "Woody", "Fern", "Moss"))) %>% 
  group_by(year, ring, co2, PFG) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  arrange(year, ring)


# co2 sum
Veg_CO2Sum <- Veg_RingSum %>% 
  group_by(year, co2, PFG) %>% 
  summarise(value = sum(value)) %>% 
  ungroup()
  
  
# approximate SE
Veg_RingSum_cst <- Veg_RingSum %>% 
  spread(PFG, value) %>% 
  mutate(Total = rowSums(.[, 4:10]))
  

ratio <- function(d, w) {
  c3SE          <- sum(d$c3_Grass * w)/sum(d$Total * w)
  c4SE          <- sum(d$c4_Grass * w)/sum(d$Total * w)
  legumeSE      <- sum(d$legume * w)/sum(d$Total * w)
  Non_legumeSE  <- sum(d$Non_legume * w)/sum(d$Total * w)
  woodSE        <- sum(d$Wood * w)/sum(d$Total * w)
  fernSE        <- sum(d$Fern * w)/sum(d$Total * w)
  mossSE        <- sum(d$Moss * w)/sum(d$Total * w)
  
  c(c3SE, c4SE, legumeSE, Non_legumeSE, woodSE, fernSE, mossSE)
}

PFG_Fraction <- ddply(Veg_RingSum_cst, .(year, co2), function(x) {
  boRes     <- summary(boot::boot(x, ratio, R = 999, stype = "w"))
  boRes$PFG <- unique(Veg_RingSum$PFG)
  return(boRes)
})

# Organise DF; get cumulative sum to generate stcked bar plot for each PFG for
# each treatment for each year
PFG_Fraction <- PFG_Fraction %>% 
  arrange(PFG) %>% 
  group_by(year, co2) %>% 
  mutate(CumSum = cumsum(original),
         ystart = cumsum(original) - bootSE, 
         yend   = cumsum(original) + bootSE) %>% 
  ungroup() %>% 
  mutate(co2 = factor(co2, labels = c("Ambient", expression(eCO[2])))) %>% 
  arrange(year, co2, PFG)


pfgLabs <- c(expression(C[3]~grass), expression(C[4]~grass), "Legume", 
             "Non_legume", "Woody", "Fern", "Moss")

p <- ggplot(PFG_Fraction,aes(x = co2, y = original, fill = PFG)) +
  facet_grid(. ~ year, labeller = label_parsed) +
  
  geom_bar(stat = "identity") + 
  geom_segment(aes(xend = co2, y = CumSum, yend = yend, col  = PFG), 
               arrow = arrow(angle = 90, length = unit(.1, "inches")), size  = .4) + 
  
  
  scale_fill_discrete(name  = "PFG", labels = pfgLabs) +
  scale_color_discrete(name = "PFG", labels = pfgLabs) +
  scale_x_discrete(labels = c("Ambient", expression(eCO[2]))) +
  
  
  science_theme +
  theme(legend.position   = "bottom",
        legend.title      = element_text(),
        legend.key.width = unit(.5, "line"),
        legend.text.align = 0) +
  labs(x = NULL, y = "Proportion")
p
StackBar_PFG <- p
StackBar_PFG
ggsavePP(plot = StackBar_PFG, filename = "output/figs/Fig_Thesis/StackBar_PFG", 
         width    = 6, height   = 3.5)




# scatter plot for PFG abundance ------------------------------------------


# prepare df
pfg_co2_mean <- Veg_RingSum %>% 
  group_by(year, co2, PFG) %>% 
  summarise_each(funs(Mean = mean(. / 4, na.rm = TRUE), 
                      SE   = ci(. / 4, na.rm = TRUE)[4], 
                      N    = sum(!is.na(.))),
                 value) %>% 
  ungroup()
  
p <- ggplot(pfg_co2_mean, aes(x = as.numeric(year), y = Mean, fill = co2, group = co2)) +
  facet_wrap( ~ PFG) +
  
  geom_errorbar(aes(x = as.numeric(year), ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(width = 0.3), width = 0) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_point(position = position_dodge(width = 0.3), shape = 21) +
  
  
  scale_fill_manual(values = c("white", "black"), labels = c("Ambient", expression(eCO[2]))) +
  scale_x_continuous(labels = 0:3) +
  
  
  labs(x = "Year", y = expression(Abundance~(Count~m^'-2'))) +
  science_theme +
  theme(legend.title = element_blank(), 
        legend.position = c(0.6, 0.1))
p

ggsavePP(plot = p,filename = "output/figs/PFG_abundance_scatter", 
         width = 6.5,  height = 6.5)




# change in evenness ------------------------------------------------------


TreatSum <- ddply(veg_FullVdf, .(year, co2, variable), summarise, value = sum(value))

EvennessPlot <- dlply(TreatSum, .(co2), function(x) {
  dd <- x
  
  # 1st year df  
  dfyear <- subsetD(dd, year == "Year1")
  
  # sum
  sumdf <- ddply(dfyear, .(variable), value = sum(value))
  
  # species order
  sporder <- with(sumdf, variable[order(value)])
  
  # reorder
  dd$variable <- factor(dd$variable, levels = sporder)
  
  p <- ggplot(dd, aes(variable, y = log10(value +1)))
  p2 <- p + 
    geom_bar(aes(fill = year), 
             stat = "identity", 
             position = "identity", 
             alpha = .5) +
    science_theme + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2, size = 5),
          legend.position = c(.1, .85)) +
    labs(x = NULL, y = expression(log[10](Frequency+1)))
  return(p2)
})


pp <- arrangeGrob(EvennessPlot[[1]] + ggtitle("Ambient"), 
                  EvennessPlot[[2]] + ggtitle(expression(eCO[2])))
grid.draw(pp)
ggsavePP(plot = pp, filename = "output/figs/EvennessChange", width = 6.5, 
         height = 7.5)
