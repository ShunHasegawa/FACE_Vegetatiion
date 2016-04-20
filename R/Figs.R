theme_set(theme_bw() + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank()))

# define graphic background
science_theme <- theme(panel.border = element_rect(color = "black"),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       legend.position = c(.91, .91),
                       # legend.text = element_text(size = 2),
                       legend.title = element_blank(),
                       legend.background = element_blank(),
                       legend.key = element_blank())

###########
# Barplot #
###########

  # number of count per m2 for each plot
  Veg_Plot <- ddply(veg, .(variable, year, ring, plot, co2, form, PFG, origin), 
                    summarise, value_m2 = sum(value)/4)
  
  # organise labels
  Veg_Plot <- within(Veg_Plot, {
    OrginalVar <- variable
    variable <- gsub("[.]", " ", as.character(variable))
    co2 <- factor(co2, labels = c("Ambient", expression(eCO[2])))
    PFG <- factor(PFG, labels = c("C[3*'\u005F'*grass]","C[4*'\u005F'*grass]", 
                                  "Legume", "Moss", "Non*-legume", "Wood"))
    origin <- factor(origin, labels = c("Native", "Naturalised"))
    })
  
  # Ring mean
  Veg_Ring <- ddply(Veg_Plot, .(variable, OrginalVar, year, ring, co2, form, PFG, origin), 
                    summarise, value = mean(value_m2), N = sum(!is.na(value_m2)))
  
  # Treatment Mean
  veg_co2 <- ddply(Veg_Ring, .(variable, OrginalVar, year, co2, form, PFG, origin), summarise, 
                   Mean = mean(value), 
                   SE = ci(value)[4],
                   N = sum(!is.na(value)))

###########
# All Spp #
###########

  ## Ring ##
    p <- ggplot(Veg_Ring, aes(x = variable, y = value, fill = year))
    p2 <- p + 
      geom_bar(
        alpha    = .4, 
        position = position_dodge(width = .4), 
        stat     = "identity"
        ) +
      theme(
        axis.text.x  = element_text(
                          angle = 90,
                          hjust = 1, 
                          vjust = .5, 
                          face  = "italic"
                          ), 
        strip.text.x = element_text(
                          size  = 6
                          )
        ) +
      expand_limits(
        x = 4
        ) +
      # set minimum size of the graphic areas of each group some of them are too 
      # small to show labels
      labs(
        x = NULL, 
        y = expression(Abundance~(Count~m^'-2'))
        ) + 
      facet_grid(
        ring ~ PFG,
        scale = "free_x",
        space = "free_x",
        labeller = label_parsed
        )
    p2
    ggsavePP(filename = "output/figs/FACE_vegetation_Ring", plot = p2, 
             width= 17, height = 11)

## CO2 ##
posdos <- .6
p <- ggplot(veg_co2, 
            aes(
              x = variable, 
              y = log10(Mean + 1), 
              fill = year)
            )
p2 <- p + 
  geom_bar(
    alpha    = .5, 
    position = position_dodge(width = posdos), 
    stat     = "identity"
    ) +
  geom_errorbar(
    aes(
      x    = variable, 
      ymin = log10(Mean-SE + 1), 
      ymax = log10(Mean+SE + 1),
      col  = year
      ),
    position = position_dodge(width = posdos), 
    width    = 0, 
    size     = .1
    ) +
  theme(
    axis.text.x     = element_text(
                        angle = 90, 
                        hjust = 1, 
                        vjust = .5, 
                        face  = "italic"
                        ), 
    strip.text.x    = element_text(size = 7),
    legend.title    = element_blank(),
    legend.position = c(.85, .87)
    ) +
  expand_limits(
    x = 4
    ) +
  labs(
    x = NULL, 
    y = expression(log[10](Abundance+1~(Count~m^'-2')))
    ) + 
  facet_grid(
    co2      ~ PFG, 
    scale    = "free_x", 
    space    = "free_x",
    labeller = label_parsed
    )
p2
ggsavePP(
  filename = "output/figs/FACE_vegetation_CO2", 
  plot     = p2, 
  width    = 9.6, 
  height   = 5.5
  )

# scatter
posdos <- 1
veg_co2$PFG <- factor(
                  veg_co2$PFG,
                  levels = c("C[3*'_'*grass]","C[4*'_'*grass]", 
                             "Legume","Non*-legume", "Wood","Moss")
                  )
veg_co2$variable2 <- factor(
                        veg_co2$variable, 
                        levels = rev(unique(veg_co2$variable))
                        )

p <- ggplot(
  veg_co2, 
  aes(
    y   = variable2, 
    x   = log10(Mean + 1), 
    col = year
    )
  )

p2 <- p + 
  geom_errorbarh(aes(xmin = log10(Mean-SE + 1), 
                    xmax = log10(Mean+SE + 1), 
                    col = year), 
                height = 0, size = .5, alpha = .6) +
  geom_point(alpha = .8) +
  theme(axis.text.y = element_text(face = "italic", size = 7), 
        strip.text.y = element_text(size = 7),
        legend.title = element_blank(),
        legend.position = c(.88, .67), 
        panel.grid.major = element_line(colour = "grey90", size = .2),
        panel.grid.major.x = element_blank(),
        panel.margin = unit(0, "lines")) +
  expand_limits(y = 4.2) +
  labs(y = NULL, x = expression(log[10](Abundance+1)~(Count~m^'-2'))) + 
  facet_grid(PFG ~ co2, scale = "free_y", space = "free", 
             labeller = label_parsed)
p2
ggsavePP(filename = "output/figs/FACE_vegetation_CO2_Scatter", plot = p2, 
         width= 6, height = 9)

# log scale
RingSummary <- ddply(
  veg, 
  .(variable, year, ring, co2, PFG), 
  summarise,
  value = sum(value, na.rm = TRUE)
  )
RingSummary$logValue <- log10(RingSummary$value + 1)
Co2Summary <- ddply(RingSummary, .(variable, year, co2, PFG), summarise, 
                    Mean = mean(logValue), SE = ci(logValue)[4], N = sum(!is.na(logValue)))
p <- ggplot(Co2Summary, aes(x = variable, y = Mean, fill = year))
p2 <- p + 
  geom_bar(alpha = 0.4, position = position_dodge(width = .4), stat = "identity") + 
  geom_errorbar(aes(x = variable, ymin = Mean + SE, ymax = Mean - SE, col = year), 
                position = position_dodge(width = .4), 
                width = 0, size = .1) +
  facet_grid(co2 ~ PFG, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = "log10(Abundance+1)")
p2
ggsavePP(filename = "output/figs/Fig_Thesis/FACE_vegetation_CO2_log", plot = p2, width= 6, height = 9)


# log scale, PFG
RingSummary_pfg <- ddply(veg, .(year, ring, co2, PFG), summarise, 
                         value = sum(value, na.rm = TRUE))
RingSummary_pfg$logValue <- log10(RingSummary_pfg$value + 1)
Co2Summary_pfg <- ddply(RingSummary_pfg, .(year, co2, PFG), summarise, 
                        Mean = mean(logValue), SE = ci(logValue)[4], N = sum(!is.na(logValue)))
p <- ggplot(Co2Summary_pfg, aes(x = PFG, y = Mean, fill = year))
p2 <- p + 
  geom_bar(alpha = 0.4, position = position_dodge(width = .4), stat = "identity") + 
  geom_errorbar(aes(x = PFG, ymin = Mean + SE, ymax = Mean - SE, col = year), 
                position = position_dodge(width = .4), 
                width = 0, size = .1) +
  facet_grid(co2 ~ ., scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = "log10(Abundance+1)")
p2
ggsavePP(filename = "output/figs/FACE_vegetation_CO2_PFG_log", plot = p2, width= 17, height = 11)


## Difference from the 1st year ##

# subset year1 and get co2 mean and grand mean
Year0DF <- subsetD(RingSummary, year == "Year0")
Year0_co2 <- ddply(Year0DF, .(variable, co2, PFG), summarise, value = sum(value))
Year0_total <- ddply(Year0DF, .(variable, PFG), summarise, value = sum(value)/sum(Year0DF$value))

# reorder according to abundance
Year0_co2$variable <- factor(Year0_co2$variable, 
                             levels = Year0_total$variable[order(Year0_total$value)])
Year0_co2 <- Year0_co2[order(as.numeric(Year0_co2$variable)), ]

# yearly difference
YearDiff <- ddply(RingSummary, .(variable, ring, co2, PFG), function(x) {
  d1 <- with(x, log10(value[year == "Year1"] + 1) - log10(value[year == "Year0"] + 1))
  d2 <- with(x, log10(value[year == "Year2"] + 1) - log10(value[year == "Year0"] + 1))
  d3 <- with(x, log10(value[year == "Year2"] + 1) - log10(value[year == "Year0"] + 1))
  data.frame(year = paste0("Year", 1:3), Dif = c(d1, d2, d3))
})
YearDiff_co2 <- ddply(YearDiff, .(variable, co2, year, PFG), summarise, 
                      Mean = mean(Dif), 
                      SE = ci(Dif)[4],
                      N = sum(!is.na(Dif)))
YearDiff_co2$variable <- factor(YearDiff_co2$variable, 
                                levels = unique(Year0_co2$variable))
YearDiff_co2 <- YearDiff_co2[order(as.numeric(YearDiff_co2$variable)), ]

Year0_total$year <- "Year0"
Year0_total$co2 <- "Year0"

p <- ggplot(YearDiff_co2, aes(x = variable, y = Mean, col = co2))
p2 <- p + 
  geom_point(position = position_dodge(.2), size = 4, alpha = .7) + 
  geom_errorbar(aes(x = variable, ymin = Mean - SE, ymax = Mean + SE, width = 0), 
                position = position_dodge(.2)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(year ~ PFG, scale = "free", space = "free_x") +
  # coord_cartesian(ylim = c(-45, 100))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
p2
ggsavePP(filename = "output/figs/FACE_vegetation_YearDif", plot = p2, width= 17, height = 11)

##################
## Dominant spp ##
##################
  DmSppBar <- subsetD(veg_co2, OrginalVar %in% DmSpp)
  
  # error bar position
  DmSppBar <- ddply(DmSppBar, 
                    .(co2, year), 
                    transform, 
                    CumSum = cumsum(Mean),
                    ystart = cumsum(Mean) - SE, 
                    yend   = cumsum(Mean) + SE
                    )
  # create sp labels
  SpLab <- gsub("[.]", " ", levels(DmSppBar$OrginalVar))
  
  p  <- ggplot(DmSppBar, aes(x = co2, y = Mean, fill = variable))
  p2 <- p + 
    geom_bar(
      size = .1, 
      stat = "identity"
      ) +
    geom_segment(
      aes(
        xend = co2, 
        y    = CumSum, 
        yend = yend, 
        col  = variable
        ), 
      arrow = arrow(
                angle  = 90, 
                length = unit(.1, "inches")
                )
      ) + 
    scale_fill_discrete(
      name   = "Dominant\nSpecies(>70%)", 
      labels = SpLab
      ) + 
    scale_colour_discrete(
      name   = "Dominant\nSpecies(>70%)", 
      labels = SpLab
      ) + 
    scale_x_discrete(
      breaks = c("Ambient", "eCO[2]"),
      labels = c("Ambient", expression(eCO[2]))
      ) +
    science_theme + 
    theme(
      legend.text       = element_text(
                            face = "italic", 
                            size = 9
                            ),
      legend.text.align = 0,
      legend.position   = "bottom",
      legend.title      = element_text(size = 9)
      ) +
    guides(
      fill = guide_legend(nrow = 2)
      ) +
    facet_grid(
      . ~ year
      ) +
    labs(
      x = NULL, 
      y = expression(Abundance~(Count~m^'-2'))
      )
  
  StackBar_DomSpp <- p2
  StackBar_DomSpp
  ggsavePP(
    plot     = StackBar_DomSpp, 
    filename = "output/figs/Fig_Thesis/RDA_DomSppBar", 
    width    = 6, 
    height   = 4
    )

######################################
# For each block for FACE manuscript #
######################################
  DmSppBar_ring <- subsetD(Veg_Ring, OrginalVar %in% DmSpp & year == "Year0")
  DmSppBar_ring$Block <- recode(DmSppBar_ring$ring, 
                              "c(1, 2) = 'Block A'; c(3, 4) = 'Block B'; c(5, 6) = 'Block C'")
  p <- ggplot(DmSppBar_ring, aes(x = ring, y = value, fill = variable))
  p2 <- p + 
    geom_bar(size = .1, stat = "identity") +
    scale_fill_discrete(name = "Dominant\nSpecies", labels = SpLab) + 
    scale_colour_discrete(name = "Dominant\nSpecies", labels = SpLab) + 
    science_theme + 
    theme(legend.text = element_text(face = "italic", size = 9),
          legend.text.align = 0,
          legend.position = "bottom",
          legend.title = element_text(size = 9)) +
    guides(fill = guide_legend(nrow = 2)) +
    facet_grid(. ~ Block, scales = "free_x", shrink = TRUE) +
    labs(x = "FACE ring number", y = expression(Abundance~(Count~m^'-2')))
  ggsavePP(plot = p2, filename = "output/figs/FACE_BlockDominantSpecies", 
         width = 6, height = 4)

#######
# PFG # 
#######
pfgLabs <- c(expression(C[3*'\u005F'*grass]), 
             expression(C[4*'\u005F'*grass]), 
             "Legume", "Moss", expression(Non*-legume), "Wood")
  
## Ring ##
Veg_PFG <- ddply(Veg_Ring, .(year, co2, ring, PFG), summarise, value =  sum(value))

p <- ggplot(Veg_PFG, aes(x = PFG, y = value, fill = year))
p2 <- p + 
  geom_bar(alpha = .4, position = position_dodge(width = .4), stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = .5), 
        strip.text.x = element_text(size = 6)) +
  scale_x_discrete(labels = c(expression(C[3*'\u005F'*grass]), 
                              expression(C[4*'\u005F'*grass]), 
                              "Legume", "Moss", expression(Non*-legume), "Wood")) +
  expand_limits(x = 4) +
  labs(x = NULL, y = expression(Abundance~(Count~m^'-2'))) + 
  facet_grid(ring ~ .,labeller = label_parsed)
p2
ggsavePP(filename = "output/figs/FACE_PFG_Ring", plot = p2, width= 8, height = 6)

## CO2 ##
Veg_PFG_co2 <- ddply(
                  Veg_PFG, 
                  .(year, co2, PFG), 
                  summarise, 
                  value = mean(value)
                  )

p <- ggplot(Veg_PFG_co2, 
            aes(
              x = PFG, 
              y = value, 
              fill = co2)
            )
p2 <- p + 
  geom_bar(
    alpha    = .4, 
    position = position_dodge(width = .4), 
    stat     = "identity"
    ) +
  theme(
    axis.text.x  = element_text(
                      angle = 90, 
                      hjust =1, 
                      vjust = .5
                      ), 
    strip.text.x = element_text(
                      size = 6
                      )
    ) +
  scale_x_discrete(
    labels = c(expression(C[3*'\u005F'*grass]), 
               expression(C[4*'\u005F'*grass]), 
               "Legume", 
               "Moss", 
               expression(Non*-legume), 
               "Wood")
    ) +
  labs(
    x = NULL, 
    y = expression(Abundance~(Count~m^'-2'))
    ) + 
  facet_grid(. ~ year,labeller = label_parsed)
p2
ggsavePP(filename = "output/figs/FACE_PFG_CO2", plot = p2, width= 8, height = 6)

####################
## Stack bar plot ##
####################
  # Ring sum
  Veg_RingSum <- ddply(veg, 
                       .(year, co2, ring, PFG), 
                       summarise, 
                       value = sum(value)
                       )
  head(Veg_RingSum)
  
  # co2 sum
  Veg_CO2Sum <- ddply(veg, 
                      .(year, co2, PFG), 
                      summarise, 
                      value = sum(value)
                      )
  
  # approximate SE
  Veg_RingSum_cst <- dcast(year + co2 + ring ~ PFG, 
                           data = Veg_RingSum
                           )
  Veg_RingSum_cst$Total <- rowSums(Veg_RingSum_cst[, 4:9])
  
  ratio <- function(d, w) {
    c3SE          <- sum(d$c3 * w)/sum(d$Total * w)
    c4SE          <- sum(d$c4 * w)/sum(d$Total * w)
    legumeSE      <- sum(d$legume * w)/sum(d$Total * w)
    mossSE        <- sum(d$moss * w)/sum(d$Total * w)
    Non_legumeSE  <- sum(d$Non_legume * w)/sum(d$Total * w)
    woodSE        <- sum(d$wood * w)/sum(d$Total * w)
    
    c(c3SE, c4SE, legumeSE, mossSE, Non_legumeSE, woodSE)
  }
  
  PFG_Fraction <- ddply(Veg_RingSum_cst, 
                        .(year, co2), 
                        function(x) {
                          boRes     <- summary(
                                          boot::boot(x,
                                                     ratio,
                                                     R     = 999, 
                                                     stype = "w")
                                          )
                          boRes$PFG <- unique(Veg_RingSum$PFG)
                          return(boRes)
                          }
                        )
  
  # Organise DF
  PFG_Fraction <- within(
                    PFG_Fraction, {
                      PFG <- factor(PFG, 
                                    levels = c("c3", "c4", "legume", 
                                               "Non_legume", "wood", "moss")
                                    )
                      co2 <- factor(co2, 
                                    labels = c("Ambient", expression(eCO[2]))
                                    )
                      }
                    )
  PFG_Fraction <- PFG_Fraction[order(as.numeric(PFG_Fraction$PFG)), ]
  PFG_Fraction <- ddply(PFG_Fraction, 
                        .(year, co2), 
                        transform, 
                        CumSum = cumsum(original),
                        ystart = cumsum(original) - bootSE, 
                        yend   = cumsum(original) + bootSE
                        )
  
  pfgLabs <- c(expression(C[3]~grass), expression(C[4]~grass), "Legume", "Forb", 
               "Woody plants", "Moss")
  
  p <- ggplot(PFG_Fraction,
              aes(
                  x = co2, 
                  y = original, 
                  fill = PFG
                  )
              )
  
  p2 <- p + 
    geom_bar(
      stat = "identity"
      ) + 
    geom_segment(
      aes(
        xend = co2, 
        y    = CumSum, 
        yend = yend, 
        col  = PFG
        ), 
      arrow = arrow(
                angle = 90, 
                length = unit(.1, "inches")
                ),
      size  = .4
      ) + 
    scale_fill_discrete(
      name   = "PFG", 
      labels = pfgLabs
      ) +
    scale_color_discrete(
      name   = "PFG", 
      labels = pfgLabs
      ) +
    scale_x_discrete(
      labels = c("Ambient", expression(eCO[2]))
      ) +
    science_theme +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.6, "line"),
      legend.title    = element_text()
      ) +
    facet_grid(. ~ year, 
               labeller = label_parsed
               ) +
    labs(
      x = NULL, 
      y = "Proportion"
      )
  p2
  StackBar_PFG <- p2
  StackBar_PFG
  ggsavePP(
    plot     = StackBar_PFG, 
    filename = "output/figs/Fig_Thesis/StackBar_PFG", 
    width    = 6, 
    height   = 3.5
    )

# ########################
# # Native or introduced #
# ########################
# Orgnplt  <- PltVeg(data = BarplDF, xval = "origin", xlab = "Orgin", size = 8) +
#   theme(strip.text.x = element_text(size = 7)) +
#   expand_limits(x = 2)
#   
# ## Ring ##
# p <- Orgnplt +
#   facet_grid(ring ~ form, scale = "free_x", space = "free_x", labeller = label_parsed, margins= "form") 
# ggsavePP(filename = "output/figs/FACE_Origin_Ring", plot = p, width= 8, height = 6)
# 
# ## CO2 ##
# p <- Orgnplt +
#   facet_grid(co2 ~ form, scale = "free_x", space = "free_x", labeller = label_parsed, margins= "form") 
# ggsavePP(filename = "output/figs/FACE_Origin_CO2", plot = p, width= 8, height = 6)

#####################
# Figure for thesis #
#####################

###############
## PFG ratio ##
###############

# C3:C4 & legume:Non_legume----
# susbset df of Grass and Form
GFdf <- subsetD(veg, form %in% c("Grass", "Forb") & PFG != "c3_4")
GFdf$prop <- factor(GFdf$form, labels = c("Legume/(Legume+Non_legume)", "C3/(C3+C4)"))

# Ring mean
SmmryPFGRing <- ddply(GFdf, .(year, co2, ring, prop), summarise, 
                    Mean = sum(value[PFG %in% c("c3", "legume")]) / sum(value)
                    )

# overall mean
SmmryPFGRAll <- ddply(GFdf, .(year, co2, prop), summarise, 
                      Mean = sum(value[PFG %in% c("c3", "legume")]) / sum(value))

# native:introduced----
# Ring mean
SmmryOrgnRing <- ddply(subset(veg, !is.na(origin)), .(year, co2, ring), summarise, 
  Mean = sum(value[origin == "native"])/sum(value),
  prop = "Native/(Native+Introduced)")

# Overall mean
SmmryOrgnAll <- ddply(subset(veg, !is.na(origin)), .(year, co2), summarise, 
  Mean = sum(value[origin == "native"])/sum(value),
  prop = "Native/(Native+Introduced)")

# merge the above data frames----
SmmryPropDfRing <- rbind.fill(SmmryPFGRing, SmmryOrgnRing)
SmmryPropDfRing$co2 <- factor(SmmryPropDfRing$co2, 
                              labels = c("Ambient", expression(eCO[2])))

SmmryPropDfAll <- rbind.fill(SmmryPFGRAll, SmmryOrgnAll)
SmmryPropDfAll$co2 <- factor(SmmryPropDfAll$co2, 
                             labels = c("Ambient", expression(eCO[2])))

# create a plot
p <- ggplot(SmmryPropDfRing, aes(x = year, y = Mean))
p2 <- p + 
  geom_line(data = SmmryPropDfAll, aes(group = co2, linetype = co2)) +
  geom_point(aes(fill = co2, group = co2), alpha = .7, shape = 21, size = 3) + 
  scale_linetype_manual(values = c(1, 2), 
                        labels = c("Ambient", expression(eCO[2])))+
  scale_fill_manual(values = c("black", "white"), 
                    labels = c("Ambient", expression(eCO[2]))) +
  facet_wrap(~prop, ncol = 2, scales = "free_y") +
  labs(y = "Proportion", x = NULL) +
  science_theme +
  theme(strip.text.x = element_text(size = 7),
        legend.position = c(.75, .25),
        legend.key.width = unit(2.5, "lines"))
ggsavePP(filename = "output//figs/FACE_CO2_PFGProportion", plot = p2,  
         width = 5, height = 4)

#######################
## Diversity indices ##
#######################
summary(DivDF)

# Mean and SE
DivDF_mlt <- melt(DivDF, id = c("year", "block", "co2", "ring", "plot", "id"))
RngSmmry_DivDF <- ddply(DivDF_mlt, .(year, co2, ring, variable), summarise, value = mean(value))
Smmry_DivDF <- ddply(RngSmmry_DivDF, .(year, co2, variable), summarise, 
                     Mean = mean(value),
                     SE = ci(value)[4],
                     N = sum(!is.na(value)))

# make a plot

# change variable names
Smmry_DivDF <- within(Smmry_DivDF, {
  variable <- factor(variable, 
                     levels = c("S", "H", "J"),
                     labels = c("(a) Species richness", "(b) Diversity", "(c) Evenness"))
})

p <- ggplot(Smmry_DivDF, aes(x = year, y = Mean, group = co2, fill = co2))
p2 <- p + 
  geom_errorbar(aes(x = year, ymin = Mean - SE, ymax = Mean + SE), 
                        position = position_dodge(.3),
                width = 0) +
  geom_line(aes(linetype = co2), position = position_dodge(.3)) + 
  geom_point(position = position_dodge(.3), size = 3, shape = 21) + 
  labs(x = NULL, y = NULL) + 
  scale_fill_manual(values = c("black", "white"), 
                    labels = c("Ambient", expression(eCO[2]))) +
  scale_linetype_manual(values = c("solid", "dashed"), 
                        labels = c("Ambient", expression(eCO[2]))) +
  facet_wrap(~variable, scales = "free_y") +
  science_theme +
  theme(legend.position = c(.18, .87),
        legend.key.width = unit(2.5, "lines"))
ggsavePP(filename = "output/figs/Fig_Thesis/FACE_CO2_DiversityIndx", 
         width = 6, height = 2.5, plot = p2)

##############################
# Fig to see evenness change #
##############################
TreatSum <- ddply(veg, .(year, co2, variable), summarise, value = sum(value))

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
ggsavePP(plot = pp, filename = "output/figs/EvennessChange", width = 6.5, 
         height = 7.5)
