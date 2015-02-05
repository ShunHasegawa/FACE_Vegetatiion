# remove rows with value of 0 as stat = bin (the number of cases in each group)
# will be used
veg <- veg[which(veg$value != 0), ]

veg <- droplevels(veg)

theme_set(theme_bw())

###########
# All Spp #
###########
pfgLabs <- c("C[3]", "C[3-4]", "C[4]", "Legume", "Lichen", "Moss", "Non_legume", "wood")
orgnLabs <- c("Native", "Introduced")

Spplt <- PltVeg(xval = "variable", size = 8) +
  theme(strip.text.x = element_text(size = 6)) +
  expand_limits(x = 4.5) 
#     set minimum size of the graphic areas of each group 
#     some of them are too small to show labels


## Ring ##
p2 <- Spplt +
  facet_grid(ring ~ form + PFG + origin, scale = "free_x", space = "free_x", labeller = label_parsed)
ggsavePP(filename = "output/figs/FACE_vegetation_Ring", plot = p2, width= 17, height = 11)

## CO2 ##
p2 <- Spplt +
  facet_grid(co2 ~ form + PFG + origin, scale = "free_x", space = "free_x", labeller = label_parsed)
ggsavePP(filename = "output/figs/FACE_vegetation_CO2", plot = p2, width= 17, height = 11)

#######
# PFG #
#######
PFGplt <- PltVeg(xval = "PFG", xlab = "PFG", size = 8) +
  theme(strip.text.x = element_text(size = 9)) +
  scale_x_discrete(breaks = pfgLabs,
                   labels=c(expression(C[3]), 
                            expression(C[3-4]),
                            expression(C[4]),
                            "Legume", "Lichen", "Moss", "Non_legume", "wood"))

## Ring ##
p <- PFGplt +
  facet_grid(ring ~ form, scale = "free_x", space = "free_x", labeller = label_parsed)
ggsavePP(filename = "output/figs/FACE_PFG_Ring", plot = p, width= 8, height = 6)

## CO2 ##
p <- PFGplt +
  facet_grid(co2 ~ form, scale = "free_x", space = "free_x", labeller = label_parsed)
ggsavePP(filename = "output/figs/FACE_PFG_CO2", plot = p, width= 8, height = 6)


########################
# Native or introduced #
########################
Orgnplt  <- PltVeg(xval = "origin", xlab = "Orgin", size = 8) +
  theme(strip.text.x = element_text(size = 7)) +
  expand_limits(x = 2)
  
## Ring ##
p <- Orgnplt +
  facet_grid(ring ~ form, scale = "free_x", space = "free_x", labeller = label_parsed, margins= "form") 
ggsavePP(filename = "output/figs/FACE_Origin_Ring", plot = p, width= 8, height = 6)

## CO2 ##
p <- Orgnplt +
  facet_grid(co2 ~ form, scale = "free_x", space = "free_x", labeller = label_parsed, margins= "form") 
ggsavePP(filename = "output/figs/FACE_Origin_CO2", plot = p, width= 8, height = 6)

#####################
# Figure for thesis #
#####################
# C4:C3 & legume:Non_legume

PfgPlotSum <- ddply(veg, .(year, ring, plot, PFG), summarise, value = sum(value))

# subset required rows
CL_PfgPlotSum <- subsetD(PfgPlotSum, PFG %in% c("c3", "c4", "legume", "Non_legume"))

pfgR <- ddply(CL_PfgPlotSum, .(year, ring, plot),
              function(x) {
               c43R <- with(x, value[PFG == "c4"]/value[PFG == "c3"])
               lgNonlgR <- with(x, value[PFG == "legume"]/value[PFG == "Non_legume"])
               # there one c4 observation which is 0, so put c4 on numerator
               return(data.frame(c43R, lgNonlgR))
             })

# native:introduced
OgnPlotSum <- ddply(veg, .(year, ring, plot, origin), summarise, value = sum(value))
orgn <- ddply(subset(OgnPlotSum, origin %in% c("naturalised", "native")), 
                     .(year, ring, plot),
                     function(x){
                       IntNatR <- with(x, value[origin == "naturalised"]/value[origin == "native"])
                       return(data.frame(IntNatR))
                     })

# biodiversity indices
summary(veg.face)

SiteName <- c("year", "ring", "plot", "position", "cell")
SppName <- names(veg.face)[!names(veg.face) %in% SiteName]

plt.veg <- ddply(veg.face, .(year, ring, plot), function(x) colSums(x[, SppName]))

vegDF <- plt.veg[, SppName]
siteDF <- plt.veg[, c("year", "ring", "plot")]

DivDF <- within(siteDF,{
  H <- diversity(vegDF) # Shannon's index
  S <- specnumber(vegDF) # number of spp
  J <- H/log(S)  # Pielou's evenness
})

# merge above data frames
BioDivDF <- Reduce(function(...) merge(..., by = c("year", "ring", "plot")), 
                   list(pfgR, orgn, DivDF))

# melt for making figure
BioDivDF_mlt <- melt(BioDivDF, id = c("year", "ring", "plot"))

# Mean and Se
RngSmmryBioDivDF <- ddply(BioDivDF_mlt, .(year, ring, variable), summarise, 
                          value = mean(value, na.rm = TRUE))
RngSmmryBioDivDF <- within(RngSmmryBioDivDF, {
  co2 = factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
  block = recode(ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
})

SmmryBioDivDF <- ddply(RngSmmryBioDivDF, .(year, co2, variable), summarise,
                       Mean = mean(value),
                       SE = ci(value)[4],
                       N = sum(!is.na(value)))

# make a plot

# define graphic background
science_theme <- theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       legend.position = c(.91, .91),
                       # legend.text = element_text(size = 2),
                       legend.title = element_blank())

# change variable names
vars <- c("C4:C3", "Legume:Non_legume",
          "Introduced:Native", "Evenness", "Species Richness", "Diversity (H')")
df <- within(SmmryBioDivDF, {
  variable <- factor(variable, 
                     levels = c("c43R", "lgNonlgR", "IntNatR", "J", "S", "H"),
                     labels = vars)
  year <- factor(year, levels = c("2012", "2014"), 
                 labels = c("2013\n(Pre-CO2)", "2014"))
})

p <- ggplot(df, aes(x = year, y = Mean, group = co2))
p2 <- p + 
  geom_errorbar(aes(x = year, ymin = Mean - SE, ymax = Mean + SE), 
                        position = position_dodge(.5),
                width = 0) +
  geom_point(aes(fill = co2, shape = co2), position = position_dodge(.5), size = 5) + 
  labs(x = "Year", y = NULL) + 
  scale_shape_manual(values = c(21, 21), labels = c("Ambient", expression(eCO[2]))) +
  scale_fill_manual(values = c("black", "white"), 
                    labels = c("Ambient", expression(eCO[2]))) +
  facet_wrap(~variable, scales = "free_y") +
  science_theme
ggsavePP(filename = "output/figs/FACE_CO2_Biodiversity", width = 6, height = 5, plot = p2)
