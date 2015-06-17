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

# remove rows with value of 0 as stat = bin (the number of cases in each group)
# will be used
BarplDF <- subsetD(veg, value != 0)

###########
# All Spp #
###########
pfgLabs <- c("C[3]", "C[3-4]", "C[4]", "Legume", "Lichen", "Moss", "Non_legume", "wood")
orgnLabs <- c("Native", "Introduced")

Spplt <- PltVeg(data = BarplDF, xval = "variable", size = 8) +
  theme(strip.text.x = element_text(size = 6)) +
  expand_limits(x = 4.5) 
#     set minimum size of the graphic areas of each group 
#     some of them are too small to show labels


## Ring ##
p2 <- Spplt +
  facet_grid(ring ~ form + PFG + origin, scale = "free_x", space = "free_x", labeller = label_parsed)
ggsavePP(filename = "output/figs/FACE_vegetation_Ring", plot = p2, 
         width= 17, height = 11)

## CO2 ##
p2 <- Spplt +
  facet_grid(co2 ~ form + PFG + origin, scale = "free_x", space = "free_x", labeller = label_parsed)
ggsavePP(filename = "output/figs/FACE_vegetation_CO2", plot = p2, width= 17, height = 11)

## Dominant spp ##

# DmSppCov <- with(subset(SppSum, variable %in% DmSpp), paste0(variable, "(", Cov, "%)"))

DmSppBar <- within(BarplDF, {
  variable <- factor(ifelse(variable %in% DmSpp, as.character(variable), "Others"))
})

# Organise df
DmSppBar <- within(DmSppBar, {
  co2 <- factor(co2, labels = c("Ambient", expression(eCO[2])))
  year <- factor(year, labels = paste0("Year", 1:3))
  variable <- factor(variable, levels = c("Others", as.character(rev(DmSpp))))})

# create sp labels
# SpLab <- gsub("[.]", "~", levels(DmSppBar$variable))
SpLab <- gsub("[.]", "~", levels(DmSppBar$variable))
SpLab <- parse(text = paste("italic(", SpLab, ")"))

p <- ggplot(DmSppBar, aes(x = year, fill = variable))
p2 <- p + geom_bar(col = "white") +
  scale_fill_discrete(name = "Species", labels = SpLab) + 
  theme(legend.text.align = 0) +
  facet_grid(. ~ co2, labeller = label_parsed)
p2

?facet_grid
?geom_bar

head(DmSppBar)
Spplt <- PltVeg(data = BarplDF, xval = "variable", size = 8) +
  theme(strip.text.x = element_text(size = 6)) +
  expand_limits(x = 4.5) 



#######
# PFG #
#######
PFGplt <- PltVeg(data = BarplDF, xval = "PFG", xlab = "PFG", size = 8) +
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
Orgnplt  <- PltVeg(data = BarplDF, xval = "origin", xlab = "Orgin", size = 8) +
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
  variable <- factor(variable, labels = c("Evenness", "Species Richness", "Diversity (H')"))
  year <- factor(year, levels = c("2013", "2014", "2015"), 
                 labels = c("2013\n(Pre-CO2)", "2014", "2015"))
})

p <- ggplot(Smmry_DivDF, aes(x = year, y = Mean, group = co2, fill = co2))
p2 <- p + 
  geom_errorbar(aes(x = year, ymin = Mean - SE, ymax = Mean + SE), 
                        position = position_dodge(.3),
                width = 0) +
  geom_line(aes(linetype = co2), position = position_dodge(.3)) + 
  geom_point(position = position_dodge(.3), size = 3, shape = 21) + 
  labs(x = "Year", y = NULL) + 
  scale_fill_manual(values = c("black", "white"), 
                    labels = c("Ambient", expression(eCO[2]))) +
  scale_linetype_manual(values = c("solid", "dashed"), 
                        labels = c("Ambient", expression(eCO[2]))) +
  facet_wrap(~variable, scales = "free_y") +
  science_theme +
  theme(legend.position = c(.51, .9),
        legend.key.width = unit(2.5, "lines"))
ggsavePP(filename = "output/figs/FACE_CO2_DiversityIndx", 
         width = 6, height = 3, plot = p2)

# #################################
# ## Dissimilarity between years ##
# #################################
# 
# # all spp ----
# 
# # Compute dissimiliraity for each plot between 2012 and 2013
# disDF_spp <- ddply(plt.veg, .(ring, plot), 
#                function(x) vegdist(x[, SppName], method = "altGower"))
# names(disDF_spp)[names(disDF_spp) == "V1"] <- "Dissim"
# disDF_spp$variable <- "AllSpp"
# 
# # PFG ----
# summary(PfgPlotSum)
# 
# # remove NA and Lichen(which is not assessed in the 2nd year)
# PFGdf <- subsetD(PfgPlotSum, !PFG %in% c(NA, "Lichen"))
# 
# # cast to create matrix-like data frame showing community compoisition
# PFGdf_cst <- dcast(year + ring + plot ~ PFG, data = PFGdf)
# 
# # Compute dissimiliraity for each plot between 2012 and 2013
# pfgName <- names(PFGdf_cst)[!names(PFGdf_cst) %in% c("year", "ring", "plot")]
# 
# disDF_pfg <- ddply(PFGdf_cst, .(ring, plot), 
#                    function(x) 
#                      vegdist(x[, pfgName], method = "altGower"))
# names(disDF_pfg)[names(disDF_pfg) == "V1"] <- "Dissim"
# disDF_pfg$variable <- "PFGs"
# 
# # make plots ----
# # merge the above data frames
# disDF <- rbind(disDF_spp, disDF_pfg)
# 
# # organise data frame
# disDF <- within(disDF, {
#   co2 = factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
# })
# 
# # ring mean
# RngDisDF <- ddply(disDF, .(co2, ring, variable), summarise, Dissim = mean(Dissim))
# 
# # co2 mean, SE
# SmmryDisDF <- ddply(RngDisDF, .(co2, variable), summarise,
#                     Mean = mean(Dissim),
#                     SE = ci(Dissim)[4],
#                     N = sum(!is.na(Dissim)))
# # plot
# 
# # change variable names
# df <- within(SmmryDisDF, {
#   co2 = factor(co2, levels = c("amb", "elev"), 
#                labels = c("Ambient", expression(eCO[2])))
# })
# p <- ggplot(data = SmmryDisDF, aes(x = co2, y = Mean))
# p2 <- p +
#   geom_errorbar(aes(x = co2, ymin = Mean - SE, ymax = Mean + SE)) +
#   geom_point(aes(shape = co2, fill = co2),size = 5) +
#   scale_x_discrete(labels = c("Ambient", expression(eCO[2]))) +
#   scale_shape_manual(values = c(21, 21)) +
#   scale_fill_manual(values = c("black", "white")) +
#   labs(x = NULL, y = "Dissimilarity between 2013 and 2014") +
#   facet_wrap(~ variable, scale = "free_y") +
#   science_theme +
#   theme(legend.position = "none")
# ggsavePP(filename = "output/figs/FACE_Vegetation_Dissimilarity", 
#          width = 4, height = 3, plot = p2)
# 
# # Dissimilarity against moisture ----
# 
# # load soil moisture data
# load("Data/FACE_TDR_ProbeDF.RData")
# summary(FACE_TDR_ProbeDF)
# soilDF <- subsetD(FACE_TDR_ProbeDF, Sample == "vegetation" & 
#                     Date <= as.Date("2013-12-31") &
#                     Date >= as.Date("2013-1-1"))
# soilDF$plot <- as.factor(soilDF$plot)
# MoistDF <- ddply(soilDF, .(ring, plot), summarise, Moist = mean(Moist, na.rm = TRUE))
# plot(Moist ~ ring, data = MoistDF)
# 
# # merge data frame
# Dis_MoistDF <- merge(disDF, MoistDF, by = c("ring", "plot"))
# 
# # regression line for all spp and produce predicted variables
# m1 <- lm(log(Dissim) ~ log(Moist), data = subsetD(Dis_MoistDF, variable == "AllSpp"))
# xv <- seq(min(Dis_MoistDF$Moist), max(Dis_MoistDF$Moist), length.out = 100)
# yv <- exp(predict(m1, list(Moist = xv)))
# PreDF <- data.frame(Moist = xv, Mean = yv, variable = "AllSpp", co2 = "amb")
# 
# # plot
# p <- ggplot(data = Dis_MoistDF, aes(x = Moist, y = Dissim, shape = co2, fill = co2))
# p2 <- p + 
#   geom_point(size = 2) + 
#   facet_grid(variable ~ ., scales = "free_y") +
#   scale_shape_manual(values = c(24, 21), labels = c("Ambient", expression(eCO[2]))) +
#   scale_fill_manual(values = c("black", "white"), labels = c("Ambient", expression(eCO[2]))) +
#   geom_line(aes(x = Moist, y = Mean), data = PreDF) +
#   labs(x = "Soil moisture contents (%)",
#        y = "Dissimilarity between 2013 and 2014") +
#   science_theme +
#   theme(legend.position = c(.8, .89))
# ggsavePP(filename = "output/figs/FACE_DissmVs.Moist", width = 4, height = 4, 
#          plot = p2)

#######
# PCoA #
#######

# ring sum
SppName
RngVegdf <- ddply(veg.face, .(year, ring, co2), function(x) colSums(x[, SppName]))
RngSppDF <- RngVegdf[, SppName]
RngSiteDF <- RngVegdf[, !names(RngVegdf) %in% SppName]

# peform PCoA

# PCoA <- cmdscale(d = vegdist(RngSppDF, method = "altGower"), eig = TRUE)

PCoA <- cmdscale(d = vegdist(log(RngSppDF + 1), method = "bray"), eig = TRUE, k = 3)
PCoA$GOF
PCoA$eig
ExpVars <- round(PCoA$eig[1:3] * PCoA$GOF[2]/(sum(PCoA$eig[1:3])) * 100, 2)

# Note: GOF is given as belows 
# PCoA$GOF[1] at k = 3
sum(PCoA$eig[1:3])/sum(abs(PCoA$eig))

# PCoA$GOF[2] at k = 2
eig2 <- ifelse(PCoA$eig < 0, 0, PCoA$eig)
sum(PCoA$eig[1:3])/sum(eig2)

# PCoA <- cmdscale(d = vegdist(log(RngSppDF + 1), method = "bray"), eig = TRUE, k = 11)


# this transformation and dissimilarity is more interpretable than above...

PCoAs <- cbind(RngSiteDF, PCoA1 = PCoA$points[, 1], 
                          PCoA2 = PCoA$points[, 2],
                          PCoA3 = PCoA$points[, 3])
PCoADF <- melt(PCoAs, id = c("year", "ring", "co2", "PCoA1"))
PCoADF$variable <- factor(PCoADF$variable, levels = c("PCoA2", "PCoA3"), 
                         labels = paste("PCoA",  c(2, 3), " (", ExpVars[2:3], " %)", sep = ""))

# make plots ----

# reorder df according to year
PCoADF <- PCoADF[order(PCoADF$year), ]

# connect the same ring in different years
p <- ggplot(data = PCoADF, aes(x = PCoA1, y = value, shape = year, col = co2))
p2 <- p + 
  geom_point(size = 3) + 
  scale_color_manual(values = c("blue", "red"), labels = c("Ambient", expression(eCO[2]))) +
  science_theme +
  labs(x = paste("PCoA1 (", ExpVars[1]," %)", sep = ""), 
       y = "PCoA axis") +
  theme(legend.position = "top",
        legend.box = "horizontal",
        legend.box.just = "left") +
  facet_grid(. ~variable)
p3 <- p2 + geom_path(aes(group = ring))
ggsavePP(filename = "output/figs/FACE_Vegetation_PCoA_Ring", width = 6, height = 5, 
         plot = p3)

# make area for each pairs of amb and elev for each year
chulDF <- ddply(PCoADF, .(year, co2, variable), 
                function(x) {chx <- chull(x[c("PCoA1", "value")]) 
                             chxDF <- data.frame(rbind(x[chx,], x[chx[1], ]))
                             return(chxDF)})
p4 <- p2 + geom_polygon(data = chulDF, alpha = .1)
ggsavePP(filename = "output/figs/FACE_Vegetation_PCoA_CO2", width = 6, height = 5, 
         plot = p4)

# different colors for each ring
p <- ggplot(data = PCoADF, aes(x = PCoA1, y = value, shape = year, col = ring))
p2 <- p + 
  geom_point(size = 3) + 
  geom_path(aes(group = ring)) +
  science_theme +
  labs(x = paste("PCoA1 (", ExpVars[1]," %)", sep = ""), 
       y = "PCoA axis") +
  theme(legend.position = "top",
        legend.box = "horizontal",
        legend.box.just = "left") +
  facet_grid(. ~variable) 
ggsavePP(filename = "output/figs/FACE_Vegetation_PCoA_EachRing", width = 6, height = 5, 
         plot = p2)
