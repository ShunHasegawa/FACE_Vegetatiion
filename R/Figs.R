theme_set(theme_bw())

# co2 factor
FACE.veg.rslt$co2 <- factor(ifelse(FACE.veg.rslt$ring %in% c(1, 4, 5), "elev", "amb"))
    
# remove unknown spp
veg <- FACE.veg.rslt[-grep("Unknown", FACE.veg.rslt$variable), ,drop = TRUE]

# natrualised(?) -> NA for the time beting
veg$origin[which(veg$origin == "naturalised(?)")] <- NA

# remove rows with value of 0 as stat = bin (the number of 
# cases in each group) will be used
veg <- veg[which(veg$value != 0), ]

veg <- droplevels(veg)

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
  facet_grid(co2 ~ form + PFG + origin, scale = "free_x", space = "free_x", labeller = label_parsed) +
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

