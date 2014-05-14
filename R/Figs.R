theme_set(theme_bw())

# co2 factor
FACE.veg.rslt$co2 <- factor(ifelse(FACE.veg.rslt$ring %in% c(1, 4, 5), "elev", "amb"))
    
# remove unknown spp
veg <- FACE.veg.rslt[-grep("Unknown", FACE.veg.rslt$variable), ,drop = TRUE]

# natrualised(?) -> NA for the time beting
veg$origin[which(veg$origin == "naturalised(?)")] <- NA

veg <- droplevels(veg)

# summary data
RngSum <- ddply(veg, .(year, ring, co2, variable, PFG, form, origin), summarise, frq = sum(value, na.rm = TRUE))
CO2Sum <- ddply(veg, .(year, co2, variable, PFG, form, origin), summarise, frq = sum(value, na.rm = TRUE))
FormSum <- ddply(veg, .(year, co2, PFG, form, origin), summarise, frq = sum(value, na.rm = TRUE))

###########
# All Spp #
###########

PltVeg <- function(data, group, ...){
  # change factor lablles for labbeling in figs
  data$co2 <- factor(data$co2, levels = c("amb", "elev"), labels = c("Ambient", "eCO[2]"))
  
  data$PFG <- factor(data$PFG, 
                     levels = c("c3", "c3_4", "c4", "legume", "Lichen", "moss", "Non_legume", "wood"),
                     labels = c("C[3]", "C[3-4]", "C[4]", "Legume", "Lichen", "Moss", "Non_legume", "wood"))
  data$origin <- factor(data$origin, 
                        levels = c("native", "naturalised"), 
                        labels = c("Native", "Naturalised"))
  data$y <- data[, group]
  p <- ggplot(data, aes(x = variable, y = frq, fill = year))
  p2 <- p + geom_bar(stat = "identity", alpha = 0.6, position = "identity") + 
    facet_grid(y ~ form + PFG + origin, scale = "free_x", space = "free_x", labeller = label_parsed) +
    theme(axis.text.y = element_text(...)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, ...)) +
    labs(x = NULL, y = "Frequency")
}

## Ring ##
p2 <- PltVeg(data = RngSum, group = "ring", size =8) +
  theme(strip.text.x = element_text(size = 6)) +
  expand_limits(x = 4.5) 
#     set minimum size of the graphic areas of each group 
#     some of them are too small to show labels
ggsavePP(filename = "output/figs/FACE_vegetation_Ring", plot = p2, width= 17, height = 11)

## CO2 ##
p2 <- PltVeg(data = CO2Sum, group = "co2", size = 8) +
  theme(strip.text.x = element_text(size = 6)) +
  expand_limits(x = 4.5) 
ggsavePP(filename = "output/figs/FACE_vegetation_CO2", plot = p2, width= 17, height = 11)




###############
# Plant forms #
###############









# plot two years on the same plot with semitransparent bars
p <- ggplot(main.spp, aes(x = factor(variable), y = frq, fill = factor(year)))
p2 <- p + geom_bar(stat = "identity", alpha=0.6, position = "identity") +
  facet_wrap(~ ring)+ 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave(filename = "output/figs/veg.summary.ring.pdf", plot = p2)

p <- ggplot(veg.sms, aes(x = factor(variable), y = frq, fill = factor(year)))
p2 <- p + geom_bar(stat = "identity", alpha=0.6, position = "identity") +
  facet_wrap(~ ring, nrow = 2) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  coord_flip()
ggsave(filename = "output/figs/allsp.ring.pdf", plot = p2, height = 20)




p <- ggplot(veg.sms, aes(x = factor(variable), y = frq, fill = factor(year)))
p2 <- p + geom_bar(stat = "identity", alpha=0.6, position = "identity") +
  facet_grid(ring ~ PFG, scale = "free_x", space = "free_x") + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  

veg.co2$co2 <- factor(veg.co2$co2, levels = c("amb", "elev"), labels= c("Ambient", "eCO[2]"))

p <- ggplot(subset(veg.co2, frq >= 20), aes(x = factor(variable), y = frq, fill = factor(year)))
p2 <- p + geom_bar(stat = "identity", alpha=0.6, position = "identity") +
  facet_grid(co2 ~ PFG, scale = "free_x", space = "free_x", labeller = label_parsed) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  labs(x = "Spp", y = "Frequency")
?labs


p <- ggplot(veg.co2, aes(x = factor(variable), y = frq, fill = factor(year)))
p2 <- p + geom_bar(stat = "identity", alpha=0.6, position = "identity") +
  facet_wrap(~ co2,) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  coord_flip()
ggsave(filename = "output/figs/allsp.co2.pdf", plot = p2, height = 20)

p <- ggplot(veg.co2, aes(x = factor(variable), y = frq, fill = factor(co2)))
p2 <- p + geom_bar(stat = "identity", alpha=0.6, position = "identity") +
  facet_grid(.~ year) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  coord_flip()
ggsave(filename = "output/figs/allsp.co2.year.pdf", plot = p2, height = 20)

# functional groups
FACE.veg.rslt$fgs <- factor(FACE.veg.rslt$form : FACE.veg.rslt$PFG)

veg.fgs <- ddply(FACE.veg.rslt, .(year, ring, fgs), summarise, frq = sum(value))
veg.fgs.co2 <- ddply(FACE.veg.rslt, .(year, co2, fgs), summarise, frq = sum(value))

# plot two years on the same plot with semitransparent bars
p <- ggplot(veg.fgs, aes(x = factor(fgs), y = frq, fill = factor(year)))
p2 <- p + geom_bar(stat = "identity", alpha=0.6, position = "identity") +
  facet_wrap(~ ring ) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  coord_flip()
ggsave(filename = "output/figs/fgs.ringxyr.pdf", plot = p2)

p <- ggplot(veg.fgs.co2, aes(x = factor(fgs), y = frq, fill = factor(year)))
p2 <- p + geom_bar(stat = "identity", alpha=0.6, position = "identity") +
  facet_wrap(~ co2) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  coord_flip()
ggsave(filename = "output/figs/fgs.co2xyr.pdf", plot = p2)

