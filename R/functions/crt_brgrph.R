#########
# Figs #
#########
library(ggplot2)

# year * ring
FACE.veg.rslt$co2 <- factor(ifelse(FACE.veg.rslt$ring %in% c("1", "4", "5"), "elev", "amb"))


veg.sms <- ddply(FACE.veg.rslt, .(year, ring, co2, variable), summarise, frq = sum(value))
veg.co2 <- ddply(FACE.veg.rslt, .(year, co2, variable), summarise, frq = sum(value))
main.spp <- subset(veg.sms, frq > 20) 

theme_set(theme_bw())
p <- ggplot(main.spp, aes(x = factor(variable), y = frq))
p2 <- p + geom_bar(stat = "identity") + facet_grid(year ~ ring) +theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
ggsave(filename = "output/figs/veg.summary.pdf", plot = p2)

# plot two years on the same plot with semitransparent bars
p <- ggplot(main.spp, aes(x = factor(variable), y = frq, fill = factor(year)))
p2 <- p + geom_bar(stat = "identity", alpha=0.6, position = "identity") +
  facet_wrap(~ ring ) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  coord_flip()
ggsave(filename = "output/figs/veg.summary.ring.pdf", plot = p2)

p <- ggplot(veg.sms, aes(x = factor(variable), y = frq, fill = factor(year)))
p2 <- p + geom_bar(stat = "identity", alpha=0.6, position = "identity") +
  facet_wrap(~ ring, nrow = 2) + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  coord_flip()
ggsave(filename = "output/figs/allsp.ring.pdf", plot = p2, height = 20)


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

