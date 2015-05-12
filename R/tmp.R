load("output/Data/FACE_Vegetation_PFG_2015.RData")
tdf <- VegRes15

DF2013 <- subset(tdf, year == 2013)
# combine sep and dec and remove duplicates
DF2013_merged <- ddply(DF2013, .(variable, year, ring, plot, position, cell), value = sum(value))
head(DF2013_merged)


tdf$y




m <- with(VegRes15, year:month)

tdf_Sum <- ddply(tdf, .(year, month, ym, variable), summarise, value = sum(value))

total_sum <- ddply(tdf, .(variable), summarise, value = sum(value))
OrderSpp <- total_sum$variable[order(total_sum$value, decreasing = TRUE)]

tdf_Sum$SpRnk <- with(tdf_Sum, as.numeric(factor(variable, levels = OrderSpp)))
tdf_Sum$variable <- factor(tdf_Sum$variable, levels = OrderSpp)
tdf_Sum <- tdf_Sum[order(tdf_Sum$SpRnk), ]
tdf_Sum$ym <- factor(tdf_Sum$ym, levels = c("2013:September", "2013:December", "2014:December", "2015:December"))

total_sum

DomSp <- total_sum$variable[total_sum$value >= 10]
p <- ggplot(subset(tdf_Sum, variable %in% DomSp), aes(x = variable, y = value))
p <- ggplot(tdf_Sum, aes(x = variable, y = value))
p2 <- p + geom_point(size = 2) + 
  facet_grid(ym ~ .) + 
  theme(axis.text.x = element_text(angle = 90))
p2
