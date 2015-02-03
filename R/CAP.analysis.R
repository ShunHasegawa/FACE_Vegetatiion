sites$co2 <- factor(ifelse(sites$ring %in% c(1, 4, 5), "elev", "amb"))
sites$YR <- sites$year:sites$ring

# transDF <- vegdist(vg.data, method = "altGower") # ln(x + 1)

capDF <- CAPdiscrim(transDF ~ YR, sites, dist = "bray", permutations = 10)

# create a plot
ldDF <- data.frame(sites, ld1 = capDF$x[, 1], ld2 = capDF$x[, 2])
chulDF <- ddply(ldDF, .(year, ring), 
                function(x) {chx <- chull(x[c("ld1", "ld2")]) 
                             chxDF <- data.frame(rbind(x[chx,], x[chx[1], ]))
                             return(chxDF)})

theme_set(theme_bw())
p <- ggplot(ldDF, aes(x = ld1, y = ld2, col = ring, shape = year))
p + geom_point(size = 5) + geom_polygon(data = chulDF, alpha = .1)

# eco2 vegeataion seems to shift towards left along ld1 in 2014
boxplot(ld1 ~ year*ring, data = ldDF)

# vegan
vegDF <- capscale(transDF ~ YR, sites, dist = "bray")
plot(vegDF)
summary(vegDF)
