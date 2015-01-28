library(BiodiversityR)

sites$co2 <- factor(ifelse(sites$ring %in% c(1, 4, 5), "elev", "amb"))
sites$Yco2 <- sites$year:sites$co2
sites$YR <- sites$year:sites$ring

?decostand
transDF <- decostand(vg.data, "log")
transDF <- vegdist(transDF, method = "altGower")

capDF <- CAPdiscrim(transDF ~ Yco2, sites, dist = "bray", permutations = 100)
capDF <- CAPdiscrim(transDF ~ YR, sites, dist = "bray", permutations = 100)

summary(capDF$manova$SS)


ldDF <- data.frame(sites, ld1 = capDF$x[, 1], ld2 = capDF$x[, 2])
chulDF <- ddply(ldDF, .(year, ring), 
                function(x) {chx <- chull(x[c("ld1", "ld2")]) 
                             chxDF <- data.frame(rbind(x[chx,], x[chx[1], ]))
                             return(chxDF)})

theme_set(theme_bw())
p <- ggplot(ldDF, aes(x = ld1, y = ld2, col = ring, shape = year))
p + geom_point(size = 5) + geom_polygon(data = chulDF, alpha = .1)

boxplot(ld1 ~ ring * year, data = ldDF)

df <- ddply(ldDF, .(ring, plot, co2), function(x) diff(x$ld1))
boxplot(V1 ~ co2, data = df)
# eco2 vegeataion shifted along ld1

ldDF$block <- factor(recode(ldDF$ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'"))
ldDF$id <- ldDF$ring:ldDF$plot

boxplot(ld1 ~ co2 * year, data = ldDF)
m1 <- lmer(ld1 ~ co2 * year + (1|block) + (1|ring) + (1|id), data = ldDF)
Anova(m1, test.statistic = "F")



#CA

ndf <- vegdist(transDF, method = "bray")

m1 <- cca(ndf)
m1.sm <- summary(m1)

summary(m1.sm)
testDF <- data.frame(sites, ca1 = m1.sm$sites[, 1], ca2 = m1.sm$sites[, 2])
chulDF2 <- ddply(testDF, .(year, ring), 
                function(x) {chx <- chull(x[c("ca1", "ca2")]) 
                             chxDF <- data.frame(rbind(x[chx,], x[chx[1], ]))
                             return(chxDF)})

theme_set(theme_bw())
p <- ggplot(testDF, aes(x = ca1, y = ca2, col = ring, shape = year))
p + geom_point(size = 5) + geom_polygon(data = chulDF2, alpha = .1)





##############
#############
library(MASS)
data(dune)
data(dune.env)

Ordination.model1 <- CAPdiscrim(dune~Management, data=dune.env, 
                                dist="bray",axes=2,m=0)
Ordination.model1
plot1 <- ordiplot(Ordination.model1)
ordisymbol(plot1,dune.env,"Management",legend=FALSE)
