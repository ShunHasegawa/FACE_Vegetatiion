library(vegan)
## 2012 & 2014 ##

# process dataset for analysis
names(veg.face)

# plot sum
plt.veg <- ddply(veg.face, .(year, ring, plot), 
                 function(x) colSums(x[, -grep("year|ring|plot|position|cell", 
                                                names(veg.face))]))

vg.data <- plt.veg[, -grep("year|ring|plot|position|cell", names(plt.veg))]
sites <- plt.veg[, grep("year|ring|plot|position|cell", names(plt.veg))]

# analysis
m1 <- cca(vg.data)
m1.sm <- summary(m1)
ca.spp <- m1.sm$species
plot(m1)
names(m1.sm)

# figs
theme_set(theme_bw())
pltCA(data = plt.veg, xv = "CA1", yv = "CA2", shape = "year")
pltCA(data = plt.veg, xv = "CA2", yv = "CA3", shape = "year")
pltCA(data = plt.veg, xv = "CA1", yv = "CA3", shape = "year")

## CA graph for each year ##
d_ply(plt.veg, .(year), plt.CA.yr)
