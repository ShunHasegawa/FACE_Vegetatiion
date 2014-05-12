library(mvabund)

# process dataset

# ring
rng.veg <- ddply(veg.face, .(year, ring), 
                 function(x) colSums(x[,  -grep("year|ring|plot|position|cell", names(veg.face))]))
vg.data <- rng.veg[, -grep("year|ring", names(rng.veg))]
sites <- rng.veg[, grep("year|ring", names(rng.veg))]
sites$co2 <- factor(ifelse(sites$ring %in% c("1", "4", "5"), "elev", "amb"))

m1 <- mvabund(vg.data)
m2 <- manyglm(m1 ~ sites$year * sites$co2 , family="negative.binomial", cor.type = "shrink")
m3 <- manyglm(m1 ~ sites$year * sites$co2 , cor.type = "shrink")
anova(m2, nBoot=500, test="wald", p.uni="adjusted", show.time=TRUE)
plot(m2)

# plot
head(plt.veg)
vg.data <- plt.veg[, -grep("year|ring|plot", names(plt.veg))]
sites <- plt.veg[, grep("year|ring|plot", names(plt.veg))]
sites$co2 <- factor(ifelse(sites$ring %in% c("1", "4", "5"), "elev", "amb"))
sites$pos <- factor(sites$ring:sites$plot)
sites$block <- factor(recode(sites$ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'"))

?recode

m1 <- mvabund(vg.data)
m2 <- manyglm(m1 ~ sites$year * sites$ring , family="negative.binomial", cor.type = "shrink")
rslt <- anova(m2, nBoot=500, test="wald", p.uni="adjusted", show.time=TRUE)
rslt.df <- data.frame(rslt["uni.p"])

# block x co2
m2 <- manyglm(m1 ~ sites$year * sites$co2 * sites$block , family="negative.binomial", cor.type = "shrink")
rslt <- anova(m2, nBoot=500, test="wald", p.uni="adjusted", show.time=TRUE)
rslt
# no three-way interaction so remove
m2 <- manyglm(m1 ~ (sites$year + sites$co2 + sites$block)^2 , family="negative.binomial", cor.type = "shrink")
rslt <- anova(m2, nBoot=500, test="wald", p.uni="adjusted", show.time=TRUE)



## functional group ##

# ring
fg.dat <- cast(FACE.veg.rslt, year + ring ~ fgs, sum)
vg.data <- fg.dat[, -grep("year|ring|notknown", names(fg.dat))]
sites <- fg.dat[, grep("year|ring", names(fg.dat))]
sites$co2 <- factor(ifelse(sites$ring %in% c("1", "4", "5"), "elev", "amb"))

m1 <- mvabund(vg.data)
m2 <- manyglm(m1 ~ sites$year * sites$co2, family="negative.binomial", cor.type = "shrink")
anova(m2, nBoot=500, test="wald", p.uni="adjusted", show.time=TRUE)
plot(m2)

# plot
fg.dat <- cast(FACE.veg.rslt, year + ring + plot ~ fgs, sum)
vg.data <- fg.dat[, -grep("year|plot|ring|notknown", names(fg.dat))]
sites <- fg.dat[, grep("year|ring|plot", names(fg.dat))]
sites$co2 <- factor(ifelse(sites$ring %in% c("1", "4", "5"), "elev", "amb"))
sites$pos <- factor(sites$ring:sites$plot)
sites$block <- factor(recode(sites$ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'"))

m1 <- mvabund(vg.data)
m2 <- manyglm(m1 ~ sites$year * sites$co2, family="negative.binomial", cor.type = "shrink")
anova(m2, nBoot=500, test="wald", p.uni="adjusted", show.time=TRUE)
plot(m2)

# co2 x block
m2 <- manyglm(m1 ~ sites$year * sites$co2 * sites$block, family="negative.binomial", cor.type = "shrink")
anova(m2, nBoot=500, test="wald", p.uni="adjusted", show.time=TRUE)
