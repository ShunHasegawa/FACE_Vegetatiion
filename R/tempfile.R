dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10)),field=gl(3,1) )
dat
Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(12)/2
Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(12)/2
total <- Agropyron + Schizachyrium
library(lattice)
dotplot(total ~ NO3, dat, jitter.x=TRUE, groups=field,
        type=c('p','a'), xlab="NO3", auto.key=list(columns=3, lines=TRUE) )

Y <- data.frame(Agropyron, Schizachyrium)
mod <- metaMDS(Y)
plot(mod)

a1 <- adonis(Y ~ NO3, data=dat, strata=dat$field, perm=1e3)
a2 <- adonis(Y ~ NO3, data=dat, perm=1e3)
summary(a1)
head(a1$f.perms)

nested.npmanova(Y ~ NO3 + field, data=dat)
Y

str(sites)
sites$id <- sites$ring:sites$plot
sites$block <- recode(sites$ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")

str(transDF)

a3 <- adonis(transDF ~ ring * year + id, data = sites, strata = sites$id, method = "bray")

a3 <- adonis(transDF ~ year * block * co2, data = sites, method = "bray")
msYC <- 0.03483
msYBC <- 0.05219
crF <- msYC/msYBC

a4 <- adonis(transDF ~ year * co2, data = sites, method = "bray")
a5 <- adonis(transDF ~ year * co2 + ring + year:ring, data = sites, method = "bray")
a6 <- adonis(transDF ~ year * co2 + ring , data = sites, method = "bray", strata = sites$ring)
a7 <- adonis(transDF ~ co2 * block * year + id, data = sites, method = "bray", strata = sites$id)

a7 <- adonis(transDF ~ co2 * year + co2:block:year, data = sites, method = "bray", 
             strata = with(sites, co2:block:year))
a7

b8 <- nested.npmanova(transDF ~ co2 + co2:block, data = sites, method = "bray")
a8 <- adonis(transDF ~ co2 + co2:block, data = sites, method = "bray", 
             strata = with(sites, co2:block))
a8
140.6/(1068.3/4)

# nested.npmanova returns correct F model but doesn't take interaction

# adonis can take interaction but desn't return correct F (P values are correct 
# though). Use adonis and specify appropriate permutable units by strata, and
# obtain P, SumsofSquares, MS and DF. then manuary calculate correct F values.


transDF <- vegdist(vg.data, method = "altGower") # ln(x + 1)
# Year
Ym <- adonis(transDF ~ co2 * block * year, data = sites, method = "bray",
             strata = with(sites, block:year))

  
Ym <- adonis(transDF ~ year + id, data = sites, method = "bray",
             strata = with(sites, id))

# CO2
Cm <- adonis(transDF ~ co2 + co2:block, data = sites, method = "bray",
             strata = with(sites, co2:block))


Cm <- adonis(transDF ~ co2 * block + id, data = sites, method = "bray",
             strata = with(sites, id))


Cm <- adonis(transDF ~ co2 + block + id, data = sites, method = "bray",
             strata = with(sites, id))


Cm2 <- nested.npmanova(transDF ~ co2 + id, data = sites, method = "bray")


140.6/(653.57/2)
140.6/(414.72/2)


# Year x CO2
YCm <- adonis(transDF ~ co2 * block * year, data = sites, method = "bray",
              strata = with(sites, co2:block:year))

MSyc = 43.542
MSycb = 174.658
Fyc = MSyc/MSycb
# it's quite small....

# ful model
fm <- adonis(transDF ~ co2 * block * year, data = sites, method = "bray",  
             strata = with(sites, id))

idm <- adonis(transDF ~ co2*year*block + id, data = sites, method = "bray")


YCm2 <- adonis(transDF ~ co2 * block * year + id, data = sites, method = "bray",
              strata = with(sites, co2:block:year))


# co2
Cm2  <- adonis(transDF ~ co2 * block * year, data = sites, method = "bray",
              strata = with(sites, co2:block))

Cm3  <- adonis(transDF ~ co2 * block, data = sites, method = "bray",
              strata = with(sites, co2:block))
# Cm3 seems correct

# year
Ym <- adonis(transDF ~ year  + id, data = sites, method = "bray",
             strata = with(sites, id))
Ym





# Year x CO2 x Block
YCBm <- adonis(transDF ~ co2 * year * block, data = sites, method = "bray")















a8 <- adonis(transDF ~ co2, data = sites, method = "bray", strata = with(sites, co2:block))
a8
140.6/(1068.3/4)


### try Anderson and Millar
df <- expand.grid(year = factor(c(1, 2)),
                  loc = paste("loc", letters[1:4], sep = ""),
                  ha = c("k", "b"),
                  sites = factor(c(1:4)),
                  trans = c(1:10))
df$sites <- df$loc:df$sites

spDF <- matrix(runif(640*5, 0, 100), 640)
m1 <- adonis(spDF ~  year * loc * ha, data = df, method = "bray", 
             strata = df$year:df$sites:df$ha)



