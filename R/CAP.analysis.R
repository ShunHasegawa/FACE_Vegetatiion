#######
# CAP #
#######


###############
# All species #
###############

# each year separately----

# create distance matrix
DisMatrix_Year <- dlply(RingSumVeg, .(year), 
                        function(x) vegdist(x[, SppName], method = "altGower")) # ln(x + 1)

# 2012 (due to the environmental setting in CAPdiscrim,, can't use llply)
Dis2012 <- DisMatrix_Year[[1]]
Cap2012 <-  CAPdiscrim(Dis2012 ~ co2, 
                       data = subset(RingSumVeg, year == 2012), 
                       permutations = 1000)
# 2014
Dis2014 <- DisMatrix_Year[[1]]
Cap2014 <-  CAPdiscrim(Dis2014 ~ co2, 
                       data = subset(RingSumVeg, year == 2014), 
                       permutations = 1000)

# merge
CapYeardf <- data.frame(CAP1 = c(Cap2012$x, Cap2014$x),
                        rbind(subset(RingSumVeg, year == 2012), 
                              subset(RingSumVeg, year == 2014))[, c("year", "co2")])
# create a plot
p <- ggplot(CapYeardf, aes(x = co2, y = CAP1))
p + geom_point(size = 4) + facet_grid(.~year)



disMatrix <- DisMatrix_Year[["2012"]]
capList <- llply(names(DisMatrix_Year),
                 function(x) {
                   disMatrix <- DisMatrix_Year[[x]]
                   environment(disMatrix) <- .GlobalEnv
                   CAPdiscrim(disMatrix ~ co2, 
                              data = subset(RingSumVeg, year == x),permutations = 1000)
})

# run cap for each dataset. Note that a matrix (or data frame) has to be 
# deifined outside of the function where CAPdiscrim is used; hence disMatrix is 
# difined out side of f() below. (it's really weird though..... it should be got
# to do with enrironment setting. but don't know how to difine environment at
# the moment. so just use this for time being)
capList <- llply(names(DisMatrix_Year),
                 function(x) {
                   disMatrix <- DisMatrix_Year[[x]]
                   f <- function(){
                     CAPdiscrim(disMatrix ~ co2, 
                                data = subset(RingSumVeg, year == x),
                                permutations = 1000)
                   }
                   return(f())})



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
