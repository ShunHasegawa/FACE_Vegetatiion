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


# Run cap for each dataset.
# Note that llply(or lapply) can't be used as a matrix object has to be difined
# out side of the function where CAPdiscrim is used (it's really weird
# though..... it should be got to do with enrironment setting. but don't know
# how to difine environment at the moment. so just use this for time being)

capList <- list()

for (i in 1:2){
  disMatrix <- DisMatrix_Year[[i]]
  capList[[i]] <- CAPdiscrim(disMatrix ~ co2, 
                                 data = subset(RingSumVeg, 
                                               year == names(DisMatrix_Year)[i]),
                                 permutations = 1000)
}

# create a plot
names(capList) <- c("2012", "2014")

CapYeardf <- ldply(names(capList), function(x) {
  lst <- capList[[x]]
  data.frame(CAP1 = lst$x[, 1], co2 = lst$group, year = x)
})
  
p <- ggplot(CapYeardf, aes(x = co2, y = CAP1))
p + geom_point(size = 4) + facet_grid(.~year)

# each co2 treatment separately----

# create distance matrix
DisMatrix_co2 <- dlply(RingSumVeg, .(co2), 
                       function(x) vegdist(log(x[, SppName] + 1), 
                                           method = "bray"))

# run CAP
capList_co2 <- list()

for (i in 1:2){
  disMatrix <- DisMatrix_co2[[i]]
  capList_co2[[i]] <- CAPdiscrim(disMatrix ~ year, 
                        data = subset(RingSumVeg, co2 == names(DisMatrix_co2)[i]),
                        permutations = 1000)
}

# create a plot
names(capList_co2) <- c("amb", "elev")

CapCo2df <- ldply(names(capList_co2), function(x) {
  lst <- capList_co2[[x]]
  data.frame(CAP1 = lst$x[, 1], year = lst$group, co2 = x)
})

p <- ggplot(CapCo2df, aes(x = year, y = CAP1))
p + geom_point(size = 4) + facet_grid(.~co2)

#######
# PFG #
#######

# each year separately----

# create distance matrix
PFG_DisMatrix_Year <- dlply(RingSumPFGMatrix, .(year), 
                            function(x) vegdist(x[, PFGName], method = "altGower")) # ln(x + 1)

PFGcapList <- list()

for (i in 1:2){
  disMatrix <- PFG_DisMatrix_Year[[i]]
  PFGcapList[[i]] <- CAPdiscrim(disMatrix ~ co2, 
                             data = subset(RingSumPFGMatrix,
                                           year == names(PFG_DisMatrix_Year)[i]),
                             permutations = 1000)
}

# create a plot
names(PFGcapList) <- c("2012", "2014")

PFGCapYeardf <- ldply(names(PFGcapList), function(x) {
  lst <- PFGcapList[[x]]
  data.frame(CAP1 = lst$x[, 1], co2 = lst$group, year = x)
})

p <- ggplot(PFGCapYeardf, aes(x = co2, y = CAP1))
p + geom_point(size = 4) + facet_grid(.~year)

# each co2 treatments separately----

# create distance matrix
PFG_DisMatrix_co2 <- dlply(RingSumPFGMatrix, .(co2),
                           function(x) vegdist(log(x[, PFGName] + 1),
                                               method = "bray")) 
  # altGower causes error when CAP is performed...

PFGcapList_co2 <- list()

for (i in 1:2){
  disMatrix <- PFG_DisMatrix_co2[[i]]
  PFGcapList_co2[[i]] <- CAPdiscrim(disMatrix ~ year, 
                                    data = subset(RingSumPFGMatrix,
                                                  co2 == names(PFG_DisMatrix_co2)[i]),
                                permutations = 1000)
}


# create a plot
names(PFGcapList_co2) <- c("amb", "elev")

PFGCapCo2df <- ldply(names(PFGcapList_co2), function(x) {
  lst <- PFGcapList_co2[[x]]
  data.frame(CAP1 = lst$x[, 1], year = lst$group, co2 = x)
})

p <- ggplot(PFGCapCo2df, aes(x = year, y = CAP1))
p + geom_point(size = 4) + 
  facet_grid(.~co2)

# vegan
vegDF <- capscale(transDF ~ YR, sites, dist = "bray")
plot(vegDF)
summary(vegDF)
