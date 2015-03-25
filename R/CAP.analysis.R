#######
# CAP #
#######

###############
# All species #
###############

# each year separately----

# create distance matrix
DisMatrix_Year <- dlply(RingSumVeg, .(year), 
                        function(x) vegdist(log(x[, SppName] + 1), method = "bray"))


# Run cap for each dataset.
# Note that llply(or lapply) can't be used as a matrix object has to be difined
# out side of the function where CAPdiscrim is used (it's really weird
# though..... it should be got to do with enrironment setting. but don't know
# how to difine environment at the moment. so just use this for time being)

capList <- list()

for (i in 1:3){
  disMatrix <- DisMatrix_Year[[i]]
  capList[[i]] <- CAPdiscrim(disMatrix ~ co2, 
                                 data = subset(RingSumVeg, 
                                               year == names(DisMatrix_Year)[i]),
                                 permutations = 1000)
}

# create a plot
names(capList) <- c("2012", "2014", "2015")

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
capSppScore_co2 <- list()
for (i in 1:2){
  disMatrix <- DisMatrix_co2[[i]]
  envDF <- subset(RingSumVeg, co2 == names(DisMatrix_co2)[i])
  capList_co2[[i]] <- CAPdiscrim(disMatrix ~ year, 
                        data = envDF,
                        permutations = 1000)
  capSppScore_co2[[i]] <- add.spec.scores(capList_co2[[i]], log(envDF[, SppName] + 1))$cproj
}

# create a plot
names(capList_co2) <- c("amb", "elev")
names(capSppScore_co2) <- c("amb", "elev")

# site scores
CapCo2df <- ldply(c("amb", "elev"), function(x) {
  lst <- capList_co2[[x]]
  data.frame(CAP1 = lst$x[, 1], CAP2 = lst$x[, 2], year = lst$group, co2 = x)
})

# species score (correlation between ld1 and species abundance (after
# transformation))
CapCo2SppScore <- ldply(c("amb", "elev"), function(x) {
  df <- capSppScore_co2[[x]]
  spdf <- data.frame(CAP1 = df[, 1], CAP2 = df[, 2], Spp = row.names(df), co2 = x, year = "Scaled species Score")
  spdf <- spdf[complete.cases(spdf), ]
  # set threshhold for correlation
  subset(spdf, abs(CAP1) > .6|abs(CAP2) > .6)
})  

p <- ggplot(CapCo2df, aes(x = CAP1, y = CAP2, fill = year, col = year))
p + geom_point(size = 4) + 
  geom_text(data = CapCo2SppScore, aes(x = CAP1 * 6, y = CAP2 * 2, label = Spp),
            size = 4) +
  facet_grid(.~co2)

#######
# PFG #
#######

# each year separately----

# create distance matrix
PFG_DisMatrix_Year <- dlply(RingSumPFGMatrix, .(year), 
                            function(x) vegdist(log(x[, PFGName] + 1), method = "bray"))

PFGcapList <- list()

for (i in 1:3){
  disMatrix <- PFG_DisMatrix_Year[[i]]
  PFGcapList[[i]] <- CAPdiscrim(disMatrix ~ co2, 
                             data = subset(RingSumPFGMatrix,
                                           year == names(PFG_DisMatrix_Year)[i]),
                             permutations = 1000)
}

# create a plot
names(PFGcapList) <- c("2012", "2014", "2015")

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
PFGcapSppScore_co2 <- list()
for (i in 1:2){
  disMatrix <- PFG_DisMatrix_co2[[i]]
  envDF <- subset(RingSumPFGMatrix, co2 == names(PFG_DisMatrix_co2)[i])
  PFGcapList_co2[[i]] <- CAPdiscrim(disMatrix ~ year, data = envDF, permutations = 1000)
  PFGcapSppScore_co2[[i]] <- add.spec.scores(PFGcapList_co2[[i]], log(envDF[, PFGName] + 1))$cproj
  rm(envDF, disMatrix)
}

# create a plot
names(PFGcapList_co2) <- c("amb", "elev")
names(PFGcapSppScore_co2) <- c("amb", "elev")

PFGCapCo2df <- ldply(names(PFGcapList_co2), function(x) {
  lst <- PFGcapList_co2[[x]]
  data.frame(CAP1 = lst$x[, 1], CAP2 = lst$x[, 2], year = lst$group, co2 = x)
})

# species score
PFGCapCo2SppScore <- ldply(c("amb", "elev"), function(x) {
  df <- PFGcapSppScore_co2[[x]]
  spdf <- data.frame(CAP1 = df[, 1],
                     CAP2 = df[, 2],
                     PFG = row.names(df), 
                     co2 = x,  
                     year = "Scaled species score")
  spdf <- spdf[complete.cases(spdf), ]
  spdf
})  

p <- ggplot(PFGCapCo2df, aes(x = CAP1, y = CAP2, col = year))
p + geom_point(size = 4) + 
  geom_text(data = PFGCapCo2SppScore, aes(x = CAP1 * 2.5, y = CAP2, label = PFG)) +
  facet_grid(.~co2)

# save all objects to create summary document later
save.image(file = "output//Data/CAP_Object.RData")

# vegan
transDF <- log(RingSumPFGMatrix[, PFGName] + 1)
sites <- RingSumPFGMatrix[, !names(RingSumPFGMatrix) %in% PFGName]
vegDF <- capscale(transDF ~ yco, sites, dist = "bray")
plot(vegDF)
summary(vegDF)
