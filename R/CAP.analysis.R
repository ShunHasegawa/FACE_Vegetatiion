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
  envDF <- subset(RingSumVeg, year == names(DisMatrix_Year)[i])
  capList[[i]] <- CAPdiscrim(disMatrix ~ co2, data = envDF, permutations = 1000)
  
  # add canonical correlation
  capList[[i]]$CanonicalCorSq <- CanonicalCor(CAPRes = capList[[i]], EnvDF = envDF, 
                                              term = "co2") 
  rm(disMatrix, envDF)
}

# create a plot
names(capList) <- c("2012", "2014", "2015")

CapYeardf <- ldply(names(capList), function(x) {
  lst <- capList[[x]]
  sqCor <- format(lst$CanonicalCorSq, digits = 3, nsmall = 3)
  data.frame(CAP1 = lst$x[, 1], co2 = lst$group, year = paste(x, sqCor, sep = "_"))
    # year: create labels for facet_grid showing year with associated canonical  
    # correlations. See notes below in more detail.
  })

p <- ggplot(CapYeardf, aes(x = co2, y = CAP1))
p2 <- p + geom_point(size = 4) + 
  facet_grid(.~year, 
             labeller = label_bquote(.(strsplit(x, "_")[[1]][1])~(sigma^2==.(strsplit(x, "_")[[1]][2])))) 
# I need to two values to create labels: year and canonical correlation, but
# label_bquote can only take one variable. Hence I concatenated those two above
# (e.g. 2012_0.597). Now I split this using strsplit and call each element
# separately within label_bquote (e.g.strsplit(x, "_")[[1]][1] = "2012")
ggsavePP(filename = "output//figs/FACE_CAPvsCO2_byYear", plot = p2,  width = 6, height = 4)

# each co2 treatment separately----

# create distance matrix
DisMatrix_co2 <- dlply(RingSumVeg, .(co2), 
                       function(x) vegdist(log(x[, SppName] + 1), method = "bray"))

# run CAP
capList_co2 <- list()
capSppScore_co2 <- list()
for (i in 1:2){
  disMatrix <- DisMatrix_co2[[i]]
  envDF <- subset(RingSumVeg, co2 == names(DisMatrix_co2)[i])
  capList_co2[[i]] <- CAPdiscrim(disMatrix ~ year, 
                        data = envDF,
                        permutations = 1000)
  # add canonical correlation
  capList_co2[[i]]$CanonicalCorSq <- CanonicalCor(CAPRes = capList_co2[[i]], EnvDF = envDF, 
                                                  term = "year") 
  # spp correlation
  capSppScore_co2[[i]] <- add.spec.scores(capList_co2[[i]], log(envDF[, SppName] + 1))$cproj
  rm(disMatrix, envDF)
}

# create a plot
names(capList_co2) <- c("amb", "elev")
names(capSppScore_co2) <- c("amb", "elev")

# site scores
CapCo2df <- ldply(c("amb", "elev"), function(x) {
  lst <- capList_co2[[x]]
  sqCor <- format(lst$CanonicalCorSq, digits = 1, nsmall = 3)
  data.frame(CAP1 = lst$x[, 1], CAP2 = lst$x[, 2], year = lst$group, 
             co2 = paste(x, sqCor[1], sqCor[2], sep = "_"))
})

# species score (correlation between ld1 and species abundance (after
# transformation))
CapCo2SppScore <- ldply(c("amb", "elev"), function(x) {
  df <- capSppScore_co2[[x]]
  spdf <- data.frame(CAP1 = df[, 1], CAP2 = df[, 2], 
                     Spp = row.names(df), 
                     co2 = x, 
                     year = "Scaled species Score")
  spdf <- spdf[complete.cases(spdf), ]
  # set threshhold for correlation
  subset(spdf, abs(CAP1) > .6|abs(CAP2) > .6)
})  
CapCo2SppScore$co2 <- factor(CapCo2SppScore$co2, labels = levels(CapCo2df$co2))

p <- ggplot(CapCo2df, aes(x = CAP1, y = CAP2, fill = year, col = year))
p2 <- p + geom_point(size = 2) + 
  geom_text(data = CapCo2SppScore, aes(x = CAP1 * 6, y = CAP2 * 2, label = Spp),
            size = 2) +
  facet_grid(.~co2,
             labeller = label_bquote(.(strsplit(x, "_")[[1]][1])~
                                       (list(sigma[1]^2==.(strsplit(x, "_")[[1]][2]),
                                             sigma[2]^2==.(strsplit(x, "_")[[1]][3])
                                             ) # list is used to insert ","
                                        )
                                     )
             ) +
  theme(legend.position = "top", 
        legend.direction = "horizontal")

ggsavePP(filename = "output//figs/FACE_CAPvsYear_byCO2", plot = p2,  width = 6, height = 6)

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
  envDF <- subset(RingSumPFGMatrix, year == names(PFG_DisMatrix_Year)[i])
  PFGcapList[[i]] <- CAPdiscrim(disMatrix ~ co2, 
                             data = envDF,
                             permutations = 1000)
  # add canonical correlation
  PFGcapList[[i]]$CanonicalCorSq <- CanonicalCor(CAPRes = PFGcapList[[i]],
                                                 EnvDF = envDF, term = "co2") 
  rm(disMatrix, envDF)
}

# create a plot
names(PFGcapList) <- c("2012", "2014", "2015")

PFGCapYeardf <- ldply(names(PFGcapList), function(x) {
  lst <- PFGcapList[[x]]
  sqCor <- format(lst$CanonicalCorSq, digits = 1, nsmall = 3)
  data.frame(CAP1 = lst$x[, 1], co2 = lst$group, year = paste(x, sqCor, sep = "_"))
})

p <- ggplot(PFGCapYeardf, aes(x = co2, y = CAP1))
p2 <- p + geom_point(size = 4) + 
  facet_grid(.~year, 
             labeller = label_bquote(.(strsplit(x, "_")[[1]][1])~
                                       (sigma^2==.(strsplit(x, "_")[[1]][2])))) 

ggsavePP(filename = "output//figs/FACE_CAPvsCO2_byYear_PFG", 
         plot = p2,  width = 6, height = 4)

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
  # add canonical correlation
  PFGcapList_co2[[i]]$CanonicalCorSq <- CanonicalCor(CAPRes = PFGcapList_co2[[i]],
                                                     EnvDF = envDF, 
                                                     term = "year") 
  # spp correlation
  PFGcapSppScore_co2[[i]] <- add.spec.scores(PFGcapList_co2[[i]], log(envDF[, PFGName] + 1))$cproj
  rm(envDF, disMatrix)
}

# create a plot
names(PFGcapList_co2) <- c("amb", "elev")
names(PFGcapSppScore_co2) <- c("amb", "elev")

PFGCapCo2df <- ldply(names(PFGcapList_co2), function(x) {
  lst <- PFGcapList_co2[[x]]
  sqCor <- format(lst$CanonicalCorSq, digits = 1, nsmall = 3)
  data.frame(CAP1 = lst$x[, 1], CAP2 = lst$x[, 2], year = lst$group, 
             co2 = paste(x, sqCor[1], sqCor[2], sep = "_"))
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
PFGCapCo2SppScore$co2 <- factor(PFGCapCo2SppScore$co2, labels = levels(PFGCapCo2df$co2))

p <- ggplot(PFGCapCo2df, aes(x = CAP1, y = CAP2, col = year))
p2 <- p + geom_point(size = 4) + 
  geom_text(data = PFGCapCo2SppScore, 
            aes(x = CAP1 * 2.5, y = CAP2, label = PFG), 
            size = 3) +
  facet_grid(.~co2,
             labeller = label_bquote(.(strsplit(x, "_")[[1]][1])~
                                       (list(sigma[1]^2==.(strsplit(x, "_")[[1]][2]),
                                             sigma[2]^2==.(strsplit(x, "_")[[1]][3])
                                       ) # list is used to insert ","
                                       )
             )
  ) +
  theme(legend.position = "top", 
        legend.direction = "horizontal")

ggsavePP(filename = "output//figs/FACE_CAPvsYear_byCO2_PFG", plot = p2,  width = 6, height = 6)

# PERMANOVA----
ambDF <- subsetD(RingSumPFGMatrix, co2 == "amb")
ambPFG <- ambDF[, PFGName]
ambSite <- ambDF[, !names(ambDF) %in% PFGName]
disMatrix <- vegdist(log(ambPFG + 1), method = "bray")

a2 <- adonis(disMatrix ~ year, data = ambSite, permutations = 9999)
a2

tdf <- log(ambPFG+1)

n1 <- nested.npmanova(tdf ~ year + ring, data = ambSite, perm=10, method = "bray")


?nested.npmanova

disMatrix <- vegdist(log(RingSumPFGMatrix[, PFGName] + 1), method = "bray")
Sitetdf <- RingSumPFGMatrix[, !names(RingSumPFGMatrix) %in% PFGName]
a2 <- adonis(disMatrix ~ year + co2, data = Sitetdf, strata = Sitetdf$ring, 
             permutations = 9999)
a2



# save all objects to create summary document later
save.image(file = "output//Data/CAP_Object.RData")






# vegan
transDF <- log(RingSumPFGMatrix[, PFGName] + 1)
sites <- RingSumPFGMatrix[, !names(RingSumPFGMatrix) %in% PFGName]
vegDF <- capscale(transDF ~ yco, sites, dist = "bray")
plot(vegDF)
summary(vegDF)
