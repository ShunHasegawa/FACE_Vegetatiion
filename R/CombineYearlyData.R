################################
# Combine 2013, 2014 2015 Data #
################################

# 2013
# source("R/veg.12.process.R")
load("output/Data/FACE_Vegetation_2013.RData")

# 2014
# source("R/prcss.veg.2014.R")
load("output/Data/FACE_Veg2014.RData")
veg.14$month <- "December"

# 2015
# source("R/prcss.veg.2015.R")
load("output//Data/FAVE_vegetation2015.RData")
veg.2015$month <- "December"

vdf <- rbind.fill(df2013, veg.14, veg.2015)

####################
# organise dataset #
####################
# turn na into 0
vdf[is.na(vdf)] <- 0

# remove spp which were not observed
SiteVec <- c("year", "month", "ring", "plot", "position", "cell")
SppVec <- names(vdf)[!names(vdf) %in% SiteVec]
all(colSums(vdf[, SppVec]) > 0)

# sort columns
vdf <- vdf[, c(SiteVec, sort(SppVec))]

# rmeove unknown spp
vdf <- vdf[!grepl("unknown", names(vdf), ignore.case = TRUE)]

# organise spp----
SpName <- names(vdf)[!names(vdf) %in% SiteVec]

# Carex.breviformis -> Carex.breviculmis
vdf <- OrgSpp(vdf, siteVec = SiteVec,  
              KeepCol = "Carex.breviculmis", 
              CombineCol = SpName[grepl("carex", SpName, ignore.case = TRUE)])

# No too sure about Cyperus, but they're small abundance so just combine
colSums(vdf[, SpName[grepl("Cyperus", SpName, ignore.case = TRUE)]])
vdf <- OrgSpp(vdf, siteVec = SiteVec,
              KeepCol = "Cyperus.flaccidus", 
              CombineCol = SpName[grepl("Cyperus", SpName, ignore.case = TRUE)])

# Eragrostis brownii and benthamii: Therse spp are treated as the same spp in 
# Flora of the Sydney Region:
# http://ausgrass2.myspecies.info/content/eragrostis-brownii
vdf <- OrgSpp(vdf, siteVec = SiteVec,
              KeepCol = "Eragrostis.brownii", 
              CombineCol = SpName[grepl("Eragrostis.b", SpName, ignore.case = TRUE)])

# Fimbristylissp -> Fimbristylis.dichotoma 
vdf <- OrgSpp(vdf, siteVec = SiteVec,
              KeepCol = "Fimbristylis.dichotoma", 
              CombineCol = SpName[grepl("Fimbristylis", SpName, ignore.case = TRUE)])

# Galium.sp -> Galium.propinquum
vdf <- OrgSpp(vdf, siteVec = SiteVec,
              KeepCol = "Galium.propinquum", 
              CombineCol = SpName[grepl("Galium", SpName, ignore.case = TRUE)])

# Glycine is really hard to identify so just combine into Glycine.sp
colSums(vdf[SpName[grepl("Glycine", SpName, ignore.case = TRUE)]])
vdf <- OrgSpp(vdf, siteVec = SiteVec,
              KeepCol = "Glycine.sp", 
              CombineCol = SpName[grepl("Glycine", SpName, ignore.case = TRUE)])

# Oplisimenus.sp -> Oplismenus.aemulus
vdf <- OrgSpp(vdf, siteVec = SiteVec,
              KeepCol = "Oplismenus.aemulus", 
              CombineCol = SpName[grepl("Oplis", SpName, ignore.case = TRUE)])

# Paspalum -> Paspalum.dilatatum
vdf <- OrgSpp(vdf, siteVec = SiteVec,
              KeepCol = "Paspalum.dilatatum", 
              CombineCol = SpName[grepl("Paspalum", SpName, ignore.case = TRUE)])

# Phyllanthussp -> Phyllanthus.sp 
vdf <- OrgSpp(vdf, siteVec = SiteVec,  
              KeepCol = "Phyllanthus.sp", 
              CombineCol = SpName[grepl("Phyllanthus", SpName, ignore.case = TRUE)])

# Poranthera.microphylla -> Poranthera.microphylla
vdf <- OrgSpp(vdf, siteVec = SiteVec,  
              KeepCol = "Poranthera.microphylla", 
              CombineCol = SpName[grepl("Poranthera", SpName, ignore.case = TRUE)])

# Rubussp -> Rubus.parvifolius
vdf <- OrgSpp(vdf, siteVec = SiteVec,  
              KeepCol = "Rubus.parvifolius", 
              CombineCol = SpName[grepl("Rubus", SpName, ignore.case = TRUE)])

# Schoenus.opogon -> Schoenus.apogon
vdf <- OrgSpp(vdf, siteVec = SiteVec,  
              KeepCol = "Schoenus.apogon", 
              CombineCol = SpName[grepl("Schoenus", SpName, ignore.case = TRUE)])

# Solanum.sp. -> Solanum.nigrum
vdf <- OrgSpp(vdf, siteVec = SiteVec,  
              KeepCol = "Solanum.nigrum", 
              CombineCol = SpName[grepl("Solanum", SpName, ignore.case = TRUE)])
names(vdf)

# check all values <2
tsp <- names(vdf)[!names(vdf) %in% SiteVec]
all(vdf[, tsp] < 2)
# FALSE
which(vdf[, tsp] > 1, arr.ind = TRUE)
names(vdf[, tsp])[42]
names(vdf[, tsp])[34]

# Glycine sp. and Eragrostis.brownii; turn those into 1
vdf[, tsp][which(vdf[, tsp] > 1, arr.ind = TRUE)] <- 1
all(vdf[, tsp] < 2)
# TRUE

# Spp list----
vdf.mlt <- melt(vdf, id = c("year", "month", "ring", "plot", "position", "cell"))
spp <- data.frame(sp = sort(levels(vdf.mlt$variable)))
write.csv(spp, file = "output/Data/spp_2015.csv", row.names = FALSE)

# Spp which were found only in 2015
Spp <- names(vdf)[!names(vdf) %in% SiteVec]
YearSum <- ddply(vdf, .(year), function(x) colSums(x[Spp]))
newSp <- names(YearSum)[apply(rbind(YearSum[1:2, ] == 0, YearSum[3, ] != 0), 2, all)]
YearSum[newSp] # very small..

# remove Lichen and Carex.breviformis for the time being----
vdf$Lichen <- NULL

# month is factor
vdf$month <- factor(vdf$month)

###################
# Data correction # 
################### 

# For 2013 data, two surveys were conducted in September and December 2012. 
# September data is used for forbs and December for grass and sedge. In 
# December, forbs that had been recorded in September were checked again in
# December if there're missing or addition.

# The issues is that some forbs below seems to have been observed in December 
# but no cell position was recorded. As a solution, compare with adjacent 
# subplots and see if they're important. If not, then just ignore If they are, 
# allocate estimated number from the adjacent subplots. not the best solution
# but shouldn't change the final result too much cause those spp were not that
# abundant.

# 1.1.C, D Commelina.cyanea
  # "Observed in numerous cells.."
InspctPlot(ringval = 1, plotval = 1, sp = "Commelina.cyanea")
# not that abundant in the adjacent plots but it still says "numerous". So add 5
# for each of C and D
vdf[vdf$year == "2013" & 
      vdf$month == "December" &
      vdf$ring == 1 & 
      vdf$plot == 1 & 
      vdf$position == "C" & 
      vdf$Commelina.cyanea == 0, 
    "Commelina.cyanea"][1:5] <- 1

vdf[vdf$year == "2013" & 
      vdf$month == "December" &
      vdf$ring == 1 & 
      vdf$plot == 1 & 
      vdf$position == "D" & 
      vdf$Commelina.cyanea == 0, 
    "Commelina.cyanea"][1:5] <- 1

# 6.2.D. Parsonsia.straminea
InspctPlot(ringval = 6, plotval = 2, sp = "Parsonsia.straminea")
  # This one is tricky.. Ovserved in September in A, B and C, but not in D. Then
  # it was found in December in D but no notes for A, B and C. Probably, it was
  # found in A, B and C in Decebmer as well, but just was not noted.
  # So copy September to December and also add 2 in D in December

# 1) Copy September to December
vdf[vdf$year == "2013" & 
      vdf$month == "December" &
      vdf$ring == 6 & 
      vdf$plot == 2,
    "Parsonsia.straminea"] <- vdf[vdf$year == "2013" & 
                                    vdf$month == "September" &
                                    vdf$ring == 6 & 
                                    vdf$plot == 2,
                                  "Parsonsia.straminea"]
InspctPlot(ringval = 6, plotval = 2, sp = "Parsonsia.straminea")

# 2) Add 2 in D in December
vdf[vdf$year == "2013" &
      vdf$month == "December" &
      vdf$ring == 6 & 
      vdf$plot == 2 & 
      vdf$position == "D" & 
      vdf$Parsonsia.straminea == 0, 
    "Parsonsia.straminea"][1:2] <- 1

# 6.3.A. Parsonsia.straminea
InspctPlot(ringval = 6, plotval = 3, sp = "Parsonsia.straminea")
  # Same problem as above. Copy September to December and add 2 in A in
  # December.

# 1) Copy September to December
vdf[vdf$year == "2013" & 
      vdf$month == "December" &
      vdf$ring == 6 & 
      vdf$plot == 3,
    "Parsonsia.straminea"] <- vdf[vdf$year == "2013" & 
                                    vdf$month == "September" &
                                    vdf$ring == 6 & 
                                    vdf$plot == 3,
                                  "Parsonsia.straminea"]
InspctPlot(ringval = 6, plotval = 3, sp = "Parsonsia.straminea")

# 2) add 2 in A in December
vdf[vdf$year == "2013" & 
      vdf$month == "December" &
      vdf$ring == 6 & 
      vdf$plot == 3 & 
      vdf$position == "A" & 
      vdf$Parsonsia.straminea == 0, 
    "Parsonsia.straminea"][1:2] <- 1

# 6.3.D. Solanum.nigrum 
InspctPlot(ringval = 6, plotval = 3, sp = "Solanum.nigrum")
  # not observed in the adjacent subplots so just ignore

# 6.4.A. Parsonsia.straminea
InspctPlot(ringval = 6, plotval = 4, sp = "Parsonsia.straminea")
  # Same problem as above. Copy September to December and add 1 in A in December.

# 1) Copy September to December
vdf[vdf$year == "2013" & 
      vdf$month == "December" &
      vdf$ring == 6 & 
      vdf$plot == 4,
    "Parsonsia.straminea"] <- vdf[vdf$year == "2013" & 
                                    vdf$month == "September" &
                                    vdf$ring == 6 & 
                                    vdf$plot == 4,
                                  "Parsonsia.straminea"]
InspctPlot(ringval = 6, plotval = 4, sp = "Parsonsia.straminea")

# add 1 in A in December
vdf[vdf$year == "2013" & 
      vdf$month == "December" &
      vdf$ring == 6 & 
      vdf$plot == 4 & 
      vdf$position == "A" & 
      vdf$Parsonsia.straminea == 0, 
    "Parsonsia.straminea"][1] <- 1

# save----
save(vdf, file = "output//Data/FACE_Vegetation_Raw_2013_2015.RData")

##############################
# Create df including plant  #
# characteristis (e.g. PFGs) #
##############################
vdf.mlt <- melt(vdf, id = c("year", "month", "ring", "plot", "position", "cell"))

# plant properties
spList <- read.csv("Data//FACE_Vegetation_sp.list.csv", na.strings = c("NA", ""))
VegRes15 <- merge(vdf.mlt, spList, by.x = "variable", by.y = "sp", all.x = TRUE)

save(VegRes15, file = "output/Data/FACE_Vegetation_PFG_2015.RData")
