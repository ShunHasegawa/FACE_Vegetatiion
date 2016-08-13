#########################
# combine all year Data #
#########################

# load data ---------------------------------------------------------------

# vector for site
SiteVec <- c("year", "month", "ring", "plot", "position", "cell")

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

# 2016
# source("R/prcss.veg.2016.R")
load("output/Data/FACE_vegetation2016.RData")
veg.2016$month <- "February"

# remove Litter
veg.2016$Litter <- NULL
# turn all values >0 into 1
Sp2016 <- names(veg.2016)[!names(veg.2016) %in% SiteVec]
veg.2016[, Sp2016][which(veg.2016[, Sp2016] > 0, arr.ind = TRUE)] <- 1

vdf <- rbind.fill(df2013, veg.14, veg.2015, veg.2016)


# organise dataset ------------------------------------------------------

# turn na into 0
vdf[is.na(vdf)] <- 0

# remove spp which were not observed
SppVec <- names(vdf)[!names(vdf) %in% SiteVec]
all(colSums(vdf[, SppVec]) > 0)

# sort columns
vdf <- vdf[, c(SiteVec, sort(SppVec))]

# rmeove unknown spp
vdf <- vdf[!grepl("unknown", names(vdf), ignore.case = TRUE)]

# organise spp
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

# check all values <2
tsp <- names(vdf)[!names(vdf) %in% SiteVec]
all(vdf[, tsp] < 2)
# FALSE
which(vdf[, tsp] > 1, arr.ind = TRUE)
names(vdf[, tsp])[36]
names(vdf[, tsp])[44]

# Glycine sp. and Eragrostis.brownii; turn those into 1
vdf[, tsp][which(vdf[, tsp] > 1, arr.ind = TRUE)] <- 1
all(vdf[, tsp] < 2)
# TRUE

# Spp list
vdf.mlt <- melt(vdf, id = c("year", "month", "ring", "plot", "position", "cell"))
spp <- data.frame(sp = sort(levels(vdf.mlt$variable)))
write.csv(spp, file = "output/Data/spp_2016.csv", row.names = FALSE)

# Spp which were found only in 2015 and in 2016
Spp <- names(vdf)[!names(vdf) %in% SiteVec]
YearSum <- ddply(vdf, .(year), function(x) colSums(x[Spp]))
newSp_2015 <- names(YearSum)[apply(rbind(YearSum[1:2, ] == 0, YearSum[3, ] != 0), 2, all)]
newSp_2016 <- names(YearSum)[apply(rbind(YearSum[1:3, ] == 0, YearSum[4, ] != 0), 2, all)]
YearSum[newSp_2015] # very small..
YearSum[newSp_2016] 

# remove Lichen and Carex.breviformis for the time being
vdf$Lichen <- NULL

# month is factor
vdf$month <- factor(vdf$month)


# data correction -------------------------------------------------------

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


# . 1.1.C, D Commelina.cyanea ---------------------------------------------

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


# . 6.2.D. Parsonsia.straminea --------------------------------------------

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


# . 6.3.A. Parsonsia.straminea --------------------------------------------

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


# . 6.3.D. Solanum.nigrum  ------------------------------------------------

InspctPlot(ringval = 6, plotval = 3, sp = "Solanum.nigrum")
  # not observed in the adjacent subplots so just ignore

# . 6.4.A. Parsonsia.straminea  ------------------------------------------------
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


# . inspect Breynia -------------------------------------------------------

# Coordination of Breynia in 2016 was not very consistent with 2015, so inspect
# this
# Creatae a df for coordination
CorDF <- expand.grid(position = LETTERS[1:4], cell = 1:25)
CorDF$yv <- as.numeric(as.character(cut(CorDF$cell, 5, labels = 5:1)))
CorDF$xv <- CorDF$cell %% 5
CorDF$xv[CorDF$xv == 0] <- 5
CorDF$xv <- ifelse(CorDF$position %in% c("A", "C"), CorDF$xv - 5.5, CorDF$xv - 0.5)
CorDF$yv <- ifelse(CorDF$position %in% c("C", "D"), CorDF$yv - 5.5, CorDF$yv - 0.5)
plot(yv ~ xv, data = CorDF, type = "n")
with(CorDF, text(x = xv, y = yv, labels = paste0(position,cell)))
veg.2015$Breynia.oblongifolia
# merge with vdf
CorVeg <- merge(vdf, CorDF, by = c("position", "cell"))
head(CorVeg)

cmdf <- subset(CorVeg, ring == 3 & Breynia.oblongifolia > 0, 
               select =c ("year", "month", "ring", "plot", "position", "cell", "xv", "yv"))

cmdf$ym <- cmdf$year:cmdf$month
p <- ggplot(cmdf, aes(x = xv, y = yv, col = ym))
p2 <- p + 
  geom_point(position= position_dodge(0.7), cex = 5, alpha = 0.5) +
  facet_wrap( ~ plot) +
  scale_x_continuous(breaks = seq(-5, 5, 1), minor_breaks = NULL, limits = c(-5, 5)) +
  scale_y_continuous(breaks = seq(-5, 5, 1), minor_breaks = NULL, limits = c(-5, 5))
p2  
# it seems fine.. suggesting that the data sheet I created may be wrong..  

# save
save(vdf, file = "output//Data/FACE_Vegetation_Raw_2013_2016_PreInspection.RData")

# process Year0 -----------------------------------------------------------

###################################################
# Handling Year0 data from September and December #
###################################################

# Three options
  # 1) Combine Sep ad Dec in Year0 (FullVdf)
  # 2) Combine Sep ad Dec in Year0 and use only spp seen in subsequent years
  # for forbs (uniqueYear0_Vdf)
  # 3) Following 2), remove species from Year0 if they are not observed in the
  # same plots in the subsequent years (uniqueYear0_plot_Vdf)

vdf.mlt <- melt(vdf, id = c("year", "month", "ring", "plot", "position", "cell"))


# . Solution 1) Combine sep and dec in Year0 ------------------------------

  # Combine Sep and Dec and remove duplicates
  DF2013 <- subset(vdf.mlt, year == 2013)
  DF2013_merged <- DF2013 %>%
    group_by(variable, year, ring, plot, position, cell) %>%
    summarise(value = sum(value, na.rm = TRUE))
  
  # Remove double count
  DF2013_merged$value[which(DF2013_merged$value == 2)] <- 1
  DF2013_merged$month <- factor("Sep_Dec")
  
  # Create a new data frame containing Dec2013, S+D2013, 2014, 2015 and 2016
  VegRes16_SD <- rbind.fill(subset(vdf.mlt, month %in% c("December", "February")), 
                            DF2013_merged)
  
  # remove Dec2013
  FullSpdf <- VegRes16_SD %>%
    filter(!(year == "2013" & month == "December")) %>% 
    select(-month) %>% 
    mutate(year = factor(year, labels = paste0("Year", 0:3)),
           cell = factor(cell)) 
  
  # create matrix
  FullVdf <- dcast(year + ring + plot + position + cell ~ variable,  
                   value.var = "value", data = FullSpdf)
  summary(FullVdf)
  
  # organise
  
  # remove Euc seedlings as it's not reliable
  FullVdf$Euc.seedling <- NULL
  
  # check unknown spp
  all(!grepl("unknown", names(FullVdf), ignore.case = TRUE))
  
  # add co2, block and id, combine sedge and grass, wood and shrub
  FullVdf <- FullVdf %>% 
    mutate(cell  = as.character(cell),
           co2   = factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb")),
           block = recode(ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'"),
           id    = ring:plot) %>% 
    select(one_of(SiteName), everything())
  summary(FullVdf)
  
# . Solution 2) - remove unique spp in Year0 for forbs--------------------------
  
  # forb and grass species
  spList <- read.csv("Data//FACE_Vegetation_sp.list.csv", na.strings = c("NA", ""))
  spList <- filter(spList, sp %in% names(FullVdf))
  
  grass_spp <- spList %>% 
    filter(form %in% c("Grass", "Rush", "sedge")) %>% 
    select(sp)
  grass_spp <- as.character(grass_spp$sp)
  
  forb_spp <-  spList %>% 
    filter(!form %in% c("Grass", "Rush", "sedge")) %>% 
    select(sp) 
  forb_spp <- as.character(forb_spp$sp)
  
  # Use full-combined Sep-Dec, but without species uniquely found in Year0 for
  # forbs
  spp_sum_by_year <- FullVdf %>% 
    group_by(year) %>%
    summarise_each_(funs(sum), forb_spp)
  
  rm_spp <- names(which(apply(spp_sum_by_year[, -1], 2, 
                              function(x) all(x[2:4] == 0))))
  length(rm_spp)
  # 19 spp are to be removed
  
  # remove those spp
  uniqueYear0_Vdf <- select(FullVdf, -one_of(rm_spp))
  
  # check if it's correct
  new_forb_species <- setdiff(names(uniqueYear0_Vdf), c(SiteName, grass_spp))
  new_spp_sum_by_year <- ddply(uniqueYear0_Vdf, .(year), 
                               function(x) colSums(x[, new_forb_species]))
  all(!apply(new_spp_sum_by_year[, new_forb_species], 
             2, function(x) all(x[2:4] == 0)))

# . Solution 3) remove spp that are not observed in the same plots-------------
  
  # some species were observed in Year0 as well as subsequent years but not in the
  # same plot. Unless they were observed in the same plot, turn their values into 
  # 0 for Year0. 
  
  # e.g. species1 was observed in plot1 in Year0 and other plots in the following
  # years. But it was not observed in the following year in the same plot. Turn
  # the value for this species in that plot in Year0 into 0.
  
  uniqueYear0_plot_Vdf <- ddply(uniqueYear0_Vdf, .(ring, plot), 
                                function(x) correct_year0(x, spp = new_forb_species))

# save
save_Rdata_csv(FullVdf, 
               file = "output//Data/EucFACE_understorey_vegetation_2012-2106_S1")
save_Rdata_csv(uniqueYear0_Vdf, 
               file = "output//Data/EucFACE_understorey_vegetation_2012-2106_S2")
save_Rdata_csv(uniqueYear0_plot_Vdf, 
               file = "output//Data/EucFACE_understorey_vegetation_2012-2106_S3")

# combine veg dfs  with PFG etc -----------------------------------------------
# Create df including plant characteristis (e.g. PFGs)

FullVdf_list <- list(FullVdf, uniqueYear0_Vdf, uniqueYear0_plot_Vdf)
FullVdf_list_mlt <- llply(FullVdf_list, function(x) melt(x, id = SiteName))

# merge with plant properties
veg_list <- llply(FullVdf_list_mlt, function(x) {
  merge(x, spList, by.x = "variable", by.y = "sp", all.x = TRUE)
})
llply(veg_list, summary)

# combine sedge and grass, wood and shrub
veg_list <- llply(veg_list, function(x)
  mutate(x, form  = factor(ifelse(form %in% c("Tree", "Shrub"), "Wood",
                                  ifelse(form %in% c("Grass", "Sedge", "Rush"), 
                                         "Grass", as.character(form))))))

# remove c3_4 as it's really small number and hard to deal with c3_4..
# (Aristida.warburgii)
llply(veg_list, function(x) sum(x$value[x$PFG == "c3_4"]))
veg_list <- llply(veg_list, function(x) filter(x, PFG != "c3_4"))

# save
veg_FullVdf <- veg_list[[1]]
veg_uniqueYear0_Vdf <- veg_list[[2]]
veg_uniqueYear0_plot_Vdf <- veg_list[[3]]


l_ply(list())


save_Rdata_csv(veg_FullVdf, 
               filename = "output/Data/EucFACE_understorey_vegetation_2012-2106_S1_PFG")
save_Rdata_csv(veg_uniqueYear0_Vdf, 
               filename = "output/Data/EucFACE_understorey_vegetation_2012-2106_S2_PFG")
save_Rdata_csv(veg_uniqueYear0_plot_Vdf, 
               filename = "output/Data/EucFACE_understorey_vegetation_2012-2106_S3_PFG")

