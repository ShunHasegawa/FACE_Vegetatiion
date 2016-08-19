rm(list=ls(all=TRUE))

## ---- LoadData
source("R/Packages.R")
source("R/functions.R")
SiteName <- c("year", "block", "ring", "co2", "plot", "id", "position", "cell")

# Process Data ------------------------------------------------------------
# source("R/CombineYearlyData.R")

# load data ---------------------------------------------------------------

# Raw data for multi variate analysis (matrix)
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S1.RData")
summary(FullVdf)

# dfs with plant functional groups (df)
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S1_PFG.RData")
summary(veg_FullVdf)

# spp
SppName <- as.character(unique(veg_FullVdf$variable))
  
# organise dfs ------------------------------------------------------------

# > all species -----------------------------------------------------------

# grass and forb spp
gfspp <- veg_FullVdf %>% 
  filter(form %in% c("Grass", "Forb")) %>% 
  select(variable, form) %>%
  mutate(variable = as.character(variable)) %>% 
  distinct()

SppName_grass <- gfspp[gfspp$form == "Grass", 1]
SppName_forb <- gfspp[gfspp$form == "Forb", 1]

# plot sum
PlotSumVeg <- ddply(FullVdf, .(year, ring, plot, block, co2, id), 
                    function(x) colSums(x[, SppName]))

# ring sum
RingSumVeg <- ddply(PlotSumVeg, .(year, ring, block, co2), 
                    function(x) colSums(x[, SppName]))


# > PFG -------------------------------------------------------------------

# plot
PlotSumPFGMatrix <- dcast(year + block + co2 + ring + plot ~ PFG, 
                          data = subset(veg_FullVdf, !is.na(PFG)), sum)
PlotSumPFGMatrix %>% 
  select(-one_of(SiteName)) %>% 
  summarise_each(funs(sum))

# add interaction term
PlotSumPFGMatrix <- mutate(PlotSumPFGMatrix, id = ring:plot, yco = year:co2)

PFGName <- c("c3", "c4", "fern", "legume", "moss", "Non_legume", "wood")

# ring
RingSumPFGMatrix <- ddply(PlotSumPFGMatrix, .(year, block, ring, co2, yco), 
                          function(x) colSums(x[, PFGName]))
  
# diversity indices -------------------------------------------------------

# Diversity & eveness
siteDF <- select(PlotSumVeg, -one_of(SppName))

vegDF_list <- llply(list(all_spp = SppName, grass_spp = SppName_grass, 
                         forb_spp = SppName_forb), 
                    function(x) PlotSumVeg[, x])

vegDF <- vegDF_list[["all_spp"]]

# compute diversity indices
DivDF_list <- llply(vegDF_list, function(x) {
  mutate(siteDF, 
         H = diversity(x),  # Shannon's index
         S = specnumber(x), # number of spp
         J = H/log(S)       # Pielou's evenness
  )})

DivDF       <- DivDF_list[["all_spp"]]
DivDF_grass <- DivDF_list[["grass_spp"]]
DivDF_forb  <- DivDF_list[["forb_spp"]]

# Identify dominant spp
SppSum <- ddply(veg_FullVdf, .(variable), summarise, value = sum(value))
SppSum <- SppSum[order(SppSum$value, decreasing = TRUE),]
SppSum <- within(SppSum, {
  Cov <- round(value * 100/sum(value), 3)
  CumSum <- cumsum(value)
  Dominant <- Cov >= 5 # species with >5 % coverage
})
DmSpp <- droplevels(SppSum$variable[SppSum$Dominant])
DmSpp
sum(SppSum$value[SppSum$Dominant])/sum(SppSum$value)


# figs --------------------------------------------------------------------
source("R//Figs.R")


# stats -------------------------------------------------------------------
source("R/Stats.R")

# save all objects. This will be used when creating a summary document all
# objects. This will be used when creating a summary document
sink("output/session_info.txt", append = TRUE)
paste0("/n/n/n",now(), "/n")
sessionInfo()
sink()

save.image(file = "output//Data/AllObj.RData")
