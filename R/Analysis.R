rm(list=ls(all=TRUE))

## ---- LoadData
source("R/Packages.R")
source("R/functions.R")
SiteName <- c("year", "block", "ring", "co2", "plot", "id", "position", "cell")

# Process Data ------------------------------------------------------------
# source("R/CombineYearlyData.R")


# load data ---------------------------------------------------------------

# Raw data for multi variate analysis (veg_matrix)
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S1.RData")
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S2.RData")
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S3.RData")
veg_matrix_list <- list(S1 = FullVdf, 
                        S2 = uniqueYear0_Vdf, 
                        S3 = uniqueYear0_plot_Vdf)
llply(veg_matrix_list, summary)

# dfs with plant functional groups (veg_df)
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S1_PFG.RData")
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S2_PFG.RData")
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S3_PFG.RData")
veg_df_list <- list(S1 = veg_FullVdf, 
                    S2 = veg_uniqueYear0_Vdf, 
                    S3 = veg_uniqueYear0_plot_Vdf)
llply(veg_df_list, summary)

# spp
SppName_list <- llply(veg_df_list, function(x) as.character(unique(x$variable)))
  
# > all species -----------------------------------------------------------


# grass and forb sp maes
gfspp <- veg %>% 
  filter(form %in% c("Grass", "Forb")) %>% 
  select(variable, form) %>%
  mutate(variable = as.character(variable)) %>% 
  distinct()

SppName_grass <- gfspp[gfspp$form == "Grass", 1]
SppName_forb <- gfspp[gfspp$form == "Forb", 1]

# plot sum
PlotSumVeg <- ddply(veg.face, .(year, ring, plot, block, co2, id), function(x) colSums(x[, SppName]))

# ring sum
RingSumVeg <- ddply(PlotSumVeg, .(year, ring, block, co2), function(x) colSums(x[, SppName]))


# > PFG -------------------------------------------------------------------

# plot
PlotSumPFGMatrix <- dcast(year + block + co2 + ring + plot ~ PFG, 
                          data = subset(veg, !is.na(PFG)), sum)
colSums(PlotSumPFGMatrix[,6:11])

# remove lichen, also add interaction term
PlotSumPFGMatrix <- within(PlotSumPFGMatrix, {
  id = ring:plot
  yco = year:co2
})

PFGName <- c("c3", "c4", "legume", "moss", "Non_legume", "wood")

# ring
RingSumPFGMatrix <- ddply(PlotSumPFGMatrix, .(year, block, ring, co2, yco), 
                          function(x) colSums(x[, PFGName]))
  

# diversity indices -------------------------------------------------------

# Diversity & eveness
siteDF <- PlotSumVeg[, !names(PlotSumVeg) %in% SppName]

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
SppSum <- ddply(veg, .(variable), summarise, value = sum(value))
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
