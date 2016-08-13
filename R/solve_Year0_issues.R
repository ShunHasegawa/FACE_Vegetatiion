# load data ---------------------------------------------------------------

# Raw data for multi variate analysis (matrix)
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S1.RData")
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S2.RData")
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S3.RData")
veg_matrix_list <- list(s1 = FullVdf, s2 = uniqueYear0_Vdf,
                        s3 = uniqueYear0_plot_Vdf)
llply(veg_matrix_list, summary)

# dfs with plant functional groups (df)
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S1_PFG.RData")
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S2_PFG.RData")
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S3_PFG.RData")
veg_df_list <- list(s1 = veg_FullVdf, s2 = veg_uniqueYear0_Vdf,
                    s3 = veg_uniqueYear0_plot_Vdf)
llply(veg_df_list, summary)

# spp
SppName_list <- llply(veg_df_list, function(x) unique(as.character((x$variable))))


# organise for further analysis -------------------------------------------

# > grass and forb spp ----------------------------------------------------

gfspp_list <- llply(veg_df_list, function(x){
  x %>% 
    filter(form %in% c("Grass", "Forb")) %>% 
    select(variable, form) %>%
    mutate(variable = as.character(variable)) %>% 
    distinct()
  })
llply(gfspp_df_list, some)

SppName_grass_list <- llply(gfspp_list, function(x) x$variable[x$form == "Grass"])
SppName_forb_list <- llply(gfspp_list, function(x) x$variable[x$form == "Forb"])


# > plot sum --------------------------------------------------------------

PlotSumVeg_list <- llply(1:3, function(x) {
  s <- SppName_list[[x]]
  d <- veg_matrix_list[[x]]
  d %>% 
    group_by(year, block, ring, plot, id, co2) %>% 
    summarise_each_(funs(sum), s) %>% 
    ungroup()
})
  
# > diversity indices -------------------------------------------------------

vegDF_list <- llply(1:3, function(x){
  d <- PlotSumVeg_list[[x]]
  d_list <- list(all_spp   = d[, SppName_list[[x]]], # all spp
                 grass_spp = d[, SppName_grass_list[[x]]], # grass spp
                 forb_spp  = d[, SppName_forb_list[[x]]]) # forb spp
  return(d_list)
})
names(vegDF_list) <- paste0("s", 1:3)
llply(vegDF_list, summary)

# compute diversity indices
siteDF <- select(PlotSumVeg_list[[1]], -one_of(SppName_list[[1]]))

DivDF_list <- llply(vegDF_list, function(x) {
  llply(x, function(y) mutate(siteDF, 
                              H = diversity(y),  # Shannon's index
                              S = specnumber(y), # number of spp
                              J = H/log(S)       # Pielou's evenness
                              ))})
llply(DivDF_list, summary)
