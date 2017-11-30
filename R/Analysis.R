rm(list=ls(all=TRUE))

source("R/Packages.R")
source("R/functions.R")
options(na.action = "na.fail")  # change na.action setting for MuMIN::dredge
SiteName <- c("year", "block", "ring", "co2", "plot", "id", "position", "cell", "RY")




# Process Data ------------------------------------------------------------
# source("R/CombineYearlyData.R")




# load data ---------------------------------------------------------------


# Raw data for multi variate analysis (matrix)
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S1.RData")

# fix species names
FullVdf <- FullVdf %>% 
  dplyr::rename(Ambrosia.artemisiifolia  = Ambrosia.sp,
                Arthropodium.minus       = Arthropodium.sp,
                Digitaria.longiflora     = Digitaria.sp,
                Drosera.auriculata       = Drosera.sp,
                Leontodon.saxatilis      = Leontodon.taraxacoides,
                Phyllanthus.gunnii       = Phyllanthus.sp,
                Sisyrinchium.iridifolium = Sisyrinchium.sp)
summary(FullVdf)

# dfs with plant functional groups (df)
load("output//Data/EucFACE_understorey_vegetation_2012-2106_S1_PFG.RData")

# fix species names
veg_FullVdf <-  veg_FullVdf %>% 
  mutate(variable   = dplyr::recode(variable,                                         
                                  "Ambrosia.sp"            = "Ambrosia.artemisiifolia",
                                  "Arthropodium.sp"        = "Arthropodium.minus",
                                  "Digitaria.sp"           = "Digitaria.longiflora",
                                  "Drosera.sp"             = "Drosera.auriculata",
                                  "Leontodon.taraxacoides" = "Leontodon.saxatilis",
                                  "Phyllanthus.sp"         = "Phyllanthus.gunnii",
                                  "Sisyrinchium.sp"        = "Sisyrinchium.iridifolium"))
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
SppName_forb  <- gfspp[gfspp$form == "Forb", 1]

# plot sum
PlotSumVeg <- ddply(FullVdf, .(year, ring, plot, block, co2, id, RY), 
                    function(x) colSums(x[, SppName])) %>% 
  arrange(ring, plot, year)


# save raw data for EucFACE vegetaiotn manuscript
forb_raw_data <- PlotSumVeg %>% 
  select(year, co2, ring, plot, one_of(SppName_forb))
gram_raw_data <- PlotSumVeg %>% 
  select(year, co2, ring, plot, one_of(SppName_grass))
gram_pfg <- veg_FullVdf %>% 
  filter(form == "Grass") %>% 
  select(variable, PFG) %>% 
  distinct()

write.csv(forb_raw_data, "output/Data/forb_data.csv", row.names = FALSE)
write.csv(gram_raw_data, "output/Data/graminoid_data.csv", row.names = FALSE)
write.csv(gram_pfg, "output/Data/graminoid_pfg.csv", row.names = FALSE)




# ring sum
RingSumVeg <- ddply(PlotSumVeg, .(year, ring, block, co2), 
                    function(x) colSums(x[, SppName])) %>% 
  arrange(ring, year)




# > PFG -------------------------------------------------------------------


# plot
PlotSumPFGMatrix <- dcast(year + block + co2 + ring + plot ~ PFG, 
                          data = subset(veg_FullVdf, !is.na(PFG)), sum) %>% 
  mutate(id = ring:plot, yco = year:co2)

PlotSumPFGMatrix %>% 
  select(-one_of(c(SiteName, "id", "yco"))) %>% 
  summarise_each(funs(sum))


PFGName <- setdiff(names(PlotSumPFGMatrix), c(SiteName, "id", "yco")) 
  

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
  siteDF %>% 
    mutate( 
      H = diversity(x),  # Shannon's index
      S = specnumber(x), # number of spp
      J = H/log(S)       # Pielou's evenness
    ) %>% 
    group_by(year, ring, co2) %>% 
    summarise_each(funs(mean(., na.rm = TRUE)), H, S, J) %>% 
    ungroup()
  })

# Identify dominant spp
SppSum <- ddply(veg_FullVdf, .(variable), summarise, value = sum(value))
SppSum <- SppSum[order(SppSum$value, decreasing = TRUE),]
SppSum <- within(SppSum, {
  Cov    <- round(value * 100 / sum(value), 3)
  CumSum <- cumsum(value)
  Dominant <- Cov >= 5 # species with >5 % coverage
})
DmSpp <- droplevels(SppSum$variable[SppSum$Dominant])
DmSpp
sum(SppSum$value[SppSum$Dominant])/sum(SppSum$value)




# species yearly summary --------------------------------------------------


# yearly summary by form
spp_year_summary <- veg_FullVdf %>% 
  group_by(year, variable, form) %>% 
  summarise(value = sum(value)) %>% 
  group_by(year, form) %>% 
  summarise(S = sum(value > 0)) %>% 
  ungroup() %>% 
  spread(key = form, value = S) %>% 
  mutate(Total = rowSums(.[, -1]))


# yearly summary by co2
spp_co2_summary <- veg_FullVdf %>% 
  group_by(year, variable, co2) %>% 
  summarise(value = sum(value)) %>% 
  group_by(year, variable) %>% 
  summarise(uni_amb   = value[co2 == "amb"] > 0  & value[co2 == "elev"] == 0,     # uniqe spp in amb
            uni_elev  = value[co2 == "amb"] == 0 & value[co2 == "elev"] > 0,      # uniqe spp in elev
            common    = value[co2 == "amb"] > 0  & value[co2 == "elev"] > 0) %>%  # common spp
  group_by(year) %>% 
  summarise_each(funs(sum), -variable) %>%                                        # number of TRUE
  mutate(tot_amb  = uni_amb + common,
         tot_elev = uni_elev + common,
         total    = uni_amb + uni_elev + common)


# total % coverage by spp
spp_cov_summary <- veg_FullVdf %>% 
  group_by(variable) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  mutate(cov = value * 100 / sum(value)) %>% 
  arrange(cov)


# total % coverage by PFG
pfg_cov_summary <- veg_FullVdf %>% 
  mutate(form    = as.character(form),
         PFG     = as.character(PFG),
         pfgform = ifelse(form == "Grass", paste0(PFG, form),
                          ifelse(form == "Forb", PFG, form))) %>% 
  group_by(pfgform) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  mutate(cov = round(value * 100 / sum(value), 1)) %>% 
  arrange(cov)
pfg_cov_summary


# process IEM (Soil N, P data)
iem <- read.xlsx2(file = "Data/newIem_smmry_tbl.xlsx", sheetIndex = 1, stringsAsFactors = FALSE)
iem <- iem %>% 
  select(-Moist, -Temp_Mean) %>% 
  mutate_each(funs(as.Date), insertion, sampling, date) %>% 
  mutate_each(funs(as.numeric), nitrate, ammonium, phosphate) %>%
  mutate_each(funs(as.factor), ring, plot) %>% 
  filter(time %in% c(6, 7, 13, 14, 16, 19)) %>% 
  mutate(year = ifelse(time %in% c(6, 7), "Year0",
                       ifelse(time %in% c(13, 14), "Year1", 
                              ifelse(time == 16, "Year2", "Year3"))))


iem %>% 
  select(insertion, sampling, time, year) %>% 
  gather(variable, value, insertion, sampling) %>% 
  mutate(day2 = yday(value)) %>% 
  distinct() %>% 
  arrange(as.numeric(as.character(time))) %>%
  mutate(day2 = ifelse(day2 < 100, 365 + day2, day2)) %>% 
  ggplot(., aes(x = day2, y = 1, group = time, col = time))+
  geom_point()+
  geom_path()+
  facet_grid(year ~ .)


iem_raw <- iem %>% 
  rename(no = nitrate, nh = ammonium, p = phosphate) %>% 
  group_by(year, time, ring) %>% 
  summarise_each(funs(mean), no, nh, p) %>% 
  group_by(year, ring) %>% 
  summarise_each(funs(mean), no, nh, p) %>% 
  ungroup() %>% 
  mutate(nitr = no + nh,
         np   = nitr / p)

# soil data for teh manuscript
soil_data <- iem_raw %>% 
  select(year, ring, no, nh, p) %>% 
  rename(nitrate = no, ammonium = nh, phosphate = p) %>% 
  filter(year != "Year0")
write.csv(soil_data, "output/table/manuscript_data/soil_data.csv",
          row.names = FALSE)


plot(nitr ~ as.numeric(factor(year)), data = iem_raw, col = ring, pch = 19)
plot(p ~ as.numeric(factor(year)), data = iem_raw, col = ring, pch = 19)




# analysis and fig --------------------------------------------------------

source("R/fig_diversity_indices.R")                   # diversity indices
source("R/fig_dominent_spp.R")                        # dominent spp
source("R/Stats_C3_4_Grass.R")                        # C3vsC4, LegvsNonLeg, NativevsIntro
source("R/PFG_Proportaion.R")                         # Grass vs Forb
source("R/analysis_c43ratio.R")                       # C4:C3 ratios
source("R/fig_pfg_propotion.R")                       # combine grassn and PFG proportino results to create a figure
# source("R/analysis_PRC_MDS.R")                        # prncipal response curve analysis and MDS analysis
source("R/analysis_deltaC34.R")                       # process HIEv data for temp and moist, and perform multivple regression for delta C3 and C4 with temp, moist and PAR
source("R/analysis_dominant_subordinate.R")           # define subordinate and doiminant spp; analyse S:D ratios
source("R/analysis_LAR_dominant_subordinate_env.R")   # LAR of subordinate and dominant C3 and C4 species against environmental variables
source("R/analysis_LAR_dominant_subordinate_soil.R")  # LAR of subordinate and dominant C3 and C4 species against soil nutrients
source("R/analysis_summary.R")                        # summarise above analyses for repeated measured anova


# Luke’s biomass harvest --------------------------------------------------

# lh <- read.csv("Data/luke_harvest.csv")
# lh %>% 
#   group_by(ring, form) %>% 
#   summarise(value = sum(Mass..g.)) %>% 
#   spread(form, value) %>% 
#   mutate(total = sum(f, g, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   mutate(f_p = sum(f, na.rm = TRUE) / sum(total),
#          g_p = sum(g)/ sum(total))
# 
# 
# 
# 
# # figs --------------------------------------------------------------------
# source("R//Figs.R")
# 
# 
# 
# 
# # stats -------------------------------------------------------------------
# source("R/Stats.R")

# save all objects. This will be used when creating a summary document all
# objects. This will be used when creating a summary document
sink("output/session_info.txt", append = TRUE)
paste0("/n/n/n",now(), "/n")
sessionInfo()
sink()

save.image(file = "output//Data/AllObj.RData")
