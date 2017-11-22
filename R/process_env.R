# download from HIEv ------------------------------------------------------


# # setToken(tokenfile = "Data/token.txt")
# # 
# # airvar_raw <- downloadTOA5(filename = "FACE_.*_AirVars_.*dat",
# #                           topath    = "Data/hievdata/raw_data/",
# #                           maxnfiles = 999)
# # 
# # save(airvar_raw, file = "output/Data/FACE_airvar_raw.RData")
# load("output/Data/FACE_airvar_raw.RData")
# 
# 
# 
# 
# # process HIEv data -------------------------------------------------------
# 
# 
# # get daily mean, min and max temparater at two layers (at a height of 2m and 15m)
# boxplot(airvar_raw$AirTC_1_Avg)
# airvar_day <- airvar_raw %>%
#   filter(Date        > as.Date("2012-09-18"),         # when co2 treatment was commenced
#          Date        < as.Date("2016-2-15"),
#          AirTC_1_Avg > -10) %>%                       # remove obviously weird value
#   distinct() %>%                                      # remove duplicates
#   mutate(ring = factor(substring(Source, 7, 7)),      # add ring number
#          year = factor(year(Date))) %>%
#   group_by(year, Date, ring) %>% 
#   summarise_each(funs(Mean = mean(., na.rm = TRUE), 
#                       Min  = min(., na.rm = TRUE),
#                       Max  = max(., na.rm = TRUE),
#                       N    = sum(!is.na(.))), 
#                  airtemp2m = AirTC_1_Avg, airtemp15m = AirTC_2_Avg) %>% 
#   group_by(year, ring) %>% 
#   mutate(annual_temp2m = mean(airtemp2m_Mean),
#          T = max(27, 1.745 * annual_temp2m  + 11.143)[1]) %>% 
#   ungroup() %>% 
#   mutate(c3growth = airtemp2m_Min >= -1 & airtemp2m_Max >= 10 & airtemp2m_Max < 24,  # c3 growing dates
#          c4growth = airtemp2m_Max >= 21 & airtemp2m_Max < T & month(Date))           # c4 growtin dates
# save(airvar_day, file = "output/Data/FACE_air_variables.RData")

load("output/Data/FACE_air_variables.RData")  # airvar_day; run the above to obtain up-to-date data from HIEv

c34growthdate <- airvar_day %>% 
  select(year, Date, ring, c3growth, c4growth, annual_temp2m)


## understorey temperature
annualtemp <- c34growthdate %>% 
  select(year, ring, annual_temp2m) %>% 
  distinct() %>% 
  filter(!year %in% c("2012", "2016")) %>% 
  mutate(year = mapvalues(year, c(2013, 2014, 2015), paste0("Year", 1:3)))


## understorey PAR
load("output/Data/FACE_FloorPAR.RData")

Lightdf$PAR <- rowMeans(Lightdf[, c("PAR_Den_1_Avg", "PAR_Den_2_Avg", "PAR_Den_3_Avg")], na.rm = TRUE)
underpar <- Lightdf %>% 
  mutate(year = year(DateTime)) %>% 
  filter(year %in% c(2013:2015)) %>% 
  group_by(year, ring) %>% 
  summarise(PAR = mean(PAR, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(year = factor(year, labels = paste0("Year", 1:3)))

## Murphy 2007. For C 3 grasses, growth was considered possible in any month 
## where the daily minimum temperature was ≥ −1 oC, and the daily maximum 
## temperature was ≥ 10 oC and < 24 oC. For C 4 grasses, growth was considered 
## possible in months where the daily maximum temperature was ≥ 21 °C and less 
## than the tem- perature define (Murphy and Bowman 2007).




# assign moisture data to each vegetation plot based on their corrdinates --------


# soil moisture for each vegetation plot assined based probe coordinates
load("Data/FACE_TDR_ProbeDF.RData")
veg_moist <- FACE_TDR_ProbeDF %>%
  filter(Date        > as.Date("2012-09-18"),         # when co2 treatment was commenced
         Date        < as.Date("2016-2-15"),
         Sample == "vegetation") %>% 
  mutate(year = factor(year(Date)), plot = factor(plot)) %>% 
  select(year, Date, ring, plot, Moist) 


# merge
summary(c34growthdate)
summary(veg_moist)
summary(Lightdf)


## survey dates: 2012-12-15 (Year0), 2014-1-15 (Year1), 2015-1-30 (Year2), 2016 (Year3)

c34growth_moist <- veg_moist %>% 
  filter(year %in% 2013:2015) %>% 
  mutate(year = factor(year, labels = paste0("Year", 1:3))) %>% 
  group_by(year, ring) %>%
  summarise(totalmoist = mean(Moist)) %>%
  ungroup() %>%
  left_join(annualtemp) %>% 
  left_join(underpar)

## save for the manuscript
env_data <- c34growth_moist %>% 
  rename(Moist = totalmoist,
         Temp = annual_temp2m)
write.csv(env_data, "output/table/manuscript_data/env_data.csv",
          row.names = FALSE)