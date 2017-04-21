

# original data -----------------------------------------------------------


PlotSumVeg2 <- rename(PlotSumVeg, ring_plot = id, ring_year = RY)

site_data      <- PlotSumVeg2[ ,c("year", "ring", "co2", "ring_plot", "ring_year")]
graminoid_data <- select(PlotSumVeg2, year, co2, ring, plot, one_of(SppName_grass))
forb_data      <- select(PlotSumVeg2, year, co2, ring, plot, one_of(SppName_forb))
env_data <- c34growth_moist %>% 
  mutate(co2 = factor(ifelse(ring %in% c(1, 4, 5), "elev", "amb"))) %>% 
  rename(Temp = annual_temp2m, Moist = totalmoist) %>% 
  select(year, co2, ring, plot, Moist, Temp, PAR)



write.csv(graminoid_data, file = "output/Data/graminoid_data.csv", row.names = FALSE)
write.csv(forb_data, file = "output/Data/forb_data.csv", row.names = FALSE)
write.csv(env_data, file = "output/Data/env_data.csv", row.names = FALSE)
                              