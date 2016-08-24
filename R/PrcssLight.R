# download from HIEv
setToken(tokenfile = "Data/token.txt")

# download files from HIEv
light_raw <- downloadTOA5(filename = "ACE_P0037_RA_GAPFRACLAI_OPEN_L2.dat",
                          topath = "Data/hievdata/raw_data/")

# organise
light_dd <- light_raw %>% 
  select(Date, Ring, Gapfraction.mean) %>% 
  filter(month(Date) >= 10) %>% # use only 26th Oct to 31st Dec
  filter(!(month(Date) == 10 & day(Date) < 26)) %>% 
  mutate(ring = substr(Ring, 2, 2), 
         year = factor(year(Date), labels = paste0("Year", 0:3)))%>% 
  select(-Ring) %>% 
  group_by(year, ring) %>% 
  summarise(gapfraction = mean(Gapfraction.mean), N = n()) %>% 
  ungroup()
save(light_dd, file = "output/Data/FACE_canopy_transmittance.RData")

plot(gapfraction ~ as.numeric(year), col = ring, pch = 19, data = light_dd)
d_ply(light_dd, .(ring), function(x) lines(gapfraction ~ as.numeric(year), 
                                           col = unique(as.numeric(ring)), data = x))
