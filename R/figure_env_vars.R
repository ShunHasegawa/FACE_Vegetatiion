
# Light -------------------------------------------------------------------

# get daily mean
load("output/Data/FACE_FloorPAR.RData")
Lightdf$PAR <- rowMeans(Lightdf[, c("PAR_Den_1_Avg", "PAR_Den_2_Avg", "PAR_Den_3_Avg")], na.rm = TRUE)

light_dayly_co2 <- Lightdf %>% 
  filter(Date >= as.Date("2012-10-15")) %>% 
  mutate(co2 = factor(ifelse(ring %in% c("1", "4", "5"), "eCO2", "Ambient"))) %>% 
  group_by(Date, co2) %>% 
  summarise(PAR = mean(PAR, na.rm = TRUE))


# add survey dates
suevey_date <- data.frame(Date = as.Date(c("2012-09-15", "2012-12-15", "2014-1-15", "2015-1-30", "2016-2-15")))


fig_light <- ggplot() +
  geom_point(data = light_dayly_co2, aes(x = Date, y = PAR), col = "grey") +
  geom_smooth(data = light_dayly_co2, aes(x = Date, y = PAR, col = co2),
              method = "loess", span = .02, n = 2648, se = FALSE, size = .6)+
  geom_point(data = suevey_date, aes(x = Date, y = 270), shape = 25, size = 3, fill = "black") +
  labs(y = expression(PAR~(mu*mol~s^'-1'~m^'-2')), x = NULL) +
  science_theme +
  theme(legend.position = "none")+
  geom_vline(xintercept = as.numeric(as.Date("2012-09-18")), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red"), labels = c("Ambient", expression(eCO[2])))

  


# temperature -----------------------------------------------------------

load("output/Data/FACE_air_variables.RData")  # airvar_day; run the above to obtain up-to-date data from HIEv

head(airvar_day)
temp_daily_co2 <- airvar_day %>%
  mutate(co2 = factor(ifelse(ring %in% c("1", "4", "5"), "eCO2", "Ambient"))) %>% 
  group_by(Date, co2) %>% 
  summarise(temp = mean(airtemp2m_Mean, na.rm = TRUE))
  


fig_temp <- ggplot()+
  geom_point(data = temp_daily_co2, aes(x = Date, y = temp), col = "grey") +
  geom_smooth(data = temp_daily_co2, aes(x = Date, y = temp, col = co2, linetype = co2),
              method = "loess", span = .02, n = 2648, se = FALSE, size = .6)+
  scale_color_manual(values = c("blue", "red"), labels = c("Ambient", expression(eCO[2]))) +
  geom_vline(xintercept = as.numeric(as.Date("2012-09-18")), linetype = "dashed") +
  labs(x = NULL, y = expression(Temperature~(degree*C))) +
  science_theme +
  theme(legend.position = "none")
fig_temp




# moisture --------------------------------------------------------------

load("Data//SoilVariables/soil.var_ring.means.RData")


# sampling period
MT_range <- data.frame(start_date = seq(as.Date("2012-8-1"), as.Date("2015-8-1"), "year"),
                       end_date   = seq(as.Date("2012-12-31"), as.Date("2015-12-31"), "year"),
                       ymin       = -Inf,
                       ymax       = Inf,
                       co2 = "Ambient")


# get daily mean
MT_dayly_co2 <- ring.means %>% 
  group_by(Date, co2) %>% 
  summarise_each(funs(mean), moist, temp) %>% 
  mutate(moist = moist * 100)



fig_moist <- ggplot()+
  geom_line(dat = MT_dayly_co2, aes(x = Date, y = moist, col = co2)) +
  labs(x = NULL, y = "Moisture (%)") +
  geom_vline(xintercept = as.numeric(as.Date("2012-09-18")), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red"), labels = c("Ambient", expression(eCO[2]))) +
  science_theme +
  theme(legend.position = "none")

fig_moist




# precipitation -----------------------------------------------------------

load("output/Data/allrain.RData")

fig_rain <- ggplot() +
  geom_bar(data = allrain, aes(x = Date, y = Rain_mm_Tot), stat = "identity") +
  geom_vline(xintercept = as.numeric(as.Date("2012-09-18")), linetype = "dashed") +
  science_theme +
  labs(x = NULL, y = "Precipitation (mm)")
fig_rain



  

# merge figures -----------------------------------------------------------

env_figs <- list(light = fig_light,
                 temp  = fig_temp,
                 moist = fig_moist,
                 rain  = fig_rain)
env_figs <- llply(env_figs, function(x) x + scale_x_date(limits = c(as.Date("2012-6-15"), as.Date("2016-2-20"))))

for (i in 1:3){
  env_figs[[i]] <- env_figs[[i]] + theme(axis.text.x = element_blank())
}

env_figs[[3]] <- env_figs[[3]] +
  theme(legend.position = c(.11, 1.1),
        legend.background = element_rect(fill = "white"))

env_figs[[4]] <- env_figs[[4]] +
  scale_x_date(breaks= date_breaks("3 month"), 
               labels = date_format("%b-%y"),
               limits = c(as.Date("2012-6-15"), as.Date("2016-2-20"))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# merge plots
env_fig_merged <- rbind( 
  ggplotGrob(env_figs[[1]]),
  ggplotGrob(env_figs[[2]]),
  ggplotGrob(env_figs[[3]]),
  ggplotGrob(env_figs[[4]]))

grid.newpage()
grid.draw(env_fig_merged)

# save
ggsavePP(filename = "output/figs/daily_env_var", plot = env_fig_merged,
         width = 6.5, height = 7)
