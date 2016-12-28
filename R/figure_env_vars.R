
# Light -------------------------------------------------------------------

# get daily mean
load("output/Data/FACE_canopy_transmittance_raw.RData")
light_dayly_co2 <- light_raw %>% 
  mutate(co2 = factor(ifelse(Ring %in% c("R1", "R4", "R5"), "eCO2", "Ambient"))) %>% 
  group_by(Date, co2) %>% 
  summarise_each(funs(mean), Gapfraction.mean)

# define sampling period
light_range <- data.frame(start_date = seq(as.Date("2012-10-26"), as.Date("2015-10-26"), "year"),
                          end_date   = seq(as.Date("2012-12-31"), as.Date("2015-12-31"), "year"),
                          ymin       = -Inf,
                          ymax       = Inf,
                          co2 = "Ambient")

# add survey dates
suevey_date <- data.frame(Date = as.Date(c("2012-09-15", "2012-12-15", "2014-1-15", "2015-1-30", "2016-2-15")))


fig_light <- ggplot() +
  geom_rect(data = light_range, 
            aes(xmin = start_date, 
                xmax = end_date, 
                ymin = -Inf, ymax = Inf), 
            fill = "gray70", col = NA) +
  geom_line(data = light_dayly_co2, aes(x = Date, y = Gapfraction.mean, col = co2)) +
  geom_point(data = suevey_date, aes(x = Date, y = .35), shape = 25, size = 3, fill = "black") +
  labs(y = "Canopy transmittance\n(Understorey light)", x = NULL) +
  science_theme +
  geom_vline(xintercept = as.numeric(as.Date("2012-09-18")), linetype = "dashed") +
  scale_color_manual(values = c("blue", "red"), labels = c("Ambient", expression(eCO[2])))

  


# Moisture & Temperature --------------------------------------------------

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



# . temperature -----------------------------------------------------------


fig_temp <- ggplot()+
  geom_rect(data = MT_range, 
            aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
            fill = "gray70", col = NA) +
  geom_line(dat = MT_dayly_co2, aes(x = Date, y = temp, col = co2, linetype = co2)) +
  scale_color_manual(values = c("blue", "red"), labels = c("Ambient", expression(eCO[2]))) +
  geom_vline(xintercept = as.numeric(as.Date("2012-09-18")), linetype = "dashed") +
  labs(x = NULL, y = expression(Soil~temperature~(degree*C))) +
  science_theme +
  theme(legend.position = "none")
fig_temp




# . moisture --------------------------------------------------------------


fig_moist <- ggplot()+
  geom_rect(data = MT_range, 
            aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
            fill = "gray70", col = NA)+
  geom_line(dat = MT_dayly_co2, aes(x = Date, y = moist, col = co2)) +
  labs(x = NULL, y = "Soil moisture (%)") +
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

env_figs[[1]] <- env_figs[[1]] +
  theme(legend.position = c(.11, .18),
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
