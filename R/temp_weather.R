site_rain <- read.csv("Data/station_weather/IDCJAC0001_067021_Data12.csv")
site_rain[, c(-1:-2)] <- apply(site_rain[, c(-1:-2)], 2, as.numeric)
site_rain <- site_rain[complete.cases(site_rain), ]
site_rain <- filter(site_rain, Year >= 1929)
plot(Annual ~ Year, data = site_rain)
summary(site_rain)
quantile(site_rain$Annual)
hist(site_rain$Annual)
site_rain %>% 
  filter(Year %in% c(2013:2015)) %>% 
  select(Year, Annual) %>% 
  t()
boxplot(site_rain$Annual)
points(site_rain$Annual[site_rain$Year == 2013], col = "red", pch = 19, cex = 1.5)
points(site_rain$Annual[site_rain$Year == 2014], col = "blue", pch = 19, cex = 1.5)
points(site_rain$Annual[site_rain$Year == 2015], col = "orange", pch = 19, cex = 1.5)

abline(h = 268)

site_rain[site_rain$Annual == 268,]


site_temp1 <- read.csv("Data/station_weather/IDCJAC0002_067033_Data12.csv")
site_temp1[, c(-1:-3)] <- apply(site_temp1[, c(-1:-3)], 2, as.numeric)
site_temp2 <- read.csv("Data/station_weather/IDCJAC0002_067105_Data12.csv")
site_temp2[, c(-1:-3)] <- apply(site_temp2[, c(-1:-3)], 2, as.numeric)

site_temp1 <- site_temp1[complete.cases(site_temp1), ]
site_temp2 <- site_temp2[complete.cases(site_temp2), ]

site_temp <- bind_rows(site_temp1, site_temp2)
plot(Annual ~ Year, data = site_temp, type = "l")
abline(coef(lm(Annual ~ Year, data = site_temp)))
quantile(site_temp$Annual, probs = c(.9, .95, .97))
?quantile
filter(site_temp, Year %in% c(2013:2015)) %>% 
  select(Year, Annual)
boxplot(site_temp$Annual)
points(site_temp$Annual[site_temp$Year == 2013], col = "red", pch = 19, cex = 1.5)
points(site_temp$Annual[site_temp$Year == 2014], col = "blue", pch = 19, cex = 1.5)
points(site_temp$Annual[site_temp$Year == 2015], col = "orange", pch = 19, cex = 1.5)

1 - pnorm(c(25.2, 24.9, 24.1), mean = mean(site_temp$Annual), sd = sd(site_temp$Annual))

?qnorm

summary(site_temp)
summary(site_temp2)
site_rain2 <- site_rain %>% 
  rename(rain = Annual) %>% 
  select(Year, rain)
site_temp2 <- site_temp %>% 
  rename(temp = Annual) %>% 
  select(Year, temp)

site_clim <- left_join(site_rain2, site_temp2)
site_clim <- site_clim[complete.cases(site_clim), ]
plot(rain ~ temp, data = site_clim)
points(rain ~ temp, data = site_clim, subset = Year %in% c(2013:2015), col = "red", pch = 19)


msummary(site_temp1)
