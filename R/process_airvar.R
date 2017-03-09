# download from HIEv
setToken(tokenfile = "Data/token.txt")

airvar_raw <- downloadTOA5(filename = "FACE_.*_AirVars_.*dat",
                          topath    = "Data/hievdata/raw_data/",
                          maxnfiles = 999)

save(airvar_raw, file = "output/Data/FACE_airvar_raw.RData")


# get daily mean, min and max temparater at two layers (at a height of 2m and 15m)
boxplot(airvar_raw$AirTC_1_Avg)
airvar_day <- airvar_raw %>%
  filter(Date        > as.Date("2012-09-18"),         # when co2 treatment was commenced
         Date        < as.Date("2016-2-15"),
         AirTC_1_Avg > -10) %>%                       # remove obviously weird value
  distinct() %>%                                      # remove duplicates
  mutate(ring = factor(substring(Source, 7, 7)),      # add ring number
         year = factor(year(Date))) %>%
  group_by(year, Date, ring) %>% 
  summarise_each(funs(Mean = mean(., na.rm = TRUE), 
                      Min  = min(., na.rm = TRUE),
                      Max  = max(., na.rm = TRUE),
                      N    = sum(!is.na(.))), 
                 airtemp2m = AirTC_1_Avg, airtemp15m = AirTC_2_Avg) %>% 
  group_by(year) %>% 
  mutate(annual_temp2m = mean(airtemp2m_Mean),
         T = max(27, 1.745 * annual_temp2m  + 11.143)[1]) %>% 
  ungroup() %>% 
  mutate(c3growth = airtemp2m_Min >= -1 & airtemp2m_Max >= 10 & airtemp2m_Max < 24, # c3 growing dates
         c4growth = airtemp2m_Max >= 21 & airtemp2m_Max < T)                        # c4 growtin dates

c34growthdate <- airvar_day %>% 
  select(year, Date, ring, c3growth, c4growth)
annualtemp <- c34growthdate %>% 
  select(year, annual_temp2m) %>% 
  distinct() %>% 
  filter(!year %in% c("2012", "2016")) %>% 
  mutate(year = mapvalues(year, c(2013, 2014, 2015), paste0("Year", 1:3)))

# For C 3 grasses, growth was considered possible in any month where the daily 
# minimum temperature was ≥ −1 oC, and the daily maximum temperature was ≥ 10 oC
# and < 24 oC. For C 4 grasses, growth was considered possible in months where 
# the daily maximum temperature was ≥ 21 °C and less than the tem- perature 
# define (Murphy and Bowman 2007).



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

c34growth_moist <- left_join(veg_moist, c34growthdate) %>% 
  filter(Date >= as.Date("2012-12-15"), Date <= as.Date("2016-2-15")) %>% 
  mutate(year = ifelse(Date < as.Date("2014-1-15"), "Year1", 
                       ifelse(Date < as.Date("2015-1-30"), "Year2", "Year3")),
         plot = factor(plot)) %>% 
  group_by(year, ring, plot) %>%
  summarise_each(funs(c3moist = sum(.[c3growth]),
                      c4moist = sum(.[c4growth]),
                      totalmoist = sum), 
                 Moist) %>% 
  ungroup() %>% 
  mutate(swa_c3 = c3moist / (c3moist + c4moist),  # seasonal water availability
         swa_c4 = c4moist / (c3moist + c4moist),
         plot = factor(plot)) %>% 
  left_join(annualtemp)


ggplot(c34growth_moist, aes(x = year, y = swa_c3))+
  geom_point() +
  facet_wrap(~ ring)


# test against C34 proportion
c34_prop <- PfgRDF$C3vsC4 %>% 
  arrange(id, desc(year)) %>% 
  group_by(id) %>% 
  mutate(dratios = c(NA, diff(ratios))) %>% 
  ungroup() %>% 
  filter(!is.na(dratios)) %>%
  left_join(c34growth_moist)
names(c34_prop)
m1 <- lmer(dratios ~ co2 * (swa_c3 + annual_temp2m) + (1|plot) + (1|id), data = c34_prop)
plot(m1)
head(c34_prop)

library(MuMIn)
options(na.action = "na.fail") 

# mm <- dredge(tm)
# mm
# # removing co2 increase more than 2 (>4) AICc units, indicating that co2 is
# # animportant term
# tm2 <- get.models(mm, 1)[[1]]
# tm3 <- update(tm2, REML = TRUE)
# Anova(tm3, test.statistic = "F")

?lag
fullm <- dredge(m1)
plot(m1)
Anova(m1, test.statistic = "F")

boxplot(dratios ~ co2 * year, data = c34_prop)
plot(logit(ratios) ~ swa_c3, data = c34_prop, col  = co2)
plot(logit(ratios) ~ c4moist, data = c34_prop, col  = co2)


plot(logit(dratios) ~ c3moist, data = c34_prop, col  = co2)
plot(dratios ~ c4moist, data = c34_prop, col  = co2)
ggplot(c34_prop, aes(x = c4moist, y = dratios, col = co2)) +
  geom_point() + 
  facet_wrap(~ co2)

top_n(c34_prop, 30)

boxplot(dratios~co2*year, data = c34_prop)




C3grassC4 <- veg_FullVdf %>% 
  filter(form == "Grass") %>% 
  mutate(yval = factor(ifelse(PFG == "c3", "p", "q")))

c34sum <- C3grassC4 %>%
  group_by(year, block, ring, plot, co2, id, PFG) %>% 
  summarise(value = sum(value)) %>% 
  spread(key = PFG, value = value) %>% 
  ungroup() %>% 
  mutate(id = ring:plot,
         c4r = c4 / (c3 + c4)) %>% 
  arrange(id, year) %>%
  group_by(id) %>%
  mutate_each(funs(delta = c(NA, diff(.))), c3, c4, c4r) %>%
  # mutate_each(funs(delta = (. + 1) / lag(. + 1) ), c3, c4) %>%
  filter(year != "Year0") %>%
  left_join(c34growth_moist)

  
summary(c34sum)
plot(c4_delta ~ c4moist, data = c34sum, col = co2)
plot(c4_delta ~ log(totalmoist), data = c34sum, col = co2)
plot(c3_delta ~ totalmoist, data = c34sum, col = co2)
plot(c3_delta ~ c3moist, data = c34sum, col = co2)


m1 <- lmer(c4_delta ~ co2  +log(totalmoist) + (1|ring) + (1|id), data = c34sum, REML = F)
m2 <- lmer(c4_delta ~ log(totalmoist) + (1|ring) + (1|id), data = c34sum, REML = F)
m3 <- lmer(c4_delta ~ co2 + (1|ring) + (1|id), data = c34sum, REML = F)
m4 <- lmer(c4_delta ~ 1 + (1|ring) + (1|id), data = c34sum, REML = F)
# m2 <- lmer(c4_delta ~ co2 + totalmoist + (1|ring) + (1|id), data = c34sum, REML = F)
# m3 <- lmer(c4_delta ~ co2 + log(c4moist) + (1|ring) + (1|id), data = c34sum, REML = F)

?layout

summary(m1)


dredge(m1)
dev.off()

pdf(file = "output/figs/c4delta.moist.pdf", width = 4, height = 6)
nf <- layout(matrix(c(1, 2)), 1, c(2, 1))
# layout.show(nf)
visreg(m1, xvar = "totalmoist", by = "co2", overlay = TRUE,
       xlab = "Cumulative moisture", ylab = "Annual change of C4")
abline(h = 0)
par(mar = c(5.1, 4.1, 0, 2.1))
boxplot(totalmoist ~ year, data = c34sum, horizontal  = TRUE,
        xlab = "cumulative moisture", cex.axis = .5)
dev.off()

# visreg(m1, xvar = "totalmoist", by = "co2", overlay = TRUE)
Anova(update(m1, REML = T), test.statistic = "F")
model.sel(m1, m2, m3, m4, extra = "r.squaredGLMM")
model.sel(m1, m2, m3, m4, extra = "deviance")

boxplot(log(totalmoist) ~ year, data = c34sum, horizontal  = TRUE)



m3 <- lmer(c4_delta ~ co2 + totalmoist + (1|ring) + (1|id), data = c34sum, REML = F)
par(mfrow = c(2, 1), mar = c(2, 4, 2, 2))
visreg(m3, xvar = "totalmoist", by = "co2", overlay = TRUE)
abline(h = 0)
boxplot(totalmoist ~ year, data = c34sum, horizontal  = TRUE)

m3 <- lmer(c3_delta ~ co2 * log(c3moist) + (1|ring) + (1|id), data = c34sum, REML = F)
m4 <- lmer(c3_delta ~ co2 * log(totalmoist) + (1|ring) + (1|id), data = c34sum, REML = F)

m5 <- lmer(c3_delta ~ log(totalmoist) + (1|ring) + (1|id), data = c34sum, REML = F)
plot(annual_temp2m ~ log(totalmoist), data = c34sum)
AICc(m3, m4)
dredge(m4, extra = "r.squaredGLMM")


dredge(m3)
dredge(m4)
dredge(m5)
Anova(update(m1, REML = TRUE), test.statistic = "F")
Anova(update(m3, REML = TRUE), test.statistic = "F")
visreg(m3, xvar = "totalmoist", by = "co2", overlay = TRUE)


visreg(m1, xvar = "totalmoist", by = co2, overlay = TRUE)
par(mfrow = c(1, 2))
visreg(m5, xvar = "meanmoist")
visreg(m5, xvar = "annual_temp2m")


anova(m1, m2)

m2 <- lmer(c3_delta ~ co2 * log(meanmoist) + (1|ring) + (1|id), data = c34sum, REML = F)
dredge(m1, extra = "r.squaredGLMM")
dredge(m2, extra = "r.squaredGLMM")


dredge(m2, extra = "r.squaredGLMM")
Anova(update(m1, REML = TRUE), test.statistic = "F")
Anova(update(m2, REML = TRUE), test.statistic = "F")

# m2 <- lmer(c4_delta ~ co2 + c3 + c4moist + (1|ring) + (1|id), data = c34sum)
visreg(m1, xvar = "c4moist", by = "co2", overlay = TRUE)
visreg(m1, xvar = "meanmoist", by = "co2", overlay = TRUE)
abline(h = 0)

g <- ggplot(aes(x = year, y = c4moist, col = co2), data = c34sum)+
  geom_boxplot()
g

fullm <- dredge(m1, extra = "r.squaredGLMM")
summary(model.avg(fullm, subset = delta <))

visreg(m1, xvar = "c4moist", by = "co2", overlay = TRUE)



boxplot(c34sum$c3)
plot(c4_delta ~ c3, data = c34sum)


# fig

figd <- left_join(veg_moist, airvar_day) %>% 
  filter(Date >= as.Date("2012-12-15"), Date <= as.Date("2016-2-15")) %>% 
  mutate(year = ifelse(Date < as.Date("2014-1-15"), "Year1", 
                       ifelse(Date < as.Date("2015-1-30"), "Year2", "Year3")),
         plot = factor(plot))

fig_m <- figd %>% 
  group_by(year, ring, plot) %>%
  summarise_each(funs(c3moist = sum(.[c3growth]),
                      c4moist = sum(.[c4growth]),
                      mean_moist = mean), 
                 Moist) %>% 
  group_by(year, ring) %>% 
  summarise_each(funs(mean), c3moist, c4moist, mean_moist)

fig_d <- figd %>% 
  select(Date, year, ring, c3growth, c4growth) %>% 
  distinct() %>% 
  group_by(year, ring) %>% 
  summarise_each(funs(sum), c3growth, c4growth) %>% 
  mutate(co2 = ifelse(ring %in% c(1, 4, 5), "elev", "amb"))

fig_d %>% 
  gather(key = "pfg", value = "days", c3growth, c4growth) %>% 
ggplot(., aes(x = year, y = days, col = co2))+
  geom_boxplot()+
  facet_grid(. ~ pfg)


  figd_t <- figd %>%
  group_by(year, ring) %>%
  summarise(mean_temp = mean(airtemp2m_Mean)) %>% 
  left_join(fig_m) %>% 
  left_join(fig_d)

summary(m1)
figd_t <- figd_t %>% 
  mutate(co2 = ifelse(ring %in% c(1, 4, 5), "elev", "amb"))
g <- ggplot(figd_t, aes(x = mean_moist, y = c4moist, size = c4growth,
                        col = co2, shape = year))+
  geom_point(alpha = .8)+
  geom_hline(yintercept = xa, linetype = "dotted", col = "red")+
  geom_hline(yintercept = xe, linetype = "dotted", col = "blue")+
  labs(x = "Mean annual moisture", y = "Seasonal water availability for C4")

g2 <- ggplot(figd_t, aes(x = mean_moist, y = c3moist, size = c3growth,
                        col = co2, shape = year))+
  geom_point(alpha = .8)+
  labs(x = "Mean annual moisture", y = "Seasonal water availability for C3")
g2  
g
visreg(m1, xvar = "c4moist", by = "co2", overlay = TRUE, ylab = "Delta C4",
       xlab = "Seasonal water availability for C4")
abline(h = 0)

ggsave(g, filename = "output/figs/SWA_c4.pdf", width = 5, height = 4)
ggsave(g2, filename = "output/figs/SWA_c3.pdf", width = 5, height = 4)
pdf(file = "output/figs/deltac4_miost.pdf", width = 4, height = 4)
visreg(m1, xvar = "c4moist", by = "co2", overlay = TRUE, ylab = "Delta C4",
       xlab = "Seasonal water availability for C4")
abline(h = 0)
dev.off()

# when delta C4 = 0
xa <- 7.4054/0.5526
xe <- 21.2866/0.5526


y = a + bx
a =  -13.8812-7.4054
b = 0.5526



