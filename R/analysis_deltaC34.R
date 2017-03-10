# Here we will:
# 1) download moisture and understorey air temperature data from HIEv
# 2) identify growing season for each of C3 and C4 species according to Murphy 2007
# 3) Fit moisture and temperature to delta C3 and C4 (annual change in abundance)



# download from HIEv ------------------------------------------------------


setToken(tokenfile = "Data/token.txt")

airvar_raw <- downloadTOA5(filename = "FACE_.*_AirVars_.*dat",
                          topath    = "Data/hievdata/raw_data/",
                          maxnfiles = 999)

save(airvar_raw, file = "output/Data/FACE_airvar_raw.RData")
# load("output/Data/FACE_airvar_raw.RData")




# process HIEv data -------------------------------------------------------


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
  select(year, Date, ring, c3growth, c4growth, annual_temp2m)
annualtemp <- c34growthdate %>% 
  select(year, annual_temp2m) %>% 
  distinct() %>% 
  filter(!year %in% c("2012", "2016")) %>% 
  mutate(year = mapvalues(year, c(2013, 2014, 2015), paste0("Year", 1:3)))

## Murphy 2007. For C 3 grasses, growth was considered possible in any month 
## where the daily minimum temperature was ≥ −1 oC, and the daily maximum 
## temperature was ≥ 10 oC and < 24 oC. For C 4 grasses, growth was considered 
## possible in months where the daily maximum temperature was ≥ 21 °C and less 
## than the tem- perature define (Murphy and Bowman 2007).




# assign moisture data to each vegetation plot based on their corr --------


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

plot(swa_c4 ~ totalmoist, data = c34growth_moist)
plot(swa_c3 ~ totalmoist, data = c34growth_moist)
plot(c4moist ~ totalmoist, data = c34growth_moist)
c34growth_moist %>% 
  gather(key = category, value = moisture, c3moist, c4moist, totalmoist) %>%
  group_by(category, year, ring) %>% 
  summarise(moisture = mean(moisture)) %>% 
ggplot(., aes(x = year, y = moisture, col = category))+
  geom_path(aes(group =  category)) +
  facet_wrap(~ ring)


# merge with C34 ratios df ---------------------------------------------

c34_prop <- PfgRDF$C3vsC4 %>% 
  arrange(id, year) %>% 
  group_by(id) %>% 
  mutate_each(funs(rat_diff = . - lag(., 1),             # Year1-Year0 and etc.
                   rat_prop = (. + 1) / lag(. + 1, 1)),  # Year1/Year0 and etc.
              ratios) 
    
# move year0 to a new column so that it can be used as a covariate
c34_prop_year0 <- c34_prop %>% 
  select(year, id, ratios) %>% 
  filter(year == "Year0") %>% 
  rename(ratios0 = ratios) %>%
  select(-year) %>% 
  right_join(c34_prop) %>% 
  left_join(c34growth_moist) %>% 
  filter(year != "Year0")



# . ratios  --------------------------------------------------------------
names(c34_prop_year0)
summary(c34_prop_year0)
ratd_m1 <- lmer(ratios ~ co2 * (swa_c4     + annual_temp2m) + ratios0 + (1|ring) + (1|RY) + (1|id), data = c34_prop_year0)
ratd_m2 <- lmer(ratios ~ co2 * (c4moist    + annual_temp2m) + ratios0 + (1|ring) + (1|RY) + (1|id), data = c34_prop_year0)
ratd_m3 <- lmer(ratios ~ co2 * (totalmoist + annual_temp2m) + ratios0 + (1|ring) + (1|RY) + (1|id), data = c34_prop_year0)


plot(ratios ~ totalmoist, data = c34_prop_year0, col = co2)
plot(ratios ~ annual_temp2m, data = c34_prop_year0, col = co2)
Anova(ratd_m3, test.statistic = "F")
ratd_m3_full <- dredge(ratd_m3, REML = F)
ratd_m3_avg <- model.avg(get.models(ratd_m3_full, subset = cumsum(weight) <= .95))
summary(ratd_m3_avg)
confint(ratd_m3_avg)

# . predicted values ------------------------------------------------------

sitedf <- c34_prop_year0 %>% 
  select(ring, id, RY, co2) %>% 
  ungroup() %>% 
  distinct()
# moistval <- seq(min(c34_prop_year0$totalmoist), max(c34_prop_year0$totalmoist), length.out = 1000)
# tempval <- quantile(c34_prop_year0$annual_temp2m)[4]

moistval <- quantile(c34_prop_year0$totalmoist)[4]
tempval  <- seq(min(c34_prop_year0$annual_temp2m), max(c34_prop_year0$annual_temp2m), length.out = 1000)

newdf    <- ldply(1:10, function(x) 
  cbind(sitedf, 
        annual_temp2m  = tempval[sample(1000, nrow(sitedf), replace = TRUE)])) %>% 
  mutate(totalmoist = moistval,
         ratios0       = mean(c34_prop_year0$ratios0))

ratd_m3_pred <- predict(ratd_m3, newdf, se.fit = TRUE, re.form = NA)
ratd_m3_pred_df <- cbind(ratd_m3_pred, newdf) %>% 
  mutate(lwr = fit - se.fit * 1.96,
         upr = fit + se.fit * 1.96) 
ggplot(ratd_m3_pred_df, aes(x = annual_temp2m, y = fit, col = co2)) +
  geom_line()+
  geom_line(aes(y = lwr), linetype = "dashed") +
  geom_line(aes(y = upr), linetype = "dashed")
names(ratd_m3_pred_df)

summary(ratd_m1)
summary(c34_prop_year0)


range(PfgRDF$C3vsC4$ratios)
test <- sample(10)
test - lag(test, 1)
test/lag(test, 1)
  

            
    rat_diff = c(NA, diff(ratios)),
         rat_prop = (. + 1) / lag(. + 1))%>% 
  
  
  ungroup() %>% 
  filter(!is.na(dratios)) %>%
  mutate(year = mapvalues(year, paste0("Year", 0:2), paste0("Year", 1:3))) %>% 
  left_join(c34growth_moist)








summary(c34_prop)

m1 <- lmer(dratios ~ co2 * (swa_c4 + annual_temp2m) + (1|ring) + (1|id) + (1|RY), data = c34_prop)
m2 <- lmer(dratios ~ co2 * (swa_c4 + annual_temp2m) + (1|ring), data = c34_prop)

plot(m1)
plot(m2)


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
fullm <- dredge(update(m2, REML = F))
plot(m1)
Anova(m2, test.statistic = "F")

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




# C3grassC4 <- veg_FullVdf %>% 
#   filter(form == "Grass") %>% 
#   mutate(yval = factor(ifelse(PFG == "c4", "p", "q")))

names(C3grassC4)
c34sum <- C3grassC4 %>%
  group_by(year, block, ring, plot, co2, id, PFG, RY) %>% 
  summarise(value = sum(value)) %>% 
  spread(key = PFG, value = value) %>% 
  ungroup() %>% 
  mutate(c4r = c4 / (c3 + c4)) %>% 
  arrange(id, year) %>%
  group_by(id) %>%
  mutate_each(funs(delta = c(NA, diff(.))), c3, c4, c4r) %>%
  # mutate_each(funs(delta = (. + 1) / lag(. + 1) ), c3, c4) %>%
  filter(year != "Year0") %>%
  left_join(c34growth_moist)

summary(c34sum)
plot(c4_delta ~ c4moist, data = c34sum, col = co2)
plot(c4_delta ~ log(c4moist), data = c34sum, col = co2)
plot(c4_delta ~ log(totalmoist), data = c34sum, col = co2)
plot(c3_delta ~ totalmoist, data = c34sum, col = co2)
plot(c3_delta ~ c3moist, data = c34sum, col = co2)

summary(c34sum)
m1 <- lmer(c4_delta ~ co2  * (log(c4moist) + annual_temp2m) + (1|ring) + (1|id) + (1|RY), data = c34sum)

c34sum <- c34sum %>% 
  mutate(logc4moist = log(c4moist))

xyplot(c4_delta ~ logc4moist | co2, group = id, data = c34sum, type=c("p","r"))
xyplot(c4_delta ~ logc4moist | co2, group = RY, data = c34sum, type=c("p","r"))
xyplot(c4_delta ~ logc4moist | co2, group = ring, data = c34sum, type=c("p","r"))


xyplot(c4_delta ~ annual_temp2m | co2, group = id, data = c34sum, type=c("p","r"))
xyplot(c4_delta ~ annual_temp2m | co2, group = RY, data = c34sum, type=c("p","r"))
xyplot(c4_delta ~ annual_temp2m | co2, group = ring, data = c34sum, type=c("p","r"))

xyplot(Claims/Holders ~ Age | Group, groups=District, 
       data = c34sum, t="l", main="Claims frequency",
       auto.key=list(space="top", columns=4, 
                     title="District", cex.title=1,
                     lines=TRUE, points=FALSE))

m1 <- lmer(c4_delta ~ co2  * log(c4moist) + (1 + c4moist|ring) + (1 + c4moist|id) + (1 + c4moist|RY), data = c34sum)
m2 <- lmer(c4_delta ~ co2  * log(c4moist) + (1|ring) + (1|id) + (1|RY), data = c34sum)



m1full <- dredge(m1, REML = F)
m1ave <- model.avg(get.models(m1full, cumsum(weight) <= .95))
confint(m1ave)
confint(m1, method = "boot")

summary(m1ave)





predm1 <- predict(m1ave, newdf)
predm1 <- cbind(newdf, predm1)
ggplot(predm1, aes(x = c4moist, y = predm1, col = co2))+
  geom_point(alpha = .2) +
  # geom_point(data = c34sum, aes(x = c4moist, y = c4_delta, col = co2)) +
  geom_hline(yintercept = 0)
  


library(merTools)
# preds2 <- predictInterval(m2, newdata = newdf, n.sims = 10)

ranfoms <- c("ring", "id", "RY")


random_l1 <- llply(1:3, function(n){
  combn(ranfoms, n, FUN = function(x) paste0("(1|", x, ")", collapse = "+"))
})
random_l2 <- llply(1:3, function(n){
  combn(ranfoms, n, FUN = function(x) paste0("(1+logc4moist|", x, ")", collapse = "+"))
})

random_l3 <- combn(ranfoms, 2, FUN = function(x) paste0("(1+logc4moist|", x[1], ")+", "(1|", x[2], ")"))
random_l4 <- combn(ranfoms, 2, FUN = function(x) paste0("(1|", x[1], ")+", "(1+logc4moist|", x[2], ")"))
random_l5 <- combn(ranfoms, 2, FUN = function(x) paste0("(1+logc4moist|", x[1], ")+", "(1+logc4moist|", x[2], ")"))


random_l6 <- combn(ranfoms, 3, FUN = function(x) paste0("(1+logc4moist|", x[1], ")+", "(1|", x[2], ")+(1|", x[3], ")"))
random_l7 <- combn(ranfoms, 3, FUN = function(x) paste0("(1|", x[1], ")+", "(1+logc4moist|", x[2], ")+(1|", x[3], ")"))
random_l8 <- combn(ranfoms, 3, FUN = function(x) paste0("(1|", x[1], ")+", "(1|", x[2], ")+(1+logc4moist|", x[3], ")"))

random_l9  <- combn(ranfoms, 3, FUN = function(x) paste0("(1+logc4moist|", x[1], ")+", "(1+logc4moist|", x[2], ")+(1|", x[3], ")"))
random_l10 <- combn(ranfoms, 3, FUN = function(x) paste0("(1|", x[1], ")+", "(1+logc4moist|", x[2], ")+(1+logc4moist|", x[3], ")"))
random_l11 <- combn(ranfoms, 3, FUN = function(x) paste0("(1+logc4moist|", x[1], ")+", "(1|", x[2], ")+(1+logc4moist|", x[3], ")"))

random_l12 <- combn(ranfoms, 3, FUN = function(x) paste0("(1+logc4moist|", x[1], ")+", "(1+logc4moist|", x[2], ")+(1+logc4moist|", x[3], ")"))

random_l <- unlist(list(random_l1, random_l2,random_l3,random_l4,random_l5,random_l6,random_l7,random_l8,random_l9,random_l10, random_l11, random_l12))





fl <- llply(random_l, function(x) {
  f  <- paste("c4_delta ~ co2  * log(c4moist) +", x)
  # f2 <- paste("c4_delta ~ co2  * c4moist +", x)
  return(formula(f))
})
# fl <- unlist(fl)
mls <- llply(fl, function(x) lmer(x, data = c34sum))
model.sel(mls)

m3 <- lmer(c4_delta ~ co2  * log(c4moist) + (1|RY) + (1|ring) + (1|id), data = c34sum)
m4 <- lmer(c4_delta ~ co2  + log(c4moist) + (1|RY), data = c34sum)
m3full <- dredge(m3, REML = F)
m3ave <- model.avg(get.models(m3full, subset = cumsum(weight) <= .95), fit = TRUE)
m3ave <- model.avg(get.models(m3full, subset = delta <= 2), fit = TRUE)

?bootMer

confint(m3ave)
summary(m3ave)

m3 <- lmer(c4_delta ~ co2  + log(c4moist) + (1|RY) + (1|ring) + (1|id), data = c34sum)
pm3 <- predictInterval(m3, newdata = c34sum, 
                       level = 0.95, n.sims = 1000,
                       stat = "median", type="linear.prediction",
                       include.resid.var = TRUE)
pm3df <- bind_cols(c34sum, pm3)
names(pm3df)
ggplot(pm3df, aes(x = c4moist, y = fit, col = co2, ymin=lwr, ymax=upr))+
  geom_point() +
  geom_linerange()+
  geom_point(aes(y = c4_delta, shape = co2), col = "black")


PI.arm.sims <- arm::sim(m3ave, 50)

PI.arm <- data.frame(
  fit=apply(fitted(PI.arm.sims, m3), 1, function(x) quantile(x, 0.500)),
  upr=apply(fitted(PI.arm.sims, m3), 1, function(x) quantile(x, 0.975)),
  lwr=apply(fitted(PI.arm.sims, m3), 1, function(x) quantile(x, 0.025))
)

newdat<-data.frame(x=seq(0,10,length=20))
mm<-model.matrix(~runif(100,0,10),newdf)
# mm%*%
  
fixef(m3)


?fixef
PIdf <- bind_cols(c34sum, PI.arm)
head(PIdf)
ggplot(PIdf, aes(x = log(c4moist), y = fit, ymin = lwr, ymax = upr, col = co2)) +
  geom_point() +
  geom_linerange() +
  geom_point(aes(y = c4_delta, shape = co2), col = "black")

predave <- predict(m3ave, newdata = newdf, se.fit = TRUE, re.form = NA)  
predave

Anova(m3, test.statistic = "F")

anova(m3)


preds2 <- predictInterval(m3, newdata = filter(c34sum, ring == "1"), n.sims = 10, which = "fixed")
preds2 <- bind_cols(filter(c34sum, ring == "1"), preds2)
ggplot(preds2, aes(x = c4moist, y = fit, col = co2))+
  geom_line(aes(x = c4moist, y = lwr), size = .5) +
  geom_line(aes(x = c4moist, y = upr), size = .5) +
  geom_point() +
  # geom_point(data = c34sum, aes(x = c4moist, y = c4_delta, col = co2)) +
  geom_hline(yintercept = 0)


?predict
?get.models
m2 <- lmer(c4_delta ~ co2  * log(c4moist) + (1 + c4moist|ring) + (1 + c4moist|id) + (1|RY), data = c34sum)
m3 <- lmer(c4_delta ~ co2  * log(c4moist) + (1 + c4moist|ring) + (1|id) + (1 + c4moist|RY), data = c34sum)
m4 <- lmer(c4_delta ~ co2  * log(c4moist) + (1 |ring) + (1 + c4moist|id) + (1 + c4moist|RY), data = c34sum)
m5 <- lmer(c4_delta ~ co2  * log(c4moist) + (1 |ring) + (1|id) + (1|RY), data = c34sum)
m6 <- lm(c4_delta ~ co2  * log(c4moist), data = c34sum)
AICc(m5, m6)

dredge(m6)

anova(m2, m5)
Anova(update(m2, REML = T), test.statistic = "F")
visreg(m2, xvar = "co2")
anova(m1, m2, m3, m4)
model.sel(m1, m2, m3, m4, m5, rank = "deviance", extra = "AICc")
summary(m1)
summary(m2)
summary(m3)
summary(m4)



AICc(m1, m2, m3, m4)
summary(m1)

Anova(update(m1, REML = T), test.statistic = "F")



m2 <- lmer(c4_delta ~ co2  * (log(c4moist) + annual_temp2m) + (1|ring) + (1|RY), data = c34sum)
m3 <- lmer(c4_delta ~ co2  * (log(c4moist) + annual_temp2m) + (1|ring) + (1|id), data = c34sum)
m4 <- lmer(c4_delta ~ co2  * (log(c4moist) + annual_temp2m) + (1|ring), data = c34sum)
# m4 <- lmer(c4_delta ~ co2  * c4moist + (1|ring), data = c34sum, REML = F)
# m5 <- lmer(c4_delta ~ co2  * totalmoist + (1|ring), data = c34sum, REML = F)


model.sel(m1, m2, m3, m4)

m4full <- dredge(update(m6, REML = F))
m4full <- dredge(m6)
m4ave <- model.avg(get.models(m4full, subset = cumsum(weight) < .95))
summary(m4ave)
confint(m4ave)
Anova(m4, test.statistic = "F")
Anova(m1, test.statistic = "F")
visreg(m4, xvar = "")
m2 <- lmer(c4_delta ~ co2  * log(c4moist) + (1|ring) + (1|id) + (1|RY), data = c34sum, REML = F)


Anova(update(m1, REML = T), test.statistic = "F")

m2 <- lmer(c4_delta ~ log(totalmoist) + (1|ring) + (1|id) + (1|RY), data = c34sum, REML = F)
m3 <- lmer(c4_delta ~ co2 + (1|ring) + (1|id) + (1|RY), data = c34sum, REML = F)
m4 <- lmer(c4_delta ~ 1 + (1|ring) + (1|id) + (1|RY), data = c34sum, REML = F)
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



