# random slope
c4d_m0_slp_1 <- lmer(c4_ddiff ~ co2 + log(c4moist) + log(c3) + (1+s_logmiost+s_c3|ring) + (1+s_logmiost+s_c3|year), data = c34sum)
c4d_m0_slp_2 <- lmer(c4_ddiff ~ 1 + (1+s_logmiost+s_c3|ring) + (1+s_logmiost+s_c3|year), data = c34sum)
c4d_m0_slp_3 <- lmer(c4_ddiff ~ co2 + log(c4moist) + log(c3) + (1+s_logmiost|ring) + (1+s_logmiost|year), data = c34sum)


c4d_m0_slp_4 <- lmer(c4_ddiff ~ co2 + log(c4moist) + (1+s_logmiost|ring) + (1+s_logmiost|year), data = c34sum)
c4d_m0_slp_4 <- lmer(c4_ddiff ~ co2 + log(c4moist) + (1+s_logmiost|id) + (1+s_logmiost|year), data = c34sum)

c4d_m0_slp_5 <- lmer(c4_ddiff ~ co2 + log(c4moist) + (1+s_logmiost|id) + (1|RY), data = c34sum)
c4d_m0_slp_7 <- lmer(c4_ddiff ~ co2 + (1+s_logmiost|ring/id) + (1|RY), data = c34sum)
c4d_m0_slp_7 <- lmer(c4_ddiff ~ 1 + (1+s_logmiost|co2/ring/plot) + (1|year:(ring:co2)), data = c34sum)


c4d_m0_slp_9   <- lmer(c4_ddiff ~ co2 * (log(c4moist)+ log(c3)) + (0+s_logmiost+s_c3|ring/id) + (0+s_logmiost + s_c3|year:ring), data = c34sum)
c4d_m0_slp_9   <- lmer(c4_ddiff ~ co2 * (log(c4moist)+ log(c3)) + 
                         (0+s_logmiost+s_c3|ring/id) + (1+s_logmiost + s_c3|year:ring), 
                       # control = lmerControl(optimizer ="Nelder_Mead"),
                       data = c34sum)

summary(c4d_m0_slp_9)

ggplot(c34sum, aes(x = s_logmiost, y = c4_ddiff, col = id))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~ring)

ggplot(c34sum, aes(x = s_c3, y = c4_ddiff, col = id))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ring)

ggplot(DivDF_list$grass_spp, aes(x = ring, y = S))+
  geom_boxplot()

# mds_all <- cmdscale(d = vegdist(prc_sp, method = "bray"), eig = TRUE, k = 2)  # MDS with bray-curtis dissimilarity
# library(BiodiversityR)
# add.spec.scores(mds_all, prc_sp)
# ?add.spec.scores
# round((mds_all$eig*100)/sum(mds_all$eig), 1)
# mds_all[, 1]
# ?cmdscale
# summary(mds_all)
# prcdf <- cbind(prc_site, mds_all$points)
# names(prcdf)[6:7] <- paste0("MDS", 1:2)
# prcdf <- prcdf %>% 
#   group_by(year, ring) %>% 
#   summarise_each(funs(mean, se), MDS1, MDS2)
# ggplot(prcdf, aes(x = MDS1_mean, y = MDS2_mean, col = ring))+
#   geom_point()+
#   geom_errorbar(aes(ymin = MDS2_mean-MDS2_se, ymax = MDS2_mean+MDS2_se))+
#   geom_errorbarh(aes(xmin = MDS1_mean-MDS1_se, xmax = MDS1_mean+MDS1_se))+
#   facet_wrap(~year)


summary(c4d_m0_slp_9)
m9full <- dredge(c4d_m0_slp_9, REML = F)
m9bs <- get.models(m9full, subset = 1)[[1]]
plot(m9bs)
qqnorm(resid(m9bs))
qqline(resid(m9bs))
Anova(m9bs, test.statistic = "F")
summary(m9bs)

summary(c4d_m0_slp_9)


# c4d_m0_slp_9   <- lmer(c4_ddiff ~ co2 * (log(c4moist)+ log(c3)) + (1+s_logmiost+s_c3|ring/id) + (1+s_logmiost + s_c3|year:ring), data = c34sum)

# c4d_m0_slp_9    <- lmer(c4_ddiff ~ co2 * (log(c4moist) + annual_temp2m) + 
#                           (0+s_logmiost|ring/id) + (1+s_logmiost|RY), data = c34sum)



# c4d_m0_slp_9    <- lmer(c4_ddiff ~ co2 * (log(c4moist)+ log(c3) + annual_temp2m) 
#                         + (0+s_logmiost+s_c3|ring/id) + (0+s_c3|year:ring) + (0+s_logmiost+s_c3|year),
#                         data = c34sum)
summary(c4d_m0_slp_9)



m9full <- dredge(c4d_m0_slp_9, REML = F, extra = "r.squaredGLMM")
m9full_nested <- subset(m9full, !nested(.))


m9bs <- get.models(m9full, subset = 1)[[1]]
m9avg <- model.avg(get.models(m9full, subset = delta <= 2))
summary(m9bs)
summary(m9avg)
summary(m1bs)
confint(m9avg, level = .9)
Anova(m9bs, test.statistic = "F")



sitedf <- c34sum %>% 
  select(ring, id, RY, co2) %>% 
  ungroup() %>% 
  distinct()
moistval <- seq(min(c34sum$c4moist), max(c34sum$c4moist), length.out = 1000)
c3val    <- median(c34sum$c3)


# c3 is median
c4d_m0_preddf <- ldply(1, function(x){
  
  newdf       <- ldply(1:10, function(y){ 
    cbind(sitedf, 
          c4moist  = moistval[sample(1000, nrow(sitedf), replace = TRUE)])
  })
  
  newdf <- newdf %>% 
    mutate(c3      = c3val[x])
  c4d_m0_pred    <- predict(m9avg, newdf, se.fit = TRUE, re.form = NA)
  c4d_m0_pred_df <- cbind(c4d_m0_pred, newdf) %>% 
    mutate(lwr = fit - se.fit * 1.96,
           upr = fit + se.fit * 1.96) 
  return(c4d_m0_pred_df)
})

c4d_p <- ggplot(c4d_m0_preddf, aes(x = log(c4moist), y = fit, col = co2)) +
  geom_line()+
  geom_line(aes(y = lwr), linetype = "dashed") +
  geom_line(aes(y = upr), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  science_theme
c4d_p



######
# against C3
moistval <- median(c34sum$c4moist)
c3val    <- seq(min(c34sum$c3), max(c34sum$c3), length.out = 1000)

c4d_m0_preddf <- ldply(1, function(x){
  
  newdf       <- ldply(1:10, function(y){ 
    cbind(sitedf, 
          c3  = c3val[sample(1000, nrow(sitedf), replace = TRUE)])
  })
  
  newdf <- newdf %>% 
    mutate(c4moist      = moistval[x])
  c4d_m0_pred    <- predict(m9avg, newdf, se.fit = TRUE, re.form = NA)
  c4d_m0_pred_df <- cbind(c4d_m0_pred, newdf) %>% 
    mutate(lwr = fit - se.fit * 1.96,
           upr = fit + se.fit * 1.96) 
  return(c4d_m0_pred_df)
})

c4d_p2 <- ggplot(c4d_m0_preddf, aes(x = c3, y = fit, col = co2)) +
  geom_line()+
  geom_line(aes(y = lwr), linetype = "dashed") +
  geom_line(aes(y = upr), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point(data = c34sum, aes(y = c4_ddiff)) +
  science_theme
c4d_p2

######









anova(m1bs, m9bs)

m1bs <- get.models(c4d_m0_full, subset = 1)[[1]]





c4d_m0_slp_10   <- lmer(c4_ddiff ~ co2 +  log(c4moist)+ log(c3) + (0+s_logmiost|ring/id) + (1+s_logmiost|year:ring), data = c34sum)
summary(c4d_m0_slp_10)

c4d_m0_slp_10  <- lmer(c4_ddiff ~ co2 + log(c4moist) + (0+s_logmiost|ring/id) + (1+s_logmiost|year:ring), data = c34sum, REML = F)
c4d_m0_slp_11  <- lmer(c4_ddiff ~ co2 + (0+s_logmiost|ring/id) + (1+s_logmiost|year:ring), data = c34sum, REML = F)
c4d_m0_slp_12  <- lmer(c4_ddiff ~ log(c4moist) + (0+s_logmiost|ring/id) + (1+s_logmiost|year:ring), data = c34sum, REML = F)
c4d_m0_slp_13  <- lmer(c4_ddiff ~ 1 + (1+s_logmiost|ring/id) + (1+s_logmiost|year:ring), data = c34sum, REML = F)
model.sel(c4d_m0_slp_9, c4d_m0_slp_10, c4d_m0_slp_11, c4d_m0_slp_12, c4d_m0_slp_13)
summary(c4d_m0_slp_9)

c4d_m0_slp_13  <- lmer(c4_ddiff ~ 1 + (1+s_logmiost|ring/id) + (1+s_logmiost|year:ring), data = c34sum)



# c4d_m0_slp_10  <- lmer(c4_ddiff ~ co2 * log(c4moist) + (1|ring/id) + (1|year:ring), data = c34sum)
# c4d_m0_slp_11  <- lmer(c4_ddiff ~ co2 * log(c4moist) + (1|ring) + (1|id) + (1|RY), data = c34sum)


summary(c4d_m0_slp_9)
summary(c4d_m0_slp_10)
summary(c4d_m0_slp_11)
Anova(c4d_m0_slp_9, test.statistic = "F")

c4d_m0_slp_10  <- lmer(c4_ddiff ~ co2 + log(c4moist) + (1+s_logmiost|ring/id) + (1|year:ring) + (1+s_logmiost|year), data = c34sum)
summary(c4d_m0_slp_10)
ll <- data.frame((VarCorr(c4d_m0_slp_10)))
round((ll$vcov * 100)/sum(ll$vcov), 2)


var(c34sum$c4_ddiff)

summary(c4d_m0_slp_10)
s10 <- summary(c4d_m0_slp_10)
names(s10)
s10$vcov
ll <- data.frame(s10$varcor) %>% 
  filter(is.na(var2))
sum(ll$vcov)
var(c34sum$c4_ddiff)



s10$vcov


c4d_m0_slp_10@pp

summary(c4d_m0_slp_7)

Anova(c4d_m0_slp_7, test.statistic = "F")


c4d_m0_slp_6 <- lmer(c4_ddiff ~ co2 + log(c4moist) + (1|ring) + (1|id) + (1|RY), data = c34sum)
c4d_m0_slp_7 <- lmer(c4_ddiff ~ co2 + log(c4moist) + (1|id) + (1|RY), data = c34sum)
summary(c4d_m0_slp_5)

confint(c4d_m0_slp_1, method = "boot")
Anova(c4d_m0_slp_1, test.statistic = "F")
summary(c4d_m0_slp_1)
summary(c4d_m0_slp_2)
summary(c4d_m0_slp_3)
summary(c4d_m0_slp_4)
summary(c4d_m0_slp_5)
ggplot(c34sum, aes(x = s_logmiost, y = c4_ddiff, col = ring))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  ggtitle("Random slopes bewteen rings")

ggplot(c34sum, aes(x = s_logmiost, y = c4_ddiff, col = year))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  ggtitle("Random slopes bewteen years")

ggplot(c34sum, aes(x = s_logmiost, y = c4_ddiff, col = id))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  ggtitle("Random slopes bewteen IDs")


?geom_smooth
