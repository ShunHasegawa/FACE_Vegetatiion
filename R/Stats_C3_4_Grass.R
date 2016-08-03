head(veg)

# Create DFs for C3total:C4, C3grass:C4, legume:non-legume and Native:introduced
C3grassC4 <- within(subsetD(veg, form %in% c("Sedge", "Grass") & PFG %in% c("c3", "c4")), {
  yval <- factor(ifelse(PFG == "c3", "p", "q"))  
})
C3totalC4 <- within(subsetD(veg, !form %in% c("Moss", "Fern")), {
  yval <- factor(ifelse(PFG == "c3", "p", "q"))  
})

legumeR <- within(subsetD(veg, form == "Forb"), {
  yval <- factor(ifelse(PFG == "legume", "p", "q"))  
})

NativeR <- within(subsetD(veg, !is.na(origin)), {
  yval <- factor(ifelse(origin == "native", "p", "q"))  
})

dfList <- list(C3grassC4 = C3grassC4, C3totalC4 = C3totalC4, legumeR = legumeR, NativeR = NativeR)

# compute ratios and total number
PfgRDF <- llply(dfList, function(x){
  x %>% 
    group_by(year, co2, block, ring, plot, id) %>%
    summarise(Total = sum(value),
              ratios = sum(value[yval == "p"]/Total)) %>%
    ungroup() # grouping informaiton is not required later
  })

# Move Year0 value to a new column to be used as covariate in the analyssis
# below
PfgRDF_year0 <- llply(PfgRDF, function(x) {
  
  year0_dd <-x %>%                   # Year0
    filter(year == "Year0") %>%
    select(id, ratios) %>%
    rename(ratios0 = ratios)

  subyear_dd <- filter(x, year != "Year0") # subsequent years

  newyear_dd <- merge(subyear_dd, year0_dd, by = "id") # merge
  
  newyear_dd$obs <- 1:nrow(newyear_dd) # add id
  return(newyear_dd)
})

length(PfgRDF_year0)
par(mfrow = c(2, 2), mar = c(5, 4, 2, 1))
l_ply(1:4, function(x) {
  d <- PfgRDF_year0[[x]]
  plot(ratios ~ ratios0, pch = 19, col = year, data = d, main = names(PfgRDF_year0)[x])})


# C3:C4 (grass + sedge) ---------------------------------------------------

c3gc4DF <- PfgRDF_year0[[1]]
head(c3gc4DF)
par(mfrow = c(1, 3))
boxplot(logit(ratios) ~ year:co2, data = c3gc4DF, main = "logit")
boxplot(logit(ratios) ~ year:ring, data = c3gc4DF, main = "logit")
plot(logit(ratios) ~ ratios0, data = c3gc4DF)


# . GLM --------------------------------------------------------------------

m1 <- glmer(ratios ~ year*co2 + ratios0 + (1|block) + (1|ring)  + (1|id), 
            family = "binomial", weights = Total, data = c3gc4DF)
# model failed to converge

# . LMM -------------------------------------------------------------------

m3 <- lmer(logit(ratios) ~ year*co2 + ratios0 + (1|block) + (1|ring)  + (1|id), data = c3gc4DF)
plot(m3)
qqnorm(resid(m3))
qqline(resid(m3))

AnvF_c3gc4 <- Anova(m3, test.statistic = "F")
AnvF_c3gc4

# . post-hoc test ---------------------------------------------------------
plot(lmerTest::lsmeans(m3))
lsmeans::lsmeans(m3, pairwise ~ co2 | year)


# C3 total:C4 -------------------------------------------------------------

###############
# C3 total:C4 #
###############

names(PfgRDF_year0)[2]
c3totalc4DF <- PfgRDF_year0[[2]]


# . GLM -------------------------------------------------------------------

m1 <- glmer(ratios ~ year*co2 + ratios0 + (1|block) + (1|ring)  + (1|id), 
            family = "binomial", weights = Total, data = c3totalc4DF)
overdisp.glmer(m1)

# overdispersed
m2 <- update(m1, ~ . + (1|obs))
  # model railed to converge

# . LMM -------------------------------------------------------------------
par(mfrow = c(2, 2))
boxplot(logit(ratios) ~ year:co2, data = c3totalc4DF, main = "logit")
boxplot(logit(ratios) ~ year:ring, data = c3totalc4DF, main = "logit")
plot(logit(ratios) ~ ratios0, data = c3totalc4DF, col = year, pch = 19)

m3 <- lmer(logit(ratios) ~ year * co2 + ratios0 + (1|block) + (1|ring)  + (1|id), 
           data = c3totalc4DF)
plot(m3)
qqnorm(resid(m3))
qqline(resid(m3))
AnvF_c3totalc4 <- Anova(m3, test.statistic = "F")
AnvF_c3totalc4


# Legume:nonlegume --------------------------------------------------------

####################
# Legume:nonlegume #
####################
names(PfgRDF_year0)[3]
legumeRDF <- PfgRDF_year0[[3]]


# . GLM -------------------------------------------------------------------

m1 <- glmer(ratios ~ year*co2 + ratios0 + (1|block) + (1|ring)  + (1|id), 
            family = "binomial", weights = Total, data = legumeRDF)
summary(m1)
overdisp.glmer(m1)
# overdispersed
m2 <- update(m1, ~ . + (1|obs))

# convergence error. so try different optimizer
m1.1 <- glmer(ratios ~ year*co2 + ratios0 + (1|block) + (1|ring)  + (1|id), 
              family = "binomial", weights = Total, data = legumeRDF, 
              control = glmerControl(optimizer = "bobyqa"))
overdisp.glmer(m1.1)
# overdispersed
m2 <- update(m1.1, ~ . + (1|obs))
Anova(m2)
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))
legumeR_CompAic <- CompAIC(m2)
legumeR_CompAic


# . LMM -------------------------------------------------------------------

par(mfrow = c(2, 2))
boxplot(logit(ratios) ~ year:co2, data = legumeRDF, main = "logit")
boxplot(logit(ratios) ~ year:ring, data = legumeRDF, main = "logit")
plot(logit(ratios) ~ ratios0, data = legumeRDF, pch = 19, col = year)

m3 <- lmer(logit(ratios) ~ year*co2 + ratios0 + (1|block) + (1|ring)  + (1|id), 
           data = legumeRDF)
plot(m3)
qqnorm(resid(m3))
qqline(resid(m3))
AnvF_legumeR <- Anova(m3, test.statistic = "F")
AnvF_legumeR


# Native:introduced -------------------------------------------------------

#####################
# Native:introduced #
#####################
names(PfgRDF_year0)[4]
NativeRDF <- PfgRDF_year0[[4]]


# . GLM -------------------------------------------------------------------

m1 <- glmer(ratios ~ year*co2 + (1|block) + (1|ring) + (1|id), 
            family = "binomial", weights = Total, data = NativeRDF)
overdisp.glmer(m1)
# overdispersed
m2 <- update(m1, ~ . + (1|obs))
overdisp.glmer(m2)
Anova(m2)
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))
NativeR_CompAic <- CompAIC(m2)
NativeR_CompAic
# indication of year effect


# . LMM -------------------------------------------------------------------

par(mfrow = c(2, 2))
boxplot(logit(ratios) ~ year:co2, data = NativeRDF, main = "logit")
boxplot(logit(ratios) ~ year:ring, data = NativeRDF, main = "logit")
plot(logit(ratios) ~ ratios0, pch = 19, col = year, data = NativeRDF)


m3 <- lmer(logit(ratios) ~ year * co2 + ratios0 + (1|block) + (1|ring) + (1|id), 
           data = NativeRDF)
plot(m3)
qqnorm(resid(m3))
qqline(resid(m3))
AnvF_NativeR <- Anova(m3, test.statistic = "F")
AnvF_NativeR


# . post-hoc test ---------------------------------------------------------

plot(lmerTest::lsmeans(m3))
lsmeans::lsmeans(m3, pairwise ~ co2 | year)


# Summary -----------------------------------------------------------------

###########
# Summary #
###########

# F from LMM
SummaryAnvF_PFG <- ldply(list(c3gc4     = AnvF_c3gc4, 
                              c3totalc4 = AnvF_c3totalc4,
                              legumeR   = AnvF_legumeR,
                              NativeR   = AnvF_NativeR ), 
                         function(x)data.frame(x, terms = row.names(x)),
                         .id = "variable")
summary(SummaryAnvF_PFG)

SummaryAnvF_PFG <- within(SummaryAnvF_PFG, {
  F = round(F, 2)
  Df.res = round(Df.res, 0)
  Pr = round(SummaryAnvF_PFG$Pr, 3)
  Pr..F. = NULL})
SummaryAnvF_PFG$terms <- factor(SummaryAnvF_PFG$terms, 
                                levels = c("ratios0","co2", "year", "year:co2"))

write.csv(SummaryAnvF_PFG, file = "output/table/FACE_EachPFG_Prop_AnvovaF.csv", 
          row.names = FALSE)


# CO2 response ratios -----------------------------------------------------

#######################
# CO2 response ratios #
#######################

# total number and count
PfgRAbundDF <- ldply(dfList, function(x) 
  ddply(x, .(year, co2, block, ring), function(y){
    Total <- sum(y$value)
    value <- sum(y$value[y$yval == "p"])
    data.frame(Total, value)}), 
  .id = "variable")

# there are some 0s so add 1
PfgRAbundDF[c("Total", "value")] <- PfgRAbundDF[c("Total", "value")] + 1

# compute SE for ratio using non-parametric bootstrap

# cast to pair each block
PfgRAbundDF_mlt <- melt(PfgRAbundDF, id = c("variable", "year", "co2", "block", "ring"), 
                        variable_name = "count")
PfgRAbundDF_cst <- dcast(PfgRAbundDF_mlt, variable + year + block ~ co2 + count)

ratio <- function(d, w) {
  AmbM  <- with(d, sum(amb_value * w) / sum(amb_Total * w))
  elevM <- with(d, sum(elev_value * w)/ sum(elev_Total * w))
  elevM/AmbM- 1
}

RatioSE <- ddply(PfgRAbundDF_cst, .(year, variable), function(x) {
  b <- boot::boot(x, ratio, R = 999, stype = "w")
  summary(b)
})
RatioSE$co2R <- RatioSE$original

# organise it to export as a table
RatioSE$co2R_se <- with(RatioSE, paste0(round(co2R, 2), "(", round(bootSE, 2), ")"))

RatioSE_cst <- dcast(RatioSE[c("variable", "year", "co2R_se")], variable ~ year, value.var = "co2R_se")
write.csv(RatioSE_cst, file = "output/table/CO2ResponseRatio_SpecificPFG.csv", row.names = FALSE)

# plot

# block mean
blockMean <- ddply(PfgRAbundDF, .(year, variable, block), 
                   function(x) {
                     am <- with(x, value[co2 == "amb"] / Total[co2 == "amb"])
                     em <- with(x, value[co2 == "elev"]/ Total[co2 == "elev"])
                     co2R <- em/am - 1
                     data.frame(co2R) 
                   })

labs <- c(expression(C3[grass]:C4), expression(C3[total]:C4), 
          "Legume:Non-legume", "Native:Introduced")
blockMean$variable <- factor(blockMean$variable, labels = labs)
RatioSE$variable <- factor(RatioSE$variable, labels = labs)

p <- ggplot(RatioSE, aes(x = year, y = co2R))
p2 <- p + 
  geom_point(aes(x = year, y = co2R), size = 4)+
  geom_errorbar(aes(x = year, ymin = co2R - bootSE, ymax = co2R + bootSE), width = 0) + 
  geom_point(data = blockMean, size = 2, col = "red", alpha = .7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "CO2 response ratio", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 7)) +
  facet_wrap( ~ variable, scale = "free_y")
p2
ggsavePP(plot = p2, filename = "output/figs/FACE_CO2ResponseRatio_SpecificPFG", width = 4, height = 4)
