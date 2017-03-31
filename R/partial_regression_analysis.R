devtools::install_github("hohenstein/remef")
library(remef)
summary(c4d_m0_bs)
p_c4d_m0_bs_moist <- keepef(c4d_m0_bs, fix = c("s_logmoist", "co2elev"))
p_c4d_m0_bs_temp  <- keepef(c4d_m0_bs, fix = c("s_temp",   "co2elev"))
p_c4d_m0_bs_par   <- keepef(c4d_m0_bs, fix = c("s_logpar", "co2elev"))
p_c4d_m0_bs_co2   <- keepef(c4d_m0_bs, fix = c("co2elev",  "co2elev"))

newd <- c34sum %>% 
  mutate(p_moist = p_c4d_m0_bs_moist, 
         p_temp  = p_c4d_m0_bs_temp, 
         p_par   = p_c4d_m0_bs_par,
         p_co2   = p_c4d_m0_bs_co2)

mcoefs <- coef(summary(c4d_m0_bs))[, "Estimate"]
par(mfrow = c(2, 2))
plot(p_moist ~ s_logmoist, newd, pch = 19, col = co2, ylim = c(-2.875544, 2.630993))
abline(mcoefs[1], mcoefs[3], col = 1)
abline(mcoefs[1]+mcoefs[2], mcoefs[3], col = 2)


plot(p_temp ~ s_temp, newd, pch = 19, col = co2, ylim = c(-2.875544, 2.630993))
abline(mcoefs[1], mcoefs[4], col = 1)
abline(mcoefs[1]+mcoefs[2], mcoefs[4], col = 2)


plot(p_par ~ s_logpar, newd, pch = 19, col = co2, ylim = c(-2.875544, 2.630993))
abline(mcoefs[1], mcoefs[5], col = 1)
abline(mcoefs[1]+mcoefs[2], mcoefs[5], col = 2)


boxplot(p_co2 ~ co2, newd, ylim = c(-2.875544, 2.630993))


plot(s_c4_ddiff ~ p_moist, data = newd)
abline(0, 1)
plot(s_c4_ddiff ~ p_temp , data = newd)
abline(0, 1)
plot(s_c4_ddiff ~ p_par  , data = newd)
abline(0, 1)
range(c34sum$s_c4_ddiff)


?keepef
