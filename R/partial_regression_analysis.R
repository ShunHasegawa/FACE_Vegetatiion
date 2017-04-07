devtools::install_github("hohenstein/remef")
library(remef)
summary(c4d_m0_bs)
?keepef
p_c4d_m0_bs_moist <- keepef(c4d_m0_bs, fix = c("s_logmoist"))
p_c4d_m0_bs_temp  <- keepef(c4d_m0_bs, fix = c("s_temp"    ))
p_c4d_m0_bs_par   <- keepef(c4d_m0_bs, fix = c("s_logpar"  ))
p_c4d_m0_bs_co2   <- keepef(c4d_m0_bs, fix = c("co2elev"   ))

newd <- c34sum %>% 
  mutate(p_moist = p_c4d_m0_bs_moist, 
         p_temp  = p_c4d_m0_bs_temp, 
         p_par   = p_c4d_m0_bs_par,
         p_co2   = p_c4d_m0_bs_co2)

mcoefs <- coef(summary(c4d_m0_bs))[, "Estimate"]
par(mfrow = c(2, 3))
plot(p_moist ~ s_logmoist, newd, pch = 19, col = ring, ylim = c(-2.875544, 2.630993))
abline(mcoefs[1], mcoefs[3], col = 1, lwd = 2)


plot(p_temp ~ s_temp, newd, pch = 19, col = ring, ylim = c(-2.875544, 2.630993))
abline(mcoefs[1], mcoefs[4], col = 1, lwd = 2)


plot(p_par ~ s_logpar, newd, pch = 19, col = ring, ylim = c(-2.875544, 2.630993))
abline(mcoefs[1], mcoefs[5], col = 1, lwd = 2)


boxplot(p_co2 ~ co2, newd, ylim = c(-2.875544, 2.630993))


plot(s_c4_ddiff ~ p_moist, data = newd)
abline(0, 1)
plot(s_c4_ddiff ~ p_temp , data = newd)
abline(0, 1)
plot(s_c4_ddiff ~ p_par  , data = newd)
abline(0, 1)
range(c34sum$s_c4_ddiff)



m4_lm <- lm(s_c4_ddiff ~ co2 + s_logmoist + s_logpar + s_temp, data = c34sum)
coef(m4_lm)
par(mfrow = c(2, 2))
avPlot(m4_lm, variable = "s_logmoist", marginal.scale = TRUE, grid = FALSE, 
       xlab = "Adjusted soil moisture", ylab = expression(Adjusted~Delta*C[4]))
avPlot(m4_lm, variable = "s_temp", marginal.scale = TRUE, grid = FALSE, 
       xlab = "Adjusted tempearture", ylab = expression(Adjusted~Delta*C[4]))
avPlot(m4_lm, variable = "s_logpar", marginal.scale = TRUE, grid = FALSE, 
       xlab = "Adjusted understorey PAR", ylab = expression(Adjusted~Delta*C[4]))


