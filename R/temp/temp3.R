ppp <- predict(m3ave, newdf, se.fit = TRUE, re.form = NA)
ppp <- cbind(newdf, ppp)
p1 <- ggplot(data = ppp, aes(x = c4moist, y = fit, col = co2)) +
  geom_line()+
  geom_line(aes(y = fit + se.fit * 1.96), linetype = "dashed")+
  geom_line(aes(y = fit -se.fit * 1.96), linetype = "dashed")+
  geom_point(data = c34sum, aes(y = c4_delta, col = co2))+
  geom_hline(yintercept = 0)


bb <- bootMer(m3, 
              FUN=function(x) predict(x, newdf, re.form = NA),
              nsim=999)
lci <- apply(bb$t, 2, quantile, 0.025)
uci <- apply(bb$t, 2, quantile, 0.975)
PredVal <- bb$t0
df <- cbind(lci, uci, PredVal, newdf)


names(df)
p2 <- ggplot(df, aes(x = c4moist, y = PredVal, col = co2))+
  geom_line(aes(y = lci), linetype = "dashed") +
  geom_line(aes(y = uci), linetype = "dashed") +
  geom_point(data = c34sum, aes(y = c4_delta)) +
  geom_hline(yintercept = 0)


pm <- cbind(ggplotGrob(p1), ggplotGrob(p2))
grid.newpage()
grid.draw(pm)
