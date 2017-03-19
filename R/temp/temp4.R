  select(ring, id, RY) %>% 
  ungroup() %>% 
  distinct()
moistval <- seq(min(c34sum$totalmoist), max(c34sum$totalmoist), length.out = 1000)
tempval1  <- quantile(c34sum$annual_temp2m)[3]


newdf_c3    <- ldply(1:10, function(x) 
  cbind(sitedf, 
        totalmoist  = moistval[sample(1000, nrow(sitedf), replace = TRUE)])) %>% 
  mutate(annual_temp2m = tempval1)

c3d_m3_pred <- predict(c3d_m3_avg, newdf_c3, se.fit = TRUE, re.form = NA)
c3d_m3_pred_df <- cbind(c3d_m3_pred, newdf_c3) %>% 
  mutate(lwr = fit - se.fit * 1.96,
         upr = fit + se.fit * 1.96) 
ggplot(c3d_m3_pred_df, aes(x = totalmoist, y = fit)) +
  geom_line()+
  geom_line(aes(y = lwr), linetype = "dashed") +
  geom_line(aes(y = upr), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")



moistval <- quantile(c34sum$totalmoist)[3]
tempval1 <- seq(min(c34sum$annual_temp2m), max(c34sum$annual_temp2m), length.out = 1000) 


newdf_c3_temp    <- ldply(1:10, function(x) 
  cbind(sitedf, 
        annual_temp2m  = tempval1[sample(1000, nrow(sitedf), replace = TRUE)])) %>% 
  mutate(totalmoist = moistval)

c3d_m3_pred <- predict(c3d_m3_avg, newdf_c3_temp, se.fit = TRUE, re.form = NA)
c3d_m3_pred_df <- cbind(c3d_m3_pred, newdf_c3_temp) %>% 
  mutate(lwr = fit - se.fit * 1.96,
         upr = fit + se.fit * 1.96) 

ggplot(c3d_m3_pred_df, aes(x = annual_temp2m, y = fit)) +
  geom_line()+
  geom_line(aes(y = lwr), linetype = "dashed") +
  geom_line(aes(y = upr), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")
