names(RingSumVeg)


# Year3, all species ------------------------------------------------------


d <- filter(RingSumVeg, year == "Year3") # spp df
e <- filter(EnvDF_3df, year == "Year3")

# remove non-observed species
sppsum <- colSums(d[, SppName])   # sp sum
sp_to_use <- names(sppsum)[sppsum != 0]  # spp to be used
names(sppsum)[sppsum == 0]
length(sp_to_use)
length(SppName)
new_d <- select(d, one_of(sp_to_use, "ring"))
new_d[, sp_to_use] <- decostand(new_d[, sp_to_use], method = "hellinger")  # transform using hellinger
new_d <- merge(new_d, e, by = "ring")

# rda
r_co2         <- rda(new_d[sp_to_use] ~ co2, data = new_d)
r_Depth_HL    <- rda(new_d[sp_to_use] ~ Depth_HL, data = new_d)
r_TotalC      <- rda(new_d[sp_to_use] ~ TotalC, data = new_d)
r_Drysoil_ph  <- rda(new_d[sp_to_use] ~ Drysoil_ph, data = new_d)
r_moist       <- rda(new_d[sp_to_use] ~ moist, data = new_d)
r_temp        <- rda(new_d[sp_to_use] ~ temp, data = new_d)
r_gapfraction <- rda(new_d[sp_to_use] ~ gapfraction, data = new_d)

# adjusted R2
llply(list(r_co2, r_TotalC, r_moist, r_Drysoil_ph, r_Depth_HL, r_gapfraction, r_temp),
      function(x) RsquareAdj(x)$adj.r.squared)

# full models
r_full <- rda(new_d[, sp_to_use] ~ TotalC + moist + Depth_HL, data = new_d)
anova(r_full, permutations = allPerms(6))
RsquareAdj(r_full)$adj.r.squared

# simplificaiton
rr2 <- rda(new_d[, sp_to_use] ~ 1, data = new_d)                                                  
rr3 <- ordiR2step(rr2, r_full, permutations = allPerms(6), direction = "forward", Pin = .1)  
anova(rr3, permutations = allPerms(6))
RsquareAdj(rr3)$adj.r.squared

# get F and P value for each term in the final model (rr3)
anova_rr3 <- anova(rr3, permutations = allPerms(6), by = "margin")

# Year0, forb species ------------------------------------------------------


d <- RingSumVeg %>% 
  filter(year == "Year0") %>% 
  select(one_of(SppName_forb, "ring"))
names(d)
e <- filter(EnvDF_3df, year == "Year0")

# remove non-observed species
sppsum <- colSums(d[, SppName_forb])   # sp sum
sp_to_use <- names(sppsum)[sppsum != 0]  # spp to be used
names(sppsum)[sppsum == 0]
length(sp_to_use)
length(SppName_forb)
new_d <- select(d, one_of(sp_to_use, "ring"))
new_d[, sp_to_use] <- decostand(new_d[, sp_to_use], method = "hellinger")  # transform using hellinger
new_d <- merge(new_d, e, by = "ring")

# rda
r_co2         <- rda(new_d[sp_to_use] ~ co2, data = new_d)
r_Depth_HL    <- rda(new_d[sp_to_use] ~ Depth_HL, data = new_d)
r_TotalC      <- rda(new_d[sp_to_use] ~ TotalC, data = new_d)
r_Drysoil_ph  <- rda(new_d[sp_to_use] ~ Drysoil_ph, data = new_d)
r_moist       <- rda(new_d[sp_to_use] ~ moist, data = new_d)
r_temp        <- rda(new_d[sp_to_use] ~ temp, data = new_d)
r_gapfraction <- rda(new_d[sp_to_use] ~ gapfraction, data = new_d)

# adjusted R2
llply(list(r_co2, r_TotalC, r_moist, r_Drysoil_ph, r_Depth_HL, r_gapfraction, r_temp),
      function(x) RsquareAdj(x)$adj.r.squared)

# full models
pos_tems <- c("TotalC", "moist", "Drysoil_ph", "Depth_HL", "gapfraction", "temp")

# 4 terms
t4 <- apply(combn(pos_tems, 4), 2, function(x) paste(x, collapse = "+"))
t4 <- paste("new_d[, sp_to_use] ~", t4)
  
r_t4 <- list()
for (i in 1:15){
  r_t4[[i]] <- rda(as.formula(t4[i]), data = new_d)
}
llply(r_t4, function(x) anova(x, permutations = allPerms(6)))
llply(r_t4, function(x) vif.cca(x))
filter(rda_summary, grepl(pattern = "^term_n4.forb", .id))

# 3 trems
t3 <- apply(combn(pos_tems, 3), 2, function(x) paste(x, collapse = "+"))
t3 <- paste("new_d[, sp_to_use] ~", t3)

r_t3 <- list()
for (i in 1:20){
  r_t3[[i]] <- rda(as.formula(t3[i]), data = new_d)
}
llply(r_t3, function(x) anova(x, permutations = allPerms(6)))
resu <- ldply(r_t3, function(x) {
  v <- vif.cca(x)
  vdfd <- data.frame(as.list(v))
  adjr <- RsquareAdj(x)$adj.r.squared
  
  data.frame(vdfd, adjr)
})
arrange(resu, -adjr)

r_full <- rda(new_d[, sp_to_use] ~ TotalC + moist + Drysoil_ph, data = new_d)
anova(r_full, permutations = allPerms(6))
RsquareAdj(r_full)$adj.r.squared

# simplificaiton
rr2 <- rda(new_d[, sp_to_use] ~ 1, data = new_d)                                                  
rr3 <- ordiR2step(rr2, r_full, permutations = allPerms(6), direction = "forward", Pin = .1)  
anova(rr3, permutations = allPerms(6))
RsquareAdj(rr3)$adj.r.squared

# get F and P value for each term in the final model (rr3)
anova(rr3, permutations = allPerms(6), by = "margin")





# Year0 grass species -----------------------------------------------------

d <- filter(RingSumVeg, year == "Year0")
e <- filter(EnvDF_3df, year == "Year0")

# remove non-observed species
sppsum <- colSums(d[, SppName_grass])   # sp sum
sp_to_use <- names(sppsum)[sppsum != 0]  # spp to be used
names(sppsum)[sppsum == 0]
length(sp_to_use)
length(SppName)
new_d <- select(d, one_of(sp_to_use, "ring"))
new_d[, sp_to_use] <- decostand(new_d[, sp_to_use], method = "hellinger")  # transform using hellinger
new_d <- merge(new_d, e, by = "ring")

# rda
r_co2         <- rda(new_d[sp_to_use] ~ co2, data = new_d)
r_Depth_HL    <- rda(new_d[sp_to_use] ~ Depth_HL, data = new_d)
r_TotalC      <- rda(new_d[sp_to_use] ~ TotalC, data = new_d)
r_Drysoil_ph  <- rda(new_d[sp_to_use] ~ Drysoil_ph, data = new_d)
r_moist       <- rda(new_d[sp_to_use] ~ moist, data = new_d)
r_temp        <- rda(new_d[sp_to_use] ~ temp, data = new_d)
r_gapfraction <- rda(new_d[sp_to_use] ~ gapfraction, data = new_d)

# adjusted R2
llply(list(r_co2, r_TotalC, r_moist, r_Drysoil_ph, r_Depth_HL, r_gapfraction, r_temp),
      function(x) RsquareAdj(x)$adj.r.squared)


# full models
r_full <- rda(new_d[, sp_to_use] ~ moist + Drysoil_ph, data = new_d)
anova(r_full, permutations = allPerms(6))

# not significant
RsquareAdj(r_Drysoil_ph)$adj.r.squared
RsquareAdj(r_moist)$adj.r.squared

anova(r_Drysoil_ph, permutations = allPerms(6))

# simplificaiton
rr2 <- rda(new_d[, sp_to_use] ~ 1, data = new_d)                                                  
rr3 <- ordiR2step(rr2, r_Drysoil_ph, permutations = allPerms(6), direction = "forward", Pin = .1)  
anova(rr3, permutations = allPerms(6))
RsquareAdj(rr3)$adj.r.squared
