## Here, I check how purmute functions works when there is not nest-plot

# 6 rings -> 6! = 720
# 4 years (time series); 4 permutations for each ring -> 4^6
# total possible permutation = 6! * 4^6


## get site df for ring-summary df
prc_sp1   <- decostand(PlotSumVeg[, SppName], method = "log")
prc_sp2   <- prc_sp1 %>% 
  bind_cols(PlotSumVeg[, c("year", "ring", "co2")]) %>% 
  group_by(year, ring, co2) %>% 
  summarise_all(mean) %>% 
  ungroup()
prc_sitet <- prc_sp2 %>%
  select(year, ring, co2) %>% 
  arrange(ring, year)


## define permutation
cntrol_year <- how(within = Within(type = "series"),                        # year within ring is a time seiries
                   plot   = Plots(strata = prc_sitet$ring, type = "free"),  # rings are freely excanged
                   nperm = 4999)  
perm_yeart   <- shuffleSet(nrow(prc_sitet), control = cntrol_year)


## give individual ids for each of 720 combinations of 6 rings
perm6df <- data.frame(cbind(rbind(1:6, allPerms(6)), id = 1:720)) %>%
  transmute(id, com =  paste(V1, V2, V3, V4, V5, V6, sep = "_"))


d2 <-  data.frame(t(perm_yeart)) %>%                                      # transpose 
  mutate_all(funs(cut(., breaks = 6, labels = paste0("ring", 1:6)))) %>%  # give ring id (ring1 = 1:6, ring2 = 7:12 etc.)
  mutate_all(funs(as.numeric(.))) %>%                                     # turn all values into numeric
  distinct() %>%                                                          # remove duplicated rows (24 rows (6 rings x 4 years) -> 6 rows) to see the used combinations of 6 rings
  t() %>%                                                                 # back-transpose 
  data.frame() %>% 
  transmute(com = paste(X1, X2, X3, X4, X5, X6, sep = "_")) %>%           # get ids for combiation of rings
  left_join(perm6df, by = "com") %>%
  group_by(id) %>% 
  summarise(N = n())                                                      # get number of rows for each id to see how many times each id shows up


range(d2$N)  # min and max freaquency of ids
hist(d2$N)
plot(N ~ id, data = d2)  # all ids are fairly randomly selected



## check the permutation we use for prc analysis
nall <- nrow(alperms_bind)
d3 <-  data.frame(t(alperms_bind[sample(nall, 9999), ])) %>%              # transpose 
  mutate_all(funs(cut(., breaks = 6, labels = paste0("ring", 1:6)))) %>%  # give ring id (ring1 = 1:6, ring2 = 7:12 etc.)
  mutate_all(funs(as.numeric(.))) %>%                                     # turn all values into numeric
  distinct() %>%                                                          # remove duplicated rows (24 rows (6 rings x 4 years) -> 6 rows) to see the used combinations of 6 rings
  t() %>%                                                                 # back-transpose 
  data.frame() %>% 
  transmute(com = paste(X1, X2, X3, X4, X5, X6, sep = "_")) %>%           # get ids for combiation of rings
  left_join(perm6df, by = "com") %>%
  group_by(id) %>% 
  summarise(N = n())                                                      # get number of rows for each id to see how many times each id shows up
hist(d3$N)                                                                # note that if you increase the number of permutation, shapes of histgram start to vary every time
plot(N ~ id, data = d3)  # all ids are fairly randomly selected


## true random sampling
randomdf <- data.frame(id = sample(720, 9999, replace = TRUE)) %>% 
  group_by(id) %>% 
  summarise(N = n())
hist(randomdf$N)
plot(N ~ id, randomdf)

