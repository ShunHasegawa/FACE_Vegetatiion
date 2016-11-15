library(plyr)
library(dplyr)
library(permute)

testd <- expand.grid(ring = paste0("ring", 1:6), 
                     plot = paste0("plot", 1:4),
                     year = paste0("year", 1:3)) %>% 
  mutate(id = ring:plot)
# %>% 
  # arrange(ring, plot, year)
testd
cc <- how(within = Within(type = "series", constant = TRUE),
          plot = Plots(strata = testd$id, type = "none"))
pp <- shuffleSet(nrow(testd), control = cc)

pp
ppl <- llply(levels(testd$ring), function(x) pp[, testd$ring == x])
perm6 <- rbind(1:6, allPerms(6))
tl <- alply(perm6, 1, function(x) do.call(cbind, ppl[x]))
tlb <- do.call(rbind, tl)

llply(split(tlb[1, ], testd$ring), matrix, ncol = 4)
llply(split(tlb[10, ], testd$ring), matrix, ncol = 4)
alply(tlb[sample(700, size = 5), ], 1, function(x) {
  llply(split(x, testd$ring), matrix, ncol = 4)
})
