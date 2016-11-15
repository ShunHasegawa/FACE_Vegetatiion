require(permute)
require(plyr)
require(data.table)

## define permutation for this experiment; year is nested within plot within
## ring.

## 1. plot is a random factor, should never be exchanged. 
## 2. year is a time series measurement. This can be exchanged, but its order 
## should be kept to take autocorrelation into account (i.e. [1, 2, 3], [2, 3, 
## 1], [3, 1, 2]). Also this order should be consistent among all plots.
## 3. ring is an experimental unit so fully exchangable

load("output/Data/prc_site.RData")

## 1) define permutation for year within each plot

cntrol_year <- how(within = Within(type = "series"),                     # year within plot is a time seiries
                   plot   = Plots(strata = prc_site$id, type = "none"),  # plot should not be exchanged
                   nperm = 49999)  
perm_year   <- shuffleSet(nrow(prc_site), control = cntrol_year)
nrow(perm_year)

## 2) the above permutation did not exchange between rings. do so as below

perm_year_or <- rbind(1:nrow(prc_site), perm_year)                               # add the original (i.e. 1:96) as this should be also permuted between rings
perm_ring_l  <- llply(levels(prc_site$ring), function(x) perm_year_or[, prc_site$ring == x])  # split each row of the permutation by ring into a list
perm6        <- rbind(1:6, allPerms(6))                                          # potential permutation for 6 rings (i.e. 6!)
perm_by_ring <- alply(perm6, 1, function(x) do.call(cbind, perm_ring_l[x]))      # permute by ring, exchanging perm_year_or between rings (e.g. origin[1, 2, 3, 4, 5, 6], 2nd[1, 2, 3, 4, 6, 5] etc.)
comp_perm    <- do.call(rbind, perm_by_ring)[-1, ]                               # row bind all results and remove the original
nrow(comp_perm)

fwrite(as.data.table(comp_perm), "output/Data/permu_matrix_prc.csv", 
       showProgress = FALSE)
