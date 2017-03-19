data("cement")
fm <- lm(y ~ X1 + X2 + X3 + X4, data = Cement, na.action = na.fail)
ms <- dredge(fm)

# filter out overly complex models according to the 
# "nesting selection rule":
subset(ms, !nested(.)) # dot represents the ms table object
subset(c4d_m6_full, !nested(.)) # dot represents the ms table object

newmod <- model.avg(get.models(subset(c4d_m6_full, !nested(.)), subset = delta < 2))
summary(newmod)
confint(newmod, level = .90, method = "boot", nsim = 50)
?scale
bm <- get.models(c4d_m6_full, subset = 1)[[1]]
confint(bm, level = .90, method = "boot", nsim = 500, boot.type = "basic")

Anova(get.models(c4d_m6_full, subset = 1)[[1]], test.statistic = "F")
?confint
!nested(c4d_m6_full)
# print model "4" and all models nested within it
nst <- nested(c4d_m6_full, indices = "row")
c4d_m6_full[c("15", nst[["15"]])]

c4d_m6_full$nested <- sapply(nst, paste, collapse = ",")

ms
