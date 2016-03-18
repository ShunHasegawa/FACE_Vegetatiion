RngSppDF <- RingSumVeg[, SppName]
RngSiteDF <- RingSumVeg[, !names(RingSumVeg) %in% SppName]


PCoA <- cmdscale(d = vegdist(log(RngSppDF + 1), method = "bray"), eig = TRUE, k = 3)

PCoA_SiteScoreDF <- cbind(RngSiteDF, PCoA$points)
names(PCoA_SiteScoreDF)[5:7] <- paste0("PCoA", 1:3)

PCoA_SiteScoreDF_mlt <- melt(PCoA_SiteScoreDF, 
                             id = c("year", "ring", "block", "co2", "PCoA1"))

p <- ggplot(PCoA_SiteScoreDF_mlt, aes(x = PCoA1, y = value, 
                                      col = ring, shape = year, group = ring))
p2 <- p + 
  geom_point(size = 3) +
  geom_path() +
  facet_grid(. ~ variable)
ggsave(p2, filename = "output/figs//FACE_Vegetation_PCoA.pdf", width = 6, height = 4)

par(mfrow = c(1, 2))
plot(PCoA2 ~ PCoA1, data = PCoA_SiteScoreDF, type = "n")
d_ply(PCoA_SiteScoreDF, .(year, ring), function(x) {
  points(PCoA2 ~ PCoA1, col = ring, pch = as.numeric(year), data = x, cex = 2)
})

plot(PCoA3 ~ PCoA1, data = PCoA_SiteScoreDF, type = "n")
d_ply(PCoA_SiteScoreDF, .(year, ring), function(x) {
  points(PCoA3 ~ PCoA1, col = ring, pch = as.numeric(year), data = x, cex = 2)
})


