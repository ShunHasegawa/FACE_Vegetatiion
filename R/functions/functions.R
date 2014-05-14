#############################
# create ellipses on ggplot #
#############################
require(proto)

StatEllipse <- proto(ggplot2:::Stat,
{
  required_aes <- c("x", "y")
  default_geom <- function(.) GeomPath
  objname <- "ellipse"
  
  calculate_groups <- function(., data, scales, ...){
    .super$calculate_groups(., data, scales,...)
  }
  calculate <- function(., data, scales, level = 0.75, segments = 51,...){
    dfn <- 2
    dfd <- length(data$x) - 1
    if (dfd < 3){
      ellipse <- rbind(c(NA,NA))	
    } else {
      require(MASS)
      v <- cov.trob(cbind(data$x, data$y))
      shape <- v$cov
      center <- v$center
      radius <- sqrt(dfn * qf(level, dfn, dfd))
      angles <- (0:segments) * 2 * pi/segments
      unit.circle <- cbind(cos(angles), sin(angles))
      ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
    }
    
    ellipse <- as.data.frame(ellipse)
    colnames(ellipse) <- c("x","y")
    return(ellipse)
  }
}
)

stat_ellipse <- function(mapping=NULL, data=NULL, geom="path", position="identity", ...) {
  StatEllipse$new(mapping=mapping, data=data, geom=geom, position=position, ...)
}

###################
# plot CA results #
###################
pltCA <- function(data, xv, yv, ...){
  # CA
  vg.data <- data[, -grep("year|ring|plot|position|cell", names(data))]
  sites <- data[, grep("year|ring|plot|position|cell", names(data))]
  
  m1 <- cca(vg.data)
  m1.sm <- summary(m1)
  ca.spp <- data.frame(m1.sm$species)
  ca.site <- data.frame(sites, m1.sm$sites)
  
  # plot spp score: higher sp score gets darker color
  ass <- (abs(ca.spp[, xv]) + abs(ca.spp[, yv])) / max((abs(ca.spp[, xv]) + abs(ca.spp[, yv])))
  
  pl.ttl <- paste(unique(ca.site$year), collapse = "&")
    
  p <- ggplot(ca.site, aes_string(x=xv, y=yv))
  # aes_string takes character vector
  p2 <- p + geom_point(data = ca.site, size = 5, alpha = 0.8, aes_string(col = "ring", ...)) +
    stat_ellipse(data = ca.site, geom = "polygon", alpha = 0.1, level=0.75, 
                 aes_string(col = "ring", fill = "ring", ...)) +
    geom_text(data = ca.spp, size = 2, aes_string(x = xv, y = yv), alpha = ass, label = rownames(ca.spp)) +
    ggtitle(pl.ttl)
  figtitle <- paste("output/figs/FACE.veg.", pl.ttl, "_", xv, "vs", yv, ".pdf", sep = "" )
  ggsave(filename = figtitle, plot = p2, width = 9, height = 6)
}

#pltCA <- function(data, xv, yv, ...){
#  p <- ggplot(data, aes_string(x=xv, y=yv, col = "ring", ...))
#  #aes_string takes character vector
#  pl.ttl <- paste(unique(data$year), collapse = "&")
#  p2 <- p + geom_point(size = 5, alpha = 0.8) + 
#    stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = ring), level=0.75) +
#    ggtitle(pl.ttl)
#  figtitle <- paste("output/figs/FACE.veg.", pl.ttl, "_", xv, "vs", yv, ".pdf", sep = "" )
#  ggsave(filename = figtitle, plot = p2)
#}

# analyse for each year separately and plot
plt.CA.yr <- function(data){
  pltCA(data, xv = "CA1", yv = "CA2")
  pltCA(data, xv = "CA2", yv = "CA3")
  pltCA(data, xv = "CA1", yv = "CA3")
}

##############################
# Save ggplot in PDF and PNG #
##############################
ggsavePP <- function(filename, plot, width, height){
  ggsave(filename = paste(filename, ".pdf", sep = ""), 
         plot = plot, 
         width = width, 
         height = height)
  
  ggsave(filename = paste(filename, ".png", sep = ""), 
         plot = plot, 
         width = width, 
         height = height, 
         dpi = 600)
}
