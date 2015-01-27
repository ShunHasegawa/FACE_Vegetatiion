#############################
# create ellipses on ggplot #
#############################
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
#       require(MASS)
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
  ass <- apply(cbind(ca.spp[, xv], ca.spp[, yv]), 1, function(x) max(abs(x))) / 
    max(c(abs(ca.spp[, xv]), abs(ca.spp[, yv])))
  
  pl.ttl <- paste(unique(ca.site$year), collapse = "&")
  
  # df to draw an outline for each site
  df <- ddply(ca.site[c("year", "ring", xv, yv)], .(year, ring),
              function(x) {chx <- chull(x[c(xv, yv)]) 
                           chxDF <- data.frame(rbind(x[chx,], x[chx[1], ]))
                           return(chxDF)})
  
  p <- ggplot(ca.site, aes_string(x=xv, y=yv))
  # aes_string takes character vector
  p2 <- p + geom_point(data = ca.site, size = 2, alpha = 0.8, aes_string(col = "ring", ...)) +
    geom_polygon(data = df, aes_string(x = xv, y = yv, fill = "ring", 
                                       col = "ring", ..., alpha = .01)) +  
    geom_text(data = ca.spp, size = 2, aes_string(x = xv, y = yv), alpha = ass, 
              label = rownames(ca.spp)) +
    ggtitle(pl.ttl)
  figtitle <- paste("output/figs/FACE.veg.", pl.ttl, "_", xv, "vs", yv, sep = "" )
  ggsavePP(filename = figtitle, plot = p2, width = 9, height = 6)
}


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

####################
# Create bargraphs #
####################
PltVeg <- function(data = veg, 
                   xval, 
                   xlab = NULL, 
                   ..., 
                   pfgLabs = c("C[3]", "C[3-4]", "C[4]", "Legume", "Lichen", "Moss", "Non_legume", "wood"),
                   orgnLabs = c("Native", "Introduced")){
  # change factor lablles for labbeling in figs
  data$co2 <- factor(data$co2, levels = c("amb", "elev"), labels = c("Ambient", "eCO[2]"))
  
  data$PFG <- factor(data$PFG, 
                     levels = c("c3", "c3_4", "c4", "legume", "Lichen", "moss", "Non_legume", "wood"),
                     labels = pfgLabs)
  data$origin <- factor(data$origin, 
                        levels = c("native", "naturalised"), 
                        labels = orgnLabs)
  data$xv <- data[, xval]
  p <- ggplot(data, aes(x = xv, fill = year))
  p2 <- p + geom_bar(alpha = 0.6, position = "identity") + 
    theme(axis.text.y = element_text(...)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, ...)) +
    labs(x = xlab, y = "Frequency")
  return(p2)
}

#########################
# subset and droplevels #
#########################
subsetD <- function(...) droplevels(subset(...))
