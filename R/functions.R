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
PltVeg <- function(data, xval, xlab = NULL, ..., 
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
  p2 <- p + geom_bar(alpha = 0.4, position = position_dodge(width = .4)) + 
    theme(axis.text.y = element_text(...)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, ...)) +
    labs(x = xlab, y = "Frequency")
  return(p2)
}

#########################
# subset and droplevels #
#########################
subsetD <- function(...) droplevels(subset(...))

###########################
# step deletion with lmer #
###########################
stepLmer <- function(model, red.rndm = FALSE, ddf = "Kenward-Roger", ...){
  update(step(model, reduce.random = red.rndm, ddf = ddf,...)$model, 
         contrasts = NULL)
}
# use "Kenward-Roger" for approximation for denominator degrees of freedom. This
# is the same as the default DF given by Anova(model, test.statistic = "F). The
# default of step gives me a warning message for IEM-NO3 for some reasons (not
# sure why.. so changed it.)

###########################################
# produce box plots with transformed data #
###########################################
# log OR sqrt OR power(1/3) OR inverse OR box-cox
bxplts <- function(value, xval, ofst = 0, data, ...){
  data$y <- data[[value]] + ofst #ofst is added to make y >0
  data$xv <- data[[xval]]
  a <- boxcox(y ~ xv * year, data = data)
  par(mfrow = c(2, 3))
  boxplot(y ~ xv*year, data, main = "raw")
  boxplot(log(y) ~ xv*year, main = "log", data)
  boxplot(sqrt(y) ~ xv*year, main = "sqrt", data)
  boxplot(y^(1/3) ~ xv*year, main = "power(1/3)", data)
  boxplot(1/y ~ xv*year, main = "inverse", data)
  BCmax <- a$x[a$y == max(a$y)]
  texcol <- ifelse(BCmax < 0, "red", "black") 
  boxplot(y^(BCmax) ~ xv*year, 
          main = "", sep = "=", 
          data = data)
  title(main = paste("Box Cox", round(BCmax, 4)), 
        col.main = texcol)
  par(mfrow = c(1,1))
}

# multiple box-cox power plot for different constant values
bxcxplts <- function(value, xval, data, sval, fval){
  par.def <- par() # current graphic conditions
  data$yval <- data[[value]]
  data$xv <- data[[xval]]
  ranges <- seq(sval, fval, (fval - sval)/9)
  
  # store parameters given from box-cox plot
  BCmax <- vector()
  for (i in 1:10){
    data$y <- data$yval + ranges[i]
    a <- boxcox(y ~ xv * year, data = data)
    BCmax[i] <- a$x[a$y == max(a$y)]
  }
  
  # plot box plot with poer given from box-box for 
  # each contstant value
  par(mfrow = c(5, 2), omi = c(0, 0, 0, 0), mai = c(0.4, 0.4, 0.4, 0))
  sapply(1:10, function(x) {
    boxplot((yval + ranges[x]) ^ BCmax[x] ~ xval * year, 
            main = "", data = data)
    texcol <- ifelse(BCmax[x] < 0, "red", "black") 
    title(main = paste("constant=", round(ranges[x], 4), 
                       ", boxcox=", round(BCmax[x], 4)),
          col.main = texcol)
  })
  par(par.def) # set the graphic conditions back
}

####################################################
# function which reads worksheet from an xcel file #
####################################################
read.veg.xlx <- function(sheetName, file) {
  a <- read.xlsx2(file, sheetName,
                  header = TRUE, startRow = 4, endRow = 29, stringsAsFactors = FALSE)
  a <- a[ ,!grepl("X.", names(a))]
  a$position <- sheetName
  a[a == ""]  <- 0 # empty cell -> 0
  xlcFreeMemory() # Frees Java Virtual Machine (JVM) memory
  # it's not normally necessary to do this every time but the file size is 
  # really huge and cannot read all the worksheets at once so free memory every time 
  return(a)
}

####################################################
# organise species in a data frame read from excel #
####################################################
OrgSpp <- function(df, KeepCol, CombineCol, siteVec = c("ring", "plot", "position", "cell")){
  df[KeepCol] <- rowSums(df[CombineCol])
  RemoveCol <- CombineCol[CombineCol != KeepCol]
  df <- df[!names(df) %in% RemoveCol]
  # sort
  SpVec <- names(df)[!names(df) %in% siteVec]
  df <- df[, c(siteVec, sort(SpVec))]
  return(df)
}

##########################################
# compute canonical correlation from CAP #
##########################################
CanonicalCor <- function(CAPRes, EnvDF, term){
  m <- CAPRes$m
  xv <- EnvDF[, term]
  ml <- lm(CAPRes$PCoA[, 1:m] ~ xv)
  if (m != 1){ # you can't run candisc when m = 1
  cdaRes <- candisc(ml, term = "xv")
  cdaRes$canrsq
  } else {
    warning("m = 1")
    return(summary(ml)$r.squared)
  }
}
# also look at CAP_CanonicalCorrelation.pdf


################################
# Return star based on P value #
################################
FormatPval <- function(Pval) {
  stars <- ifelse(Pval > .1, "",
                  ifelse(Pval > .05, "scriptstyle('\u2020')",
                         ifelse(Pval > .01, "*",
                                ifelse(Pval > .001, "**",
                                       c("***")))))
  
  p <- as.character(ifelse(Pval > .1, round(Pval, 3),
                           ifelse(Pval < .001, "bold('<0.001')", 
                                  # shown with bold font. Note that inside of
                                  # bold needs to be in ''
                                  paste("bold(", round(Pval, 3), ")", sep = "'"))))
  return(data.frame(stars, p))
} 

####################################
# create table of contrast results #
####################################
cntrstTbl <- function(cntrstRes, data, variable, ...){
  Df <- data.frame(
    year = unique(data[, "year"]),
    contrast  =  cntrstRes$Contrast,
    SE = cntrstRes$SE,
    t = cntrstRes$testStat,
    df = cntrstRes$df,
    P.value = cntrstRes$Pvalue,
    FormatPval(cntrstRes$Pvalue),
    variable = variable)
  return(Df)
}

#######################
# Compute R2 for GLMM #
#######################
source("R/rsquaredglmm.R")


########################
# Yearly dissimilarity #
########################
# Compute dissimilarity for each plot between 2013 & 2014 and 2014 & 2015
YearDssmlrty <- function(x, spp) {
  df1 <- subset(x, year %in% c(2013, 2014))
  dis1 <- vegdist(log(df1[, spp] + 1), method = "bray")
  df2 <- subset(x, year %in% c(2014, 2015))
  dis2 <- vegdist(log(df2[, spp] + 1), method = "bray")
  dfs <- data.frame(year = c("Year1", "Year2"), BC = c(dis1, dis2))
  return(dfs)
}

#######################################
# Use special symbols with facet_wrap #
#######################################
facet_wrap_labeller <- function(gg.plot,labels=NULL) {
  #works with R 3.0.1 and ggplot2 0.9.3.1
  require(gridExtra)
  
  g <- ggplotGrob(gg.plot)
  gg <- g$grobs      
  strips <- grep("strip_t", names(gg))
  
  for(ii in seq_along(labels))  {
    modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
  }
  
  g$grobs <- gg
  class(g) = c("arrange", "ggplot",class(g)) 
  g
}

###################################
# df to create a circle on ggplot #
###################################
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}


#####################################################
# make df to set ranges of x and y axies for biplot #
#####################################################
rangeDF <- function(df, xv = "CAP1", yv = "CAP2", year = "2013"){
  ddply(df, .(co2), function(x) {
    xyrange <- c(range(x[, xv]), range(x[, yv]))
    MaxAbs <- max(abs(xyrange)) # highest absolute value
    xydf <- data.frame(CAP1 = c(MaxAbs, -MaxAbs), 
                       CAP2 = c(MaxAbs, -MaxAbs), 
                       year = year) 
    return(xydf)}
  )
}

#############################################
# Compute Species correlation with CAP axes #
#############################################
SpScorCorFun <- function(CapSpScorDF, CorVal = .7, co2) {
  # CorVal is threshhold for correlation
  # CapSpScorDF is species score from CAP analysis
  spdf <- data.frame(CAP1 = CapSpScorDF[, 1],
                     CAP2 = CapSpScorDF[, 2], 
                     Spp = row.names(CapSpScorDF), 
                     co2 = co2)
  spdf <- spdf[complete.cases(spdf), ]
  # set threshhold for correlation
  spdf <- subset(spdf, abs(CAP1) > CorVal|abs(CAP2) > CorVal)
  return(spdf)
}


###################################
# Create species correlation plot #
###################################
SpCorpPlot <- function(df, xv = "CAP1", yv = "CAP2", textpos = 1.07) {
  # textpos adjusts species names position in the grap
  df$xval <- df[, xv]
  df$yval <- df[, yv]
  
  df <- within(df, {
    # position for species names
    xval2 <- xval * textpos
    yval2 <- yval * textpos
  })
  
  p <- ggplot(data = df, aes(xval, yval))
  CorPl <- p + 
    geom_segment(aes(x = 0, y = 0, xend = xval, yend = yval), 
                 arrow = arrow(length = unit(.1, "cm")),
                 alpha = .7) +
    geom_path(data = circleFun(diameter = 2), aes(x, y), alpha = .7) +
    geom_text(data = df, 
              aes(x = xval2, y = yval2, label = Spp),
              size = 2, fontface = "bold.italic", 
              colour = "red", alpha = .7, lineheight = .7) +
    coord_fixed() +
    facet_wrap(~co2, scales = "free") +
    labs(x = "Correlation with CAP1", y = "Correlation with CAP2") +
    theme(plot.margin=unit(c(0, 0.5, 0, 0), "lines"),
          strip.background = element_blank(),
          strip.text.x = element_blank(), 
          axis.title.y = element_text(lineheight = 1.2))
  return(CorPl)
}

##########################
# Create CAP result plot #
##########################
CapPlot <- function(df){
  p <- ggplot(df, aes(x = CAP1, y = CAP2, fill = year, col = year))
  p2 <- p + geom_point(size = 3.5) + 
    geom_blank(aes(x = CAP1, y = CAP2), data = rangeDF(df)) +
    facet_wrap(~co2, scales = "free") +
    theme(legend.position = "top", 
          legend.direction = "horizontal",
          plot.margin=unit(c(0, 0.5, 0, 0.5), "lines"), 
          axis.title.y = element_text(lineheight = 1.2)) +
    coord_fixed()
  return(p2)
}

###############################
# Barplot for data inspection #
###############################
InspctPlot <- function(df = vdf, ringval, plotval, sp){
  cmdf <- subset(df, ring == ringval & plot == plotval, 
                 select =c ("year", "month", "ring", "plot", "position", "cell", sp))
  cmdf$ym <- cmdf$year:cmdf$month
  p <- ggplot(cmdf, aes_string(x = "position", y = sp))
  p2 <- p + geom_bar(stat = "identity") + facet_grid(. ~ ym)
  print(p2)
}

#############################################
# Remove year, ring and co2 columns from DF #
#############################################
Rm_ymc <- function(x) x[, !names(x) %in% c("year", "ring", "co2")]

###########
# Triplot #
###########
TriPlot <- function(MultValRes, env, yaxis, axispos, EnvNumeric = TRUE, lowx = .5, lowy = .5,
                    spcons = 2.5, biplcons = 3, centcons = 1){
  #   lowx, lowy are minimum sp correlation with axes to plot
  ResultScore <- summary(MultValRes)
  
  # Fitted or sp-weighted site score
  # sites <- if(FittedSite) data.frame(ResultScore$sites, env) else data.frame(ResultScore$constraints, env)
  
  # extract scores
  Rlist <- list(sites = data.frame(ResultScore$sites, env),
                species = data.frame(spp = gsub("[.]", "\n", row.names(ResultScore$species)), 
                                     ResultScore$species * spcons, 
                                     row.names = NULL), 
                biplot = data.frame(ResultScore$biplot * biplcons, 
                                    predictor = row.names(ResultScore$biplot),
                                    row.names = NULL))
  if(!is.na(ResultScore$centroids)) {
    Rlist$centroids <- data.frame(treatment = row.names(ResultScore$centroids),
                                 ResultScore$centroids * centcons,
                                 row.names = NULL)
  } 
  
  # remove categorical variables from biplot when there's numerc
  if(EnvNumeric){
    tdf <- Rlist$biplot
    Rlist$biplot <- tdf[!grepl("co2|year", as.character(tdf$predictor)), ]
  }
  
  # add ring and year columns for later
  for (i in 2:length(Rlist)) Rlist[[i]] <- data.frame(Rlist[[i]], ring = "1", year = "2013")
  axisnames <- colnames(ResultScore$cont$importance)[axispos] # axes to be used
  Rlist_mlt <- llply(Rlist, function(x) {
    mltid <- names(x)[!names(x) %in% axisnames[-1]]
    melt(x, id = mltid)
  })
  
  # % variance explained by axes
  VarProp <- ResultScore$cont$importance["Eigenvalue",]/ResultScore$tot.chi
  axislabs <- paste0(axisnames, "(", round(VarProp[axispos] * 100, 2), "%)")
  Rlist_mlt <- llply(Rlist_mlt, function(x) {
    names(x)[names(x) == axisnames[1]] <- "xval"
    x$variable <- factor(x$variable, levels = axisnames[-1], labels = axislabs[-1])
    return(x)
  })
  
  # make a plot
  theme_set(theme_bw())
  p <- ggplot(data = Rlist_mlt$sites, aes(x = xval, y = value, col = ring, shape = year))
  p2 <- p + geom_point(size = 3) + facet_grid(variable ~ .)
  
  # treatment
  if(!is.na(ResultScore$centroids)){
  p2 <- p2 + 
    geom_segment(data = Rlist_mlt$centroids, 
                 aes(x = 0, y = 0, xend = xval, yend = value), 
                 arrow = arrow(length = unit(.1, "cm")), 
                 alpha = .6,
                 color = "black") + 
    geom_text(data = Rlist_mlt$centroids, 
              aes(x = xval * 1.1 , y = value * 1.1, label = treatment), 
              alpha = .6, lineheight = .7, 
              color = "black", size = 2, 
              fontface = "bold") }
  
  # species
  spdf <- subset(Rlist_mlt$species, abs(xval) > lowx * 2.5| abs(value) > lowy * 2.5)
  p4 <- p2 + 
    geom_segment(data = spdf, 
                 aes(x = 0, y = 0, xend = xval, yend = value), 
                 arrow = arrow(length = unit(.1, "cm")),
                 alpha = .6) +
    geom_text(data = spdf, 
              aes(x = xval * 1.2, y = value * 1.2, label = spp),
              size = 2, fontface = "bold.italic", 
              colour = "red", alpha = .6, lineheight = .7) +
    labs(x = axislabs[1], y = yaxis)
  
  # environmental variables
  if(EnvNumeric){
    p4 <- p4 + 
      geom_segment(data = Rlist_mlt$biplot, 
                   aes(x = 0, y = 0, xend = xval, yend = value), 
                   arrow = arrow(length = unit(.1, "cm")), 
                   alpha = .6,
                   color = "blue") + 
      geom_text(data = Rlist_mlt$biplot, 
                aes(x = xval * 1.1, y = value * 1.1, label = predictor), 
                alpha = .6, lineheight = .7, 
                color = "blue", size = 2, 
                fontface = "bold")
  }
  p5 <- p4 + 
    scale_color_hue(labels = paste0(1:6, c("e", "a", "a", "e", "e", "a"))) +
    geom_hline(aes(yintercept = 0), linetype = "dotted") +
    geom_vline(aes(xintercept = 0), linetype = "dotted")
  p5
}
