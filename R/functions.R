# #############################
# # create ellipses on ggplot #
# #############################
# StatEllipse <- proto(ggplot2:::Stat,
# {
#   required_aes <- c("x", "y")
#   default_geom <- function(.) GeomPath
#   objname <- "ellipse"
#   
#   calculate_groups <- function(., data, scales, ...){
#     .super$calculate_groups(., data, scales,...)
#   }
#   calculate <- function(., data, scales, level = 0.75, segments = 51,...){
#     dfn <- 2
#     dfd <- length(data$x) - 1
#     if (dfd < 3){
#       ellipse <- rbind(c(NA,NA))	
#     } else {
# #       require(MASS)
#       v <- cov.trob(cbind(data$x, data$y))
#       shape <- v$cov
#       center <- v$center
#       radius <- sqrt(dfn * qf(level, dfn, dfd))
#       angles <- (0:segments) * 2 * pi/segments
#       unit.circle <- cbind(cos(angles), sin(angles))
#       ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
#     }
#     
#     ellipse <- as.data.frame(ellipse)
#     colnames(ellipse) <- c("x","y")
#     return(ellipse)
#   }
# }
# )
# 
# stat_ellipse <- function(mapping=NULL, data=NULL, geom="path", position="identity", ...) {
#   StatEllipse$new(mapping=mapping, data=data, geom=geom, position=position, ...)
# }

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
                   pfgLabs = c("C[3]", "C[4]", "Legume", "Moss", "Non_legume", "wood"),
                   orgnLabs = c("Native", "Introduced")){
  # change factor lablles for labbeling in figs
  data$co2 <- factor(data$co2, levels = c("amb", "elev"), labels = c("Ambient", "eCO[2]"))
  
  data$PFG <- factor(data$PFG, 
                     levels = c("c3", "c4", "legume", "moss", "Non_legume", "wood"),
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
  data$y <- data[[value]] + ofst # ofst is added to make y >0
  data$xv <- data[[xval]]
  a <- boxcox(y ~ xv * year, data = data, plotit = FALSE)
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
  df1 <- subset(x, year %in% c("Year0", "Year1"))
  dis1_bc <- vegdist(log(df1[, spp] + 1), method = "bray")
  dis1_eu <- vegdist(log(df1[, spp] + 1), method = "euclidean")
  df2 <- subset(x, year %in% c("Year1", "Year2"))
  dis2_bc <- vegdist(log(df2[, spp] + 1), method = "bray")
  dis2_eu <- vegdist(log(df2[, spp] + 1), method = "euclidean")
  df3 <- subset(x, year %in% c("Year2", "Year3"))
  dis3_bc <- vegdist(log(df3[, spp] + 1), method = "bray")
  dis3_eu <- vegdist(log(df3[, spp] + 1), method = "euclidean")
  dfs <- data.frame(year = c("Year0-1", "Year1-2", "Year2-3"), 
                    BC = c(dis1_bc, dis2_bc, dis3_bc), 
                    EU  = c(dis1_eu, dis2_eu, dis3_eu))
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
  for (i in 2:length(Rlist)) Rlist[[i]] <- data.frame(Rlist[[i]], ring = "1", year = "Year1")
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

#######################################
# Plot RDA result against year by co2 #
#######################################
PlotRDA_Year <- function(rdaResLst, spscore = .3, env){
  names(rdaResLst) <- c("amb", "elev")
  
  # % variance for RDA1
  Rda1Prop <- laply(rdaResLst, function(x) {
    ss <- round(summary(x)$cont$importance["Proportion Explained", "RDA1"] * 100, 2)
    paste0(ss, "%")}
  )
  
  # Spp score
  RDAsppDF <- ldply(c("amb", "elev"), function(x) {
    ll <- rdaResLst[[x]]
    spdf <- vegan::scores(ll)$species
    data.frame(spdf, 
               co2 = x,
               year = "Species score",
               sp = row.names(spdf))
  })
  
  RDAsppDF$co2 <- factor(RDAsppDF$co2, labels = paste0(c("Ambient (", "eCO2 ("), Rda1Prop, ")"))
  RDAsppDF <- subset(RDAsppDF, abs(RDA1) > spscore)
  
  RDAsiteDF <- ldply(c("amb", "elev"), function(x) {
    ll <- rdaResLst[[x]]
    if(x == "amb") dd <- env[[1]] else dd <- env[[2]]
    data.frame(vegan::scores(ll)$sites, co2 = x, dd)
  }) 
  RDAsiteDF$co2 <- factor(RDAsiteDF$co2, labels = paste0(c("Ambient (", "eCO2 ("), Rda1Prop, ")"))
  
  p <- ggplot(RDAsiteDF, aes(x = year, y = RDA1))
  p2 <- p + 
    geom_point(size = 4, alpha = .7) + 
    geom_text(data = RDAsppDF, aes(x = year, y = RDA1 * 3, label = sp), size = 2) + 
    geom_hline(xintercept = 0, linetype = "dashed") + 
    facet_grid(. ~ co2)
  p2
}

###############################
# Compare AIC and R2 for glmm #
###############################
CompAIC <- function(model){
  m1 <- update(model, ~ . - year:co2)
  m2 <- update(m1, ~ . - co2)
  m3 <- update(m1, ~ . - year)
  res <- ldply(list(model, m2, m3, m1), r.squared)
  # delta AIC, for co2 and year, calculated from the 2nd model.
  res$dAIC <- with(res, c(0, AIC[c(2, 3)] - AIC[4], AIC[4] - AIC[1]))
  row.names(res) <- c("Full", "co2", "year", "year:co2")
  return(res)
}

############################
# Add labels in facet_wrap #
############################
facet_wrap_labeller <- function(gg.plot,labels=NULL) {
  #works with R 3.0.1 and ggplot2 0.9.3.1
  # copied from http://stackoverflow.com/questions/19282897/
  # how-to-add-expressions-to-labels-in-facet-wrap
  # require(gridExtra)
  
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

# save png file with 600 dpi 

save_png600 <- function(...) png(..., res = 600, units = "in")


# transform P values to star maks
get_star <- function(pval, dagger = TRUE){
  
  # dagger mark for p < 0.1
  dg <- ifelse(dagger, expression("\u2020"), ".")
  
  cut(pval, right = FALSE,
      breaks = c(0, .1, .05, .01, .001, 1),  
      labels = c("***", "**", "*", dg, ""))
}


# ggplot theme setting ----------------------------------------------------

theme_set(theme_bw() + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank()))

# define graphic background
science_theme <- theme(panel.border      = element_rect(color = "black"),
                       panel.grid.major  = element_blank(), 
                       panel.grid.minor  = element_blank(), 
                       legend.position   = c(.91, .91),
                       legend.title      = element_blank(),
                       legend.background = element_blank(),
                       legend.key        = element_blank(),
                       legend.key.width  = unit(2.5, "lines"),
                       legend.key.height = unit(.8, "lines"),
                       axis.ticks.length = unit(-.2, "lines"),
                       axis.text.x       = element_text(margin = margin(5)),
                       axis.text.y       = element_text(margin = margin(0, 5)),
                       axis.title.y      = element_text(margin = margin(0, 10)))

# correct year0 vaalue; this remove species that were observed only in Year0 yet
# not in the subseqent years
correct_year0 <- function(x, spp){
  d <- x
  
  # column sum for each of spp
  d_sum <- d %>% 
    group_by(year) %>% 
    summarise_each_(funs(sum), spp)
  
  # identify species which was observed only in Year0 in this plot
  a <- names(which(apply(d_sum[, -1], 2, function(x) all(x[2:4] == 0) & x[1] != 0)))
  
  # replace their values with 0
  d[, a] <- 0
 
  return(d)
}

# save as .RData and .csv files
save_Rdata_csv <- function(filename, ...){
  save(..., file = paste0(filename, ".RData"))
  write.csv(..., file = paste0(filename, ".csv"), row.names = FALSE)
}


# functions for RDA -------------------------------------------------------


is_exist_spdd <- function(ignoredd = FALSE){
  
  # this function checks if dd and spdd exsit as they will be overwritten when
  # performing rda using the functions below
  
  if(ignoredd) {
    
    
    # return warning message if they exist but to be ignored
    if(any(exists("dd"), exists("spdd")))
      warning("dd or spdd have been overwritten...")  
    
    
  } else {
    
    
    # stop if they aren't ignored
    if(any(exists("dd"), exists("spdd")))
      stop("dd and spdd exist...\nuse ignoredd = TRUE to overwrite")   
    
    
  }
}




get_adjR_singl <- function(x,                # df containing environmental variables and site 
                           ignoredd = FALSE, # ignore dd and spdd
                           expl,             # encironmental variables
                           SiteName_rda      # site
) {
  
  
  # check object names
  is_exist_spdd(ignoredd = ignoredd)
  
  
  # prepare data frames (dd and spdd)
  dd <<-  x                                       # df for environmental variables
  spdd <<- select(dd, -one_of(expl, SiteName_rda)) # df for spp
  
  # formula for each variable
  singl_fmls        <- llply(paste("spdd ~", expl), as.formula)
  names(singl_fmls) <- expl
  
  # run rda and compute R2
  r2 <- ldply(singl_fmls, 
              function(y) {
                adjR <- RsquareAdj(rda(y, data = dd))$adj.r.squared
                return(data.frame(adjR))
              },
              .id = "variable")
  
  rm(dd, spdd, envir = .GlobalEnv) # remove dd and spdd from global env
  return(r2)
}




get_full_formula <- function(x, 
                             n_term = 4 # number of terms to make combination
                             ){
  
  if(length(x) >= n_term) { # when there are >n_term(4 as default) terms
    comb_exp <- matrix(combn(x, n_term), nrow = n_term) 
  } else {             # when there are <n_term
    comb_exp <- matrix(combn(x, length(x)), nrow = length(x))
  }
  
  # create formula
  expl_fml <- apply(comb_exp, 2, function(y) paste(y, collapse = "+"))
  
  return(expl_fml)
}




get_rda_summary <- function(formula_list, df, expl, SiteName_rda){
  
  # check object names
  is_exist_spdd(ignoredd = TRUE)
  
  # create formulas
  f <- llply(paste("spdd ~", formula_list), as.formula)
  names(f) <- 1:length(f)
  
  
  # prepare df for rda
  dd <<-  df # df needs to be assigned to global environment
  spdd <<- select(dd, -one_of(expl, SiteName_rda))
  r <- llply(f, function(x) rda(x, data = dd))
  
  
  # generate summary results of rda
  r_summary <- ldply(r, function(x){
    
    # vif
    v <- vif.cca(x)
    vdfd <- data.frame(as.list(v))
    
    # adjusted R2
    adjr <- RsquareAdj(x)$adj.r.squared
    
    # permutation anova
    ar      <- anova(x, permutations = allPerms(6))
    p_value <- tidy(ar)$p.value[1]
    
    return(cbind(vdfd, adjr, vif_less10 = !any(v > 10), p_value))
  }, 
  .id = "f_id") %>%  # organise df
    select(-adjr, -vif_less10, -p_value, everything()) %>% 
    mutate(f_id = as.numeric(f_id)) %>% 
    arrange(-vif_less10, -adjr)
  
  rm(dd, spdd, envir = .GlobalEnv)
  return(r_summary)
  
}




get_simple_rda <- function(f, df, expl, SiteName_rda){
  
  
  # check object names
  is_exist_spdd(ignoredd = TRUE)
  
  
  # prepare df for rda
  dd <<-  df # df needs to be assigned to global environment
  spdd <<- select(dd, -one_of(expl, SiteName_rda))
  
  
  # make rda models
  rr1 <- rda(as.formula(paste("spdd ~", f)), data = dd)                            # full model determined by adjusted R2 and P value
  rr2 <- rda(spdd ~ 1, data = dd)                                                  # null model
  rr3 <- ordiR2step(rr2, rr1, permutations = allPerms(6), direction = "forward", Pin = .1)  # forward model simplification
  
  
  # get F and P value for each term in the final model (rr3)
  anova_rr3 <- anova(rr3, permutations = allPerms(6), by = "margin")
  
  rm(dd, spdd, envir = .GlobalEnv) # remove dd and spdd from the global environment
  return(list(full_mod = rr1, final_mod = rr3, anova_final_mod = anova_rr3))
  
}




get_rda_model_summary <- function(x){
  mod_adjr <- RsquareAdj(x)$adj.r.squared
  if(length(mod_adjr) == 0) mod_adjr <- NA
  
  mod_p    <- anova(x, permutations = allPerms(6))$'Pr(>F)'[1]
  return(data.frame(mod_adjr, mod_p))
}




get_rda_scores <- function(x){ # x:  RDA results
  
  RdaAllRes <- summary(x)                                                         # rda summary
  sitedd   <- data.frame(RdaAllRes$site, all_4y_d[, SiteName_rda])                # site score + site variables
  bipldd   <- data.frame(RdaAllRes$biplot, co2 = "amb", year = "Year0",           # scores for numeric predictor
                         variable = row.names(RdaAllRes$biplot)) %>%            
    mutate(variable = ifelse(variable == "as.numeric(year)", "Year", 
                             as.character(variable)))                             # re-label
  VarProp  <- RdaAllRes$cont$importance["Eigenvalue",] / RdaAllRes$tot.chi             
  # proportion of variablce explained by RDA1 and RDA2
  axislabs <- paste0(c("RDA1", "RDA2"), "(", round(VarProp[c(1, 2)] * 100, 2), "%)")
  
  return(list(sitedd = sitedd, bipldd = bipldd, VarProp = VarProp, axislabs = axislabs))
  
}


create_rda_plots <- function(sitedd,    # site score
                             bipldd,    # biplot score (numeric predictor) 
                             VarProp,   # variance propotion explained by RDA1 and 2
                             axislabs,  # axis labels
                             b_cons    # constant values for rescaling biplot
                             ){
  # Arguments are inhereted from get_rda_scores
  
  p <- ggplot(data = sitedd, aes(x = RDA1, y = RDA2, shape = year)) +
    
    # base lines and faceting (only for making top label)
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~ Form) +
    
    
    geom_path(aes(group = ring), col = "black") +
    geom_point(aes(fill = co2), size = 2, alpha = .7) +
    
    
    # numeric predictor
    geom_segment(data  = bipldd, 
                 aes(x = 0, y = 0, xend = RDA1 * b_cons, yend = RDA2 * b_cons), 
                 arrow = arrow(length = unit(.2, "cm")),  color = "red", alpha = .7) +
    geom_text(data = bipldd, 
              aes(x = RDA1 * b_cons, y = RDA2 * b_cons, label = variable), 
              color = "red", size = 3, 
              fontface = "bold") +
    
    # ring id
    geom_text(data = filter(sitedd, year == "Year0"), 
              aes(label = ring), size = 3, hjust = 1, 
              fontface = "italic") +
    
    
    # scales
    scale_fill_manual(name   = expression(CO[2]),
                      values = c("grey40", "white"), 
                      labels = c("Ambient", expression(eCO[2])),
                      guide  = guide_legend(override.aes = list(shape = 21))) +
    scale_shape_manual(name   = "Year",
                       values = c(21, 22, 23, 24)) +
    science_theme +
    theme(legend.position = "none") +
    labs(x = axislabs[1], y = axislabs[2])
  
  
  return(p)
}
