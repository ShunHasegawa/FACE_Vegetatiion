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
OrgSpp <- function(df, KeepCol, CombineCol){
  df[KeepCol] <- rowSums(df[CombineCol])
  RemoveCol <- CombineCol[CombineCol != KeepCol]
  df <- df[!names(df) %in% RemoveCol]
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
# Compute dissimilarity for each plot between 2012 & 2014 and 2014 & 2015
YearDssmlrty <- function(x, spp) {
  df1 <- subset(x, year %in% c(2012, 2014))
  dis1 <- vegdist(log(df1[, spp] + 1), method = "bray")
  df2 <- subset(x, year %in% c(2014, 2015))
  dis2 <- vegdist(log(df2[, spp] + 1), method = "bray")
  dfs <- data.frame(year = c("Year1", "Year2"), BC = c(dis1, dis2))
  return(dfs)
}
