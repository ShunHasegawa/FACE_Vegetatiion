###############
# All species #
###############


# Each year separately ----------------------------------------------------

  # combine environment and spp df
  seDF <- merge(RingSumVeg, EnvDF_3df, by = c("year", "ring", "block", "co2"))
  
  
# . Single term -----------------------------------------------------------

  # R2adj for each single term
  
  adjR_singl_Lst <- list()
  for (i in 1:4) {
    dd                  <- subsetD(seDF, year == levels(seDF$year)[i])
    spdd                <- dd[ , SppName]
    # formula for each variable
    singl_fmls          <- llply(paste("log(spdd + 1) ~", expl), as.formula)
    names(singl_fmls)   <- expl
    adjR_singl          <- ldply(singl_fmls, 
                                 function(y) {
                                   RsquareAdj(rda(y, data = dd))$adj.r.squared
                                   },
                                 .id = "variable"
                                 )
    adjR_singl_Lst[[i]] <- adjR_singl
    rm(dd, spdd, singl_fmls, adjR_singl)
  }
  
  names(adjR_singl_Lst) <- paste0("Year", 0:3)
  
  # Get variables with positive R2adj
  PosAdjR <- llply(adjR_singl_Lst, function(x) as.character(x$variable[x$V1 > 0]))
  
  # Formula for full models
  FullFormula <- llply(PosAdjR, 
                       function(x) {
                          comb_exp <- combn(x, 4) # combination of four
                          expl_fml <- apply(comb_exp, 
                                            2, 
                                            function(x) paste(x, collapse = "+")
                                            )
                          return(expl_fml)
                          }
                       )
  
  # Y matric (lefthand part)
  LH   <- list("log(df2013[ , SppName] + 1) ~",
               "log(df2014[ , SppName] + 1) ~",
               "log(df2015[ , SppName] + 1) ~", 
               "log(df2016[ , SppName] + 1) ~" )
  fmls <- llply(list(Year0 = 1, Year1 = 2, Year2 = 3, Year3 = 4), 
                function(x) llply(paste(LH[[x]], FullFormula[[x]]), as.formula))


# . 1st year --------------------------------------------------------------

  df2013 <- subsetD(seDF, year == "Year0")
  
  # There are too many environmental variables to fit. so choose four which showed
  # highest R2adj 
  # adjusted R2
  adjR <- ldply(fmls$Year0, function(x) RsquareAdj(rda(x, data = df2013))$adj.r.squared)
  
  # highest R2
  rr <- rda(fmls$Year0[[which(max(adjR) == adjR)]], df2013)
  
  # check multicollinearity
  vif.cca(rr)
  anova(rr, permutations = allPerms(6))
  rr2 <- rda(log(df2013[ , SppName] + 1) ~ 1, df2013)
  rr3 <- ordiR2step(rr2, rr, permutations = allPerms(6), direction = "forward", Pin = .1)
  summary(rr3)
  
  # summary result
  rda2013 <- list(IniRda = rr, FinRda = rr3)

  
# . 2nd year --------------------------------------------------------------

  df2014 <- subsetD(seDF, year == "Year1")
  
  # adjusted R2
  adjR <- laply(fmls$Year1, function(x) RsquareAdj(rda(x, data = df2014))$adj.r.squared)
  
  # highest R2
  rr <- rda(fmls$Year1[[which(max(adjR) == adjR)]], df2014)
  
  # check multicolliniarity using vif
  vif.cca(rr)
  anova(rr, permutations = allPerms(6))
  # not significant
  
  # choose only three terms
  comb_exp <- combn(PosAdjR$Year1, 3)
  expl_fml <-apply(comb_exp, 2, function(x) paste(x, collapse = "+"))
  fmls_3 <- llply(paste("log(df2014[ , SppName] + 1) ~", expl_fml), as.formula)
  
  adjR <- laply(fmls_3, function(x) RsquareAdj(rda(x, data = df2014))$adj.r.squared)
  rr <- rda(fmls_3[[which(max(adjR) == adjR)]], df2014)
  vif.cca(rr)
  anova(rr, permutations = allPerms(6))
  # good
  rr2 <- rda(log(df2014[ , SppName] + 1) ~ 1, df2014)
  rr3 <- ordiR2step(rr2, rr, permutations = allPerms(6), direction = "forward", Pin = .1)
  
  # summary result
  rda2014 <- list(IniRda = rr, FinRda = rr3)


# . 3rd year --------------------------------------------------------------

  df2015 <- subsetD(seDF, year == "Year2")
  
  # adjusted R2
  adjR <- laply(fmls$Year2, function(x) RsquareAdj(rda(x, data = df2015))$adj.r.squared)
  
  # highest R2
  rr <- rda(fmls$Year2[[which(max(adjR) == adjR)]], df2015)
  
  # check multicollinearity
  vif.cca(rr)
  # TotalC >10. 
  
  # There are only four terms for this year, so use three terms
  comb_exp <- combn(PosAdjR$Year2, 3)
  expl_fml <-apply(comb_exp, 2, function(x) paste(x, collapse = "+"))
  fmls_3 <- llply(paste("log(df2015[ , SppName] + 1) ~", expl_fml), as.formula)
  
  adjR <- laply(fmls_3, function(x) RsquareAdj(rda(x, data = df2015))$adj.r.squared)
  rr <- rda(fmls_3[[which(max(adjR, na.rm = TRUE) == adjR)]], df2015)
  
  vif.cca(rr)
  anova(rr, permutations = allPerms(6))
  # good
  
  rr2 <- rda(log(df2015[ , SppName] + 1) ~ 1, df2015)
  rr3 <- ordiR2step(rr2, rr, permutations = allPerms(6), direction = "forward", Pin = .1)
  
  # summary result
  rda2015 <- list(IniRda = rr, FinRda = rr3)

  

# . 4th year --------------------------------------------------------------

  df2016 <- subsetD(seDF, year == "Year3")
  
  # adjusted R2
  adjR <- laply(fmls$Year3, function(x) RsquareAdj(rda(x, data = df2016))$adj.r.squared)
  
  # highest R2
  rr <- rda(fmls$Year3[[which(max(adjR) == adjR)]], df2016)
  
  # check multicollinearity
  vif.cca(rr)
  anova(rr, permutations = allPerms(6))
    # not significant
  
  # choose only three terms
  comb_exp <- combn(PosAdjR$Year3, 3)
  expl_fml <- apply(comb_exp, 2, function(x) paste(x, collapse = "+"))
  fmls_3   <- llply(paste("log(df2016[ , SppName] + 1) ~", expl_fml), as.formula)
  adjR     <- laply(fmls_3, function(x)
                      RsquareAdj(rda(x, data = df2016))$adj.r.squared)
  rr       <- rda(fmls_3[[which(max(adjR, na.rm = TRUE) == adjR)]], df2016)
  vif.cca(rr)
  anova(rr, permutations = allPerms(6))
  # good
  
  rr2 <- rda(log(df2016[ , SppName] + 1) ~ 1, df2016)
  rr3 <- ordiR2step(rr2, rr, permutations = allPerms(6), direction = "forward", Pin = .1)
  
  # summary result
  rda2016 <- list(IniRda = rr, FinRda = rr3)
  

# Summary ---------------------------------------------------------------

  # Adjusted R2 for each term
  AdjTbl        <- dcast(variable ~ .id, 
                         data      = ldply(adjR_singl_Lst), 
                         value.var = "V1")
  AdjTbl[, 2:5] <- round(AdjTbl[, 2:5], 3)
  # replace negative values with <0
  AdjTbl[AdjTbl < 0] <- "<0" 
    
  ## Results of anova for full/parsimonious models ##
  ## R2adj for initial full model ##
  RdaLst     <- list(Year0 = rda2013, 
                     Year1 = rda2014, 
                     Year2 = rda2015, 
                     Year3 = rda2016)
  
  # FinRDA in Year1 doesn't have any explanatory variable so remove
  RdaLst$Year1$FinRda <- NULL
  
  # AdjustedR2 and P values for each model
  Extract_Adj_P <- function(x) {
                   data.frame(AdjR = RsquareAdj(x)$adj.r.squared, 
                              Pr   = anova(x, 
                                           permutations = allPerms(6))$Pr[1])
                   }



  FuladjR_pv <- ldply(RdaLst, 
                      function(x) ldply(x, Extract_Adj_P, .id = "Model"),
                      .id = "year")
  
  FuladjR_pv_tbl <- dcast(variable ~ year + Model, 
                          data = melt(FuladjR_pv, id = c("year", "Model")))
  FuladjR_pv_tbl[ , 2:8] <- round(FuladjR_pv_tbl[ , 2:8], 4)
  
  # F and P values for each term in parsimonious models
  rda_anova <- ldply(RdaLst[-2], 
                     function(x) {
                     a <- anova(x$FinRda, 
                                permutations = allPerms(6), 
                                by           = "margin")
                     ad <- data.frame(term = row.names(a), 
                                      a[c(1, 3, 4)])
                     ad[, 2:4] <- round(ad[, 2:4], 3)
                     return(ad)
                     },
                     .id = "year")
  rda_anova$term <- factor(rda_anova$term, 
                           levels = c("TotalC", "moist", "Residual"))
  rda_anova_tbl <- dcast(term ~ year + variable, 
                         data = melt(rda_anova, id = c("year", "term")))

  # save
  wb <- createWorkbook()
  sheet1 <- createSheet(wb, sheetName = "adjustedR2")
  addDataFrame(AdjTbl, 
               sheet1, 
               showNA      = TRUE, 
               row.names   = FALSE,
               characterNA = "NA")
  sheet2 <- createSheet(wb, sheetName = "summary_models")
  addDataFrame(FuladjR_pv_tbl, 
               sheet2, 
               showNA      = TRUE, 
               row.names   = FALSE, 
               characterNA = "NA")
  sheet3 <- createSheet(wb, sheetName = "F_P_values_parsimonious_model")
  addDataFrame(rda_anova_tbl, 
               sheet3, 
               showNA      = TRUE, 
               row.names   = FALSE, 
               characterNA = "NA")
  saveWorkbook(wb, "output/table/RDA_Restuls_AllSpp.xlsx")
  


# 4-year data set ---------------------------------------------------------

  # From the above analysis, moist and TotalC are determied to be imporatnt
  # driver
  
  rda_all <- rda(log(seDF[, SppName] + 1) ~ TotalC + moist + year, seDF)
  rda_all <- rda(log(seDF[, SppName] + 1) ~ TotalC + moist + as.numeric(year), seDF)
  
  # can't run anova as it is. cause different perumutation units need to be
  # defined for year and co2. so anyway create a triplot and see the pattern.
  rda_all
  
  # plot
  p <- TriPlot(MultValRes = rda_all, env = seDF, yaxis = "RDA axis", axispos = c(1, 2, 3), centcons = 2)
  ggsavePP(filename = "output//figs/FACE_RDA_EnvVar_Year0_3", plot = p, width = 6, height = 6)


# Fig for thesis ----------------------------------------------------------

  RdaAllRes <- summary(rda_all)
  seDF$year <- factor(seDF$year, labels = paste0("Year", 0:3))
  sitedd <- data.frame(RdaAllRes$site, seDF)
  
  bipldd <- data.frame(RdaAllRes$biplot, co2 = "amb", year = "Year0", 
                       variable = row.names(RdaAllRes$biplot))
  bipldd <- subsetD(bipldd, variable %in% c("TotalC", "moist", "as.numeric(year)"))
  bipldd$variable <- factor(bipldd$variable,
                            levels = c("TotalC", "moist", "as.numeric(year)"),
                            labels = c("Moist", "Total C", "Year"))
  
  VarProp <- RdaAllRes$cont$importance["Eigenvalue",] / RdaAllRes$tot.chi
  axislabs <- paste0(c("RDA1", "RDA2"), "(", round(VarProp[c(1, 2)] * 100, 2), "%)")
  
  # make a plot
  p <- ggplot(data = sitedd, aes(x = RDA1, y = RDA2, shape = year))
  p2 <- p + 
    geom_path(aes(group = ring), col = "black") +
    geom_point(aes(fill = co2), size = 4) + 
    scale_fill_manual(values = c("black", "white"), 
                      labels = c("Ambient", expression(eCO[2])),
                      guide = guide_legend(override.aes = list(shape = 21))) +
    # need to add overrisde here to make white circle in the legend with shape of
    # 21
    scale_shape_manual(values = c(21, 22, 23, 24)) + 
    geom_segment(data = bipldd,
                 aes(x = 0, y = 0, xend = RDA1 * 2, yend = RDA2 * 2), 
                 arrow = arrow(length = unit(.2, "cm")), 
                 color = "red") +
    geom_text(data = bipldd, 
              aes(x = RDA1 * 2.3 , y = RDA2 * 2.3, label = variable), 
              lineheight = .7, 
              color = "red", size = 4, 
              fontface = "bold") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    science_theme +
    theme(legend.position = c(.17, .15), 
          legend.box = "horizontal", 
          legend.box.just = "top") +
    labs(x = axislabs[1], y = axislabs[2])
  RDA_Plot_AllSpp <- p2
  RDA_Plot_AllSpp
  ggsavePP(plot     = RDA_Plot_AllSpp, 
           filename = "output/figs/Fig_Thesis/RDA_3yr_AllSpp", 
           width    = 6, 
           height   = 4)
