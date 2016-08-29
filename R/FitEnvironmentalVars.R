
# OrganiseEnvVars ---------------------------------------------------------

# Environmental variables
  load("output//Data/FACE_EnvironmenVars.RData")
  
  names(EnvVarDF)

# remove all minor metals
  rmMetal  <- c("Silver", "Arsenic", "Lead", "Cadmium", "Chromium", "Copper", 
                 "Manganese","Nickel", "Selenium", "Zinc", "Mercury", "Iron", 
                 "Aluminium", "Boron", "Silicon", "Vanadium", "Cobalt", 
                 "Molybdenum", "Barium", "Calcium", "Magnesium", "Potassium", 
                 "Sodium")

  EnvVarDF <- EnvVarDF[, !names(EnvVarDF) %in% rmMetal]

# correlations
  CorMatrix <- cor(Rm_ymc(EnvVarDF), use = "pairwise.complete.obs")
  p = CorMatrix
  p[is.na(CorMatrix)]=0.2 
  p[is.na(CorMatrix)==F]=0
  CorMatrix[is.na(CorMatrix)]=0
  corrplot.mixed(CorMatrix, p.mat=p, tl.cex = .5)

# some values are really redundant so remove. Also some of them are hard to
# interpret so remove. e.g. EC, Theta75
  RmVar     <- c("sand", "silt", "clay", "OrganicMatter",
                 "OrganicC", "OrganicC", "TotalP_CM", 
                 "Theta5", "Theta30", "Theta75", "EC", 
                 "T5", "T10", "T20", "T30", "T50", "T100", 
                 "TotalN", "Phosphorus", "nitrification",
                 "ThetaHL")
  EnvVarDF  <- EnvVarDF[, !names(EnvVarDF) %in% RmVar]
  CorMatrix <- cor(Rm_ymc(EnvVarDF), use = "pairwise.complete.obs")
  
  corrplot.mixed(CorMatrix, tl.cex = .5, mar = par()$mar)

# some of the variables are not complete for four years

# Variables measured only in the 1 and 2nd years----
  naCol_2ndYear <-cbind(EnvVarDF[, c("year", "ring")], 
                        EnvVarDF[, apply(EnvVarDF, 2, function(x) any(is.na(x)))])

###########################
# Analyse 3-year data set #
###########################
  
  naCol_2ndYear
  # these variables may vary year by year. so remove them for the time being
  EnvDF_3df <- EnvVarDF[, apply(EnvVarDF, 2, function(x) !any(is.na(x)))]
  
  EnvDF_3df$block <- recode(EnvDF_3df$ring, "1:2 = 'A'; 3:4 = 'B'; 5:6 = 'C'")
  save(EnvDF_3df, file = "output/Data/EucFACE_understorey_env_vars_2012-2016.RData")
  write.csv(EnvDF_3df, file = "output/Data/EucFACE_understorey_env_vars_2012-2016.csv",
            row.names = FALSE)
  
  
# RDA ---------------------------------------------------------------------
  
  # possible explanatory variables
  expl <- c("co2",  "TotalC", "moist", "Drysoil_ph", "Depth_HL", "gapfraction", "temp")
  
  # All species
  source("R/Stats_RDA_AllSpp.R")
  
  # PFG
  source("R/Stats_RDA_PFG.R")


