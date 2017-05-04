
# packages ----------------------------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(car)
library(lmerTest)
library(vegan)
library(tidyr)
library(lsmeans)
library(visreg)
library(MuMIn)
library(LMERConvenienceFunctions)



# load data ---------------------------------------------------------------

graminoid_data <- read.csv("output/Data/graminoid_data.csv")  # data frame for graminoids
forb_data      <- read.csv("output/Data/forb_data.csv")       # data frame for forbs
site_var       <- c("year", "ring", "co2", "plot")            # vector for site variables
site_data      <- graminoid_data[ ,site_var]                  # data frame for site
SppName_gram   <- setdiff(names(graminoid_data), site_var)    # graminoid species names
SppName_forb   <- setdiff(names(forb_data), site_var)         # forb species names
env_data       <- read.csv("output/Data/env_data.csv")        # environmental variables (moisture, temperature and PAR)
sp_pfg         <- read.csv("output/Data/graminoid_pfg.csv")        # graminoid species and and their corresponding plant functional groups




# analysis ----------------------------------------------------------------

source("R/analysis_diversity.R")  # analysis on diversity indices
source("R/analysis_abundance.R")  # analysis on abundance of Microlaena, Cynodon, and total C3 and C4
source("R/analysis_LAR.R")        # multiple regression analysis on log annual change ratios
source("R/analysis_PRC.R")        # multivariate analysis with principal response curve
