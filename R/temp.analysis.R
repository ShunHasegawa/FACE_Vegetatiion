
# load data ---------------------------------------------------------------

graminoid_data <- read.csv("output/Data/graminoid_data.csv")          # data frame for graminoids
forb_data      <- read.csv("output/Data/forb_data.csv")               # data frame for forbs
site_var       <- c("year", "ring", "co2", "plot")                    # vector for site variables
site_data      <- graminoid_data[ ,site_var]                          # data frame for site
SppName_gram   <- setdiff(names(graminoid_data), site_var)            # graminoid species names
SppName_forb   <- setdiff(names(forb_data), site_var)                 # forb species names




