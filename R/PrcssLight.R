#########
# Light #
#########
# download from HIEv
setToken(tokenfile = "Data/token.txt")

# download files from HIEv
Light_Hiev <- downloadTOA5("FACE.*AirVars.*dat", 
                           cachefile = "Data/hievdata/tmp.RData",
                           topath = "Data/hievdata/raw_data",
                           maxnfiles = 999)
  # LI190SB_PAR_Den_Avg at 23 m
  # PAR_Den_1-3_Avg at 0.15 m
  # Look at
  # https://sites.google.com/site/hievuws/facilities/eucface/collection-codes for
  # the description of each variables

# save(Light_Hiev, file = "output/Data/FACE_AirVars_Raw.RData")
# load("output/Data/FACE_AirVars_Raw.RData")

# organise data frame. Note that this is really huge so use dplyr----

# Add ring number and remove Source
Light_Hiev <- select(mutate(Light_Hiev, 
                            ring = factor(substr(Light_Hiev$Source, 7, 7))), 
                     -Source)
# remove duplicates
Light_Hiev <- distinct(Light_Hiev)

# save
save(Light_Hiev, file = "output/Data/FACE_AirVars_Processed.RData")

# select only required columns as it's too large
Lightdf <- select(Light_Hiev, DateTime, Date, ring, PAR_Den_1_Avg, PAR_Den_2_Avg, PAR_Den_3_Avg)
save(Lightdf, file = "output/Data/FACE_FloorPAR.RData")



summary(Light_Hiev)

Light_Hiev <- Light_Hiev[, 1:6]
Light_Hiev <- within(Light_Hiev, {
  Ring <- as.character(Ring)
  ring <- as.factor(gsub("R", "", Ring))
  Ring <- NULL
  DateTime <- ymd_hms(Light_Hiev$DateTime)
  date <- as.Date(DateTime)
})

# remove NA in ring column
Light_Hiev <- Light_Hiev[!is.na(Light_Hiev$ring),]

# Daily mean
SiteVec <- c("DateTime", "date", "ring")
Vars <- names(Light_Hiev)[!names(Light_Hiev) %in% SiteVec]

Light_Hiev_Day <- ddply(Light_Hiev, .(date, ring),
                        function(x) colMeans(x[, Vars], na.rm = TRUE))
Light_Hiev_Day_mlt <- melt(Light_Hiev_Day, id = c("date", "ring"))

Light_Hiev_Day_mlt$Type <- factor(ifelse(Light_Hiev_Day_mlt$variable == "LI190SB_PAR_Den_Avg", 
                                         "CanopyPar", "FloorPar"))

p <- ggplot(Light_Hiev_Day_mlt, aes(x = date, y = value, col = variable))
p2 <- p + geom_point() + facet_grid(ring ~. )
p2

Light_type <- ddply(Light_Hiev_Day_mlt, .(ring, date, Type), summarise, 
                    Mean = mean(value, na.rm = TRUE))

boxplot(Mean ~ ring, data = subset(Light_type, Type == "FloorPar"))




