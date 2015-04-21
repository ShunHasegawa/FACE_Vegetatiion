# function which reads worksheet from an xcel file
read.veg.xlx <- function(sheetName, file, ...) {
  a <- read.xlsx2(file, sheetName, ...)
  a <- a[ ,!grepl("X.", names(a))]
  a[a == ""]  <- 0 # empty cell -> 0
  xlcFreeMemory() # Frees Java Virtual Machine (JVM) memory
  # it's not normally necessary to do this every time but the file size is 
  # really huge and cannot read all the worksheets at once so free memory every time 
  return(a)
}

##################
# September 2012 #
##################
# get sheetnames
wb <- loadWorkbook("Data/FACE_Vegetation_SEP2012.xlsx")
sheets <- names(getSheets(wb))
options(java.parameters = "-Xmx100m") # increase java memory

# read files
fls <- ldply(sheets, function(x) 
  read.veg.xlx(file = "Data/FACE_Vegetation_SEP2012.xlsx", sheetName = x ,
               header = TRUE, startRow = 1, endRow = 401, stringsAsFactors = FALSE),
             .progress = "text")

# only forbs are used from this data so subset
forbs <- fls[,which(substring(names(fls),1,5)!="Grass" & substring(names(fls),1,5)!="Sedge")]

# change values into numeric
forbs[, !names(forbs) %in% c("ring", "plot", "position", "cell")] <- 
  apply(forbs[, !names(forbs) %in% c("ring", "plot", "position", "cell")], 2, as.numeric)

# few of them were not observed in December 2012 so turn them into 0----
# 2.1.D 
forbs[forbs$ring == 2 & forbs$plot == 1 & forbs$position == "D", "Forb.Senecio.madagascariensis"] <- 0
# 2.4.D 
forbs[forbs$ring == 2 & forbs$plot == 4 & forbs$position == "D", "Forb.Desmodium.rhytidophyllum"] <- 0
# 6.1.D 
forbs[forbs$ring == 6 & forbs$plot == 1 & forbs$position == "D", "Shrub.Bossiaea.prostrata"] <- 0
# 6.4.D
forbs[forbs$ring == 6 & forbs$plot == 4 & forbs$position == "A", "Forb.Unknown.herb"] <- 0

# 5.4.B Hypochaeris.radicata -> Leontodon taraxacoides
CellIndx <- which(forbs[forbs$ring == 5 & forbs$plot == 4 & forbs$position == "B", "Forb.Hypochaeris.radicata"] == 1)
## Turn Hypochaeris.radicata into 0
forbs[forbs$ring == 5 & forbs$plot == 4 & forbs$position == "B", "Forb.Hypochaeris.radicata"] <- 0
## Add Leontodon taraxacoides
forbs$Forb.Leontodon.taraxacoides <- 0
forbs[forbs$ring == 5 & forbs$plot == 4 & forbs$position == "B", "Forb.Leontodon.taraxacoides"][CellIndx] <- 1

# correct sp names----

# Shrub.Opercularia.sp -> Shrub.Opercularia.diphylla
names(forbs)[grep("Opercularia", names(forbs))] <- "Shrub.Opercularia.diphylla"

# Wahlenbergia.sp
names(forbs)[grep("Wahlenbergia.sp", names(forbs))] <- "Forb.Wahlenbergia.gracilis"

# remove PFG
names(forbs) <- gsub("Forb.|Shrub.|Moss.|Tree.|Fern.|Lichen.", "", names(forbs))

#################
# December 2012 #
#################

# get sheetnames
wb <- loadWorkbook("Data/FACE_Vegetation_DEC2012.xlsx")
sheets <- names(getSheets(wb))

# read files
fls <- ldply(sheets, function(x) {
  a <- read.veg.xlx(file = "Data/FACE_Vegetation_DEC2012.xlsx", sheetName = x ,
                    header = TRUE, startRow = 1, endRow = 26, stringsAsFactors = FALSE)
  site <- ldply(strsplit(x, split = "[.]"))
  names(site) <- c("ring", "plot", "position")
  a <- data.frame(a, site)
  return(a)},
  .progress = "text")

# orgaise data frame----
dec2012 <- fls

# re-oder columns
SiteVec <- c("ring", "plot", "position", "Cell")
decName <- names(dec2012)[!names(dec2012) %in% SiteVec]
dec2012 <- dec2012[, c("ring", "plot", "position", "Cell", sort(decName))]

# Cell -> cell
names(dec2012)[names(dec2012) == "Cell"] <- "cell"

# turn "..." or ".." to .
names(dec2012) <- gsub("[.][.][.]|[.][.]", ".", names(dec2012))

#######################
# Combine Sep and Dec #
#######################
names(forbs)
names(dec2012)

df2013 <- rbind.fill(forbs, dec2012)
# re-organise data----

# turn NA into 0
df2013[is.na(df2013)] <- 0 

# turn values into numeric
SiteVec <- c("ring", "plot", "position", "cell")
Spps <- names(df2013)[!names(df2013) %in% SiteVec]
df2013[, Spps] <- apply(df2013[, Spps], 2, as.numeric)

# two data frames are row-binded, so each cell is duplicated. Take sum of each cell
nrow(df2013)
df2013 <- ddply(df2013, .(ring, plot, position, cell), function(x) colSums(x[, Spps]))

# remove spp with no count
df2013 <- df2013[, c(SiteVec, Spps[!colSums(df2013[, Spps]) == 0])]

# rmeove unknown spp
df2013 <- df2013[!grepl("unknown|unkown|Perennial.Grass.Tuft", names(df2013), ignore.case = TRUE)]

# remove and combien some spp----
SpName <- names(df2013)[!names(df2013) %in% SiteVec]

# reorder columns
df2013 <- df2013[, c(SiteVec, sort(SpName))]

# Axonopsis and Axonopus.sp -> Axonopus.fissifolius
df2013 <- OrgSpp(df2013, KeepCol = "Axonopus.fissifolius", 
                 CombineCol = SpName[grepl("Axonop", SpName, ignore.case = TRUE)])

# Bursaria -> Bursaria.spinosa
df2013 <- OrgSpp(df2013, KeepCol = "Bursaria.spinosa", 
                 CombineCol = SpName[grepl("Bursaria", SpName, ignore.case = TRUE)])

# Carex.breviformis -> Carex.beviculmis
df2013 <- OrgSpp(df2013, KeepCol = "Carex.breviculmis", 
                 CombineCol = SpName[grepl("carex", SpName, ignore.case = TRUE)])

# Commelina -> Commelina.cyanea
df2013 <- OrgSpp(df2013, KeepCol = "Commelina.cyanea", 
                 CombineCol = SpName[grepl("Commelina", SpName, ignore.case = TRUE)])

# Conyza. -> Conyza.sumatrensis
df2013 <- OrgSpp(df2013, KeepCol = "Conyza.sumatrensis", 
                 CombineCol = SpName[grepl("Conyza", SpName, ignore.case = TRUE)])

# Cyperus.flacidus, Cyperus.sp -> Cyperus.flaccidus
# Cyperus.sp is not too sure, but its really small abundance, so just combine
sum(df2013[, "Cyperus.sp"])
df2013 <- OrgSpp(df2013, KeepCol = "Cyperus.flaccidus", 
                 CombineCol = SpName[grepl("Cyperus", SpName, ignore.case = TRUE)])

# Glycine.clandestina, Glycine.sp, Glycine.sp.1 -> Glycine.sp These ones are
# really hard to identfy. Also other than Glycine.clandestina abundance is so
# small so just combine
df2013 <- OrgSpp(df2013, KeepCol = "Glycine.sp", 
                 CombineCol = SpName[grepl("Glycine", SpName, ignore.case = TRUE)])

# Hypercium -> Hypericum.gramineum
df2013 <- OrgSpp(df2013, KeepCol = "Hypericum.gramineum", 
                 CombineCol = SpName[grepl("Hyper", SpName, ignore.case = TRUE)])

# Juncus.sp, Juncus.continuus. -> Juncus.continuus
df2013 <- OrgSpp(df2013, KeepCol = "Juncus.continuus", 
                 CombineCol = SpName[grepl("Juncus", SpName, ignore.case = TRUE)])

# Lachnagrostis -> Lachnagrostis.filiformis
names(df2013)[names(df2013) == "Lachnagrostis"] <- "Lachnagrostis.filiformis"

# Laxmannia, Laxmannia.gracilis. -> Laxmannia.gracilis
df2013 <- OrgSpp(df2013, KeepCol = "Laxmannia.gracilis", 
                 CombineCol = SpName[grepl("Laxmannia", SpName, ignore.case = TRUE)])

# Opercularia. -> Opercularia.diphylla
df2013 <- OrgSpp(df2013, KeepCol = "Opercularia.diphylla", 
                 CombineCol = SpName[grepl("Opercularia", SpName, ignore.case = TRUE)])
# assume that paspalidium was paspalidium distans at this time --> needs to be
# corrected later. Also Paspalidium.distans = Paspalidium.radiatum
df2013 <- OrgSpp(df2013, KeepCol = "Paspalidium.distans", 
                 CombineCol = SpName[grepl("paspalidium", SpName, ignore.case = TRUE)])

# Poranthera -> Poranthera.microphylla 
df2013 <- OrgSpp(df2013, KeepCol = "Poranthera.microphylla ", 
                 CombineCol = SpName[grepl("Poranthera", SpName, ignore.case = TRUE)])

# Schoenus -> Schoenus.opogon
df2013 <- OrgSpp(df2013, KeepCol = "Schoenus.opogon", 
                 CombineCol = SpName[grepl("Schoenus", SpName, ignore.case = TRUE)])

# Sisyrinchium. -> Sisyrinchium.Iridaceae
df2013 <- OrgSpp(df2013, KeepCol = "Sisyrinchium.Iridaceae", 
                 CombineCol = SpName[grepl("Sisyrinchium", SpName, ignore.case = TRUE)])

# sp. -> sp
names(df2013) <- gsub("sp[.]", "sp", names(df2013))

# add .sp at the end of sp name that were not identified
## end with "."
names(df2013)[grepl("[.]$", names(df2013))] <- paste0(names(df2013)[grepl("[.]$", names(df2013))], "sp")

## no species
names(df2013)[!grepl("[.]", names(df2013))][5:7] <- paste0(names(df2013)[!grepl("[.]", names(df2013))][5:7], "sp") 

# check if there're valuse > 2----
SpName <- names(df2013)[!names(df2013) %in% SiteVec]
all(df2013[SpName] < 2)
# false; turn them into 1
df2013[SpName][which(df2013[SpName] > 1, arr.ind = TRUE)] <- 1

# check agian
all(df2013[SpName] < 2)

# organise site vectors
df2013 <- within(df2013, {
  ring <- factor(ring)
  plot <- factor(plot)
  position <- factor(position)
  cell <- factor(cell)
  year <- "2013"
  })

# save
save(df2013, file = "output//Data/FACE_Vegetation_2013.RData")
