# produce sheetname
a <- as.vector(outer(1:6, 1:4, paste, sep = "."))
shts <- as.vector(outer(a, LETTERS[1:4], paste, sep = "."))

# raed all files
options(java.parameters = "-Xmx100m") # increase java memory

fls <- llply(shts, function(x) 
  read.veg.xlx(file = "Data/Result_FACE_Vegetation_Datasheet_2014.xlsx",
               sheetName = x),
             .progress = "text")

# combine
veg.2014.raw <- rbind.fill(fls) 

# ring, plot, pos, year
splt <- strsplit(veg.2014.raw$position, "[.]")
veg.2014.raw$ring <- factor(sapply(splt, "[", 1))
veg.2014.raw$plot <- factor(sapply(splt, "[", 2))
veg.2014.raw$pos <- factor(sapply(splt, "[", 3))
veg.2014.raw$year <- factor("2014")

# sort colmuns
NotSpp <- c("year", "ring", "plot", "pos", "cell", "position")
Spp <- sort(names(veg.2014.raw)[which(!(names(veg.2014.raw) %in% NotSpp))])

veg.2014.raw <- veg.2014.raw[c(NotSpp, Spp)]
names(veg.2014.raw)

# Turn Spp into numeric
veg.2014.raw[Spp] <- apply(veg.2014.raw[Spp], 2, as.numeric)

# turn na into 0
veg.2014.raw[is.na(veg.2014.raw)] <- 0


# remove the spp which were not found
kpt.sp<- names(which(colSums(veg.2014.raw[Spp]) != 0))
veg.14 <- veg.2014.raw[c(NotSpp, kpt.sp)]
veg.14$position <- NULL
names(veg.14)[4] <- "position"
save(veg.14, file = "output/Data/FACE_Veg2014.RData")
