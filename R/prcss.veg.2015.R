# produce sheetname
a <- as.vector(outer(1:6, 1:4, paste, sep = "."))
shts <- as.vector(outer(a, LETTERS[1:4], paste, sep = "."))

# raed all files
options(java.parameters = "-Xmx100m") # increase java memory
veg.2015.raw <- ldply(shts, function(x) 
  read.veg.xlx(file = "Data/Result_FACE_Vegetation_Datasheet_2015.xlsx",
               sheetName = x), .progress = "text")

Spp <- names(veg.2015.raw)[!names(veg.2015.raw) %in% c("cell", "position")]

# turn P/p into 1 and others into 0
veg.2015 <- veg.2015.raw
veg.2015[Spp] <- apply(veg.2015[Spp], 2, 
                           function(x) recode(x, "c('p', 'P') = 1; else = 0"))

# remove spp which were not observed
veg.2015 <- cbind(veg.2015[c("cell", "position")], 
                      veg.2015[Spp][colSums(veg.2015[Spp]) != 0])

# organise species
ns <- names(veg.2015)
sort(ns)

# Combine anagallis
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Anagallis.arvensis", 
                       CombineCol = ns[grepl("Anagallis", ns)])
# Arthropodium
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Arthropodium.sp", 
                       CombineCol = ns[grepl("Arthropodium", ns)])

# Axonopus
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Axonopus.fissifolius", 
                       CombineCol = ns[grepl("Axonopus", ns)])

# Bidens.pilosa
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Bidens.pilosa", 
                       CombineCol = ns[grepl("Bidens.pilosa", ns)])

# Bidens.alternans = Bidens.subalternans (?)
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Bidens.subalternans", 
                       CombineCol = ns[grepl("alternans", ns)])

# Breynia.oblongifolia
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Breynia.oblongifolia", 
                       CombineCol = ns[grepl("Breynia.oblongifolia", ns)])

# Conyza.sumatrensis
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Conyza.sumatrensis", 
                       CombineCol = ns[grepl("Conyza", ns)])

# Echinopogon.caespitosus
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Echinopogon.caespitosus", 
                       CombineCol = ns[grepl("Echinopogon", ns)])

# Eragrostis.curvula
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Eragrostis.curvula", 
                       CombineCol = ns[grepl("curvula", ns)])

# Eragrostis.leptostachya
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Eragrostis.leptostachya", 
                       CombineCol = ns[grepl("leptostachya", ns)])

# Oplismenus.aemulus
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Oplismenus.aemulus", 
                       CombineCol = ns[grepl("Oplismenus|Oplismemus", ns)])

# Paspalidium.distans
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Paspalidium.distans", 
                       CombineCol = ns[grepl("Paspalidium", ns)])

# Phyllanthus.sp
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Phyllanthus.sp", 
                       CombineCol = ns[grepl("Phyllanthus|Phyllathus", ns)])

# Pratia.purpurascens
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Pratia.purpurascens", 
                       CombineCol = ns[grepl("Pratia", ns)])

# Senecio.madagascariensis
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Senecio.madagascariensis", 
                       CombineCol = ns[grepl("Senecio.madagascariensis", ns)])

# Vernonia.cinerea
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Vernonia.cinerea", 
                       CombineCol = ns[grepl("Vernonia", ns)])

# Viola.betonicifolia
veg.2015 <- OrgSpp(df = veg.2015, 
                       KeepCol = "Viola.betonicifolia", 
                       CombineCol = ns[grepl("Viola", ns)])

# add ring, plot, pos, year
splt <- strsplit(veg.2015$position, "[.]")
veg.2015$ring <- factor(sapply(splt, "[", 1))
veg.2015$plot <- factor(sapply(splt, "[", 2))
veg.2015$position <- factor(sapply(splt, "[", 3))
veg.2015$year <- factor("2015")

# sort colmuns
NotSpp <- c("year", "ring", "plot", "position", "cell")
Spp <- sort(names(veg.2015)[which(!(names(veg.2015) %in% NotSpp))])
veg.2015 <- veg.2015[c(NotSpp, Spp)]

# save
save(veg.2015, file = "output//Data/FAVE_vegetation2015.RData")

######################################
# combine with previous year dataset #
######################################
load("output/Data/FACE_Vegetation_Raw.RData")
vdf <- rbind.fill(veg.face, veg.2015)

# turn na into 0
vdf[is.na(vdf)] <- 0

# organise spp----

# Solanum.nigrum
vdf <- OrgSpp(df = vdf,
              KeepCol = "Solanum.nigrum",
              CombineCol = names(vdf)[grepl("Solanum", names(vdf))])

# Paspalum.dilatatum
vdf <- OrgSpp(df = vdf,
              KeepCol = "Paspalum.dilatatum",
              CombineCol = names(vdf)[grepl("Paspalum", names(vdf))])

## Tricoryne smplix -> Tricoryne elatior?
# http://www.environment.nsw.gov.au/determinations/cumber landplainpd.htm
vdf <- OrgSpp(df = vdf,
              KeepCol = "Tricoryne.elatior",
              CombineCol = names(vdf)[grepl("Tricoryne", names(vdf))])

## uknonw spp ##
unknownSp <- vdf[, grepl("unknown", names(vdf), ignore.case = TRUE)]
colSums(unknownSp)

# very samll proportaion, so remove
vdf <- vdf[, !grepl("unknown", names(vdf), ignore.case = TRUE)]

# sort colmuns----
Spp <- sort(names(vdf)[which(!(names(vdf) %in% NotSpp))])
vdf <- vdf[c(NotSpp, Spp)]

all(vdf[Spp] <= 1)

# save
save(vdf, file = "output//Data/FACE_Vegetation_Raw_2013_2015.RData")

# Spp list
vdf.mlt <- melt(vdf, id = c("year", "ring", "plot", "position", "cell"))

spp <- data.frame(sp = sort(levels(vdf.mlt$variable)))
write.csv(spp, file = "output/Data/spp_2015.csv", row.names = FALSE)

# Spp which were found only in 2015
YearSum <- ddply(vdf, .(year), function(x) colSums(x[Spp]))
newSp <- names(YearSum)[apply(YearSum, 2, function(x) all(c(x[1:2] == 0), x[3] != 0))]
YearSum[newSp] # very small..
