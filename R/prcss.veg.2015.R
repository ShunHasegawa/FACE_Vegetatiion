## ---- Process2015Data

######################
# Organise 2015 data #
######################

# produce sheetname
a <- as.vector(outer(1:6, 1:4, paste, sep = "."))
shts <- as.vector(outer(a, LETTERS[1:4], paste, sep = "."))

# raed all files
options(java.parameters = "-Xmx100m") # increase java memory
veg.2015.raw <- ldply(shts, function(x){ 
  read.veg.xlx(file = "Data/Result_FACE_Vegetation_Datasheet_2015.xlsx", 
               sheetName = x)}, 
  .progress = "text")

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
