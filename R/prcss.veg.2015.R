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
veg.2015.raw[Spp] <- apply(veg.2015.raw[Spp], 2, 
                           function(x) recode(x, "c('p', 'P') = 1; else = 0"))

# remove spp which were not observed
veg.2015.raw <- cbind(veg.2015.raw[c("cell", "position")], 
                      veg.2015.raw[Spp][colSums(veg.2015.raw[Spp]) != 0])

# organise species
ns <- names(veg.2015.raw)
sort(ns)

# Combine anagallis
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Anagallis.arvensis", 
                       CombineCol = ns[grepl("Anagallis", ns)])
# Arthropodium
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Arthropodium.sp", 
                       CombineCol = ns[grepl("Arthropodium", ns)])

# Axonopus
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Axonopus.fissifolius", 
                       CombineCol = ns[grepl("Axonopus", ns)])

# Bidens.pilosa
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Bidens.pilosa", 
                       CombineCol = ns[grepl("Bidens.pilosa", ns)])

# Breynia.oblongifolia
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Breynia.oblongifolia", 
                       CombineCol = ns[grepl("Breynia.oblongifolia", ns)])

# Conyza.sumatrensis
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Conyza.sumatrensis", 
                       CombineCol = ns[grepl("Conyza", ns)])

# Echinopogon.caespitosus
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Echinopogon.caespitosus", 
                       CombineCol = ns[grepl("Echinopogon", ns)])

# Eragrostis.curvula
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Eragrostis.curvula", 
                       CombineCol = ns[grepl("curvula", ns)])

# Eragrostis.leptostachya
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Eragrostis.leptostachya", 
                       CombineCol = ns[grepl("leptostachya", ns)])

# Oplismenus.aemulus
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Oplismenus.aemulus", 
                       CombineCol = ns[grepl("Oplismenus|Oplismemus", ns)])

# Paspalidium.distans
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Paspalidium.distans", 
                       CombineCol = ns[grepl("Paspalidium", ns)])

# Phyllanthus.sp
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Phyllanthus.sp", 
                       CombineCol = ns[grepl("Phyllanthus|Phyllathus", ns)])

# Pratia.purpurascens
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Pratia.purpurascens", 
                       CombineCol = ns[grepl("Pratia", ns)])

# Senecio.madagascariensis
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Senecio.madagascariensis", 
                       CombineCol = ns[grepl("Senecio.madagascariensis", ns)])

# Senecio.madagascariensis
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Vernonia.cinerea", 
                       CombineCol = ns[grepl("Vernonia", ns)])

# Viola.betonicifolia
veg.2015.raw <- OrgSpp(df = veg.2015.raw, 
                       KeepCol = "Viola.betonicifolia", 
                       CombineCol = ns[grepl("Viola", ns)])

# add ring, plot, pos, year
splt <- strsplit(veg.2015.raw$position, "[.]")
veg.2015.raw$ring <- factor(sapply(splt, "[", 1))
veg.2015.raw$plot <- factor(sapply(splt, "[", 2))
veg.2015.raw$pos <- factor(sapply(splt, "[", 3))
veg.2015.raw$year <- factor("2015")

# sort colmuns
NotSpp <- c("year", "ring", "plot", "pos", "cell", "position")
Spp <- sort(names(veg.2015.raw)[which(!(names(veg.2015.raw) %in% NotSpp))])
veg.2015.raw <- veg.2015.raw[c(NotSpp, Spp)]

# save
save(veg.2015.raw, file = "output//Data/FAVE_vegetation2015.RData")


