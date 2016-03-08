## ---- Process2016Data

######################
# Organise 2016 data #
######################

# produce sheetname
a <- as.vector(outer(1:6, 1:4, paste, sep = "."))
shts <- as.vector(outer(a, LETTERS[1:4], paste, sep = "."))

# raed all files
options(java.parameters = "-Xmx100m") # increase java memory
veg.2016.raw <- ldply(shts, function(x){ 
  read.veg.xlx(file = "Data/Result_FACE_Vegetation_Datasheet_2016.xlsx", 
               sheetName = x)}, 
  .progress = "text")

Spp <- names(veg.2016.raw)[!names(veg.2016.raw) %in% c("cell", "position")]

# Turn Spp into numeric
veg.2016 <- veg.2016.raw
veg.2016[, Spp] <- apply(veg.2016[, Spp], 2, as.numeric)

# turn NA into into 0
veg.2016[is.na(veg.2016)] <- 0

# remove spp which were not observed
veg.2016 <- cbind(veg.2016[c("cell", "position")], 
                      veg.2016[Spp][colSums(veg.2016[Spp]) != 0])

# add ring, plot, pos, year
splt <- strsplit(veg.2016$position, "[.]")
veg.2016$ring <- factor(sapply(splt, "[", 1))
veg.2016$plot <- factor(sapply(splt, "[", 2))
veg.2016$position <- factor(sapply(splt, "[", 3))
veg.2016$year <- factor("2016")

# organise species
ns <- names(veg.2016)
sort(ns)

# Asteraceace
summary(veg.2016$Asteraceace)
veg.2016$Asteraceace[which(veg.2016$Asteraceace > 0)]
  # only one observation, so remove
veg.2016$Asteraceace <- NULL

# remove bare soil
veg.2016$Bare.soil <- NULL

# Centella.asiatica
veg.2016 <- OrgSpp(df = veg.2016, 
                       KeepCol = "Centella.asiatica", 
                       CombineCol = ns[grepl("Centella", ns)])

# cinnamomum.camphora
veg.2016 <- OrgSpp(df = veg.2016, 
                       KeepCol = "Cinnamomum.camphora", 
                       CombineCol = ns[grepl("^cin|^Cin", ns)])

# Cynodon.dactylon
veg.2016 <- OrgSpp(df = veg.2016, 
                       KeepCol = "Cynodon.dactylon", 
                       CombineCol = ns[grepl("Cynodon.dactylon", ns)])

# Dichondra.repens
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Dichondra.repens", 
                   CombineCol = ns[grepl("Dichondra|dirhondria", ns)])

# Echinopogon.caespitosus
veg.2016 <- OrgSpp(df = veg.2016, 
                       KeepCol = "Echinopogon.caespitosus", 
                       CombineCol = ns[grepl("Echinopogon", ns)])

# Eragrostis.benthamii
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Eragrostis.benthamii", 
                   CombineCol = ns[grepl("benthamii", ns)])

# Glycine.like
summary(veg.2016$Glycine.like)
veg.2016$Glycine.like[which(veg.2016$Glycine.like > 0)]
 # only one observation, so remove
veg.2016$Glycine.like <- NULL

# Hydrocotyle.peduncularis
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Hydrocotyle.peduncularis", 
                   CombineCol = ns[grepl("Hydrocotyle", ns)])

# Hypoxis.hygrometrica
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Hypoxis.hygrometrica", 
                   CombineCol = ns[grepl("Hypox", ns, ignore.case = TRUE)])

# Juncus.continuus
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Juncus.continuus", 
                   CombineCol = ns[grepl("Juncus", ns, ignore.case = TRUE)])

# Lantana.camara
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Lantana.camara", 
                   CombineCol = ns[grepl("^Lant", ns, ignore.case = TRUE)])

# Leontodon.taraxacoides
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Leontodon.taraxacoides", 
                   CombineCol = ns[grepl("Leontodon", ns, ignore.case = TRUE)])

# Litter
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Litter", 
                   CombineCol = ns[grepl("litter", ns, ignore.case = TRUE)])

# Lomandra.sp
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Lomandra.sp", 
                   CombineCol = ns[grepl("Lom.ndra", ns, ignore.case = TRUE)])

# Parsonsia.straminea
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Parsonsia.straminea", 
                   CombineCol = ns[grepl("Parsonsia.straminea", ns, ignore.case = TRUE)])

# Phyllanthus.sp
veg.2016 <- OrgSpp(df = veg.2016, 
                       KeepCol = "Phyllanthus.sp", 
                       CombineCol = ns[grepl("^Phyl", ns)])

# Poranthera.microphylla
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Poranthera.microphylla", 
                   CombineCol = ns[grepl("Poranthera.microphylla", ns)])

# Senecio.madagascariensis
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Senecio.madagascariensis", 
                   CombineCol = ns[grepl("Senecio.madagascariensis", ns)])

# sida.like
summary(veg.2016$sida.like)
veg.2016$sida.like[which(veg.2016$sida.like > 0)]
 # Only one observation so remove
veg.2016$sida.like <- NULL

# Solanum.nigrum
veg.2016 <- OrgSpp(df = veg.2016, 
                   KeepCol = "Solanum.nigrum", 
                   CombineCol = ns[grepl("^Sol.*um", ns)])

# Wahlenbergia.gracilis <- Wahlenbergia
names(veg.2016)[which(names(veg.2016) == "Wahlenbergia")] <- "Wahlenbergia.gracilis"

# remove unknown
summary(veg.2016[, ns[grepl("^unkn", ns, ignore.case = TRUE)]])
# number of observations
apply(veg.2016[, ns[grepl("^unkn", ns, ignore.case = TRUE)]], 2, 
      function(x) length(which(x > 0)))
veg.2016 <- veg.2016[, !grepl("^unkn", names(veg.2016), ignore.case = TRUE)]

# sort colmuns
NotSpp <- c("year", "ring", "plot", "position", "cell")
Spp <- sort(names(veg.2016)[which(!(names(veg.2016) %in% NotSpp))])
veg.2016 <- veg.2016[c(NotSpp, Spp)]

# save
save(veg.2016, file = "output//Data/FACE_vegetation2016.RData")
