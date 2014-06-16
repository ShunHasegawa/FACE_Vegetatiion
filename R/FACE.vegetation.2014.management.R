################
# process 2014 #
################
# source("R/prcss.veg.2014.R")

#############
# 2012 data #
#############
load("output/veg.12.Rdata")

# remove plant type
sp.12 <- names(veg.12)[-1:-4]
splt.12 <- strsplit(sp.12, "[.]")
new.sp.12 <- paste(sapply(splt.12, "[", 2),
                sapply(splt.12, "[", 3), sep = ".")
names(veg.12)[-1:-4] <- new.sp.12

# year
veg.12$year <- factor("2012")


############### 
# 2012 & 2014 #
###############
load("output/Data/FACE_Veg2014.RData")

veg.face <- rbind.fill(veg.12, veg.14)

######################
## organise dataset ##
######################
# turn na into 0
veg.face[is.na(veg.face)] <- 0

## Glycine ##
# Glycine: could not really identify spp so combine
GlySum <- rowSums(veg.face[, grep("Glycine", names(veg.face))], na.rm = TRUE)

# if there are values larger than 1, trun them into 1
GlySum <- ifelse(GlySum > 1, 1, GlySum)

# remove Glycine
veg.face <- veg.face[,-grep("Glycine", names(veg.face))]

# add Sum of Glycine
veg.face$Glycine.sp <- GlySum
  
## Galium ##
# Galium sp = Galium propinquum
veg.face$Galium.propinquum <- rowSums(veg.face[, c("Galium.propinquum", "Galium.sp")], na.rm = TRUE)
veg.face$Galium.sp <- NULL

## Oplismenus.aemulus ##
# Oplismenus.aemulus = Oplismenus.sp
veg.face$Oplismenus.aemulus <- rowSums(veg.face[, c("Oplismenus.aemulus", "Oplismenus.sp")], na.rm = TRUE)
veg.face$Oplismenus.sp <- NULL

## Eragrostis brownii and benthamii ##
# Therse spp are treated as the same spp in 
# Flora of the Sydney Region 
# http://ausgrass2.myspecies.info/content/eragrostis-brownii 

# Eragrostis brownii = benthamii
veg.face$Eragrostis.brownii <- rowSums(veg.face[, c("Eragrostis.brownii", "Eragrostis.benthamii")], na.rm = TRUE)
veg.face$Eragrostis.benthamii <- NULL

# if there are values larger than 1, trun them into 1
veg.face$Eragrostis.brownii <- ifelse(veg.face$Eragrostis.brownii > 1, 1, veg.face$Eragrostis.brownii)

## sort column order ##
NotSpp <- c("year", "ring", "plot", "position", "cell")
Spp <- sort(names(veg.face)[which(!(names(veg.face) %in% NotSpp))])

veg.face <- veg.face[c(NotSpp, Spp)]

save(veg.face, file = "output/Data/FACE_Vegetation_Raw.RData")

## Lachnagrostis may be Lachnagrostis filiformis ##
# http://www.environment.nsw.gov.au/determinations/cumber
# landplainpd.htm
names(veg.face)[grep("Lachnagrostis", names(veg.face))] <- "Lachnagrostis.filiformis"


## Tricoryne smplix -> Tricoryne elatior?
# http://www.environment.nsw.gov.au/determinations/cumber
# landplainpd.htm
veg.face$Tricoryne.simplex
names(veg.face)[grep("Tricoryne.simplex", names(veg.face))] <- "Tricoryne.elatior"

# Spp list
veg.face.mlt <- melt(veg.face, id = c("year", "ring", "plot", "position", "cell"))

spp <- data.frame(sp = sort(levels(veg.face.mlt$variable)))
write.csv(spp, file = "output/Data/spp.csv",row.names = FALSE)

# plant properties
spList <- read.csv("Data//FACE_Vegetation_sp.list.csv")

FACE.veg.rslt <- merge(veg.face.mlt, spList, by.x = "variable", by.y = "sp", all = TRUE)

save(FACE.veg.rslt, file = "output/Data/FACE_Vegetation_PFG.RData")

