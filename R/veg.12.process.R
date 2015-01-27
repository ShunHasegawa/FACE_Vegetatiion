########
# Forb #
########

# time1(sep'12) NOT not grass/sedges 
time1 <- read.table("Data/time1.veg.txt",header=T)
forbs<-time1[,which(substring(names(time1),1,5)!="Grass"&substring(names(time1),1,5)!="Sedge")]

# check maximum value
apply(forbs, 2, max, na.rm = TRUE)

# Forb.Parsonsia.straminea max is 2
forbs$Forb.Parsonsia.straminea[which(forbs$Forb.Parsonsia.straminea == 2)] <- 1

# Shrub.Opercularia.sp -> Shrub.Opercularia.diphylla
names(forbs)[grep("Opercularia", names(forbs))] <- "Shrub.Opercularia.diphylla"

# Wahlenbergia.sp
names(forbs)[grep("Wahlenbergia.sp", names(forbs))] <- "Forb.Wahlenbergia.gracilis"

###############
# Grass/Sedge #
###############
# time2(dec'12) for grass/sedge
time2 <- read.table("Data/time2.veg.txt",header=T)
grs<-time2[,c(1:4,which(substring(names(time2),1,5)=="Grass"|substring(names(time2),1,5)=="Sedge"))]
grs<-grs[,c(1:34,36:38)] #exclude unknown spp
length(names(grs))
grs<-grs[,order(names(grs))]

# assume that paspalidium sp was paspalidium distans at this time --> needs to be corrected later
# Paspalidium.distans = Paspalidium.radiatum
grs[,"Grass.Paspalidium.distans"]<- rowSums(grs[, grep("Paspalidium", names(grs))])
summary(grs[,"Grass.Paspalidium.distans"])  
names(grs)

# remover the combined spp
rm.sp <- c("Grass.Paspalidium", "Grass.Paspalidium.sp.","Grass.Paspalidium.radiatum", "Grass.Paspalum.dilatatum")
for (i in rm.sp) {grs[i] <- NULL}
names(grs)[grep("Paspalidium", names(grs))]

# check if there is duble count
all(grs[,"Grass.Paspalidium.distans"] <= 1)

# Cyperus.flacidus,cyperus.sp-->Cyperus.flaccidus
names(grs[grep("Cyperus", names(grs))])
grs[,"Sedge.Cyperus.flaccidus" ] <- rowSums(grs[, grep("Cyperus", names(grs))])
grs[,"Sedge.Cyperus.flacidus" ]<-NULL
grs[,"Sedge.Cyperus.sp" ]<-NULL
all(grs[,"Sedge.Cyperus.flaccidus"] <= 1)

# schoenus-->Schoenus.apogon
grs[,"Sedge.Schoenus.apogon"]<-grs[,"Sedge.Schoenus.apogon"]+grs[,"Sedge.Schoenus"]
grs[,"Sedge.Schoenus"]<-NULL
all(grs[,"Sedge.Schoenus.apogon" ] <= 1)
names(grs)

# Sedge.Carex.breviformis -> Sedge.Carex.beviculmis
grs[, "Sedge.Carex.beviculmis"] <- rowSums(grs[, grep("Carex", names(grs))])
grs[, "Sedge.Carex.breviformis"] <- NULL
summary(grs$Sedge.Carex.beviculmis)

# Sedge.Carex.beviculmis -> Sedge.Carex.breviculmis
names(grs)[grep("Sedge.Carex.beviculmis", names(grs))] <- "Sedge.Carex.breviculmis"
all(grs[,"Sedge.Carex.breviculmis"] <= 1)

# Axonopsis -> Axonopus.fissifolius
names(grs)[grep("Axonopsis", names(grs))] <- "Grass.Axonopus.fissifolius"

# Sedge.Fimbristylis -> Fimbristylis.dichotoma
names(grs)[grep("Fimbristylis", names(grs))] <- "Sedge.Fimbristylis.dichotoma"

# Oplisimenus.sp -> Oplismenus.sp
names(grs)[grep("Oplisimenus.sp", names(grs))] <- "Grass.Oplismenus.sp"

#########################
# combine time1 & time2 #
#########################
veg.12 <- merge(forbs,grs,by=c("ring","plot","position","cell"))

# remove unknown spp
veg.12 <- veg.12[, -grep("Unknown", names(veg.12))]

#remove columns which is zero
names(veg.12)
veg.12 <- veg.12[ ,!names(veg.12) %in% (names(which(colSums(veg.12[ ,5:ncol(veg.12)]) == 0))),
                 drop=FALSE]
veg.12$ring <- as.factor(veg.12$ring)
veg.12$plot<-as.factor(veg.12$plot)
veg.12$cell<-as.factor(veg.12$cell)
save(veg.12,file="output/veg.12.Rdata")#save
#write.table(veg.12,"output/vegetation.2012.csv",sep=",") ##save the summary table as csv
