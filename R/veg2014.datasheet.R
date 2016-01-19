
################################################################################
################################################################################
#                                                                              #
# Note that this script doesn't work anymore becuase some of the file and data #
# frames were modified so needs to be updated                                  #
#                                                                              #
################################################################################
################################################################################


load("output/veg.12.Rdata")
names(veg.12)
# Papalidium distans = Paspalidium radiatum, so combine them
veg.12$Grass.Paspalidium.distans <- veg.12$Grass.Paspalidium.distans+veg.12$Grass.Paspalidium.radiatum

#check if there is a value hiegher than 1 
any(veg.12$Grass.Paspalidium.distans>1)
#all(veg.12$Grass.Paspalidium.distans %in% c(0,1))

#remove Grass.Paspalidium.radiatum
veg.12  <- veg.12[,-grep("Grass.Paspalidium.radiatum",names(veg.12))]
names(veg.12)

#summary table
#plot
names(veg.12)
plot.sum <- with(veg.12,aggregate(veg.12[,5:ncol(veg.12)],list(ring=ring,plot=plot),sum))
summary(plot.sum)
head(plot.sum)

#ring
ring.mean <- with(plot.sum,aggregate(plot.sum[,3:ncol(plot.sum)],list(ring=ring),mean))
head(ring.mean)

#save as escel
library(xlsx)

#create a work book
wb <- createWorkbook()

#create a sheet
sheet1 <- createSheet(wb,sheetName="row_data")
sheet2 <- createSheet(wb,sheetName="FACE_Vegetation_2012_Plot.sum")
sheet3 <- createSheet(wb,sheetName="FACE_Vegetation_2012_Ring.mean")

#add data to the sheets
addDataFrame(veg.12,sheet1,showNA=TRUE,row.names=FALSE,startRow=1)
addDataFrame(plot.sum,sheet2,showNA=TRUE,row.names=FALSE,startRow=1)
addDataFrame(ring.mean,sheet3,showNA=TRUE,row.names=FALSE,startRow=1)

#save
saveWorkbook(wb,"table/FACE_Vegetation.xlsx")


#modify sp names

spp <- names(veg.12)[c(-4:-1)]
?strsplit
spp.sp <- strsplit(spp,"[.]") #split spp by ".", 
summary(spp.sp)
#NOTE: split is regular expression. in regular expression . means any single character, so use ["."] instead
sp1 <- sapply(spp.sp,"[",2) 
#"[" is a function which subsets an index from a dataset. the index is specified in the following argument (i.e. 2 this time) 
sp2 <- sapply(spp.sp,"[",3) 
sp2[is.na(sp2)] <- "sp" #replace NA with sp
new.spp <- paste(sp1,sp2,sep=" ") #combine and make new spp
new.spp[grep("Axon",new.spp)] <- "Axonopus fissifolius" #correct sp name

#rename spp
names(veg.12)[c(-4:-1)] <- new.spp 
names(veg.12)

#reorder columns in the alphabetical order
veg.12 <- cbind(veg.12[,1:4],veg.12[,4+order(names(veg.12)[c(-1:-4)])])
names(veg.12)

#list of spp for each ring
ringsum <- with(veg.12,aggregate(veg.12[,5:ncol(veg.12)],list(ring=ring),sum))#sum for each ring
library(reshape)
ringsum.ml <- melt(ringsum,id="ring")

ring.sp <-list() #make an empty list
for (i in unique(ringsum.ml$ring)){#select spp whose sum is not 0 for each ring
  ring.sp[[i]] <- as.character(ringsum.ml$variable[ringsum.ml$value!=0&ringsum.ml$ring==i])#treated as character rather than factor
}
ring.sp

#list of spp for each 2m x 2m plot
plotsum <- with(veg.12,aggregate(veg.12[,5:ncol(veg.12)],list(ring=ring,plot=plot),sum))#sum for each plot
plotsum.ml <- melt(plotsum,id=c("ring","plot"))
head(plotsum.ml)

plot.sp <-list() #make an empty list
for (i in unique(plotsum.ml$ring)){#make an empty list for each ring
  plot.sp[[i]] <- list()
} 

for (i in unique(plotsum.ml$ring)){ #select spp whose sum is not 0 for each plot
  for (j in unique(plotsum.ml$plot)) {
    plot.sp[[i]][[j]] <- as.character(plotsum.ml$variable[plotsum.ml$value!=0&plotsum.ml$ring==i&plotsum.ml$plot==j])
    #treated as character rather than factor
  } 
}
plot.sp

#datasheet for each ring which contains only spp which occurred in that ring
summary(ring.sp)

# Sp list for Luke


data.frame(spp =  ring.sp[[1]], a = rep("y", length(ring.sp[[1]])))

# fucntion which prodeces datafreme with sp names as row names
SpLstRng <- function(x, ring.num) {
  a <- data.frame(spp =  x, presence = rep("y", length(x)))
  names(a)[2] <- paste("Ring_", ring.num, sep = "")
  return(a)
}

# apply this to all rings
ring.sp.ls <- lapply(c(1:6), function(x) SpLstRng(ring.sp[[x]], ring.num = x))
ring.sp.ls
# merge all rings
ring.sp.mrg <- ring.sp.ls[[1]]
for (i in 2:6){
  ring.sp.mrg <- merge(ring.sp.mrg, ring.sp.ls[[i]], all = TRUE)
}

# save in excel
#create a work book
wb <- createWorkbook()

#create a sheet
sheet1 <- createSheet(wb,sheetName="Ring_Sp_List")

#add data to the sheets
addDataFrame(ring.sp.mrg,sheet1,showNA=TRUE,row.names=FALSE,startRow=1)

#save
saveWorkbook(wb,"table/FACE_Sp_List.xlsx")
#####


#split dataset individual rings and plots
ring.veg <- split(veg.12,list(veg.12$ring,veg.12$plot)) #for some reasons "drop" argument doesn't work
ring.veg <- lapply(ring.veg,droplevels) #droplevels
summary(ring.veg[[1]])

#function which extract only spp which were present in a ring
ringSpList <- function(dataset){
  ringnum <- unique(as.numeric(as.character(dataset$ring)))#numeric factor for ring number
  ringSpData <- cbind(dataset[,c(1:4)],dataset[ring.sp[[ringnum]]])#ring.sp[[ringnum]]] specifies ring spp to be subsetted
  return(ringSpData)
}

#using the function above(ringSpData), make a ring sp list
ringdata <- lapply(ring.veg,ringSpList)
lapply(ringdata,ncol)
lapply(ring.sp,function(x) length(x)+4) #num of spp are the same. function worked fine

#re-order colomuns: spp which were present in a 2m x 2m plot come first (left), 
#and the others (spp which were in the same ring but not in that 2m x 2m plot) come after (right)
summary(plot.sp)

# make a function for that
  # function which extracts the first letter of spp
  firstL <- function(x){
    h <- substr(names(x),1,1)
    h[duplicated(h)] <- "" #make duplicated empty
    return (h)
  }

plotSpList <- function(dataset){
  ringnum <- as.numeric(as.character(unique(dataset$ring)))
  plotnum <- as.numeric(as.character(unique(dataset$plot)))
  plotspname <- plot.sp[[ringnum]][[plotnum]]
  notPlotsp <- setdiff(names(dataset)[-1:-4],plotspname)#setdiff(x,y) gives components that are not in y from x.
  loc <- rbind(dataset[,1:4],NA) 
    #add another row NA to make it length same as the following datasets
    #these colmuns are factor, not character, so empty cell "" can't be added
  plotSpData <- dataset[plotspname]
  plotSpData <- rbind(plotSpData,firstL(plotSpData))#add another row showing the first letter of spp
  NotPlotSpData <- dataset[notPlotsp]
  NotPlotSpData <- rbind(NotPlotSpData,firstL(NotPlotSpData))
  #combine
  plotSpData <- cbind(loc,plotSpData,NA,NotPlotSpData)#NA is just an empty column
  names(plotSpData)[grep("NA",names(plotSpData))] <- "" 
  #deleate "blank" NOTE: the column name for "blank" is removed, but the column itself remains
  return(plotSpData)
}


#replace 1->P, 0->"", if replace is done after making a plot sp list below, 
#for some reasons, the blank column has got a column name, "var.X", so replace first
ringdataP <- lapply(ringdata,function(x) cbind(x[,1:4],replace(x[5:ncol(x)],x[5:ncol(x)]==1,"P"))) #replace 1 with P
ringdataP <- lapply(ringdataP,function(x) cbind(x[,1:4],replace(x[5:ncol(x)],x[5:ncol(x)]==0,""))) #replace 0 with "" (blank)
#reorder accoding to cell $ position
ringdataP <- lapply(ringdataP,function(x) x[order(x$position,x$cell),])

#using the funciton above(plotSpList), make a plots sp list
plotdata <- lapply(ringdataP,plotSpList)
summary(plotdata)
# save in an excel file, ring-plot-position per sheet
library(xlsx)

# create a workbook
wb <- createWorkbook()

#specifiy the order of worksheet in "ival"
ival <- c(seq(1,24,6),seq(2,24,6),seq(3,24,6),seq(4,24,6),seq(5,24,6),seq(6,24,6))

for (i in ival){
  for (j in c("A","B","C","D")){
    d <- subset(plotdata[[i]],position==j|is.na(ring))
    ringnum <- as.character(unique(d$ring))[1]#there's NA as well, so just use the first itme
    plotnum <- as.character(unique(d$plot))[1]
    
    # create a sheet in the workbook
    sheet <- createSheet(wb,sheetName=paste(ringnum,plotnum,j,sep="."))
    
    # add data to the sheet
    addDataFrame(d[-26,c(-1:-3)],sheet,showNA=FALSE,row.names=FALSE,startRow=4)
    
    # labels for the datasheet
    labs <- vector()
    labs[c(1,16,20,25,31,41)] <- c("EucFace understory vegetation survey Jan'14",
                                   "Date:",
                                   paste("Ring:",ringnum,sep=" "),
                                   paste("2m^2 plot:",plotnum,sep=" "),
                                   paste("1m^2 position:",j,sep=" "),
                                   "People&Roles:")
    addDataFrame(t(labs),sheet,showNA=FALSE,row.names=FALSE,col.names=FALSE,
                 startRow=1)
    #t:transpors cause it's a row
    
    # add the first letters of spp
    #fls <- firstL(d[,c(-1:-4)])
    addDataFrame(d[is.na(d$ring),c(-1:-3)],sheet,showNA=FALSE,row.names=FALSE,col.names=FALSE,
                 startRow=3)
  }
}

#save
saveWorkbook(wb, "datasheet.2014/FACE_Vegetation_Datasheet_2014.xlsx")

