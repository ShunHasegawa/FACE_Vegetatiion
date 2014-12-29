# Here, I'm going to create a datasheet for vegetation survey at FACE in 2015
# based on 2014 result

# load previous data
load("output/Data//FACE_Vegetation_Raw.RData")

# subset 2014
veg14 <- subset(veg.face, year == 2014)
subset(veg.face, ring == 3 & plot == 1 & position == "A")[, 1:6]

# turn cell into numeric
veg14$cell <- as.numeric(veg14$cell) 

# remove unknown
veg14 <- veg14[, -grep("Unknown", names(veg14))]

NotPlntCol <- names(veg14)[1:5]
PlntCol <- names(veg14)[-1:-5]

#################
# summary table #
#################
# plot sum
plot.sum <- ddply(veg14, .(ring, plot), function(x) colSums(x[,PlntCol]))

# ring mean
ring.mean <- ddply(plot.sum, .(ring), function(x) colMeans(x[,PlntCol]))

#############################
# list of spp for each ring #
#############################
# ring sum
ringsum <- ddply(veg14, .(ring), function(x) colSums(x[,PlntCol]))

# melt
ringsum.ml <- melt(ringsum,id = "ring")

# remove spp which was not observed for each ring
RingSpp <- subsetD(ringsum.ml, value != 0)

# turn variable into character
RingSpp$variable <- as.character(RingSpp$variable)

# list of spp for each 2m x 2m plot
plotsum.ml <- melt(plot.sum,id=c("ring","plot"))

# remove spp which was not observed for each plot
PlotSpp <- subsetD(plotsum.ml, value != 0)
PlotSpp$variable <- as.character(PlotSpp$variable)

# split dataset to for each plot wihtin each ring
ring.veg <- split(veg14,list(veg14$ring, veg14$plot))

# function which extract only spp which were present in a ring
ringSpList <- function(dataset){
  ringnum <- unique(dataset$ring) # ring number
  ringSpData <- cbind(dataset[,c("ring", "plot", "position", "cell")],
                      dataset[,RingSpp[RingSpp$ring == ringnum, "variable"]])
                      # RingSpp specifies ring spp to be subsetted
  return(ringSpData)
}

# using the function above(ringSpData), make a ring sp list
ringdata <- lapply(ring.veg, ringSpList)
head(ringdata[[1]])
sapply(ringdata, ncol)
ddply(RingSpp, .(ring), function(x) nrow(x) + 4)
  # num of spp are the same. function worked fine

# replace 1->P, 0->"", if replace is done after making a plot sp list below, 
# replace 1 with P
ringdataP <- lapply(ringdata, function(x) cbind(x[,1:4],replace(x[5:ncol(x)],x[5:ncol(x)]==1,"P"))) 

# replace 0 with "" (blank)
ringdataP <- lapply(ringdataP, function(x) cbind(x[,1:4],replace(x[5:ncol(x)],x[5:ncol(x)]==0,""))) 

# reorder accoding to cell & position
ringdataP <- lapply(ringdataP, function(x) x[order(x$position,x$cell),])


# re-order colomuns: spp which were present in a 2m x 2m plot come first (left),
# and the others (spp which were in the same ring but not in that 2m x 2m plot)
# come after (right)

# make a function to do so
  # function which extracts the first letter of spp
  firstL <- function(x){
    h <- substr(names(x), 1, 1)
    h[duplicated(h)] <- "" # make duplicated characters blank
    return (h)
  }

plotSpList <- function(dataset){
  ringnum <- unique(dataset$ring)
  plotnum <- unique(dataset$plot)
  plotspname <- PlotSpp[PlotSpp$ring == ringnum & PlotSpp$plot == plotnum, "variable"]
  notPlotsp <- setdiff(names(dataset)[-1:-4], plotspname)
    # setdiff(x,y) gives components that are not in y from x.
  loc <- rbind(dataset[,1:4],NA) 
    # add another row NA to make its row length same as the following datasets. 
    # these colmuns are factors, not character, so empty cell "" can't be added
  plotSpData <- dataset[plotspname]
  plotSpData <- rbind(plotSpData, firstL(plotSpData))
    # add another row showing the first letter of spp
  NotPlotSpData <- dataset[notPlotsp]
  NotPlotSpData <- rbind(NotPlotSpData, firstL(NotPlotSpData))
  # combine
  plotSpData <- cbind(loc,plotSpData,NA,NotPlotSpData)#NA is just an empty column
  names(plotSpData)[grep("NA",names(plotSpData))] <- "" 
  # deleate "blank" NOTE: the column name for "blank" is removed, but the column
  # itself remains
  return(plotSpData)
}

# using the funciton above(plotSpList), make a plots sp list
plotdata <- lapply(ringdataP,plotSpList)

#######################################################
# save in an excel file, ring-plot-position per sheet #
#######################################################

# create a workbook
wb <- createWorkbook()

# specifiy the order of worksheet in "ival"
summary(plotdata)
ival <- c(seq(1,24,6),seq(2,24,6),seq(3,24,6),seq(4,24,6),seq(5,24,6),seq(6,24,6))

for (i in ival){
  for (j in c("A","B","C","D")){
    d <- subset(plotdata[[i]],position==j|is.na(ring))
    ringnum <- as.character(unique(d$ring))[1]#there's NA as well, so just use the first itme
    plotnum <- as.character(unique(d$plot))[1]
    
    # create a sheet in the workbook
    sheet <- xlsx::createSheet(wb, sheetName = paste(ringnum, plotnum, j, sep="."))
      # createSheet is also difined by XLConnect so specify which package to use
    
    # add data to the sheet
    addDataFrame(d[-26,c(-1:-3)],sheet,showNA=FALSE,row.names=FALSE,startRow=4)
    
    # labels for the datasheet
    labs <- vector()
    labs[c(1,16,20,25,31,41)] <- c("EucFace understory vegetation survey 2015",
                                   "Date:",
                                   paste("Ring:",ringnum,sep=" "),
                                   paste("2m^2 plot:",plotnum,sep=" "),
                                   paste("1m^2 position:",j,sep=" "),
                                   "People&Roles:")
    addDataFrame(t(labs),sheet,showNA=FALSE,row.names=FALSE,col.names=FALSE,
                 startRow=1)
    #t:transpors cause it's a row
    
    # add the first letters of spp
    addDataFrame(d[is.na(d$ring),c(-1:-3)],sheet,showNA=FALSE,row.names=FALSE,col.names=FALSE,
                 startRow=3)
  }
}

# save
xlsx::saveWorkbook(wb, "output/FACE_Vegetation_Datasheet_2015.xlsx")
  # saveWorkbook is also difined by XLConnect so specify which package to use
