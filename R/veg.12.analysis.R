###vegetation survey
###combine Forbs(and others except grass/sedge) from Sep'12 and Grass/Sedge from Dec'12
rm(list=ls(all=TRUE))

#library
source("functions/list_library.R")

#melt
veg12_ml <-melt(veg.12,id=c("ring","plot","position","cell")) 
head(veg12_ml)
#sp sum
veg12_sum <- aggregate(value~variable,data=veg12_ml,FUN=sum)

#re-order accoding to freaquency
veg12_sum <- veg12_sum[order(veg12_sum$value,decreasing=TRUE),]
save(veg12_sum,file="output/veg12_sum.R")
#load("output/veg12_sum.R")
summary(veg12_sum)
barplot(veg12_sum$value)
abline(h=sum(veg12_sum$value)*0.05)
summary(veg12_sum)
head(veg12_sum)

ac.val <- vector()
for (i in 1:length(veg12_sum$value)){
  ac.val[i] <- sum(veg12_sum$value[1:i])
}
head(ac.val)
plot(ac.val)
abline(h=sum(veg12_sum$value)*0.95)
veg12_sum$variable[ac.val<sum(veg12_sum$value)*0.97]

subset(veg12_sum,ac.val>sum(veg12_sum$value)*0.97)





sum(veg12_sum$value[veg12_sum$value > sum(veg12_sum$value)*0.05])/sum(veg12_sum$value)

#================================
#summerise the data in plot level
#####
#load("output/veg.12.R") #call file if necessary 
summary(veg.12)
names(veg.12)
length(names(veg.12))
plots <-with(veg.12,aggregate(veg.12[,5:length(names(veg.12))],list("ring"=ring,"plot"=plot),sum))
save(plots,file="output/plots.R")
names(plots)

#####
##summerise the data in ring level
names(plots)
rings<-with(plots,aggregate(plots[,3:length(names(plots))],list("ring"=ring),mean))
names(rings)
save(rings,file="output/rings.R")

#extract veg data
veg1<-plots[,c(3:72)]
veg1<-veg1[,order(names(veg1))]
names(veg1)
sites<-plots[,c("ring","plot")]
sites$ring<-as.factor(sites$ring)
sites$plot<-as.factor(sites$plot)

#####
#pca with vegan package
library(vegan)
#####
model2 <-rda(veg1,scale=TRUE)
summary(model2)
plot(model2,scaling=3)
plot(model2,scaling=1)

xs <-barplot(summary(model2)$sp[,1],horiz=TRUE,axisnames=F)
axis(2,at=xs,lab=names(veg1),las=1,line=-27,cex.axis=0.5)
xs <-barplot(summary(model2)$sp[,2],horiz=TRUE,axisnames=F)
axis(2,at=xs,lab=names(veg1),las=1,line=-27,cex.axis=0.5)

str(veg1)

par(mfrow=c(1,2))
#plot pc1 vs pc2
pc1 <-summary(model2)$site[,1]
pc2 <-summary(model2)$site[,2]
data.pca <-data.frame(plots$ring,plots$plot,pc1,pc2)
data.pca
names(data.pca)[c(1,2)] <-c("ring","plot")
with(data.pca,plot(pc1,pc2,type="n"))
with(data.pca,points(pc1[ring=="1"],pc2[ring=="1"],col="red",cex=1,pch="1"))
with(data.pca,points(pc1[ring=="2"],pc2[ring=="2"],cex=1,pch="2"))
with(data.pca,points(pc1[ring=="3"],pc2[ring=="1"],cex=1,pch="3"))
with(data.pca,points(pc1[ring=="4"],pc2[ring=="4"],col="red",cex=1,pch="4"))
with(data.pca,points(pc1[ring=="5"],pc2[ring=="5"],col="red",cex=1,pch="5"))
with(data.pca,points(pc1[ring=="6"],pc2[ring=="6"],cex=1,pch="6"))

#plot pc1 vs pc3
pc1 <-summary(model2)$site[,1]
pc3 <-summary(model2)$site[,3]
data.pca <-data.frame(plots$ring,plots$plot,pc1,pc2)
names(data.pca)[c(1,2)] <-c("ring","plot")
with(data.pca,plot(pc1,pc3,type="n"))
with(data.pca,points(pc1[ring=="1"],pc3[ring=="1"],col="red",cex=1,pch="1"))
with(data.pca,points(pc1[ring=="2"],pc3[ring=="2"],cex=1,pch="2"))
with(data.pca,points(pc1[ring=="3"],pc3[ring=="3"],cex=1,pch="3"))
with(data.pca,points(pc1[ring=="4"],pc3[ring=="4"],col="red",cex=1,pch="4"))
with(data.pca,points(pc1[ring=="5"],pc3[ring=="5"],col="red",cex=1,pch="5"))
with(data.pca,points(pc1[ring=="6"],pc3[ring=="6"],cex=1,pch="6"))

#boxplot of pc1
boxplot(pc1~sites$ring)
range(pc1)
boxplot(log(pc1+2)~sites$ring)
#####

#corespondence analysis
#####
library(vegan)
mod3<-cca(veg1)
summary(mod3)
plot(mod3)
summary(mod3)
names(summary(mod3))

#####
#fit environmental variables
#####
env<-read.table("env.txt",header=T)
env$ring<-as.factor(env$ring)
nrow(env)
summary(env)
env.r<-apply(env,2,function(x) rep(x,each=4)) #use the same env variables for plots within each ring
#it is analysed suing permutation so there's no need to worry about psudo-replication
summary(env.r)
str(env.r)
head(env.r)

ef<-envfit(mod3,env.r[,2:ncol(env)],permu=999)
ef
#####

#extract sp score (ca1,2,3) and save the table
ca.sp<-summary(mod3)$species[,c(1:3)]
write.table(ca.sp,"vegetation.2012.ca.species_score.csv",sep=",")

#plot results
ca1 <-summary(mod3)$site[,1]
ca2 <-summary(mod3)$site[,2]
data.ca <-data.frame(plots$ring,plots$plot,ca1,ca2)
names(data.ca)[c(1,2)] <-c("ring","plot")
with(data.ca,plot(ca1,ca2,type="n"))
with(data.ca,points(ca1[ring=="1"],ca2[ring=="1"],col="red",cex=1,pch="1"))
with(data.ca,points(ca1[ring=="2"],ca2[ring=="2"],cex=1,pch="2"))
with(data.ca,points(ca1[ring=="3"],ca2[ring=="3"],cex=1,pch="3"))
with(data.ca,points(ca1[ring=="4"],ca2[ring=="4"],col="red",cex=1,pch="4"))
with(data.ca,points(ca1[ring=="5"],ca2[ring=="5"],col="red",cex=1,pch="5"))
with(data.ca,points(ca1[ring=="6"],ca2[ring=="6"],cex=1,pch="6"))

##biplot spp
names(summary(mod3))

sptp<-c("Fern","Forb","Gras","Lich","Moss","Sedg","Shru","Tree")
sp.lgth<-NULL
for (i in 1:length(sptp)) {
  sp.lgth<-c(sp.lgth,length(which(substring(names(veg1),1,4)==sptp[i])))
}
sp.lgth

sptp.ful<-c("Fern","Forb","Grass","Lichen","Moss","Sedge","Shrub","Tree")
sp.type<-NULL
for (i in 1:length(sptp.ful)){
  sp.type<-c(sp.type,rep(sptp.ful[i],sp.lgth[i]))
}
stype <-as.factor(levels(as.factor(sp.type)))

sp.ca<-data.frame(as.vector(names(veg1)),sp.type,summary(mod3)$species)
names(sp.ca)
names(sp.ca)[1]<-"spp"

png(file="Veg2012.ca2.with_env.png",height=600,width=1200)
##spp
par(mfrow=c(1,2),mar=c(5,5,0.5,0.5), oma=c(0,0,0,0))
with(sp.ca,(plot(CA1,CA2,type="n",xlab="",ylab="")))
for (i in 1:length(stype))
{
  with(sp.ca,points(CA1[sp.type==stype[i]],
                    CA2[sp.type==stype[i]],
                    cex=2,pch=as.numeric(stype[i]),
                    lwd=2,col=i))
}
mtext(1,text="CA1",cex=1.5,line=3)
mtext(2,text="CA2",cex=1.5,line=3)
legend("topright",pch=c(1:8),legend=stype,col=c(1:8),cex=2,box.lty=2,pt.lwd=2,ncol=3)
##ring
ca1 <-summary(mod3)$site[,1]
ca2 <-summary(mod3)$site[,2]
data.ca <-data.frame(plots$ring,plots$plot,ca1,ca2)
names(data.ca)[c(1,2)] <-c("ring","plot")
with(data.ca,(plot(ca1,ca2,type="n",xlab="",ylab="")))
for (i in 1:6)
{
  with(data.ca,points(ca1[ring==i],ca2[ring==i],cex=2,pch=as.character(i)))
}
mtext(1,text="CA1",cex=1.5,line=3)
mtext(2,text="CA2",cex=1.5,line=3)


#fit environmental variables to CA model
env<-read.table("env.txt",header=T)
ncol(env)
mod<-cca(rings[,2:ncol(rings)])
ef<-envfit(mod,env[,2:ncol(env)],permu=999)
ef
names(env)
env2<-env[,c(2,3,4,5,6,7,9,21:ncol(env))]
names(env2)<-c("pH","Depth_HL","Moist_5cm","MOist_HL","Moist_30*","Moist_Total","Temp","Sand",
               "Silt","Clay","EC","Al","Cu","Fe","Ma","Zn","TP","TC","TN","OM","TOC")
names(env2)
ef<-envfit(mod,env2,permu=999)
ef
plot(ef,p.max=0.5,arrow.mul=1.8,cex=1.5)

dev.off()

#summary table of env fit
ef[[1]]$pvals
ef.table<-data.frame(scores(ef,"vectors"),ef[[1]]$r,ef[[1]]$pvals)
write.table(ef.table,"ca.env.table.csv",sep=",")

ordiarrows(ef,display="vectors", label = TRUE)
ef


names(ef)

?envfit



####
####plot key spp
with(sp.ca,text(CA1[spp=="Grass.Paspalidium.distans"],
                CA2[spp=="Grass.Paspalidium.distans"],
                "Pas.di",cex=1.5,col="white"))
with(sp.ca,text(CA1[spp=="Forb.Hydrocotyle.peduncularis"],
                CA2[spp=="Forb.Hydrocotyle.peduncularis"],
                "Hyd.pe",cex=1.5,col="white"))

###key spp
spps<-summary(mod3)$species


ca1.odr<-spps[order(spps[,1]),]
ca1.odr
ca2.odr<-spps[order(spps[,2]),]
ca2.odr
###ca1vsca3

with(sp.ca,(plot(CA1,CA3,type="n",xlab="",ylab="",xlim=c(-2,2))))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black") #background of plotting region
col.label<-c("white","2","3","4","5","6","7","8")
for (i in 1:length(stype))
{
  with(sp.ca,points(CA1[sp.type==stype[i]],
                    CA3[sp.type==stype[i]],
                    cex=2,pch=as.numeric(stype[i]),col=col.label[i]))
}
mtext(1,text="CA1",cex=1.5,line=3)
mtext(2,text="CA3",cex=1.5,line=3)
legend("bottomleft",pch=c(1:8),legend=stype,col=col.label,cex=2,text.col="white")


####what about ring level
names(rings)
mod<-cca(rings2[2:71])
mod
plot(mod)
summary(mod)


?envfit

############################################
############################################
#fit ev
f.dat<-apply(veg1,1,function(x) tapply(x,list(sp.type),sum)) # mean for each form
f.dat<-t(f.dat) # transpose table for multivariate analysis
f.dat

f.cca<-cca(f.dat)
summary(f.cca)
plot(f.cca)
ca1 <-summary(f.cca)$site[,1]
ca2 <-summary(f.cca)$site[,2]
data.ca <-data.frame(plots$ring,plots$plot,ca1,ca2)
names(data.ca)[c(1,2)] <-c("ring","plot")
with(data.ca,plot(ca1,ca2,type="n"))
with(data.ca,points(ca1[ring=="1"],ca2[ring=="1"],col="red",cex=1,pch="1"))
with(data.ca,points(ca1[ring=="2"],ca2[ring=="2"],cex=1,pch="2"))
with(data.ca,points(ca1[ring=="3"],ca2[ring=="1"],cex=1,pch="3"))
with(data.ca,points(ca1[ring=="4"],ca2[ring=="4"],col="red",cex=1,pch="4"))
with(data.ca,points(ca1[ring=="5"],ca2[ring=="5"],col="red",cex=1,pch="5"))
with(data.ca,points(ca1[ring=="6"],ca2[ring=="6"],cex=1,pch="6"))

f.data<-cbind(sites,f.dat)
f.data
f.ring<-with(f.data,aggregate(f.data[,3:ncol(f.data)],list(ring=ring),mean))
f.ring
f.ring.cca<-cca(f.ring[,2:ncol(f.ring)])
ef<-envfit(f.ring.cca,env[,2:34],permu=999)
ef
plot(f.ring.cca,display="sites")
plot(ef,p.max=0.1)
plot(f.ring.cca)


ef<-envfit(mod,env[,c("an.mo_30","sum.mo_30")],permu=999)
ef

env$ring<-as.factor(env$ring)
ef<-envfit(mod,env,permu=999)



ef2<-envfit(mod,env[,2:34],permutations=999)
ef2

####dummry
env2<-data.frame(apply(env,2,function (x) rep(x,each=4)))
length(names(env2))
env2[,2:34]<-apply(env2[,2:34],2,as.numeric)
summary(env2)
mod<-cca(veg1)
ef<-envfit(mod,env2[,2:34],permu=999)
ef


####
#


##
names(veg1)
mod<-cca(veg1)
is.factor(plots$ring)
ef<-envfit(mod,plots$ring,permu=999)
ef


###plot level with dummy data
env2<-apply(env,2,function (x) rep(x,each=4))
names(plots)
mod<-cca(plots[3:72])
ef<-envfit(mod,env2[,2:length(names(env2))],permu=999)
ef




#####################no standerdise
pca1 <-prcomp(veg1)
summary(pca1)
pc1 <-pca1$x[,1]
pc2 <-pca1$x[,2]
data.pca<-data.frame(sites,pc1,pc2)

names(data.pca)[c(1,2)]<-c("ring","plot")

jpeg(file="Time1.Veg.pca.jpg",quality=100,height=600,width=1200)
par(mfrow=c(1,2),mar=c(5,5,0.5,0.5), oma=c(0,0,0,0))
##plot functional groups pc1 vs pc2
names(pca1)
pca1$rotation

names(veg1)

sptp<-c("Fern","Forb","Gras","Lich","Moss","Sedg","Shru","Tree")
sp.lgth<-NULL
for (i in 1:length(sptp)) {
  sp.lgth<-c(sp.lgth,length(which(substring(names(veg1),1,4)==sptp[i])))
}
sp.lgth

sptp.ful<-c("Fern","Forb","Grass","Lichen","Moss","Sedge","Shrub","Tree")
sp.type<-NULL
for (i in 1:length(sptp.ful)){
  sp.type<-c(sp.type,rep(sptp.ful[i],sp.lgth[i]))
}
stype <-as.factor(levels(as.factor(sp.type)))
sp.pca<-data.frame(as.vector(names(veg1)),sp.type,pca1$rotation)
names(sp.pca)
names(sp.pca)[1]<-"spp"


with(sp.pca,(plot(PC1,PC2,type="n",xlab="",ylab="",xlim=c(-0.8,0.8))))
with(sp.pca,(plot(PC1,PC2,type="n",xlab="",ylab="")))
for (i in 1:length(stype))
{
  with(sp.pca,points(PC1[sp.type==stype[i]],
                     PC2[sp.type==stype[i]],cex=2,pch=as.numeric(stype[i])))
}
mtext(1,text="PC1",cex=1.5,line=3)
mtext(2,text="PC2",cex=1.5,line=3)
legend("bottomleft",pch=c(1:8),legend=stype,cex=2)
with(sp.pca,text(PC1[spp=="Grass.Paspalidium.distans"],
                 PC2[spp=="Grass.Paspalidium.distans"],"Pas.di",cex=1.5))
with(sp.pca,text(PC1[spp=="Grass.Microlaena.stipoides"],
                 PC2[spp=="Grass.Microlaena.stipoides"],"Mic.st",cex=1.5))
with(sp.pca,text(PC1[spp=="Forb.Hydrocotyle.peduncularis"],
                 PC2[spp=="Forb.Hydrocotyle.peduncularis"],"Hyd.pe",cex=1.5))
with(sp.pca,text(PC1[spp=="Sedge.Schoenus.apogon"],
                 PC2[spp=="Sedge.Schoenus.apogon"],"Sch.ap",cex=1.5))
with(sp.pca,text(PC1[spp=="Forb.Pratia.purpurascens"],
                 PC2[spp=="Forb.Pratia.purpurascens"],"Pra.pu",cex=1.5))
with(sp.pca,text(PC1[spp=="Grass.Cynodon.dactylon"],
                 PC2[spp=="Grass.Cynodon.dactylon"],"Cyn.da",cex=1.5))

###plot ring pc1 vs pc2  
with(data.pca,(plot(pc1,pc2,type="n",xlab="",ylab="")))
for (i in 1:6)
{
  with(data.pca,points(pc1[ring==i],pc2[ring==i],cex=2,pch=as.character(i)))
}
mtext(1,text="PC1",cex=1.5,line=3)
mtext(2,text="PC2",cex=1.5,line=3)
dev.off()

round(pca1$rotation[,1:3],4)

#####
#mvabund
#####
library(mvabund)
names(plots)
veg.dat <-mvabund(veg1)
is.factor(plots$ring)
##ring
plot(veg.dat ~ plots$ring, col=as.numeric(plots$ring))

model3 <- manyglm(veg.dat ~ plots$ring, family="negative.binomial",cor.type="shrink")
result<-anova(model3,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
summary(result)
plot(model3)
result
result["uni.p"]
write.table(result["uni.p"],"vegetation.2012.mvabund.p.csv",sep=",") ##save the summary table as csv
write.table(result["uni.test"],"vegetation.2012.mvabund.test.csv",sep=",") ##save the summary table as csv



#plant form
length(sp.type)
names(veg1)
f.dat<-apply(veg1,1,function(x) tapply(x,list(sp.type),sum)) # mean for each form
f.dat<-t(f.dat) # transpose table for multivariate analysis
veg.dat<-mvabund(f.dat)
plot(veg.dat ~ plots$ring, col=as.numeric(plots$ring))
model1 <- manyglm(veg.dat ~ plots$ring, family="negative.binomial",cor.type="shrink")
result<-anova(model1,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
result
plot(model1)

#model2 <- manyglm(veg.dat ~ plots$ring,cor.type="shrink")
#result<-anova(model2,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
#plot(model2)
#model2 is not really different than model1


##use soil moisture
setwd("C:\\Users\\sh3410\\SkyDrive\\Documents\\PhD.HIE\\R\\soil.variable\\Data")
so.var<-read.csv("Daily.soil_vars.csv",header=T)
names(so.var)
so.var$TIMESTAMP<-as.Date(so.var$TIMESTAMP,"%Y-%m-%d")
an.mo.mean<-with(so.var,tapply(moist,ring,mean,na.rm=TRUE))
an.total.mo<-with(so.var,tapply(moist,ring,sum,na.rm=TRUE))
an.total.mo
an.temp.mean<-with(so.var,tapply(temp,ring,mean,na.rm=TRUE))
data.smmr<-subset(so.var,so.var$TIMESTAMP<as.Date("2013-02-28")&so.var$TIMESTAMP>as.Date("2012-12-01"))
mo.summer<-with(data.smmr,tapply(moist,ring,sum,na.rm=TRUE))
mo.summer/an.total.mo

mo.summer<-as.vector(mo.summer)
mo.sm.an<-as.vector(mo.summer/an.total.mo)
mo.an<-as.vector(an.total.mo)
temp.an<-as.vector(an.temp.mean)

###perform mvabund with environmental variables
rings2<-with(plots,aggregate(plots[,3:length(names(plots))],list("ring"=ring),sum))#use sum rather than mean to make values integlar for mvabund
names(rings2)
veg.dat <-mvabund(rings2[,2:length(names(rings))])
model1 <- manyglm(veg.dat ~ mo.summer, family="negative.binomial",cor.type="shrink")
plot(model1)
model2 <- manyglm(veg.dat ~ mo.sm.an, family="negative.binomial",cor.type="shrink")
plot(model2)
model3 <- manyglm(veg.dat ~ mo.an, family="negative.binomial",cor.type="shrink")
plot(model3)
model4 <- manyglm(veg.dat ~ temp.an, family="negative.binomial",cor.type="shrink")
plot(model4)
model5 <- manyglm(veg.dat ~ ph, family="negative.binomial",cor.type="shrink")
plot(model5)

result1<-anova(model1,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
result2<-anova(model2,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
result3<-anova(model3,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
result4<-anova(model4,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
anova(model5,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
result1
result2
result3
result4


###try with dummy data
env<-cbind(mo.summer,mo.sm.an,mo.an,temp.an,sand,silt,clay,ph)
env2<-apply(env,2,function (x) rep(x,each=4))
veg.dat <-mvabund(veg1)
env2
model1 <- manyglm(veg.dat ~ env2[,4], family="negative.binomial",cor.type="shrink")
plot(model1)

####


summary(env)
env
length(names(rings))
veg.dat <-mvabund(rings[,2:71])
mod1<-manyglm(veg.dat~nenv,family="nagativ.binomial")


summary(mod1)
mod1
summary(mod1,resamp="residual")


ef<-mvformula(veg.dat~env)
best.r.sq( ef, n.xvars= 5)
env[,c(4,5)]
## Fit a multivariate linear model:
length(names(env))
nenv<-env[,2:34]

veg.dat <-mvabund(rings2[,2:length(names(rings2))])
model1<-manyglm(veg.dat~nenv,family="negative.binomial",cor.type="shrink") ##only two environmental variable can be included in the model probably due to small sample size
anova(model1,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
summary(model1,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)


model1<-manylm(mvformula(veg.dat~nenv))
summary(model1,n.boot=500)
anova(model1,n.Boot=500)
plot(model1)


model1
result1<-
  anova(model1,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
summary(model1,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
result1
result1$uni.p[2,]

anova(model4, nBoot=500)

summary(result)





###plot main spp
names(rings)
sps<-rings[,c("ring","Grass.Paspalidium.distans","Forb.Hydrocotyle.peduncularis")]
is.factor(sps$ring)
sps
table(sps[,2:3])
max(sps[,2]+sps[,3])


jpeg(file="veg2012.main.spp.form.jpg",quality=100,height=600,width=1200)

sps<-rings[,c("ring","Grass.Paspalidium.distans","Forb.Hydrocotyle.peduncularis")]
par(mfrow=c(1,2),mar=c(4,5,0.5,0.5),oma=c(0,0,0,0))
xs<-barplot(t(sps[,2:3]),ylim=c(0,150),axes=F,axisnames=F,col=c(0,"black"))
axis(2,las=1,cex.axis=1)
axis(1,at=xs,lab=c(1:6),cex.axis=1)
labs<-c(expression(Forb~italic(Hydrocotyle~peduncularis)),
        expression(Grass~italic(Paspalidium~distans)))
legend("topright",leg=labs,
       fill=c("black",0),
       col=c("black",0),bty="n",cex=1.5)
mtext(2,text="Abundance",line=3.5,cex=1.5)
mtext(2,text="(% coverage)",line=2.2,cex=1.5)
mtext(1,text="Ring",cex=1.5,line=2.5)
box(bty="o")

#pdf
windows(4,4)
pdf(file="veg2012.main.spp.form.pdf",height=4,width=4)
sps<-rings[,c("ring","Grass.Paspalidium.distans","Forb.Hydrocotyle.peduncularis")]
par(mfrow=c(1,1),mar=c(4,6,0.5,0.5),oma=c(0,0,0,0))
xs<-barplot(t(sps[,2:3]),axes=F,axisnames=F,col=c(0,"black"),ylim=c(0,130))
axis(2,las=1,cex.axis=1)
axis(1,at=xs,lab=c(1:6),cex.axis=1)
mtext(2,text="Abundance",line=4.5,cex=1.5)
mtext(2,text="(Accumulative % coverage)",line=3.2,cex=1.5)
mtext(1,text="Ring",cex=1.5,line=2.5)
box(bty="o")
dev.off()


##barplot plant form composition
f.data<-data.frame(sites,f.dat)
names(f.data)
sites
sps<-f.data[,c("ring","Forb","Grass","Sedge","Shrub")]
is.factor(sps$ring)
means<-with(sps,aggregate(sps[,2:5],list("ring"=ring),mean))
xs<-barplot(t(means[,2:5]),ylim=c(0,350),axes=F,axisnames=F,col=c(0,"black","gray80","gray50"))
axis(2,las=1,cex.axis=1)
axis(1,at=xs,lab=c(1:6),cex.axis=1)
labs<-c("Forb P= 0.04","Grass P < 0.01","Sedge P = 0.04","Shrub P = 0.04")
legend("topright",leg=labs,
       fill=c(0,"black","gray80","gray50"),
       col=c(0,"black","gray80","gray50"),bty="n",cex=1.5)
mtext(1,text="Ring",cex=1.5,line=2.5)
box(bty="o")
dev.off()

####
sd(temp.an)
mean(temp.an)
mean(mo.an)
sd(mo.an)


####mvabund with enviromental variables
names(rings)
rings2<-4*rings[,2:length(names(rings))] #make value integer

library(mvabund)
env<-read.table("env.txt",header=T)
names(env)

env2<-as.matrix(env[,2:length(names(env))])
rings2
veg.dat <-mvabund(rings2)
model <- mvformula(veg.dat~env2 ) 


plot(mo.lm,which=1:2,col.main="red",cex=3,overlay=FALSE)




a<-NULL
for (i in 2:ncol(env)){
  model <- mvformula(veg.dat~env[,i]) 
  mo.lm<- manylm(model)
  a<-list(a,summary(mo.lm)$coefficients)
}
a

names(env)[11]

rings2/400
plot(veg.dat ~ plots$ring, col=as.numeric(plots$ring))



names(env)[11]
model1<-manyglm(veg.dat~env[,11], family="negative.binomial")
summary(model1)
anova(model1, nBoot=500)
result<-anova(model1,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
summary(result)
pval<-result$uni.p[2,]
subset(pval,pval<0.1)


model1<-manyglm(veg.dat~env[,6], family="negative.binomial")
smr<-summary(model1)
smr$coefficients[2,2]

b<-NULL
for (i in 2:ncol(env)){
  model1<-manyglm(veg.dat~env[,i], family="negative.binomial")
  smr<-summary(model1)
  b<-c(b,smr$coefficients[2,2])
}



result<-anova(model1,nBoot=500,test="wald",p.uni="adjusted",show.time=TRUE)
summary(result)
pval<-result$uni.p[2,]
subset(pval,pval<0.8)





model1<-manyglm(veg.dat~env2, family="binomial")
model1<-manyglm(veg.dat~env2, family="poisson")

apply(rings2,2,sum)

a[[1]]

names(env)
###
env

env3<-apply(env[,2:ncol(env)],2,function(x) rep(x,each=4))
env3<-as.matrix(env3)
veg.dat<-mvabund(plots[,3:length(names(plots))])

model <- mvformula(veg.dat~env3) 
mo.lm<- manylm(model)
summary(mo.lm)
plot(mo.lm,which=1:2,col.main="red",cex=3,overlay=FALSE)



model1<-manyglm(veg.dat~env3, family="negative.binomial")
plot(model1)
summary(model1)


plot(model1)
model1
summary(model1)
rings2

mo.lm<- manylm(model)
plot(model)
env2
summary(mo.lm, nBoot=500)
anova(lm.spider, nBoot=500)
plot(mo.lm)


veg.dat <-mvabund(veg1)
