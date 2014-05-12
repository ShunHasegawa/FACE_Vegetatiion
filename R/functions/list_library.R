#library
#list of packages
pcs<-c("nlme","lme4","AICcmodavg","gtools","gdata","gmodels","contrast","doBy","Hmisc",
       "lubridate","car","reshape","xlsx","effects")
#function to call packages. If one hasn't been installed yet, it will install it and call
call.pcs <- function(x) if(any(.packages(all.available=TRUE)==x)) 
  library(x,character.only=TRUE) else {
    install.packages(x)
    library(package(x,character.only=TRUE))}          
for (i in 1:length(pcs)) call.pcs(pcs[i])


