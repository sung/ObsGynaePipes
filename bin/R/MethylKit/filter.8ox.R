#!/usr/bin/Rscript --vanilla
#http://code.google.com/p/methylkit

source ("~/Pipelines/config/methylkit.R") # load config
###############################################
# 2. Filtering samples based on read coverage #
###############################################
load (paste0(methylkit_dir,"/",my_prefix,".myobj.RData"))

filtered.myobj2=filterByCoverage(myobj2,lo.count=min.count,lo.perc=NULL, hi.count=NULL,hi.perc=NULL) # is a 'methylRawList' 
filtered.myobj3=filterByCoverage(myobj3,lo.count=min.count,lo.perc=NULL, hi.count=NULL,hi.perc=NULL) # is a 'methylRawList' 
filtered.myobj5=filterByCoverage(myobj5,lo.count=min.count,lo.perc=NULL, hi.count=NULL,hi.perc=NULL) # is a 'methylRawList' 
filtered.myobj8=filterByCoverage(myobj8,lo.count=min.count,lo.perc=NULL, hi.count=NULL,hi.perc=NULL) # is a 'methylRawList' 

filtered.myobj=filterByCoverage(myobj,lo.count=min.count,lo.perc=NULL, hi.count=NULL,hi.perc=NULL) # is a 'methylRawList' 

save(filtered.myobj2,filtered.myobj3,filtered.myobj5,filtered.myobj8,filtered.myobj,file=paste0(methylkit_dir,"/",my_prefix,".filtered.myobj.RData")) # 8oxBS.myobj.RData
