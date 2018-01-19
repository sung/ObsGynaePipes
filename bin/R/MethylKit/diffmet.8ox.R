#!/usr/bin/Rscript --vanilla
#http://code.google.com/p/methylkit

source ("~/Pipelines/config/methylkit.R") # load config

load (paste0(methylkit_dir,"/",my_prefix,".meth.RData"))

#########################################
# 5. Calculate differential methylation #
#########################################
# calculate differential methylation p-values and q-values
myDiff.destr2=calculateDiffMeth(meth.destr2,num.cores=np) # is a 'methylDiff'
myDiff.destr3=calculateDiffMeth(meth.destr3,num.cores=np) # is a 'methylDiff'
myDiff.destr5=calculateDiffMeth(meth.destr5,num.cores=np) # is a 'methylDiff'
myDiff.destr8=calculateDiffMeth(meth.destr8,num.cores=np) # is a 'methylDiff'

myDiff.destr=calculateDiffMeth(meth.destr,num.cores=np) # is a 'methylDiff'

save(myDiff.destr2,myDiff.destr3,myDiff.destr5,myDiff.destr8,myDiff.destr,file=paste0(methylkit_dir,"/",my_prefix,".diff.RData")) # 8oxBS.meth.RData
