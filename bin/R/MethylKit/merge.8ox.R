#!/usr/bin/Rscript --vanilla
#http://code.google.com/p/methylkit

source ("~/Pipelines/config/methylkit.R") # load config

load (paste0(methylkit_dir,"/",my_prefix,".filtered.myobj.RData"))

# merge all samples to one table by using base-pair locations that are covered in all samples
# setting destrand=TRUE, will merge reads on both strans of a CpG dinucleotide. This provides better 
# coverage, but only advised when looking at CpG methylation
# 
meth.destr2=unite(filtered.myobj2,destrand=TRUE) # is a 'methylBase'. Merge on both strand
meth.destr3=unite(filtered.myobj3,destrand=TRUE) # is a 'methylBase'. Merge on both strand
meth.destr5=unite(filtered.myobj5,destrand=TRUE) # is a 'methylBase'. Merge on both strand
meth.destr8=unite(filtered.myobj8,destrand=TRUE) # is a 'methylBase'. Merge on both strand

meth.destr=unite(filtered.myobj,destrand=TRUE) # is a 'methylBase'. Merge on both strand

save(meth.destr2,meth.destr3,meth.destr5,meth.destr8,meth.destr,file=paste0(methylkit_dir,"/",my_prefix,".meth.RData")) # 8oxBS.meth.RData
