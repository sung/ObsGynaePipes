#!/usr/bin/Rscript --vanilla
#http://code.google.com/p/methylkit

source ("~/Pipelines/config/methylkit.R") # load config

load (paste0(methylkit_dir,"/",my_prefix,".filtered.myobj.RData"))

# merge all samples to one table by using base-pair locations that are covered in all samples
# setting destrand=TRUE, will merge reads on both strans of a CpG dinucleotide. This provides better 
# coverage, but only advised when looking at CpG methylation
# 
meth2=unite(filtered.myobj2,destrand=FALSE) # is a 'methylBase'
meth3=unite(filtered.myobj3,destrand=FALSE) # is a 'methylBase'
meth5=unite(filtered.myobj5,destrand=FALSE) # is a 'methylBase'
meth8=unite(filtered.myobj8,destrand=FALSE) # is a 'methylBase'

meth=unite(filtered.myobj,destrand=FALSE) # is a 'methylBase'. 

# ~/scratch/results/Methyl-Seq/SGA.AGA/methylKit/8oxBS.meth.strand.RData
save(meth2,meth3,meth5,meth8,meth,file=paste0(methylkit_dir,"/",my_prefix,".meth.strand.RData")) #
