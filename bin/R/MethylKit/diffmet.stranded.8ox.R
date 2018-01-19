#!/usr/bin/Rscript --vanilla
#http://code.google.com/p/methylkit

source ("~/Pipelines/config/methylkit.R") # load config

# ~/scratch/results/Methyl-Seq/SGA.AGA/methylKit/8oxBS.meth.strand.RData
cat("loading meth obj...\n")
load (paste0(methylkit_dir,"/",my_prefix,".meth.strand.RData")) # strand-specific (do *not* merge CpG on both strand)

#########################################
# 5. Calculate differential methylation #
#########################################
# calculate differential methylation p-values and q-values
cat("diff for meth2...\n")
myDiff2=calculateDiffMeth(meth2,num.cores=np) # is a 'methylDiff'
cat("diff for meth3...\n")
myDiff3=calculateDiffMeth(meth3,num.cores=np) # is a 'methylDiff'
cat("diff for meth5...\n")
myDiff5=calculateDiffMeth(meth5,num.cores=np) # is a 'methylDiff'
cat("diff for meth8...\n")
myDiff8=calculateDiffMeth(meth8,num.cores=np) # is a 'methylDiff'
cat("diff for meth...\n")
myDiff=calculateDiffMeth(meth,num.cores=np) # is a 'methylDiff'
cat("saving to file...\n")
save(myDiff2,myDiff3,myDiff5,myDiff8,myDiff,file=paste0(methylkit_dir,"/",my_prefix,".diff.strand.RData")) # 8oxBS.meth.RData
