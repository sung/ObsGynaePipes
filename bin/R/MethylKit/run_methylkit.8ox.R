#!/usr/bin/R/MethylKitscript --vanilla
#http://code.google.com/p/methylkit

#options(width=200) # text width 
					# default 80 (option()$width)

source ("~/Pipelines/config/methylkit.R") # load config

###################################
## Save MethylKit Data Structure ##
###################################
#cat("reading files...\n");
#source ("~/Pipelines/bin/R/MethylKit/read.8ox.R") # this makes RData image (8oxBS.myobj.RData)
#cat("filtering ...\n");
#source ("~/Pipelines/bin/R/MethylKit/filter.8ox.R") # this makes RData image (8oxBS.filtered.myobj.RData)
#cat("merging...\n");
##source ("~/Pipelines/bin/R/MethylKit/merge.8ox.R") # this makes RData image (8oxBS.meth.RData)
#source ("~/Pipelines/bin/R/MethylKit/merge.stranded.8ox.R") # this makes RData image (8oxBS.meth.strand.RData)
cat("diffmet...\n");
#source ("~/Pipelines/bin/R/MethylKit/diffmet.8ox.R") # this makes RData image (8oxBS.diff.RData)
source ("~/Pipelines/bin/R/MethylKit/diffmet.stranded.8ox.R") # this makes RData image (8oxBS.diff.RData)

if(FALSE){

library(biomaRt)
library(GenomicRanges)
#library(GenomicFeatures)
################
## Load RData ##
################
#load (paste0(methylkit_dir,"/",my_prefix,".myobj.RData")) # set 'myobj'
load (paste0(methylkit_dir,"/",my_prefix,".filtered.myobj.RData")) # set 'filtered.myobj'
#load (paste0(methylkit_dir,"/",my_prefix,".meth.RData")) # set 'meth.destr'
load (paste0(methylkit_dir,"/",my_prefix,".meth.strand.RData")) # set 'meth'
#load (paste0(methylkit_dir,"/",my_prefix,".diff.RData")) # set 'myDiff'
load (paste0(methylkit_dir,"/",my_prefix,".diff.strand.RData")) # set 'myDiff'

######################################
# Calculate differential methylation #
######################################
# get differentially methylated regions with 25% difference and qvalue<0.01
# difference: cutoff for absolute value of methylation change between test and control (default:25)
# qvalue: cutoff for qvalue of differential methylation statistic (default:0.01)
# type: one of the "hyper","hypo" or "all" strings.  Specifies what
#	type of differentially menthylated bases/regions should be
#	returned.  For retrieving Hyper-methylated regions/bases
#	type="hyper", for hypo-methylated type="hypo" (default:"all")
top.dmc=get.methylDiff(myDiff,difference=myDiffvalue,qvalue=myQvalue) # isa methylDiff (resulution: base)
# get differentially hypo methylated regions with 25% difference and qvalue<0.01
#top.dmc.Hypo =get.methylDiff(myDiff,difference=myDiffvalue,qvalue=myQvalue,type="hypo") 
# get differentially hyper methylated regions with 25% difference and qvalue<0.01
#top.dmc.Hyper=get.methylDiff(myDiff,difference=myDiffvalue,qvalue=myQvalue,type="hyper")

gene.obj=read.transcript.features(location=gene.file, up.flank=promoter.up, down.flank=promoter.down) # a ‘GRangesList’ containing locations of exons/introns/promoters/TSSes

###################
## Set up Biomart #
###################
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#listFilters(ensembl)  # check which filters are available
#listAttributes(ensembl) # check attributes are available to select.More information on ensembl data base

#source ("~/Pipelines/bin/R/MethylKit/8ox.genome.wide.R") # genome-wide methylated-C analysis 
source ("~/Pipelines/bin/R/MethylKit/8ox.deg.R") # genome-wide methylated-C analysis 

}
