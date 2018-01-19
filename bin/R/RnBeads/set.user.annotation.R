#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

library(RnBeads)

##########
# Config #
##########
enst.tiling500.file <- "~/scratch/data/Annotation/dummy.rnb.user.annotation.txt" # a path 

enst.tiling500=read.table(enst.tiling500.file, header=TRUE, stringsAsFactors=TRUE ) # a dataframe
# Create RnBeads annotation by providing a data.frame
rnb.set.annotation("enstTiling500", enst.tiling500, assembly="hg19")

# Set the options to include the new annotations
rnb.options(region.types=c(rnb.getOption("region.types"),"enstTiling500"))

# check the status
rnb.get.annotation(type="enstTiling500")

enst.tiling500.RData="~/scratch/data/Annotation/dummy.rnb.user.annotation.RData" # a path 

# Save the annotation to disk
rnb.save.annotation(enst.tiling500.RData, "enstTiling500",assembly="hg19")

# Load the annotation 
rnb.load.annotation(enst.tiling500.RData, "enstTiling500")

# Check that the annotation has been successfully loaded
rnb.region.types()

warnings() # print out the warnings
logger.completed()
logger.close()
cat("All done\n")
