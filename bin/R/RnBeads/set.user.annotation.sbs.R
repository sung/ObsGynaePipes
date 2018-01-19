#!/usr/bin/Rscript --vanilla
#Sung Gong <sung@bio.cc>

library(RnBeads)

set.user.annotation<-function(my.annotation.file, my.annotation.RData, my.region){
	if(file.exists(my.annotation.RData)){
		# Load the annotation 
		cat(paste0("Loading User-defined Annotation ", my.region, "\n"))
		rnb.load.annotation(my.annotation.RData, my.region)
		cat(paste0("User-defined Annotation ", my.region, " Loaded\n"))
		cat("rnb.region.types:\n")
		print(rnb.region.types())
	}else{
		if(file.exists(my.annotation.file)){ # data/Annotation/hg19.enst.tiling500.csv.gz 
			cat(paste0("before adding ",my.region, ". rnb.region.types:\n"))
			print(rnb.region.types())

			cat(paste0("Reading input file: ", my.annotation.file, " \n"))
			anno.dataframe<-read.csv(my.annotation.file, header=TRUE, stringsAsFactors=TRUE ) # a dataframe

			cat(paste0("Setting RnBeads Annotation for ", my.region, " \n"))
			# Create RnBeads annotation by providing a data.frame
			rnb.set.annotation(my.region, anno.dataframe, assembly="hg19")

			# Set the options to include the new annotations
			# rnb.getOption("region.types")=NULL (by default)
			rnb.options(region.types=c(rnb.getOption("region.types"),my.region))

			cat(paste0("after adding ",my.region, ". rnb.region.types:\n"))
			print(rnb.region.types())

			# Save the annotation to disk
			cat(paste0("Saving RnBeads Annotation for ", my.region, " \n"))
			rnb.save.annotation(my.annotation.RData, my.region, assembly="hg19")
			cat(paste0("User-defined Annotation ", my.region, " Saved\n"))
		}else{
			stop(paste0(my.annotation.file, " not found. Stopped!"))
		}
	}
}

if(my.region=='all'){
	for(i in row.names(my.annotations)){
		set.user.annotation(my.annotations[i,]$csv.file, my.annotations[i,]$rdata.file, i)
	}
	rnb.options(region.types=NULL)  #‘NULL’ (default value)
									#signifies that all available region annotations (as returned
									#by ‘rnb.region.types’) are summarized upon loading and
									#normalization, and the other modules analyze all regions
									#summarized in the dataset. If this option is set to an empty
									#vector, analysis on the region level is skipped.
	rnb.options(export.types=c("sites",rnb.region.types(assembly = "hg19"))) 
	#rnb.options(export.types=NULL) # No exporting in Track
}else{
	set.user.annotation(my.annotations[my.region,]$csv.file, my.annotations[my.region,]$rdata.file, my.region)
	rnb.options(region.types=my.region)  
	rnb.options(export.types=my.region)
	rnb.options(analyze.sites=FALSE)
}
