#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

TR_PREFIX='GRCh37' # GRCh37|GRCh38

library(RnBeads)
library(parallel)
##########
# Config #
##########
source("~/Pipelines/config/RnBeads.sbs.R")
cat("config loaded\n")

##########################
# 1. Set RnBeads Options #
##########################
if (!logger.isinitialized()) logger.start(paste0(log.file, " started"), fname = log.file)

parallel.setup(num.cores)
if (parallel.isEnabled()) parallel.getNumWorkers()

###########################
# User-defined annotation #
###########################
# function 'set.user.annotation' defined
if(grepl("FLD",my.project)){
	rnb.options(region.types='')    # disable regional analysis; or it would fail
}else{
	source("~/Pipelines/bin/R/RnBeads/set.user.annotation.sbs.R") 
}
cat("user annotation loaded\n")

#################################
# Check user-defined annotation #
#################################
#rnb.get.annotation(type=my.region)

# assembly='hg19' (default)
# filtering.coverage.threshold=5 (default)
# rnb.region.types(assembly = "hg19"): "tiling","genes","promoters","cpgislands"
# import.default.data.type="bs.bed.dir" # failed 7/Jan/2015
rnb.options(
	analysis.name=my.title, identifiers.column="SampleName", 

	import=TRUE, # make data import report
	qc=TRUE, # QC
	preprocessing=TRUE, # TRUE if importing bed files for the first time
	differential=TRUE, # differential methylation analysis

	region.aggregation="coverage.weighted", # 'mean' by default
	gz.large.files=TRUE, # compress large file
	
	qc.coverage.plots=TRUE, qc.coverage.histograms=TRUE, qc.coverage.violins=TRUE,
	filtering.high.coverage.outliers=TRUE,
	filtering.snp='no', 
	filtering.greedycut=FALSE,
	filtering.context.removal=NULL, 
	filtering.coverage.threshold=min.doc, 
	filtering.low.coverage.masking=TRUE, # methylation values for low coverage sites should be set to missing
	export.to.bed=FALSE,
	export.to.csv=FALSE,
	export.to.trackhub=NULL, # c("bigWig"), # default: c("bigBed","bigWig")

	strand.specific=TRUE, 

	disk.dump.big.matrices=TRUE, enforce.destroy.disk.dumps=TRUE,
	logging.disk=TRUE, enforce.memory.management=TRUE
	#colors.gradient= c("#56B1F7", "#132B43") # 1=~reddish (#56B1F7) 0=~bluish (#132B43)
)

if(grepl("^SGA.AGA",my.project)){
	rnb.options(
		exploratory=TRUE, # PCA/Clustering
		differential.enrichment=FALSE,  # GO analysis
		differential.comparison.columns=c("Condition"), columns.pairing=c("Condition"="Pair"),
		filtering.missing.value.quantile=0  # 0 restricts sites identified by all samples only
											# 1 (default) retains all sites
											# 0.25 allows upto two samples (out of 8; hence 0.25) missing their %met 
	)
}else if(grepl("Boy.Girl",my.project)){
	rnb.options(
		exploratory=TRUE, 
		differential.enrichment=TRUE,  # GO analysis
		differential.comparison.columns=c("FetalSex"),
		filtering.missing.value.quantile=0  # 0 restricts sites identified by all samples only
	)
	if(my.meta.file=="~/results/RnBeads/Meta/SureSelect.Boy.Girl/Set1.cpg.sample.annotation.csv"){
		rnb.options(
			min.group.size=1, # one sample only with the group 
			filtering.missing.value.quantile=0
		)
	}else if(my.meta.file=="~/results/RnBeads/Meta/SureSelect.Boy.Girl/Set2.cpg.sample.annotation.csv"){
		rnb.options(
			filtering.missing.value.quantile=0.3
		)
	}
}else if(my.project=="BS.oxBS"){
	rnb.options(
		exploratory=FALSE, 
		differential.enrichment=FALSE,  # GO analysis
		differential.comparison.columns=c("SeqType"),
		filtering.missing.value.quantile=0,  # 0 restricts sites identified by all samples only
		qc.coverage.threshold.plot=seq(5,50,5), # 5,10,15...50
		min.group.size=1
	)
}else if(my.project=="5hmC.SGA.AGA"){
	rnb.options(
		exploratory=TRUE, 
		differential.enrichment=FALSE,  # GO analysis
		differential.comparison.columns=c("Condition"), 
		filtering.missing.value.quantile=0,  # 0 restricts sites identified by all samples only
		min.group.size=1
	)
}else if(grepl("^FLD",my.project)){
	rnb.options(
		exploratory=FALSE, 
		differential.enrichment=FALSE,  # GO analysis
		differential.comparison.columns=c("Condition"), columns.pairing=c("Condition"="Pair"),
		filtering.missing.value.quantile=0 # allow this % of missing values out of the total number of samples 
	)
}else{
	stop(paste0(my.project, " not supported. Stopped!"))
}

# Check that the annotation has been successfully loaded
cat("final rnb.region.types:\n")
rnb.region.types()

#######################
# 2. import bed files #
#######################
# load previously processed rnb.set for the same region types
# always check rnb.region.types.for.analysis(rnb.set) returns the region(s) that you intend to analyse
if(exists("preprocessed.rnb.set.dir")){ # results/RnBeads/SGA.AGA/CpG/enstTiling500/reports-2015-01-16_12PM/rnbSet_preprocessed 
	cat(paste0("Importing (previously) preprocessed RnBeads data ",preprocessed.rnb.set.dir,"...\n"))
	if(file.exists(paste0(preprocessed.rnb.set.dir,".zip"))){ 
		rnb.set<-load.rnb.set(path=paste0(preprocessed.rnb.set.dir,".zip"))
	}else{
		rnb.set<-load.rnb.set(path=preprocessed.rnb.set.dir)
	}

	#"~/results/RnBeads/SGA.AGA/CpG/all/reports-2015-01-21_04_30_PM/rnbSet_unnormalized" # 4 SGA(2F+2M) / 4 AGA(2F+2M)
	if(grepl("reports-2015-01-21_04_30_PM",preprocessed.rnb.set.dir)){
		if(my.project=="AGA.Boy.Girl"){
			rnb.set <- remove.samples(rnb.set, as.character(samples(rnb.set)[grep("SGA", samples(rnb.set))]) ) # remove SGA samples
		}else if(my.project=="SGA.Boy.Girl"){
			rnb.set <- remove.samples(rnb.set, as.character(samples(rnb.set)[rnb.sample.groups(rnb.set)$Condition$AGA])) # remove AGA samples
		}else if(my.project=="SGA.AGA.Boy"){
			rnb.set <- remove.samples(rnb.set, as.character(samples(rnb.set)[rnb.sample.groups(rnb.set)$FetalSex$F]) ) # remove Girl samples
		}else if(my.project=="SGA.AGA.Girl"){
			rnb.set <- remove.samples(rnb.set, as.character(samples(rnb.set)[rnb.sample.groups(rnb.set)$FetalSex$M]) ) # remove Boy samples
		}
	}
}else{
	# This is unprocessed data
	if(file.exists(paste0(rnb.save.dir, ".zip"))){ # RnBeads/SGA.AGA/CpG/enstTiling500/CpG.rnb.set.zip
		cat("Importing rnb.set.zip file...\n")
		rnb.set <- load.rnb.set(path=paste0(rnb.save.dir, ".zip"))
	# Read original bed files
	}else{
		cat("Importing bed files...\n")
		rnb.set <- rnb.execute.import(data.source=data.source, data.type="bs.bed.dir") # isa "RnBSet"
		cat(paste0("Saving rnb.set from bed files to ",rnb.save.dir,"...\n"))
		# 'save.rdata=TRUE' within 'rnb.run.analysis' will automatically make report_dir/rnbSet_unnormalized/*.RData
		# therefore, skip below if you don't want .zip archived file
		#save.rnb.set(rnb.set, path=rnb.save.dir, archive=TRUE) 
	}
}

#############
# 3. Run it #
#############
#cat("Preprocessing...\n")
#rnb.set<-rnb.run.preprocessing(rnb.set, dir.reports=report.dir)$rnb.set

#rnb.run.exploratory(rnb.set, dir.reports=report.dir, close.report = TRUE, show.report = FALSE )# for AGA.Boy.Girl
#Error in ifelse(is.null(selected), "Total v", "V") :
#  object 'selected' not found


#cat("Running rnb.run.differential...\n")
#rnb.run.differential(rnb.set, report.dir, close.report = TRUE, show.report = FALSE )

cat("Running the main analysis...\n")
# make sure 'report.dir' does not exist
rnb.run.analysis(dir.reports=report.dir, data.source=rnb.set, data.type="rnb.set", build.index=TRUE, save.rdata=TRUE)
destroy(rnb.set) #The ff files behind an RnBeads object can be deleted completely from the hard disk by executing the destructor method:

warnings() # print out the warnings
logger.completed()
logger.close()
cat(paste0("All done for ",my.title))
