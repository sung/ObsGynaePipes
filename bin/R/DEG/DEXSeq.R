#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

myCaller='DEXSeq'
source ("~/Pipelines/config/DEG.R") # load config
suppressPackageStartupMessages( library( "DEXSeq" ) )
library("BiocParallel")
np<-4

# RData file
dexseq.RData <- paste0(dexseq.dir, "/dexseq.",sampleType,'.RData')
if(file.exists(dexseq.RData)){
	cat("loading DESeq2 RData...\n")
	load (dexseq.RData) # 
	cat("dxd and dxr loaded\n")
}else{
	#####################
	# read DEXSeq count #
	#####################
	# dex is a 'DEXSeqDataSet' object
	cat("DEXSeqDataSetFromHTSeq...\n")
	dxd = DEXSeqDataSetFromHTSeq(
		countfiles=samples$HtseqFile,
		sampleData=samples[,-2],
		design= dexseq.design,
		flattenedfile
	)

	BPPARAM = BiocParallel::MulticoreParam(workers=np)
	cat("estimateSizeFactors...\n")
	dxd = estimateSizeFactors( dxd )
	cat("estimateDispersions...\n")
	dxd = estimateDispersions( dxd, BPPARAM=BPPARAM) # formula=design(object)
	cat("testForDEU...\n")
	if(exists("formulaReducedModel")){
		dxd = testForDEU( dxd, reducedModel = formulaReducedModel, BPPARAM=BPPARAM ) # fullModel = design(object)
																					 # Internally, it calls the DESeq2 function ‘nbinomLRT’
	}else{
		dxd = testForDEU( dxd, BPPARAM=BPPARAM )                                     # fullModel = design(object), reducedModel = ~sample + exon (by default)
	}
	cat("estimateExonFoldChanges...\n")
	dxd = estimateExonFoldChanges(dxd, fitExpToVar = my.contrast, PBPARAM=BPPARAM)  # The expression values will be fitted to 
																					# this variable using the the formula ~ sample + fitExpToVar * exon
	cat("DEXSeqResults...\n")
	dxr = DEXSeqResults(dxd)

	#geneIDs(dxd)
	#sampleAnnotation(dxd)
	#colData(dxd)
	#head(rowData(dxd))
	#head(counts(dxd))
	#head(featureCounts(dxd))
	#plotDispEsts(dxd)
	cat("saving DEXSeq RData\n")
	save(dxd, dxr,file=dexseq.RData) #
}
cat("All is done\n")
