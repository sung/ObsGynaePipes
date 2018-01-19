#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# For detial, see http://www.nature.com/nprot/journal/v8/n9/full/nprot.2013.099.html
# complex design 14.B based on the publication above

myCaller='edgeR'
TR_PREFIX='GRCh38' # GRCh37 or GRCh38 (IA samples)
source ("~/Pipelines/config/DEG.R") # load config

# pdf output filename 
pdf_file_name<- paste0(edgeR.dir,'/edgeR.',sampleType)
pdf (file=paste(pdf_file_name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

# RData file
edgeR.RData<- paste0(edgeR.dir, "/edgeR.",sampleType,".RData")

if(file.exists(edgeR.RData)){
	cat("loading edgeR RData...\n")
	load (edgeR.RData) # d2, de, edgeR.deg,
}else{
	# 1. Load the edgeR package and use the utility function, readDGE, to read in the COUNT files created from htseq-count
	counts = readDGE(samples$HtseqFile, header=FALSE, sep="")$counts # isa 'matrix'
	colnames(counts) = samples$SampleName
	#> length(grep("^ENSG", rownames(counts)))
	#[1] 63676

	# 2. Filter weakly expressed and noninformative (e.g., non-aligned) features using a command like:
	noint = rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
	counts = counts[!noint,]
	cpms = cpm(counts)
	if(sampleType=="PAPPA.SGA"){
		keep = rowSums(cpms >1) >= nrow(samples[samples$PAPPA==1,])
	}else{
		keep = rowSums(cpms >1) >= nrow(samples)/2 
	}
	#> table(keep)
	#keep
	#FALSE  TRUE 
	#50659 13022
	counts = counts[keep,]

	# 3. Visualize and inspect the count table
	#head( counts[,order(samples$Condition)], 5 ) # print first 5 case (SGA) records
	#summary(counts)

	# 4. Create a DGEList object (edgeR's container for RNA-seq count data)
	cat("creating d...\n")

	if(sampleType=="PAPPA.SGA"){
		d = DGEList(counts=counts, group=samples$PAPPA) # isa 'DGEList'
	}else{
		d = DGEList(counts=counts, group=samples$Condition) # isa 'DGEList'
	}

	# 5. Estimate normalization factors
	# The calcNormFactors function normalizes for RNA composition by finding a set of scaling
	# factors for the library sizes that minimize the log-fold changes between the samples for most genes
	# The default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples 
	d = calcNormFactors(d)

	# 6. Inspect the relationships between samples using a multidimensional scaling (MDS) plot
	#plotMDS(d, labels=samples$SampleName, col=c("darkgreen","blue")[factor(samples$Condition)])

	# 7. Create a design matrix to specify the factors that are expected to affect expression levels:
	# see the config file (DEG.R)

	# 8. Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR)-adjusted likelihood
	cat("creating d2...\n")
	d2 = estimateGLMCommonDisp(d, design, verbose=TRUE) # is a DGEList
	d2 = estimateGLMTrendedDisp(d2, design, verbose=T)
	d2 = estimateGLMTagwiseDisp(d2, design)

	# d2$counts is same with counts
	# cpm(d2) is different from cpm(counts)
	# cpm(d2)["ENSG00000174697",]
	# cpm(counts)["ENSG00000174697",]

	#plotMeanVar(d2, show.tagwise.vars=TRUE, NBline=TRUE)
	#plotBCV(d2)

	# 9. Given the design matrix and dispersion estimates, fit a GLM to each feature:
	cat("creating f...\n")
	f = glmFit(d2, design) # is a DGEGLM

	# 10. Perform a likelihood ratio test, specifying the difference of interest (e.g. case vs. control)
	# glmLRT(glmfit, coef=ncol(glmfit$design), contrast=NULL, test="chisq")
	# test: which test (distribution) to use in calculating the p-values. Possible values are ‘"F"’ or ‘"chisq"’.
	# If ‘coef’ is used, the null hypothesis is that all the coefficients indicated by ‘coef’ are equal to zero.
	# The data frame ‘table’ contains the following columns:
	#   logFC: log2-fold change of expression between conditions being tested.
	#  logCPM: average log2-counts per million, the average taken over all libraries in ‘y’.
	#      LR: likelihood ratio statistics (only for ‘glmLRT’).
	#       F: F-statistics (only for ‘glmQFTest’).
	#  PValue: p-values.
	cat("creating de...\n")
	de = glmLRT(f)  # the column number for Condition, which is ConditionCTR in our case, in design data.frame
					# isa "DGELRT"
					# check out de$table (logFC, logCPM, LR, PValue)
					# check out de$samples
					# check out de$comparison

	#o = order(de$table$PValue)
	#cpm(d2)[o[1:10],]

	# 11. Use the topTags function to present a tabular summary of the differential expression statistics 
	# topTags(object, n=10, adjust.method="BH", sort.by="PValue")
	# topTags extracts the top DE tags in a data frame for a given pair of
	# groups, ranked by p-value or absolute log-fold change.
	# (note that topTags operates on the output of exactTest or glmLRT)
	tt = topTags(de, n=nrow(d2))    # default sort by pvalue
									# is a TopTags
									# BH is a default adjust.method
									# check out tt$comparison (should be the same with de$comparison)
	#tt2 = topTags(de, n=nrow(d2), adjust.method="bonferroni")
	#tt2 = topTags(de, n=nrow(d2), sort.by="logFC", adjust.method="bonferroni")
	edgeR.deg=tt$table

	################
	## Save RData ##
	################
	cat("saving edgeR RData","\n")
	save(d2,de,edgeR.deg, file=edgeR.RData) #  save 'de, edgeR.deg'
}

cat("Bottom 3 low count\n")
d2$counts[rownames(head(edgeR.deg[order(edgeR.deg$logCPM),],n=3)),]

# 12. Inspect the depth-adjusted reads per million for some of the top differentially expressed genes:
#nc = cpm(d2, normalized.lib.sizes=TRUE) # isa 'matrix'
#rn = rownames(edgeR.deg)
#head(nc[rn,order(samples$Condition)],5) # print first 5 entries

# 13. Create a graphical summary, such as an M (log-fold change) versus A (log-average expression) plot, 
# here showing the genes selected as differentially expressed (with a 5% false discovery rate
#table(edgeR.deg$logFC<0)
plotSmear(de, col=myDotCol, main=paste0("edgeR (",sampleType,"): ",nrow(edgeR.deg), " genes (", table(edgeR.deg$logFC<0)["TRUE"], "<0, ", table(edgeR.deg$logFC<0)["FALSE"], ">=0)" ))
abline(h=c(0), col="blue")

################################################
# conventional filtering |logFC|>1 & FDR < 0.1 #
################################################
deg = rownames(edgeR.deg[(edgeR.deg$logFC <= -1 | edgeR.deg$logFC >= 1) & edgeR.deg$FDR <= 0.1,])
plotSmear(de, de.tags=deg, col=myDotCol, main=paste0("edgeR (",sampleType,"): conventional (FC>=1 and FDR<=0.1)"))
abline(h=c(-1, 1), col="blue")

# P value distribution
hist(edgeR.deg$PValue, breaks=100, xlab='P-value', main=paste0("edgeR (",sampleType,"): P-value distribution of ", nrow(edgeR.deg), " genes"))

#############
# by Pvalue #
#############
#table(edgeR.deg$PValue <= myPValue)
deg = rownames(edgeR.deg[edgeR.deg$PValue <= myPValue, ])
minlogFC=round(min(abs(edgeR.deg[deg,]$logFC)),3)
plotSmear(de, de.tags=deg, col=myDotCol, main=paste0("edgeR (",sampleType,"): ", length(deg)," genes (", table(edgeR.deg[deg,]$logFC<0)["TRUE"], "<0, ", table(edgeR.deg[deg,]$logFC<0)["FALSE"], ">=0) having P-value <=", myPValue))
abline(h=c(-minlogFC, minlogFC), col="blue")

# FDR distribution of genes having P-value <= myPValue
hist(edgeR.deg[deg,]$FDR, breaks=10, xlab='FDR', main=paste0("edgeR (",sampleType,"): FDR distribution of ", length(deg), "genes having P-value <=", myPValue))


####################
# by Pvalue && FDR #
####################
#table(Pvalue=(edgeR.deg$PValue <= myPValue), FDR=(edgeR.deg$FDR <= myFDR))
deg = rownames(edgeR.deg[edgeR.deg$PValue <= myPValue & edgeR.deg$FDR <= myFDR, ])
# get the list of genes having P-value <= myPvalue and FDR <= myFDR
#edgeR.deg[order(edgeR.deg[deg,]$logFC),]
# get the min FC from above
minlogFC=round(min(abs(edgeR.deg[deg,]$logFC)),3)
plotSmear(de, de.tags=deg, col=myDotCol, main=paste0("edgeR (",sampleType,"): ",length(deg), " genes (", table(edgeR.deg[deg,]$logFC<0)["TRUE"], "<0, ", table(edgeR.deg[deg,]$logFC <0 )["FALSE"], ">=0) of FDR <= ", myFDR, " (logFC>=",minlogFC,")"))
abline(h=c(-minlogFC, minlogFC), col="blue")

################
## DEG q<=0.4 ##
################
cat("FDR<=0.4\n")
edgeR.top.deg=edgeR.deg[edgeR.deg$FDR <= 0.4,] # if you'd like to cutoff by top 100 genes
getTopEdgeR(edgeR.top.deg,"FDR.4")
################
## DEG q<=0.3 ##
################
cat("FDR<=0.3\n")
edgeR.top.deg=edgeR.deg[edgeR.deg$FDR <= 0.3,] # if you'd like to cutoff by top 100 genes
getTopEdgeR(edgeR.top.deg,"FDR.3")
################
## DEG q<=0.2 ##
################
cat("FDR<=0.2\n")
edgeR.top.deg=edgeR.deg[edgeR.deg$FDR <= 0.2,] # if you'd like to cutoff by top 100 genes
getTopEdgeR(edgeR.top.deg,"FDR.2")
################
## DEG q<=0.1 ##
################
cat("FDR<=0.1\n")
edgeR.top.deg=edgeR.deg[edgeR.deg$FDR <= 0.1,] # if you'd like to cutoff by top 100 genes
getTopEdgeR(edgeR.top.deg,"FDR.1")
################
## DEG q<=0.05 ##
################
cat("FDR<=0.05\n")
edgeR.top.deg=edgeR.deg[edgeR.deg$FDR <= 0.05,] # if you'd like to cutoff by top 100 genes
getTopEdgeR(edgeR.top.deg,"FDR.05")
###################
## Top 100 Genes ##
###################
cat("Top100 genes (meanBase>=20) by p-value\n")
select=rownames( head(edgeR.deg[order(edgeR.deg$PValue),],n=100) ) # if you'd like to cutoff by top 100 genes
edgeR.top.deg=edgeR.deg[select,]
getTopEdgeR(edgeR.top.deg,"top100")
###################
## Top 200 Genes ##
###################
cat("Top200 genes (meanBase>=20) by p-value\n")
select=rownames( head(edgeR.deg[order(edgeR.deg$PValue),],n=200) ) # if you'd like to cutoff by top 200 genes
edgeR.top.deg=edgeR.deg[select,]
getTopEdgeR(edgeR.top.deg,"top200")

###########
# save it #
###########
write.csv(edgeR.deg, file=paste0(edgeR.dir,"/toptags_all_edgeR.",sampleType,".csv"))

# check the count for a specific gene
#mean(counts[top.deg.anno[top.deg.anno$hgnc_symbol=="COL12A1",]$ensg, samples$PAPPA=='Normal'])
#mean(counts[top.deg.anno[top.deg.anno$hgnc_symbol=="COL12A1",]$ensg, samples$PAPPA=='Low'])


#################
# scatterplot3D #
#################
#scatterplot3d(y=top.deg.anno$FDR, z=top.deg.anno$logFC, x=top.deg.anno$logCPM)

##########
# plot3D #
##########
#library(rgl)
#plot3d(y=top.deg.anno$FDR, z=top.deg.anno$logFC, x=top.deg.anno$logCPM, ylab='FDR', zlab='logFC', xlab='logCPM')
#text3d(y=top.deg.anno$FDR, z=top.deg.anno$logFC, x=top.deg.anno$logCPM, cex=2, top.deg.anno$hgnc_symbol, col=top.deg.anno$col)

# Pheno
#top.deg.pheno.anno=merge(top.deg, top.deg.pheno, by.x="ensg", by.y="ensembl_gene_id")
#write.csv(top.deg.pheno.anno[order(top.deg.pheno.anno$logFC,decreasing=TRUE),], file=paste0(edgeR.dir,'/filtered_pheno_toptags_pair_edgeR.',sampleType,'.csv'))

# GO
#top.deg.go.anno=merge(top.deg, top.deg.go, by.x="ensg", by.y="ensembl_gene_id")
#write.csv(top.deg.go.anno[order(top.deg.go.anno$logFC,decreasing=TRUE),], file=paste0(edgeR.dir,'/filtered_go_toptags_pair_edgeR.',sampleType,'.csv'))

dev.off()
cat("All is done\n")
