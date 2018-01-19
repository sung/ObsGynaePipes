#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
source ("~/Pipelines/config/DEG.R") # load config

library(gplots)
library(DESeq2)
library(genefilter)
library(RColorBrewer)

if(!file.exists(cluster.dir)){dir.create(cluster.dir)}
# pdf output filename 
pdf_file_name<- paste0(cluster.dir,'/cluster.',sampleType)
pdf(file=paste(pdf_file_name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

###############################
# Phenotype correlation & PCA #
###############################
if(sampleType=="ALL.PHEN"){
	#################################
	# 2. based on real measurements #
	#################################
	my.phen<-samples.phen[samples.phen$SampleName %in% samples$SampleName,c('SampleName','Pair','Condition','BW','PAPPA','GVel','UBD','UAD')]
	my.phen<-merge(my.phen, samples[,c("SampleName","Code")], by="SampleName")
	rownames(my.phen)<-my.phen$SampleName

	rn=rownames(my.phen)
	rn.sga=rownames(my.phen[my.phen$Condition==1,]) # SGA only
	rn.aga=rownames(my.phen[my.phen$Condition==0,]) # AGA only

	################
	# 2.1. SGA-AGA #
	################
	cat("correlation\n")
	mat<-my.phen[rn.sga,c('BW','PAPPA','GVel','UBD','UAD')]-my.phen[rn.aga,c('BW','PAPPA','GVel','UBD','UAD')] # row: Sample, col: Phenotypic measurements)
	samples.cor<-cor(mat) 
	hmcol <- colorRampPalette( brewer.pal(9, "RdBu"))(255)
	heatmap.2(samples.cor, trace="none", col = hmcol )

	cat("PCA for binary sub-type correlation matirx\n") # 1. row:samples, column:sub-types
	pca<-prcomp(mat, scale=TRUE)
	mat$SampleName<-rownames(mat); mat<-merge(mat, samples[,c("SampleName","Code")], by="SampleName"); rownames(mat)<-mat$SampleName; 
	plot(pca$x, main="PCA(SGA-AGA)", cex=0.2)
	#text(pca$x[,1], pca$x[,2], labels=paste(mat[rownames(pca$x),]$SampleName, substr(mat[rownames(pca$x),]$Code,2,5), sep=":"))
	text(pca$x[,1], pca$x[,2], labels=substr(mat[rownames(pca$x),]$Code,2,6))
	biplot(pca)
	#sum((pca$sdev)^2) # equal to the number of standardised variables
					# should be ncol(x)
	#sum((pca$rotation)^2)

	##################
	# 2.2. SGA & AGA #
	##################
	mat<-my.phen[rn,c('BW','PAPPA','GVel','UBD','UAD')] 

	cat("correlation\n")
	samples.cor<-cor(mat)
	hmcol <- colorRampPalette( brewer.pal(9, "RdBu"))(255)
	heatmap.2(samples.cor, trace="none", col = hmcol )
	cat("PCA for binary sub-type correlation matirx\n")
	# 1. row:samples, column:sub-types
	pca<-prcomp(mat, scale=TRUE)
	mat$SampleName<-rownames(mat); mat<-merge(mat, samples[,c("SampleName","Code")], by="SampleName"); rownames(mat)<-mat$SampleName; 
	plot(pca$x, main="PCA(SGA & AGA)", cex=0.2)
	text(pca$x[,1], pca$x[,2], labels=mat[rownames(pca$x),]$Code)
	biplot(pca)

	# 2. row:sub-types, column:samples
	#pca<-prcomp(t(samples) )
	#plot(pca$x, main='row(sub-types) column(samples)', cex=0.2)
	#text(pca$x[,1], pca$x[,2], labels=rownames(pca$x))
	#biplot(pca)
	#sum((pca$sdev)^2)
	#sum((pca$rotation)^2)
}

# DESeq2 RData file
deseq.RData <- paste0(deseq.dir, "/deseq.",sampleType,'.RData')
if(!file.exists(deseq.RData)){
	stop("RData file not exists\n")
}else{
	###################
	## Loading RData ##
	###################
	cat("loading DESeq2 RData...\n")
	load (deseq.RData) # dds,ddsFlt,res,resFlt,rld,vsd

	#rld.nb <- rlog(dds,blind=FALSE)
	#rlogMat.nb <- assay(rld.nb)
	rlogMat <- assay(rld) #
	vstMat <- assay(vsd)

	rn=samples[,c("SampleName")]
	rn.sga=samples[samples$Condition==1,c("SampleName")] # SGA only
	rn.aga=samples[samples$Condition==0,c("SampleName")] # AGA only

	if(sampleType=="ALL.PHEN"){
		deseq.result = read.csv("~/scratch/results/RNA-Seq/SGA.AGA/DESeq2/LOW.PAPPA/filtered_toptags_deseq.padj.1.LOW.PAPPA.csv", stringsAsFactors=FALSE)
	}else{
		deseq.result = read.csv(paste0(deseq.dir,"/filtered_toptags_deseq.padj.1.",sampleType,".csv"), stringsAsFactors=FALSE)
	}
	deseq.result<-deseq.result[order(deseq.result$pvalue),]

	#select<-deseq.result$ensg # Top DES names 
	#select <- order(rowMeans(counts(dds,normalized=T)),decreasing=TRUE)[1:30] # top 30 highly expressed genes
	pct<-c(1, 5, 10, 20, 50, 100) 
	select<- head( order(rowVars(rlogMat), decreasing=TRUE ), n=round(nrow(rlogMat)*i*0.01) ) # requires 'genefilter' library

	###########################
	# PCA plot of the samples #
	###########################
	cat("Plot PCA\n")
	print(plotPCA(rld, intgroup=c("Condition")))
	print(plotPCA(rld, intgroup=c("Pair")))
	print(plotPCA(rld, intgroup=c("Condition","Pair")))
	#plotPCA(rld, intgroup=c("PAPPA"))
	#plotPCA(rld, intgroup=c("PAPPA", "Condition"))
	#plotPCA(rld, intgroup=c("PET", "Condition"))

	#####################
	## Sample Outliers ##
	#####################
	#cat("Sample Outlier\n")
	#rlogMat[topVarGenes,c("88","75", "96")]
	#assay(dds)[topVarGenes,c("88","75", "96")]
	#samples.all[samples.all$SampleName %in% c("88","75", "96"),]

	#par(mfrow=c(1,3))
	################################
	## Heatmap of the count table ##
	################################
	if(FALSE){
	cat("Heatmap of the count table...","\n")
	hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
	# Heatmaps showing the expression data of the 30 most highly expressed genes. The data is 
	# of raw counts (1st), from regularized log transformation (2nd) and from variance stabilizing transformation (3rd)
	heatmap.2(counts(dds,normalized=T)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10,6)) # require 'gplots'
	heatmap.2(rlogMat[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
	heatmap.2(vstMat[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
	}

	#par(mfrow=c(1,1))
	################################
	# Heatmap with gene clustering #
	################################
	cat("Heatmap with gene clustering (SGA * AGA)","\n")
	hmcol <- colorRampPalette(brewer.pal(9, "RdBu"))(255)
	# Individual Clustering (SGA & AGA)
	mat<-rlogMat[select,rn] # row: ENSG, col: Pair(SGA-AGA)
	if(sampleType=="ALL.PHEN"){
		colnames(mat) <- paste(samples$SampleName, samples$Code, sep=":")
	}else{
		colnames(mat) <- paste(samples$SampleName, samples$Condition, sep=":")
	}
	rownames(mat) <- deseq.result$hgnc_symbol
	heatmap.2(mat, col = rev(hmcol), scale="row", trace="none", dendrogram="column")
	# Paired Clustering (SGA-AGA)
	cat("Heatmap with gene clustering (Pair: SGA-AGA)","\n")
	hmcol <- colorRampPalette(brewer.pal(9, "RdBu"))(255)
	mat<-rlogMat[select,rn.sga]-rlogMat[select,rn.aga] # row: ENSG, col: Pair(SGA-AGA)
	if(sampleType=="ALL.PHEN"){
		colnames(mat) <- samples[samples$Condition==1,c("Code")]
	}
	#else{
	#	colnames(mat) <- rownames(mat)<-paste0('P',samples[samples$Condition==1,c("Pair")])
	#}
	rownames(mat) <- deseq.result$hgnc_symbol
	heatmap.2(mat, col = rev(hmcol), scale="row", trace="none", dendrogram="column")

	#############################################
	# Heatmap of the sample-to-sample distances #
	#############################################
	cat("Heatmap of the sample-to-sample distances\n")
	#hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
	hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
	# 1. Individual Clustering (SGA & AGA)
	sampleDists <- dist(t(rlogMat[select,rn]))
	mat <- as.matrix(sampleDists) # row: samples, col: samples
	if(sampleType=="ALL.PHEN"){
		colnames(mat) <- rownames(mat)<-paste(samples$SampleName, samples$Code, sep=":")
	}else{
		colnames(mat) <- rownames(mat)<-paste(samples$SampleName, samples$Condition, sep=":")
	}
	heatmap.2(mat, trace="none", col = hmcol, margin=c(13, 13))

	# 2. Paired Clustering (SGA-AGA)
	mat<-rlogMat[select,rn.sga]-rlogMat[select,rn.aga] # row: ENSG, col: Pair(SGA-AGA)
	sampleDists <- dist(t( mat ))
	mat <- as.matrix(sampleDists)
	if(sampleType=="ALL.PHEN"){
		colnames(mat) <- rownames(mat)<- samples[samples$Condition==1,c("Code")]
	}
	#else{
	#	colnames(mat) <- rownames(mat)<-paste0('P',samples[samples$Condition==1,c("Pair")])
	#}
	heatmap.2(mat, trace="none", col = hmcol, margin=c(13, 13))

	#######
	# PCA #
	#######
	cat("PCA based on DEG profile\n")
	# 1. Individual Clustering (SGA & AGA)
	mat<-t(rlogMat[select,rn]) # row: samples, col: ENSG
	if(sampleType=="ALL.PHEN"){
		rownames(mat)<-paste(samples$SampleName, samples$Code, sep=":")
	}else{
		rownames(mat)<-paste(samples$SampleName, samples$Condition, sep=":")
	}
	pca<-prcomp(mat, scale=TRUE)
	plot(pca$x, main=paste0("PCA of ",sampleType), cex=0.1)
	text(pca$x[,1], pca$x[,2], labels=rownames(pca$x), cex=0.8)
	# 2. Paired Clustering (SGA-AGA)
	mat<-t(rlogMat[select,rn.sga]-rlogMat[select,rn.aga]) # row:  Pair(SGA-AGA), col: ENSG
	if(sampleType=="ALL.PHEN"){
		rownames(mat)<-samples[samples$Condition==1,c("Code")]
	}
	#else{
	#	rownames(mat)<-paste0('P',samples[samples$Condition==1,c("Pair")])
	#}
	pca<-prcomp(mat, scale=TRUE)
	plot(pca$x, main=paste0("PCA of ",sampleType), cex=0.1)
	text(pca$x[,1], pca$x[,2], labels=rownames(pca$x), cex=0.8)
	#sum((pca$sdev)^2) # should be ncol(x)
	#sum((pca$rotation)^2) # should be ncol(x)
}

#for Ulla
if(FALSE){
	rlogMat<-as.data.frame(rlogMat)
	fg.rld.diff<-subset(rlogMat,select=samples[samples$Condition==1,c("SampleName")])-subset(rlogMat,select=samples[samples$Condition==0,c("SampleName")])
	#fg.rld.diff<-t(fg.rld.diff) # sample-wise
	write.csv(fg.rld.diff, "~/scratch/results/RNA-Seq/SGA.AGA/Ulla/fg.rld.diff.csv") # gene-wise

	read.counts<-as.data.frame(counts(dds))
	fg.count.diff<-subset(read.counts,select=samples[samples$Condition==1,c("SampleName")])-subset(read.counts,select=samples[samples$Condition==0,c("SampleName")])
	#fg.count.diff<-t(fg.count.diff) # sample-wise
	write.csv(fg.count.diff, "~/scratch/results/RNA-Seq/SGA.AGA/Ulla/fg.count.diff.csv") # gene-wise

	pct<-c(1, 5, 10, 20, 50) 
	makeFile<-function(i){
		select<- head( order(rowVars(rlogMat), decreasing=TRUE ), n=round(nrow(rlogMat)*i*0.01) ) # requires 'genefilter' library

		#select<- head( order( rowVars(fg.count.diff), decreasing=TRUE ), round(nrow(fg.count.diff)*i*0.01)) # requires 'genefilter' library
		write.csv(t(fg.count.diff[select,]), paste0("~/scratch/results/RNA-Seq/SGA.AGA/Ulla/fg.count.diff.top.",i,".csv")) # sample-wise

		#select<- head( order( rowVars(fg.rld.diff), decreasing=TRUE ), round(nrow(fg.rld.diff)*i*0.01)) # requires 'genefilter' library
		write.csv(t(fg.rld.diff[select,]), paste0("~/scratch/results/RNA-Seq/SGA.AGA/Ulla/fg.rld.diff.top.",i,".csv")) # sample-wise
	}
	sapply(pct, makeFile)
}

dev.off()
cat("All is done\n")
