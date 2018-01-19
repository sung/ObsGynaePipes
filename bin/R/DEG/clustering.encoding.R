#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

myCaller="Encoding"
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
if(FALSE){
	#####################################
	# 1. based on binary classification #
	#####################################
	dummy<-samples.all; rownames(dummy)<-samples.all$SampleName;
	dummy=dummy[!is.na(dummy$PAPPA) & !is.na(dummy$GVel) & !is.na(dummy$UAD),c('Condition','PAPPA','GVel','PET','HT','UBD','UAD')]
	samples=sapply(dummy, as.numeric)
	rownames(samples)=rownames(dummy)

	cat("correlation\n")
	samples.cor<-cor(samples)
	#cor.test(samples[,c('Condition')],samples[,c('PAPPA')])
	hmcol <- colorRampPalette( brewer.pal(9, "RdBu"))(255)
	heatmap.2(samples.cor, trace="none", col = hmcol )

	cat("PCA for binary sub-type correlation matirx\n")
	# 1. row:samples, column:sub-types
	#pca<-prcomp(samples.cor, scale=TRUE)
	pca<-prcomp(samples, scale=TRUE)
	plot(pca$x, main='row(samples) column(sub-types)', cex=0.2)
	text(pca$x[,1], pca$x[,2], labels=rownames(pca$x))
	biplot(pca)
	sum((pca$sdev)^2) # equal to the number of standardised variables
					# should be ncol(x)
	sum((pca$rotation)^2)

	# 2. row:sub-types, column:samples
	pca<-prcomp(t(samples) )
	plot(pca$x, main='row(sub-types) column(samples)', cex=0.2)
	text(pca$x[,1], pca$x[,2], labels=rownames(pca$x))
	biplot(pca)
	sum((pca$sdev)^2)
	sum((pca$rotation)^2)

	#################################
	# 2. based on real measurements #
	#################################
	dummy<-samples.phen; rownames(dummy)<-samples.phen$SampleName;
	dummy=dummy[!is.na(dummy$PAPPA) & !is.na(dummy$GVel) & !is.na(dummy$UAD),c('BW','PAPPA','GVel','UBD','UAD')]
	samples<-dummy

	cat("correlation\n")
	samples.cor<-cor(samples)
	hmcol <- colorRampPalette( brewer.pal(9, "RdBu"))(255)
	heatmap.2(samples.cor, trace="none", col = hmcol )

	cat("PCA for binary sub-type correlation matirx\n")
	# 1. row:samples, column:sub-types
	#pca<-prcomp(samples.cor, scale=TRUE)
	pca<-prcomp(samples, scale=TRUE)
	plot(pca$x, main='row(samples) column(sub-types)', cex=0.2)
	text(pca$x[,1], pca$x[,2], labels=rownames(pca$x))
	biplot(pca)
	sum((pca$sdev)^2) # equal to the number of standardised variables
					# should be ncol(x)
	sum((pca$rotation)^2)

	# 2. row:sub-types, column:samples
	pca<-prcomp(t(samples) )
	plot(pca$x, main='row(sub-types) column(samples)', cex=0.2)
	text(pca$x[,1], pca$x[,2], labels=rownames(pca$x))
	biplot(pca)
	sum((pca$sdev)^2)
	sum((pca$rotation)^2)
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

	#~/scratch/results/RNA-Seq/SGA.AGA/DESeq2/ALL/filtered_toptags_deseq.padj.4.ALL.csv
	deseq.result = read.csv(paste0(deseq.dir,"/filtered_toptags_deseq.padj.4.",sampleType1,".csv"), stringsAsFactors=FALSE)
	deseq.result<-deseq.result[order(deseq.result$pvalue),]
	select<-deseq.result[deseq.result$padj<=0.1,]$ensg # Top DES (padj<=0.1)

	###########################
	# PCA plot of the samples #
	###########################
	cat("Plot PCA\n")
	plotPCA(rld, intgroup=c("Condition"))
	#plotPCA(rld, intgroup=c("PAPPA", "Condition"))
	#plotPCA(rld, intgroup=c("PET", "Condition"))
	#####################
	## Sample Outliers ##
	#####################
	#cat("Sample Outlier\n")
	#rlogMat[topVarGenes,c("88","75", "96")]
	#assay(dds)[topVarGenes,c("88","75", "96")]
	#samples.all[samples.all$SampleName %in% c("88","75", "96"),]

	##############################################################
	## 10 Paired Samples having 
	## 0) Paired, 1)PAPPA, 2)GVel, 3)HT, 4)UBD, and 5)UAD        #
	## check 'sampleType'<-'ALL'
	##############################################################
	cat("10 Paired Sample clustering...\n")
	#my.samples<-as.data.frame(colData(rld)) # all samples
	#my.samples<-as.data.frame(colData(rld)[!is.na(colData(rld)$PAPPA) & !is.na(colData(rld)$GVel) & !is.na(colData(rld)$UAD),])
	my.samples<-samples.all[!is.na(samples.all$PAPPA) & !is.na(samples.all$UAD),]
	# adjust levels
	for(i in colnames(my.samples[1:ncol(my.samples)]))(my.samples[i]<-droplevels(my.samples[i]))
	dummy<-table(my.samples$Pair)
	my.samples<-my.samples[my.samples$Pair %in% names(dummy[dummy==2]),] # Pair Number having both SGA/AGA

	rn=rownames(my.samples)
	rn.sga=rownames(my.samples[my.samples$Condition==1,]) # SGA only
	rn.aga=rownames(my.samples[my.samples$Condition==0,]) # AGA only

	##############
	## Encoding ##
	##############
	#my.samples$Code<-''
	# SGA/AGA
	my.samples$Code<-'A'
	my.samples$Code[my.samples$Condition==1]='S'
	# PAPPA
	my.samples$Code[my.samples$PAPPA==1]=paste0(my.samples$Code[my.samples$PAPPA==1],'p')
	my.samples$Code[my.samples$PAPPA==0]=paste0(my.samples$Code[my.samples$PAPPA==0],'P')
	# HT 
	my.samples$Code[my.samples$HT==1]=paste0(my.samples$Code[my.samples$HT==1],'h')
	my.samples$Code[my.samples$HT==0]=paste0(my.samples$Code[my.samples$HT==0],'H')
	# GVel
	my.samples$Code[my.samples$GVel==1]=paste0(my.samples$Code[my.samples$GVel==1],'g')
	my.samples$Code[my.samples$GVel==0]=paste0(my.samples$Code[my.samples$GVel==0],'G')
	# UBD
	my.samples$Code[my.samples$UBD==1]=paste0(my.samples$Code[my.samples$UBD==1],'b')
	my.samples$Code[my.samples$UBD==0]=paste0(my.samples$Code[my.samples$UBD==0],'B')
	# UAD
	my.samples$Code[my.samples$UAD==1]=paste0(my.samples$Code[my.samples$UAD==1],'u')
	my.samples$Code[my.samples$UAD==0]=paste0(my.samples$Code[my.samples$UAD==0],'U')

	################################
	## Heatmap of the count table ##
	################################
if(FALSE){
	cat("Heatmap of the count table...","\n")
	# Heatmaps showing the expression data of the 30 most highly expressed genes. The data is 
	# of raw counts (1st), from regularized log transformation (2nd) and from variance stabilizing transformation (3rd)
	hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
	#select <- order(rowMeans(counts(dds,normalized=T)),decreasing=TRUE)[1:30] # top 30 highly expressed genes
	heatmap.2(counts(dds,normalized=T)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10,6)) # require 'gplots'
	heatmap.2(rlogMat[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
	heatmap.2(vstMat[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
}

	################################
	# Heatmap with gene clustering #
	################################
	cat("Heatmap with gene clustering","\n")
	rn=rownames(my.samples)
	#select<- head( order( rowVars( rlogMat ), decreasing=TRUE ), 30 ) # requires 'genefilter' library
	hmcol <- colorRampPalette(brewer.pal(9, "RdBu"))(255)
	mat<-rlogMat[select,rn.sga]-rlogMat[select,rn.aga] # row: ENSG, col: Pair(SGA-AGA)
	colnames(mat) <- paste0('P',my.samples[my.samples$Condition==1,]$Pair)
	heatmap.2(mat, col = rev(hmcol), scale="row", trace="none", dendrogram="column")

	#############################################
	# Heatmap of the sample-to-sample distances #
	#############################################
	cat("Heatmap of the sample-to-sample distances\n")
	#hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
	hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
	# Individual Clustering (10 SGA & 10 AGA)
	rn=rownames(my.samples)
	sampleDists <- dist(t(rlogMat[,rn]))
	mat <- as.matrix(sampleDists) # row: samples, col: samples
	colnames(mat) <- rownames(mat)<- paste(my.samples$Code, rownames(mat), sep=" : ")
	heatmap.2(mat, trace="none", col = hmcol, margin=c(13, 13))
	# Paired Clustering (SGA-AGA)
	mat<-rlogMat[select,rn.sga]-rlogMat[select,rn.aga] # row: ENSG, col: Pair(SGA-AGA)
	sampleDists <- dist(t( mat ))
	mat <- as.matrix(sampleDists)
	colnames(mat) <- rownames(mat)<- paste(paste0(my.samples[my.samples$Condition==1,]$Code,'-',my.samples[my.samples$Condition==0,]$Code),paste0('P',my.samples[my.samples$Condition==1,]$Pair) , sep=" : ")
	heatmap.2(mat, trace="none", col = hmcol, margin=c(13, 13))

	#######
	# PCA #
	#######
	mat<-t(rlogMat[select,rn.sga]-rlogMat[select,rn.aga]) # row:  Pair(SGA-AGA), col: ENSG
	rownames(mat) <- paste0('P',my.samples[my.samples$Condition==1,]$Pair)
	pca<-prcomp(mat, scale=TRUE)
	plot(pca$x, main=paste0("PCA of ",sampleType), cex=0.1)
	text(pca$x[,1], pca$x[,2], labels=rownames(pca$x), cex=3)
	#sum((pca$sdev)^2) # should be ncol(x)
	#sum((pca$rotation)^2) # should be ncol(x)

}

dev.off()
cat("All is done\n")
