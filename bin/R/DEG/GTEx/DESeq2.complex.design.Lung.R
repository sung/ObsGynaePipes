#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# For detial, see http://www.nature.com/nprot/journal/v8/n9/full/nprot.2013.099.html
# complex design 14.D based on the publication above
# make sure to change my.local.db from config/Annotation.R

myCaller='DESeq2'
TR_PREFIX='GRCh37' # GRCh37 | GRCh38 

myProject="GTEx"
#sampleType='Adipose_Tissue'      # Boy vs. Girl 
#sampleType='Adrenal_Gland'      # Boy vs. Girl 
#sampleType='Blood'      # Boy vs. Girl 
#sampleType='Blood_Vessel'      # Boy vs. Girl 
#sampleType='Brain'      # Boy vs. Girl 
#sampleType='Colon'      # Boy vs. Girl 
#sampleType='Esophagus'      # Boy vs. Girl 
#sampleType='Heart'      # Boy vs. Girl 
#sampleType='Liver'      # Boy vs. Girl 
sampleType='Lung'      # Boy vs. Girl 
#sampleType='Pancreas'      # Boy vs. Girl 
#sampleType='Pituitary'      # Boy vs. Girl 
#sampleType='Small_Intestine'      # Boy vs. Girl 
#sampleType='Spleen'      # Boy vs. Girl 
#sampleType='Thyroid'      # Boy vs. Girl 

source ("~/Pipelines/config/DEG.R") # load config

# Set RData file
deseq.RData <- paste0(deseq.dir, "/deseq.",sampleType,'.RData')
# Load RData if exists
if(file.exists(deseq.RData)){
	cat("loading DESeq2 RData...\n")
	load (deseq.RData) # 
	cat("dds, res, rld?, ddsExonicGene? loaded\n")
# create RData for the first time 
}else{
	####################
	# read HTSeq count #
	####################
	print(system.time(ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = samples, directory="", design= deseq.design ))) # isa 'DESeqDataSet'
	##########################
	# Filter genes from chrY #
	##########################
	if(grepl("^Boy.Girl",myProject) | grepl("^GTEx",myProject)){
		if(grepl("^Boy.Girl",myProject)){
			cat("1. collapsing replicates...\n")
			# Error: sum(assay(object)) == sum(assay(collapsed)) is not TRUE (R version 3.1.1 (2014-07-10), DESeq2_1.6.3)
			# try with R 3.2.1 (DESeq2_1.10.1)
			# or collapseReplicatesBugFix
			#print(system.time(ddsCollapsed <- collapseReplicates(ddsHTSeq, groupby=colData(ddsHTSeq)$CRN, run=colData(ddsHTSeq)$Source))) # isa 'DESeqDataSet'
			print(system.time(ddsCollapsed <- collapseReplicatesBugFix(ddsHTSeq, groupby=colData(ddsHTSeq)$CRN, run=colData(ddsHTSeq)$Source))) # isa 'DESeqDat# 'groupby' will be rownames of colData(dds)
																											# 'groupby' will be rownames of colData(dds)
																											# 'run' will be join by ','
		}else{
			ddsCollapsed<-ddsHTSeq
		}
		# Boy.Girl.exon.GRCh38 (or GRCh37)
		if(grepl("exon",myProject)|grepl("iRNA",myProject)){ 
			ddsKeep <- rep(TRUE, nrow(ddsCollapsed))
		}else{
			cat("2. removing chrY...\n")
			my.anno <- getBM(attributes = c("chromosome_name", "ensembl_gene_id"), filters = my.filter, values = rownames(ddsCollapsed), mart = myMart)
			ddsKeep <- rownames(ddsCollapsed) %in% my.anno[my.anno$chromosome_name!='Y',c(my.filter)] # remove ensg of chrY
		}

		cat("\nDESeq from filtered count data...\n")
		newColData<- colData(ddsCollapsed); rownames(newColData)<-newColData$SampleId # double-check SampleId is available?
		newCounts <- counts(ddsCollapsed)[ddsKeep,]
		ddsMatrix <- DESeqDataSetFromMatrix(newCounts, newColData, deseq.design) # isa 'DESeqDataSet'

		# every gene contains at least one zero (mainly for exon-count analysis)
		# https://support.bioconductor.org/p/63229/
		if(grepl("exon",myProject)){
			if(any(rowSums(counts(ddsMatrix))==0)){
				geoMeans = apply(counts(ddsMatrix), 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
				ddsMatrix = estimateSizeFactors(ddsMatrix, geoMeans=geoMeans)
			}
		}
		print(system.time(dds<- DESeq(ddsMatrix, parallel=TRUE))) # isa 'DESeqDataSet'
	}else if(grepl("^RoadMap.exon",myProject)){ 
		cat("1. collapsing replicates...\n")
		print(system.time(ddsCollapsed <- collapseReplicatesBugFix(ddsHTSeq, groupby=colData(ddsHTSeq)$Tissue))) # isa 'DESeqDat# 'groupby' will be rownames of colData(dds)
		ddsKeep <- rep(TRUE, nrow(ddsCollapsed))
		cat("\nDESeq from filtered count data...\n")
		newColData<- colData(ddsCollapsed); rownames(newColData)<-newColData$SampleId # double-check SampleId is available?
		newCounts <- counts(ddsCollapsed)[ddsKeep,]
		ddsMatrix <- DESeqDataSetFromMatrix(newCounts, newColData, deseq.design) # isa 'DESeqDataSet'
		# every gene contains at least one zero (mainly for exon-count analysis)
		# https://support.bioconductor.org/p/63229/
		if(any(rowSums(counts(ddsMatrix))==0)){
			geoMeans = apply(counts(ddsMatrix), 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
			ddsMatrix = estimateSizeFactors(ddsMatrix, geoMeans=geoMeans)
		}
		print(system.time(dds<- DESeq(ddsMatrix))) # isa 'DESeqDataSet'
	}else{
		#make sure the base level
		#ddsHTSeq$PAPPA <- relevel(ddsHTSeq$PAPPA,"Normal")
		cat("main DESeq2 pipeline: calculating dds and res...\n")
		print(system.time(dds <- DESeq(ddsHTSeq, parallel=TRUE))) # isa 'DESeqDataSet'
	}
	cat("\ncalculating results...\n")
	print(system.time(res <- results(dds))) # by default indepedant filtering at FDR (alpha) 0.1 (10%) 
											# independant filtering will set baseMean to filter genes automatically
											# by maximizing the number of genes satisfying < 0.1 FDR
											# isa 'DESeqResults'
	#res <- results(dds, lfcThreshold=2)
	cat("\ncalculating rld and vsd...\n")
	# Extracting transformed values
	#print(system.time(rld <- rlog(dds)))# blind=TRUE (default)
										# takes quite log time!	
	#rld.nb <- rlog(ddsHTSeq, blind=FALSE) # use thie design matrix 
	#vsd <- varianceStabilizingTransformation(ddsHTSeq)
	###########################################
	# Filter low count genes from HTSeq count #
	###########################################
	if(FALSE){
		# filter genes by cpm (or fpm)
		#ddsCpms<-cpm(counts(ddsHTSeq)) # require 'edgeR'
		ddsCpms<-fpm(ddsHTSeq) 
		cat("applying filters...","\n")
		if(sampleType=="PAPPA.SGA"){
			ddsKeep = rowSums(ddsCpms >1) >= nrow(samples[samples$PAPPA==1,])
		}else if(sampleType=="4AGA.FSEX"){ # 2M vs. 2F
			cat("removing chrY...\n")
			fields <- c("chromosome_name", "ensembl_gene_id")
			my.anno <- getBM(attributes = my.fields, filters = my.filter, values = rownames(ddsHTSeq), mart = myMart)
			ddsKeep <- rownames(ddsHTSeq) %in% my.anno[my.anno$chromosome_name!='Y',c(my.filter)] # remove ensg of chrY
		}else{
			ddsKeep = rowSums(ddsCpms >1) >= nrow(samples)/2
		}
		cat("\nDESeq from filtered count data ...\n")
		newCounts <- counts(ddsHTSeq)[ddsKeep,]
		ddsMatrix <- DESeqDataSetFromMatrix(newCounts, colData(ddsHTSeq), deseq.design) # isa 'DESeqDataSet'
		ddsFlt <- DESeq(ddsMatrix, parallel=TRUE) # isa 'DESeqDataSet'
		resFlt <- results(ddsFlt) # by default indepedant filtering at FDR (alpha) 0.1 (10%) 
	}#end of FALSE

	####################################
	# set exon size per gene for FPKM  #
	####################################
	keep <- rownames(dds) %in% names(my.target.list) # remove entries of no GRange info
	# if all of rownames(dds) are in my.target (e.g. piRNA or mRNA-GRCh38)
	cat("\nsetting rowData(dds)...\n")
	cat("\nsaving DESeq2 RData\n")
	if(table(keep)["TRUE"]==length(keep)){
		dummy<-mcols(dds)
		rowData(dds) <- my.target.list[rownames(dds)[keep]] # sorted by the original rownames
		#rowRanges(dds) <- my.target.list[rownames(dds)[keep]] # sorted by the original rownames (DESeq2_1.10.1)
		mcols(rowData(dds))<-dummy
		save(dds,res, file=deseq.RData) #  save 'dds'
		#save(dds,res,rld, file=deseq.RData) #  save 'dds'
	# if not all of rownames(dds) are in my.target (e.g. miRNA and mRNA)
	}else{
		newCounts <- counts(dds)[keep,]
		ddsMatrix <- DESeqDataSetFromMatrix(newCounts, colData(dds)[1:ncol(colData(dds))-1], deseq.design) # isa 'DESeqDataSet'
		ddsExonicGene <- DESeq(ddsMatrix, parallel=TRUE) # isa 'DESeqDataSet'
		dummy<-mcols(ddsExonicGene)
		rowData(ddsExonicGene) <- my.target.list[rownames(ddsExonicGene)] # set the length of this object via rowData
		mcols(rowData(ddsExonicGene))<-dummy
		save(dds,res,ddsExonicGene, file=deseq.RData) #  save 'dds'
		#save(dds,res,rld,ddsExonicGene, file=deseq.RData) #  save 'dds'
	}
}# end of loading RData

my.res<-res
rn=rownames(colData(dds)) #as.character(samples[["SampleName"]]) # samples$SampleName (or samples[,c("SampleName")])
rn.control=rownames(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1==0,]) 
			#as.character(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1==0,c("SampleName")]) #as.character(samples[as.numeric(samples[[my.contrast]])-1==0,c("SampleName")]) # control only (AGA or Male), which is background
rn.case=rownames(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1!=0,])
			#as.character(samples[as.numeric(samples[[my.contrast]])-1!=0,c("SampleName")]) # case only (SGA or female)
cat("setting FPKM...\n")
keep <- rownames(dds) %in% names(my.target.list) # remove entries of no GRange info
if(table(keep)["TRUE"]==length(keep)){
	ddsFpkm <-fpkm(dds) # isa 'matrix'
	ddsFpm <-fpm(dds) # isa 'matrix'
}else{
	ddsFpkm <-fpkm(ddsExonicGene) # isa 'matrix'
	ddsFpm <-fpm(ddsExonicGene) # isa 'matrix'
}

######################
## Data Exploration ##
######################
#if(TRUE){
if(sampleType=="ALL"){
	if(!exists("rld")){
		rld <- rlog(dds) # isa 'SummarizedExperiment'
	}
	#############################################
	# Heatmap of the sample-to-sample distances #
	# Sample Clustering Based on rlog           #
	#############################################
	rlogMat <- assay(rld) # isa 'matirx' (row: genes, col: samples); blind=TRUE by default

	hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
	cat("making sampleDists of case & control...\n")
	sampleDists<- dist(t(rlogMat)) # isa 'dist'
	mat <- as.matrix(sampleDists) # NxN distance matrix (row: samples, col: samples)
	cat("making heatmap of sampleDists...\n")
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".sample.clustering.heatmap"))
	tiff(filename=paste0(my.filename,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
	print(heatmap.2(mat, trace="none", col = hmcol, key.xlab="Distance", key.title="", main="Sample Distance"))
	dev.off()

	if(!grepl("^Retosiban",myProject)){
		cat("making sampleDists of case...\n")
		sampleDists<- dist(t(rlogMat[,rn.case])) # isa 'dist'
		mat <- as.matrix(sampleDists) # NxN distance matrix (row: samples, col: samples)
		cat("making heatmap of sampleDists...\n")
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".sample.clustering.heatmap.case"))
		tiff(filename=paste0(my.filename,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
		print(heatmap.2(mat, trace="none", col = hmcol, key.xlab="Distance", key.title="", main="Sample Distance of case"))
		dev.off()

		cat("making sampleDists of control...\n")
		sampleDists<- dist(t(rlogMat[,rn.control])) # isa 'dist'
		mat <- as.matrix(sampleDists) # NxN distance matrix (row: samples, col: samples)
		cat("making heatmap of sampleDists...\n")
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".sample.clustering.heatmap.control"))
		tiff(filename=paste0(my.filename,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
		print(heatmap.2(mat, trace="none", col = hmcol, key.xlab="Distance", key.title="", main="Sample Distance of control"))
		dev.off()

		# summary of boxplot of rlog diff by pair
		if(any(grepl("Pair", deseq.design))){
			diff.rlogMat<-rlogMat[,rn.case]-rlogMat[,rn.control]
			#summary(rlogMat[,rn.case]-rlogMat[,rn.control])
			my.filename <- file.path(cluster.dir, paste0(sampleType, ".boxplot.rlog.case-control"))
			tiff(filename=paste0(my.filename,".tiff"),width=16,height=9,units="in",res=300, compression = 'lzw')
			original.parameters<-par()
			par(xaxt="n") 
			boxplot(diff.rlogMat, main=paste0("Boxplot for rlog(case)-rlog(control)"), xlab='case', ylab='rlog(case)-rlog(control)')
			axis(1, at=seq(1, ncol(diff.rlogMat), by=1), labels = FALSE)
			text(seq(1, ncol(diff.rlogMat), by=1), par("usr")[3] - 0.2, labels = colnames(diff.rlogMat), srt = 45, pos = 1, xpd = TRUE)
			dev.off()
			par(original.parameters)
		}
	}

	#########################################
	## PCA Based on Sample-Sample Distance ##
	#########################################
	select<-!mcols(dds)[["allZero"]] # genes having counts from at least one sample (not allZero)
	cat("making pca based on rlogMat...\n")
	pca<-prcomp(t(rlogMat[select,]), scale.=TRUE) # input (row:samples, col:genes)
	percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
	#ggplot way
	p<- ggplot(as.data.frame(pca$x), aes(PC1,PC2)) + geom_text(label=rownames(pca$x)) 
	p<- p + xlab(paste0("PC1 (",round(percentVar[1]*100,2)," % variance)")) + ylab(paste0("PC2 (",round(percentVar[2]*100,2)," % variance)"))
	p<- p + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
	cat("making pca plot...\n")
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".sample.pca"))
	tiff(filename=paste0(my.filename,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
	print(plot)
	dev.off()

	##############################
	# Heatmap of Gene Clustering # 
	# Based on top variable gene #
	##############################
	hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
	cat("making heatmap of gene clustering ...\n")
	select<- head( order(rowVars(rlogMat), decreasing=TRUE ), n=20 ) # top 20 variable genes (requires 'genefilter' library)
	my.table<-rlogMat[select,]
	my.table<-as.data.frame(my.table[order(rownames(my.table)),]) # order by rownames (ensembl_gene_id or mirbase_id)

	if(!grepl("iRNA",myProject)){ # total RNA
		fields <- c(my.filter, "hgnc_symbol","description")
		anno<-getBM(attributes = fields, filters = my.filter, values = rownames(rlogMat[select,]), mart = myMart)
		anno.collapsed <- plyr::ddply(anno, my.filter, function(df) paste(df[[my.id]], collapse=';')) 
		colnames(anno.collapsed) <- c(my.filter,my.id)
		my.table[,my.filter]<-rownames(my.table) # add a new column id
		my.table<-merge(as.data.frame(my.table), anno.collapsed, by=my.filter,  all.x=TRUE) # left join to add column 'hgnc_symbol'
		rownames(my.table) <- ifelse(my.table[[my.id]]=="" | is.na(my.table[[my.id]]), my.table[[my.filter]], my.table[[my.id]]) # replace with ensembl id if no gene name
		my.table[[my.filter]]<-NULL; my.table[[my.id]]<-NULL # drop columns
		write.csv(my.table, file=gzfile(file.path(cluster.dir, paste0(sampleType, ".top.var.gene.csv.gz"))), row.names = T)
	}else{
		write.csv(my.table, file=gzfile(file.path(cluster.dir, paste0(sampleType, ".top.var.gene.csv.gz"))), row.names = TRUE)
	}

	my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.var.gene.clustering.heatmap"))
	tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
	print(heatmap.2(as.matrix(my.table), scale="row", trace="none", dendrogram="column", col = hmcol, keysize=1, key.title="", margin=c(5,7), offsetRow=0.1, main="Top Variable Genes"))
	dev.off()

	if(!grepl("^Retosiban",myProject)){
		cat("making heatmap of gene clustering from case...\n")
		select<- head( order(rowVars(rlogMat[,rn.case]), decreasing=TRUE ), n=20 ) # top 20 variable genes (requires 'genefilter' library)
		my.table<-rlogMat[select,rn.case]
		my.table<-as.data.frame(my.table[order(rownames(my.table)),]) # order by rownames (ensembl_gene_id or mirbase_id)

		if(!grepl("iRNA",myProject)){ # total RNA
			fields <- c(my.filter, "hgnc_symbol","description")
			anno<-getBM(attributes = fields, filters = my.filter, values = rownames(rlogMat[select,rn.case]), mart = myMart)
			anno.collapsed <- plyr::ddply(anno, my.filter, function(df) paste(df[[my.id]], collapse=';')) 
			colnames(anno.collapsed) <- c(my.filter,my.id)
			my.table[,my.filter]<-rownames(my.table) # add a new column id
			my.table<-merge(as.data.frame(my.table), anno.collapsed, by=my.filter,  all.x=TRUE) # left join to add column 'hgnc_symbol'
			rownames(my.table) <- ifelse(my.table[[my.id]]=="" | is.na(my.table[[my.id]]), my.table[[my.filter]], my.table[[my.id]]) # replace with ensembl id if no gene name
			my.table[[my.filter]]<-NULL; my.table[[my.id]]<-NULL # drop columns
			write.csv(my.table, file=gzfile(file.path(cluster.dir, paste0(sampleType, ".top.var.gene.case.csv.gz"))), row.names = FALSE)
		}else{
			write.csv(my.table, file=gzfile(file.path(cluster.dir, paste0(sampleType, ".top.var.gene.case.csv.gz"))), row.names = TRUE)
		}

		my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.var.gene.clustering.heatmap.case"))
		tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
		print(heatmap.2(as.matrix(my.table), scale="row", trace="none", dendrogram="column", col = hmcol, keysize=1, key.title="", margin=c(5,7), offsetRow=0.1, main="Top Variable Genes of case"))
		dev.off()

		cat("making heatmap of gene clustering from control...\n")
		#select<- head( order(rowVars(rlogMat[,rn.control]), decreasing=TRUE ), n=round(nrow(rlogMat[,rn.control])*0.1*0.01) ) # top 0.1% variable genes (requires 'genefilter' library)
		select<- head( order(rowVars(rlogMat[,rn.control]), decreasing=TRUE ), n=20 ) # top 20 variable genes (requires 'genefilter' library)
		my.table<-rlogMat[select,rn.control]
		my.table<-as.data.frame(my.table[order(rownames(my.table)),]) # order by rownames (ensembl_gene_id or mirbase_id)

		if(!grepl("iRNA",myProject)){ # total RNA
			fields <- c(my.filter, "hgnc_symbol","description")
			anno<-getBM(attributes = fields, filters = my.filter, values = rownames(rlogMat[select,rn.case]), mart = myMart)
			anno.collapsed <- plyr::ddply(anno, my.filter, function(df) paste(df[[my.id]], collapse=';')) 
			colnames(anno.collapsed) <- c(my.filter,my.id)
			my.table[,my.filter]<-rownames(my.table) # add a new column id
			my.table<-merge(as.data.frame(my.table), anno.collapsed, by=my.filter,  all.x=TRUE) # left join to add column 'hgnc_symbol'
			rownames(my.table) <- ifelse(my.table[[my.id]]=="" | is.na(my.table[[my.id]]), my.table[[my.filter]], my.table[[my.id]]) # replace with ensembl id if no gene name
			my.table[[my.filter]]<-NULL; my.table[[my.id]]<-NULL # drop columns
			write.csv(my.table, file=gzfile(file.path(cluster.dir, paste0(sampleType, ".top.var.gene.control.csv.gz"))), row.names = FALSE)
		}else{
			write.csv(my.table, file=gzfile(file.path(cluster.dir, paste0(sampleType, ".top.var.gene.control.csv.gz"))), row.names = TRUE)
		}

		my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.var.gene.clustering.heatmap.control"))
		tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
		print(heatmap.2(as.matrix(my.table), scale="row", trace="none", dendrogram="column", col = hmcol, keysize=1, key.title="", margin=c(5,7), offsetRow=0.1, main="Top Variable Genes of control"))
		dev.off()
	}

	###############################################################
	## Box plot of FPKM of genes having > 1 CPM for all samples  ##
	## All samples should have more than 1 CPM/FPM               ##
	###############################################################
	cat("making ddsFpkm from dds(ExonicGene)...\n")
	keep = rowSums(ddsFpm >1) >= ncol(ddsFpm)
	if(FALSE){
		if(grepl("piRNA",myProject)){ # If piRNA 
			keep = rowSums(ddsFpm >1) >= ncol(dds)
		}else{ 
			keep = rowSums(ddsFpm >1) >= ncol(ddsExonicGene)
		}
	}

	#cat("making heatmap of boxplot of FPKM in log2 scale...\n")
	#boxplot(log2(1+ddsFpkm[keep,]), xlab="Samples", ylab="log2(1+FPKM)", main=paste0("FPKM boxplots of ",nrow(ddsFpkm[keep,]), " genes having > 1 CPM across all samples"))

	cat("making heatmap of boxplot of FPKM...\n")
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".fpkm.boxplot"))
	tiff(filename=paste0(my.filename,".tiff"),width=10,height=9,units="in",res=300, compression = 'lzw')
	print(boxplot(ddsFpkm[keep,], xlab="Samples", ylab="FPKM", main=paste0("FPKM boxplots of ",nrow(ddsFpkm[keep,]), " genes having > 1 CPM across samples")))
	dev.off()
	if(!grepl("^Retosiban",myProject)){
		cat("making heatmap of boxplot of case FPKM...\n")
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".fpkm.boxplot.case"))
		tiff(filename=paste0(my.filename,".tiff"),width=10,height=9,units="in",res=300, compression = 'lzw')
		print(boxplot(ddsFpkm[keep,rn.case], xlab="case Samples", ylab="FPKM", main=paste0("FPKM boxplots of ",nrow(ddsFpkm[keep,rn.case]), " genes having > 1 CPM across case samples")))
		dev.off()

		cat("making heatmap of boxplot of control FPKM...\n")
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".fpkm.boxplot.control"))
		tiff(filename=paste0(my.filename,".tiff"),width=10,height=9,units="in",res=300, compression = 'lzw')
		print(boxplot(ddsFpkm[keep,rn.control], xlab="control Samples", ylab="FPKM", main=paste0("FPKM boxplots of ",nrow(ddsFpkm[keep,rn.control]), " genes having > 1 CPM across control samples")))
		dev.off()
	}

	##################
	## FPKM  of IQR ##
	##################
	hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
	if(!grepl("^Retosiban",myProject)){
		# Top 200 FPKM Genes
		cat("making heatmap of Top FPKM based on IQR mean...\n")
		myTopExp <- getTopExp(ddsFpkm[keep,], iqr=FALSE, 200) # isa 'list'
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.iqr.fpkm.gene.heatmap.all"))
		tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
		print(heatmap.2(head(as.matrix(myTopExp$fpkm),n=round(nrow(ddsFpkm)*0.1*0.01)), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="FPKM", main="Top FPKM Genes"))
		dev.off()

		myTopExp$fpkm$mean<-rowMeans(myTopExp$fpkm,na.rm=T) # mean over all the samples
		myTopExp$fpkm$meanIQR<-apply(myTopExp$fpkm, 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]])) # mean over trimmed samples
		if(grepl("piRNA",myProject)){ # If piRNA
			my.fpkm<-myTopExp$fpkm
		}else{
			my.fpkm<-merge(cbind(myTopExp$fpkm, hgnc_symbol=rownames(myTopExp$fpkm)), myTopExp$gene, all.x=TRUE)
		}
		write.csv(my.fpkm, file=gzfile(file.path(cluster.dir, paste0(sampleType, ".top.200.iqr.fpkm.csv.gz"))))

		# Top 200 case FPKM
		cat("making heatmap of Top FPKM based on IQR mean of case...\n")
		myTopExp.case <- getTopExp(ddsFpkm[keep,rn.case], iqr=TRUE, 200) # isa 'list'
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.iqr.fpkm.gene.heatmap.case"))
		tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
		print(heatmap.2(head(as.matrix(myTopExp.case$fpkm),n=round(nrow(ddsFpkm)*0.1*0.01)), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="FPKM", main="Top IQR(FPKM) of case"))
		dev.off()

		myTopExp.case$fpkm$mean<-rowMeans(myTopExp.case$fpkm,na.rm=T) # mean over all the samples
		myTopExp.case$fpkm$meanIQR<-apply(myTopExp.case$fpkm, 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]])) # mean over trimmed samples
		if(grepl("piRNA",myProject)){ # If piRNA
			my.fpkm<-myTopExp.case$fpkm
		}else{
			my.fpkm<-merge(cbind(myTopExp.case$fpkm, hgnc_symbol=rownames(myTopExp.case$fpkm)), myTopExp.case$gene, all.x=TRUE)
		}
		write.csv(my.fpkm, file=gzfile(file.path(cluster.dir, paste0(sampleType, ".top.200.iqr.fpkm.case.csv.gz"))))

		# Top 200 control FPKM
		cat("making heatmap of Top FPKM based on IQR mean of control...\n")
		myTopExp.ctl<-getTopExp(ddsFpkm[keep,rn.control], iqr=TRUE, 200) # isa 'list' 
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.iqr.fpkm.gene.heatmap.control"))
		tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
		print(heatmap.2(head(as.matrix(myTopExp.ctl$fpkm),n=round(nrow(ddsFpkm)*0.1*0.01)), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="FPKM", main="Top IQR(FPKM) of control"))
		dev.off()

		myTopExp.ctl$fpkm$mean<-rowMeans(myTopExp.ctl$fpkm,na.rm=T) # mean over all the samples
		myTopExp.ctl$fpkm$meanIQR<-apply(myTopExp.ctl$fpkm, 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]])) # mean over trimmed samples
		if(grepl("piRNA",myProject)){ # If piRNA
			my.fpkm<-myTopExp.ctl$fpkm
		}else{
			my.fpkm<-merge(cbind(myTopExp.ctl$fpkm, hgnc_symbol=rownames(myTopExp.ctl$fpkm)), myTopExp.ctl$gene, all.x=TRUE)
		}
		write.csv(my.fpkm, file=gzfile(file.path(cluster.dir, paste0(sampleType, ".top.200.iqr.fpkm.control.csv.gz"))))
	}

	#####################
	## All Counts/FPKM ##
	#####################
	write.csv(ddsFpkm, file=gzfile(file.path(cluster.dir, paste0(sampleType, ".fpkm.csv.gz"))))
	write.csv(counts(dds), file=gzfile(file.path(cluster.dir, paste0(sampleType, ".count.csv.gz"))))

	if(any(grepl("Pair", deseq.design))){
		write.csv(ddsFpkm[,rn.case]-ddsFpkm[,rn.control], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".fpkm.diff.csv.gz"))))
		write.csv(ddsFpkm[,rn.case]/ddsFpkm[,rn.control], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".fpkm.ratio.csv.gz"))))

		write.csv(counts(dds)[,rn.case]-counts(dds)[,rn.control], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".count.diff.csv.gz"))))
		write.csv(counts(dds)[,rn.case]/counts(dds)[,rn.control], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".count.ratio.csv.gz"))))
	}
	###############################
	## SD-Mean plot of IQR(FPKM) ##
	###############################
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".iqr.fpkm.log2.one.plus.mean.log2sd"))
	tiff(filename=paste0(my.filename,".tiff"),width=9,height=10,units="in",res=300, compression = 'lzw')
	plot(log2(1+apply(ddsFpkm[keep,], 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]]))), log2(1+apply(ddsFpkm[keep,], 1, function(i) sd(i[i>=quantile(i)[2] & i<=quantile(i)[4]]))), xlab="log2(1+mean(IQR(FPKM)))", ylab="log2(1+SD(IQR(FPKM)))", main="log2(mean)-log2(SD) plot") # SD-Mean plot
	dev.off()

	if(!grepl("^Retosiban",myProject)){
		cat("making Sdplot based on ddsFpkm of case...\n")
		mean.iqr.fpkm<-apply(ddsFpkm[keep,rn.case], 1, function(i) ifelse(length(i)>=4,mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]]),mean(i))) # isa numeric vector
		sd.iqr.fpkm<-apply(ddsFpkm[keep,rn.case], 1, function(i) ifelse(length(i)>=4,sd(i[i>=quantile(i)[2] & i<=quantile(i)[4]]),sd(i))) # isa numeric vector
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".iqr.fpkm.mean.sd.case"))
		tiff(filename=paste0(my.filename,".tiff"),width=9,height=10,units="in",res=300, compression = 'lzw')
		plot(mean.iqr.fpkm, sd.iqr.fpkm, xlab="mean(IQR(FPKM))", ylab="SD(IQR(FPKM))", main="Mean-SD plot of case") # SD-Mean plot
		text(cbind(mean.iqr.fpkm, sd.iqr.fpkm)[myTopExp.case$select[1:5],], labels=rownames(myTopExp.case$fpkm)[1:5], pos=2) # label to the left (pos=2)
		dev.off()

		cat("making Sdplot based on ddsFpkm of control...\n")
		mean.iqr.fpkm<-apply(ddsFpkm[keep,rn.control], 1, function(i) ifelse(length(i)>=4,mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]]),mean(i))) # isa numeric vector
		sd.iqr.fpkm<-apply(ddsFpkm[keep,rn.control], 1, function(i) ifelse(length(i)>=4,sd(i[i>=quantile(i)[2] & i<=quantile(i)[4]]),sd(i))) # isa numeric vector
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".iqr.fpkm.mean.sd.control"))
		tiff(filename=paste0(my.filename,".tiff"),width=9,height=10,units="in",res=300, compression = 'lzw')
		plot(mean.iqr.fpkm, sd.iqr.fpkm, xlab="mean(IQR(FPKM))", ylab="SD(IQR(FPKM))", main="Mean-SD plot of control") # SD-Mean plot
		text(cbind(mean.iqr.fpkm, sd.iqr.fpkm)[myTopExp.ctl$select[1:5],], labels=rownames(myTopExp.ctl$fpkm)[1:5], pos=2) # label to the left (pos=2)
		dev.off()
	}
} #end of sampleType=="ALL"

###############################
# prepare pdf output filename #
###############################
my.file.name<- paste0(deseq.dir,'/deseq.',sampleType)
pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title=paste0(myProject,":",sampleType)) # A4 size

# Boy.Girl chrX only
#if(myProject=="Boy.Girl"){
if(FALSE){
	my.anno <- getBM(attributes = my.fields, filters = my.filter, values = rownames(dds), mart = myMart)

	##############
	## Outliers ##
	##############
	df.cooks<-assays(dds)$cooks
	top.cooks.genes<-names(head(sort(rowSums(df.cooks), decreasing=T), n=10)) # top cooks genes
	ddsFpkm[top.cooks.genes,] # top cooks genes
	top.cooks.samples<-sort(colSums(ddsFpkm[top.cooks.genes,]), decreasing=T) # top cooks samples
	round(ddsFpkm[top.cooks.genes, head(names(top.cooks.samples),n=10)], 2)
	round(sort(ddsFpkm[top.cooks.genes[1],], decreasing=T))

	my.anno[my.anno$ensembl_gene_id%in% top.cooks.genes,]

	#IGFBP1: 84C, 88, 80C, 06C
	#SNORD3C: 02C, 02C
	dummy<-rbind(
		merge(cbind(reshape2::melt(ddsFpm[,rn.case]),`Sex`=c("F")), unique(my.anno[,c("ensembl_gene_id","chromosome_name")]), by.x="Var1", by.y="ensembl_gene_id"), # case only (SGA or feale)
		merge(cbind(reshape2::melt(ddsFpm[,rn.control]),`Sex`=c("M")), unique(my.anno[,c("ensembl_gene_id","chromosome_name")]), by.x="Var1", by.y="ensembl_gene_id") # control only (AGA or male)
	)
	#> head(dummy)
	#             Var1 Var2    value Sex chromosome_name
	#1 ENSG00000000003   76 25.35166   M               X
	#2 ENSG00000000003   66 16.27970   M               X
	#3 ENSG00000000003  152 33.64118   M               X
	#4 ENSG00000000003  100 26.35651   M               X
	#5 ENSG00000000003  102 13.24184   M               X
	#6 ENSG00000000003   82 23.23487   M               X
	#p<-ggplot(dummy, aes(Sex,value)) + geom_boxplot(aes(fill=Sex))
	#print(p)

	p<-ggplot(dummy[dummy$chromosome_name=="X",], aes(log10(value))) + geom_density(aes(fill=Sex),alpha=.7)
	print(p)

	#########################################
	# density plot of log2foldchange by sex #
	#########################################
	foo<-merge(cbind(as.data.frame(res), `ensembl_gene_id`=rownames(res)), unique(my.anno[,c("ensembl_gene_id","chromosome_name")]))
	#> head(foo)
	#  ensembl_gene_id   baseMean log2FoldChange      lfcSE        stat    pvalue      padj chromosome_name
	#1 ENSG00000000003 1385.57760    0.129245753 0.09597355  1.34668096 0.1780830 0.7565668               X
	#2 ENSG00000000005   49.46438   -0.093313517 0.13485629 -0.69194782 0.4889701 0.9154923               X
	select<-foo$baseMean>=minRead & !is.na(foo$padj) & foo$padj<0.1

	p<-ggplot(foo[select,], aes(chromosome_name, log2FoldChange)) + 
		geom_boxplot() 
	print(p)

	p<-ggplot(foo[select,], aes(log2FoldChange)) +
		geom_density(aes(fill=chromosome_name),alpha=.5) + 
		xlim(0.2,5) + 
		labs(x="Log2FoldChange(Girl - Boy)")
	print(p)

	p<-ggplot(foo[select,], aes(log2FoldChange,..density..)) + 
		geom_freqpoly(aes(color=chromosome_name),alpha=.7) + 
		labs(x="Log2FoldChange(Girl - Boy)")
	print(p)

	###########
	## CSMD1 ##
	###########
	my.entry="ENSG00000183117" # CSMD1
	foo<-data.table(`SampleName`=colnames(dds),`FPKM`=ddsFpkm[my.entry,], `FPM`=ddsFpm[my.entry,], `RawCount`=counts(dds)[my.entry,])
	write.csv(foo, file=file.path(deseq.dir, paste0(my.entry, ".csv")))
	# p-value
	my.res[my.entry,]
	ddsFpm[my.entry,rn.case]
	ddsFpm[my.entry,rn.control]
	wilcox.test(ddsFpm[my.entry,rn.case], ddsFpm[my.entry,rn.control])$p.value
	# up/down 1Mb
	my.upflank=10^6 # 1Mb up from target TSS
	my.downflank=10^6 # 1Mb down from target TSS 
	my.gr.ext<-promoters(gr.ensg[my.entry], upstream=my.upflank, downstream=end(gr.ensg[my.entry])-start(gr.ensg[my.entry])+1+my.downflank) # includes +- flanking from TSS and TES
	my.overlap=as.matrix( findOverlaps(gr.ensg, my.gr.ext, ignore.strand=TRUE, minoverlap=10L) ) #By default, any overlap is accepted
	my.target.ensg<-names(gr.ensg[my.overlap[,"queryHits"]])
	my.anno[my.anno$ensembl_gene_id %in% my.target.ensg,]

	ggplot(dt.csmd.fpkm, aes(tissue, FPKM, fill=Gender)) + geom_bar(stat="identity", position="dodge") + scale_fill_manual(name="Sex",values=my.col[["Gender"]]) + theme_Publication()


	##
	## chrX
	##
	dt.res<-data.table(mirbase_id=rownames(res),data.frame(res))
	dt.chr<-as.data.table(getBM(attributes =c("mirbase_id","ensembl_gene_id","chromosome_name"), filters = "mirbase_id", values = dt.res[,unique(mirbase_id)], mart = myMart)) # is a data.frame
	dt.res.chrX=merge(dt.res, dt.chr, by="mirbase_id", all.x=TRUE, allow.cartesian=TRUE)[chromosome_name=="X"]
	dt.res.chrX<-dt.res.chrX[,`new.padj` := p.adjust(pvalue, method="BH")]

} #myProject=="Boy.Girl"


# CSMD1 exon-level expression analysis
# Boy.Girl.exon.GRCh37 
if(grepl("exon",myProject)){ 
	library(data.table)
	source("~/lib/ggplot_dual_axis.R") # only R>=3.2

	my.ensg="ENSG00000183117" # CSMD1
	my.first.exon_id="ENSE00002127078" # first exon of placenta CSMD1 transcript

	##############################
	# import Ensembl exon region #
	##############################
	#my.target; my.target.list # defined within DEG.R
	subject.exon.keys=c("seqnames","strand","start","end")
	subject.exon<-as.data.frame(my.target)
	cat("\tconvert to data.table...\n")
	dt.ens.exon<-as.data.table(subject.exon) # isa data.table
	setkeyv(dt.ens.exon,subject.exon.keys) # regions of interests 
	cat("Ensembl Exon loaded\n")

	######################################
	# import reduced Ensembl exon region #
	######################################
	subject.rdx.exon<-as.data.frame(GenomicRanges::reduce(my.target))
	dt.rdx.exon<-as.data.table(subject.rdx.exon) # isa data.table
	dt.rdx.exon[,feature_id:=paste0("exon.",.N-.I+1)]
	dt.rdx.exon[,feature_type:="exon"]
	cat("Ensembl reduced Exon loaded\n")

	########################
	# import intron region #
	########################
	subject.intron<-as.data.frame(gaps(reduce(my.target))[2:length(reduce(my.target))])
	dt.intron<-as.data.table(subject.intron) # isa data.table
	dt.intron[,feature_id:=paste0("intron.",.N-.I+1)]
	dt.intron[,feature_type:="intron"]
	cat("Ensembl reduced Intron loaded\n")
	#rm(subject.exon, subject.rdx.exon, subject.intron)

	#################################
	# merge reduced exon and intron #
	#################################
	dt.rdx.feature=rbind(dt.rdx.exon, dt.intron)[order(start)]
	setkeyv(dt.rdx.feature,subject.exon.keys) # regions of interests 

	###########################
	# import methylation data #
	###########################
	source("~/Pipelines/bin/R/RoadMap/local.R")
	dt.query=load.my.tissue.dt.merged(my.tissue="PT",my.cpg.type="CG")

	cat("\tFinding CpGs overlap...\n")
	query.key=c("V1","V3","V2","End")
	# merge methylation level by reduced feature of CSMD1 
	#####################################
	# methylation level at Ensembl exon #
	#####################################
	system.time(dt.overlap<-foverlaps(dt.query, dt.ens.exon, by.x=query.key, type="any", nomatch=0L))
	system.time(dt.dmr.exon<-dt.overlap[,list("num.sites"=.N,c.f=sum(V5.x), t.f=sum(V6.x)-sum(V5.x), c.m=sum(V5.y), t.m=sum(V6.y)-sum(V5.y),
										met=round(sum(V5.x,V5.y)/sum(V6.x,V6.y)*100,2), 
										met.f=round(sum(V5.x)/sum(V6.x)*100,2), 
										met.m=round(sum(V5.y)/sum(V6.y)*100,2), 
										p.value=fisher.test(matrix(c(sum(V5.x),sum(V6.x)-sum(V5.x),sum(V5.y),sum(V6.y)-sum(V5.y)),nrow=2))$p.value),
										by="exon_id"])
	########################################################
	# methylation level at redued features (exon & intron) #
	########################################################
	system.time(dt.overlap<-foverlaps(dt.query, dt.rdx.feature, by.x=query.key, type="any", nomatch=0L))
	system.time(dt.dmr.feature<-dt.overlap[,list("num.sites"=.N,c.f=sum(V5.x), t.f=sum(V6.x)-sum(V5.x), c.m=sum(V5.y), t.m=sum(V6.y)-sum(V5.y),
										met=round(sum(V5.x,V5.y)/sum(V6.x,V6.y)*100,2), 
										met.f=round(sum(V5.x)/sum(V6.x)*100,2), 
										met.m=round(sum(V5.y)/sum(V6.y)*100,2), 
										p.value=fisher.test(matrix(c(sum(V5.x),sum(V6.x)-sum(V5.x),sum(V5.y),sum(V6.y)-sum(V5.y)),nrow=2))$p.value),
										by="feature_id"])
	dt.dmr.feature<-merge(dt.dmr.feature, dt.rdx.feature, by="feature_id")[order(start)]

	# methylation level across reduced exons and intron
	p.feature.meth<-ggplot(dt.dmr.feature, aes(feature_id, met, group=1)) + 
				geom_point(aes(col=feature_type), size=rel(1)) +  
				geom_smooth(aes(group=feature_type,col=feature_type),se=FALSE) +
				scale_x_discrete(limits=dt.dmr.feature[,feature_id]) + 
				theme_bw() +
				theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.title = element_text(face="italic",size=rel(.9)))

	# methylation level boxplot at exon and intron
	p.feature.meth.box<-ggplot(dt.dmr.feature[num.sites>5]) + 
						geom_boxplot(aes(feature_type, met)) + 
						theme_Publication()

	###########################
	# process expression data #
	###########################
	my.pval=0.05
	exon.pval<-sapply(rownames(my.res), function(i) wilcox.test(ddsFpm[i,rn.case], ddsFpm[i,rn.control])$p.value)
	dt.fpkm.exon<-data.table(`exon_id`=rownames(ddsFpkm), 
							 `FPKM.F`=rowMeans(ddsFpkm[,rn.case]), `FPKM.M`=rowMeans(ddsFpkm[,rn.control]), 
							 `FPM.F`=rowMeans(ddsFpm[,rn.case]), `FPM.M`=rowMeans(ddsFpm[,rn.control]),
							 `CNT.F`=rowMeans(counts(dds,normalized=T)[,rn.case]), `CNT.M`=rowMeans(counts(dds,normalized=T)[,rn.control]),
							 `pvalue`=exon.pval)

	# write Differentially Expressed Exon (DEE) in GTF
	my.file=file.path(cluster.dir, paste0(sampleType, ".CSMD1.DEE.",TR_PREFIX,".gtf"))
	export.gff(my.target[mcols(my.target)$exon_id  %in% names(exon.pval[!is.nan(exon.pval) & exon.pval< my.pval]),], my.file)

	##########################################
	# merge expression & methylation by exon #
	##########################################
	dt.exp.meth.exon<-merge(dt.fpkm.exon, dt.dmr.exon, by="exon_id", all.x=TRUE)[exon_id %in% dt.fpkm.exon[,exon_id]]
	dt.cnt.melt<-reshape2::melt(dt.exp.meth.exon, id=c("exon_id","pvalue"), measure=c("CNT.F","CNT.M"), variable.name="Gender", value.name="Count"); dt.cnt.melt[,Gender:=ifelse(Gender=="CNT.F","Female","Male")] # melt: col to row
	dt.cnt.melt[,significant:=ifelse(pvalue<my.pval,"Yes","No")]
	dt.meth.melt<-reshape2::melt(dt.exp.meth.exon, id=c("exon_id","num.sites","p.value"), measure=c("met.f","met.m"), variable.name="Gender", value.name="Meth"); dt.meth.melt[,Gender:=ifelse(Gender=="met.f","Female","Male")] # melt: col to row
	#merge(dt.cnt.melt, dt.meth.melt, by=c("exon_id","Gender"))

	dt.cnt.meth.diff<-rbind(dt.exp.meth.exon[,.(exon_id, `measure`="Expression", `diff`=CNT.F-CNT.M, `size`=(CNT.F+CNT.M)/2, `pval`=pvalue)],
							dt.exp.meth.exon[,.(exon_id, `measure`="Methylation",`diff`=met.f-met.m, `size`=num.sites, `pval`=p.value)])

	write.csv(dt.exp.meth.exon[,.(exon_id,CNT.F,CNT.M,`p.value.exp`=pvalue,num.sites,met.f,met.m,`p.value.met`=p.value)]
			  , file=gzfile(file.path(cluster.dir, paste0(sampleType, ".CSMD1.cnt.meth.by.exon",TR_PREFIX,".csv.gz"))))
	#####################
	## P-value plot    ##
	## with Read Count ##
	## with % Meth     ##
	#####################
	p.cnt.pval<-ggplot(dt.exp.meth.exon, aes(exon_id, pvalue, group=1)) + 
		geom_point(na.rm=T,size=rel(.5)) + 
		geom_line(aes(linetype="Read-Count Difference"),na.rm=T) + 
		geom_hline(aes(yintercept=my.pval)) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_y_continuous(breaks=c(0, my.pval, 0.25, .5, .75, 1)) +
		scale_linetype_manual(name="P-value", values=c(3)) + # dotted
		labs(x="Exon ID", y="P-value") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.position="top" ,legend.title = element_text(face="italic",size=rel(.9)))

	p.meth.pval<-ggplot(dt.exp.meth.exon, aes(exon_id, p.value, group=1)) + 
		geom_point(na.rm=T,size=rel(.5)) + 
		geom_line(linetype="dotted",na.rm=T) + 
		geom_hline(aes(yintercept=my.pval)) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_y_continuous(breaks=c(0, my.pval, 0.25, .5, .75, 1)) +
		labs(x="Exon ID", y="P-value") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.title = element_text(face="italic",size=rel(.9)))

	# p-value both from cnt & meth
	p.cnt.meth.pval<-ggplot(dt.cnt.meth.diff, aes(exon_id, pval, group=measure)) + 
		geom_point(aes(shape=measure), na.rm=T) + 
		geom_line(aes(linetype=measure),na.rm=T) +
		geom_hline(aes(yintercept=my.pval)) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_y_continuous(breaks=c(0, my.pval, 0.25, .5, .75, 1)) +
		scale_linetype_manual(name="Source of P-value",values=c(3,1)) +
		labs(x="Exon ID", y="P-value") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.title = element_text(face="italic",size=rel(.9)))

	# p-value both from cnt & meth with diff and size
	p.cnt.meth.pval2<-ggplot(dt.cnt.meth.diff, aes(exon_id, pval, group=measure)) + 
		geom_point(aes(size=size,col=diff,shape=measure), na.rm=T) + 
		geom_line(aes(linetype=measure),col='gray70',na.rm=T) +
		geom_hline(aes(yintercept=my.pval)) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_y_continuous(breaks=c(0, my.pval, 0.25, .5, .75, 1)) +
		scale_linetype_manual(name="Source of P-value",values=c(1,2)) +
		scale_colour_gradient2(name="Difference\n(Female-Male)",mid="gray48", high="red", low="blue", midpoint=0) +
		labs(x="Exon ID", y="P-value") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.title = element_text(face="italic",size=rel(.9)))

	# methylation level without no line and with No. CpG by shape
	p.meth<-ggplot(dt.meth.melt, aes(exon_id, Meth, group=Gender, colour=Gender)) +
		geom_point(aes(size=num.sites),shape=8,na.rm=T, alpha=0.9) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_size_continuous(name="No. CpG") +
		scale_colour_manual(name="% Methylation\n(by sex)",values=my.col[["Gender"]]) +
		labs(x="Exon ID", y="% Methylation") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.position="top", legend.title = element_text(face="italic",size=rel(.7)))
	# dotted-line with No. CpG by shape
	p.meth2<-ggplot(dt.meth.melt, aes(exon_id, Meth, group=Gender, colour=Gender)) +
		geom_point(aes(size=num.sites),shape=8,na.rm=T, alpha=0.9) + 
		geom_line(linetype="dotted",na.rm=T) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_size_continuous(name="No. CpG") +
		scale_colour_manual(name="% Methylation\n(by sex)",values=my.col[["Gender"]]) +
		labs(x="Exon ID", y="% Methylation") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.justification="top", legend.position=c(0.8,1), legend.title = element_text(face="italic",size=rel(.7)))

	# read-count without line with read-count by shape
	p.cnt<-ggplot(dt.cnt.melt, aes(exon_id, Count, colour=Gender)) +
		geom_point(aes(size=Count),na.rm=T, alpha=0.9) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_size_continuous(name="Read Count") +
		scale_colour_manual(name="Read Count\n(by sex)",values=my.col[["Gender"]]) +
		labs(x="Exon ID", y="Read Count") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.position="top", legend.title = element_text(face="italic",size=rel(.7)))

	# read-count solid-line with different shape by p-value threshold
	p.cnt2<-ggplot(dt.cnt.melt, aes(exon_id, Count, colour=Gender)) +
		geom_point(aes(shape=factor(significant)),size=3,na.rm=T, alpha=0.9) + 
		geom_line(aes(group=Gender),na.rm=T) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_colour_manual(name="Read Count\n(by sex)",values=my.col[["Gender"]]) +
		scale_shape(name=paste("P-value<",my.pval)) +
		labs(x="Exon ID", y="Read Count") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.justification="top", legend.position=c(0.8,1), legend.title = element_text(face="italic",size=rel(.7)))

	# log read-count without line with read-count by shape
	p.cnt3<-ggplot(dt.cnt.melt, aes(exon_id, log(Count+0.1), colour=Gender)) +
		geom_point(aes(size=Count),na.rm=T, alpha=0.9) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_size_continuous(name="Read Count") +
		scale_colour_manual(name="Sex",values=my.col[["Gender"]]) +
		labs(x="Exon ID", y="Log(Read-Count + 0.1)") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.position="top", legend.title = element_text(face="italic",size=rel(.7)))

	# 1. count & p-val
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.exp.exon.pval.by.sex.dual.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=8,units="in",res=300, compression = 'lzw')
	ggplot_dual_axis(p.cnt,p.cnt.pval,"y")
	dev.off()

	# 2. log-count & p-val
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.log.cnt.pval.per.exon.by.sex.dual3.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=8,units="in",res=300, compression = 'lzw')
	ggplot_dual_axis(p.cnt3,p.cnt.meth.pval,"y") # p-val (exp) & p-val (meth)
	#ggplot_dual_axis(p.cnt3,p.cnt.pval,"y")
	#ggplot_dual_axis(p.cnt.pval,p.cnt3,"y")
	dev.off()
		
	# 3. meth & p-val
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.exon.pval.by.sex.dual.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=8,units="in",res=300, compression = 'lzw')
	ggplot_dual_axis(p.meth,p.cnt.pval,"y")
	dev.off()

	# 4. meth & p-val (meth)
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.exon.pval2.by.sex.dual.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=8,units="in",res=300, compression = 'lzw')
	ggplot_dual_axis(p.meth,p.meth.pval,"y")
	dev.off()

	# 5-1. cnt & meth
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.exp.exon.by.sex.dual.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=7,units="in",res=300, compression = 'lzw')
	ggplot_dual_axis(p.cnt2,p.meth2,"y")
	dev.off()

	# 5-2. meth & cnt
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.exp2.exon.by.sex.dual.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=7,units="in",res=300, compression = 'lzw')
	ggplot_dual_axis(p.meth2,p.cnt2,"y")
	dev.off()

	# 6. diff meth & diff exp
	hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
	#p.diff<-ggplot(dt.exp.meth.exon, aes(CNT.F-CNT.M, met.f-met.m)) + 
	p.diff<-ggplot(dt.exp.meth.exon[(CNT.F+CNT.M)/2>1 & num.sites>1], aes(CNT.F-CNT.M, met.f-met.m)) + 
		geom_point(aes(size=(CNT.F+CNT.M)/2, col=met)) +
		scale_size_continuous(name="Average\nRead Count") +
		scale_colour_gradientn(name="Average\n% Methylation",colors=hmcol) +
		labs(x="Read Count Difference (Female-Male)", y="% Methylation Difference (Female-Male)") +
		geom_hline(aes(yintercept=0)) +
		geom_vline(aes(xintercept=0)) +
		theme_Publication()
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.diff.exp.by.exon.by.sex.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=8.5,units="in",res=300, compression = 'lzw')
	print(p.diff)
	dev.off()

	##################
	## P-value plot ##
	##################
	hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)

	p1<-ggplot(dt.fpkm.exon, aes(exon_id, pvalue, group=1)) + 
		geom_point(aes(size=(CNT.F+CNT.M)/2, colour=log(FPKM.F/FPKM.M)),na.rm=T) + 
		geom_line(linetype="dotted",na.rm=T) + 
		geom_hline(aes(yintercept=my.pval)) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_size_continuous(name="Mean\nRead Count") +
		#scale_colour_gradient2(mid="grey", high="red", low="blue", midpoint=0) +
		scale_colour_gradientn(name="logFC",colors=hmcol) +
		scale_y_continuous(breaks=c(0, my.pval, 0.25, .5, .75, 1)) +
		labs(x="Exon ID", y="P-value") +
		#theme_Publication() + 
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.title = element_text(face="italic",size=rel(.9)))
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.exp.exon.pval.by.sex.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=4,units="in",res=300, compression = 'lzw')
	print(p1)
	dev.off()

	p2<-ggplot(dt.fpkm.exon, aes(exon_id, pvalue, group=1)) + 
		geom_point(aes(size=(CNT.F+CNT.M)/2, colour=(CNT.F-CNT.M)),na.rm=T) + 
		geom_line(linetype="dotted",na.rm=T) + 
		geom_hline(aes(yintercept=my.pval)) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_size_continuous(name="Mean\nRead Count") +
		scale_colour_gradient2(name="Count Difference\n(Female-Male)",mid="gray48", high="red", low="blue", midpoint=0) +
		scale_y_continuous(breaks=c(0, my.pval, 0.25, .5, .75, 1)) +
		labs(x="Exon ID", y="P-value") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.title = element_text(face="italic",size=rel(.5)))
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.exp.cnt.exon.pval.by.sex.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=4,units="in",res=300, compression = 'lzw')
	print(p2)
	dev.off()

	########################
	# group-level pheatmap #
	########################
	library(pheatmap)
	merged.fpkm<-cbind(`Female`=rowMeans(ddsFpkm[,rn.case]), `Male`=rowMeans(ddsFpkm[,rn.control]))
	rownames(merged.fpkm)=rownames(ddsFpkm)

	pval<-data.frame(`p.value`=ifelse(rownames(ddsFpkm) %in% names(exon.pval[!is.nan(exon.pval) & exon.pval<my.pval]), paste('p<',my.pval), paste('p>=',my.pval), row.names=rownames(ddsFpkm)))
	ann_colors<-list(`p.value`=c(`p< 0.05`='black',`p>= 0.05`='grey'))

	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.domain.fpkm.pheatmap.by.sex.",TR_PREFIX,".tiff"))
	pheatmap(t(merged.fpkm), cluster_rows=F, cluster_cols=F, annotation_col=pval, annotation_colors=ann_colors, main=paste("FPKMs of CSMD1 Exons by Sex (", TR_PREFIX,")"), filename=my.filename, width=10, height=2.5, fontsize=7, cellheight=10)
	pheatmap(t(merged.fpkm), cluster_rows=F, cluster_cols=F, annotation_col=pval, annotation_colors=ann_colors, main=paste("FPKMs of CSMD1 Exons by Sex (", TR_PREFIX,")"), width=10, height=2.5, fontsize=7, cellheight=10)

	#########################
	# sample-level pheatmap #
	#########################
	library(pheatmap)
	# Schultz.ddsFpkm is 'ddsFpkm' of RoadMap Schultz tissues (no placenta)
	# RoadMap.ddsFpkm is 'ddsFpkm' of RoadMap Consortium placenta data
	# POPS.ddsFpkm is 'ddsFpkm' of Boy.Girl.exon.GRCh37 AGA samples
	merged.fpkm<-cbind(Schultz.ddsFpkm, RoadMap.ddsFpkm, `PT`=rowMeans(POPS.ddsFpkm))

	my.filename<-"~/results/RNA-Seq/RoadMap.exon.GRCh37/Cluster/Consortium/RoadMap.POPS.CSMD1.domain.fpkm.pheatmap.by.tissue.GRCh37.tiff"
	pheatmap(t(log(merged.fpkm/10^3+1)), cluster_rows=T, cluster_cols=F, width=10, height=5, fontsize_col=7, fontsize_row=10, cellheight=9, main="Log(FPKM/10^3+1)", filename=my.filename)
	pheatmap(t(log(merged.fpkm/10^3+1)), cluster_rows=T, cluster_cols=F, width=10, height=5, fontsize_col=7, fontsize_row=10, cellheight=9, main="Log(FPKM/10^3+1)")

	#merged.cnt<-cbind(counts(RoadMap.dds), `PT`=rowMeans(counts(dds)))
	#pheatmap(t(merged.cnt), cluster_rows=F, cluster_cols=F, width=10, height=2.5, fontsize=7, cellheight=10)

	if(FALSE){
		##########################
		# sample-level heatmap.2 #
		#######################3##
		# exons of P<my.pvalmy.pval
		exon.col<-ifelse(rownames(ddsFpkm) %in% names(exon.pval[!is.nan(exon.pval) & exon.pval<my.pval]), 'red', 'black')
		names(exon.col)<-rownames(ddsFpkm)

		sample.col<-c(rep(my.col[[my.contrast]][names(my.col[[my.contrast]])=="F"], length(rn.case)),
					rep(my.col[[my.contrast]][names(my.col[[my.contrast]])=="M"], length(rn.control))
					)
		names(sample.col)<-names(c(rn.case, rn.control))
		#heatmap.2(log(t(ddsFpkm[,c(rn.case,rn.control)])+0.01), Rowv=FALSE, Colv=FALSE, scale="none", trace="none", dendrogram="none", col=hmcol, keysize=1,  key.xlab="log(FPKM)", main="FPKMs of CSMD1 domains", margins=c(8,4), srtCol=45, colRow=sample.col, colCol=exon.col)
		heatmap.2(t(ddsFpkm[,c(rn.case,rn.control)]), Rowv=FALSE, Colv=FALSE, scale="none", trace="none", dendrogram="none", col=hmcol, keysize=1,  key.xlab="FPKM", main="FPKMs of CSMD1 domains", margins=c(8,4), srtCol=45, colRow=sample.col, colCol=exon.col)

		#########################
		# group-level heatmap.2 #
		#########################
		sample.col<-cbind(`Female`=my.col[[my.contrast]][names(my.col[[my.contrast]])=="F"], `Male`=my.col[[my.contrast]][names(my.col[[my.contrast]])=="M"])
		merged.fpkm<-cbind(`Female`=rowMeans(ddsFpkm[,rn.case]), `Male`=rowMeans(ddsFpkm[,rn.control]))
		#merged.fpkm<-cbind(`F`=rowMeans(counts(dds,normalized=T)[,rn.case]), `M`=rowMeans(counts(dds,normalized=T)[,rn.control]))
		rownames(merged.fpkm)=rownames(ddsFpkm)

		my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.domain.fpkm.heatmap.by.sex.",TR_PREFIX,".v2.tiff"))
		tiff(filename=my.filename,width=10,height=3,units="in",res=300, compression = 'lzw')
		heatmap.2(t(merged.fpkm), Rowv=FALSE, Colv=FALSE, scale="none", trace="none", dendrogram="none", col=hmcol, keysize=1,  key.xlab="FPKM", main="Mean FPKMs of CSMD1 Domains", srtCol=45, colRow=sample.col, colCol=exon.col, offsetCol=0.1)
		dev.off()
	}

	################################
	## CSMD1 expression by tissue ##
	## PT + RoadMap
	################################
	my.exp.RData="~/results/RoadMap/X-inactivation/RData/dt.exp.roadmap.RData"
	cat(paste0("loading ",my.exp.RData,"...\n"))
	load (my.exp.RData) 
	cat("dt.tissue.exp.list loaded\n")

	rbindlist(dt.tissue.exp.list)[ensembl_gene_id==my.ensg,.(tissue,female.exp,male.exp)]
	dt.csmd.fpkm<-reshape2::melt(rbindlist(dt.tissue.exp.list)[ensembl_gene_id==my.ensg,.(tissue,female.exp,male.exp)], id="tissue", measure=c("female.exp", "male.exp"), variable.name="Gender", value.name="FPKM")
	dt.csmd.fpkm[,Gender:=ifelse(Gender=="female.exp","Female","Male")]

	my.filename <- file.path("~/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/Cluster/AGA/CSMD1.PT.RoadMap.by.sex.tiff")
	tiff(filename=my.filename,width=12,height=10,units="in",res=300, compression = 'lzw')
	ggplot(dt.csmd.fpkm, aes(tissue, FPKM, fill=Gender)) + geom_bar(stat="identity", position="dodge") + scale_fill_manual(name="Sex",values=my.col[["Gender"]]) + theme_Publication()
	dev.off()

	###############
	## Via ggbio ##
	###############
	if(FALSE){
		library(ggbio)
		####################
		# transcript track # 
		####################
		p.tr<-autoplot(hg19ensGene, which=GRanges(seqnames="chr8", range=ranges(range(my.target))), fill = "orange", color = "grey") + theme_Publication() # hg19ensGene defined in 'Annotation.R'

		my.target$order=seq(5,1000,by=3)[1:nrow(mcols(my.target))] # amplicon left-most order 
		names(my.target)<-my.target$exon_id
		p1<-ggplot(my.target) +  geom_rect(stat="identity", rect.height = 2, aes(y=order), col="black", ylab=NULL) 
			#scale_fill_gradient2(name="coveage(%)",low=hmcol[1], mid=hmcol[round(length(hmcol)/2)], high=hmcol[length(hmcol)], midpoint=50) + 
			#theme(axis.ticks.y = element_blank(), axis.text.y = element_blank() ) +
			#geom_text(x=start(gr.target)+250, y=seq(5,200,by=3)[1:nrow(mcols(gr.target))], label=names(gr.target), size=2.6)
		gr.reduced<-reduce(my.target)
		gr.reduced$order=max(my.target$order)+3 
		p2<-ggplot(gr.reduced) +  geom_rect(stat="identity", rect.height = 2, aes(y=order), col="blue", ylab=NULL)
	}


}

# FG plasma samples 
if(myProject=="Plasma"){
	library(data.table)

	# get chr and gene name
	cat("getting gene annotations...\n")
	my.anno <- getBM(attributes = my.fields, filters = my.filter, values = rownames(dds), mart = myMart)
	
	# get 12wk and 36wk plasma samples
	rn.12wk=rownames(colData(dds)[colData(dds)$GA=="12",])
	rn.36wk=rownames(colData(dds)[colData(dds)$GA=="36",])

	# set up data.table for FPKM by group
	#dt.meanFpkm<-data.table(`ensembl_gene_id`=rownames(ddsFpkm), `Plasma.12wk`=rowMeans(ddsFpkm[,rn.12wk]), `Plasma.36wk`=rowMeans(ddsFpkm[,rn.36wk]), `Plasma.all`=rowMeans(ddsFpkm[,rn.case]), `Placenta`=rowMeans(ddsFpkm[,rn.control]))
	dt.meanFpkm<-data.table(`ensembl_gene_id`=rownames(ddsFpkm), `Plasma.2008.12wk`=ddsFpkm[,"85PL"], `Plasma.2008.36wk`=ddsFpkm[,"87PL"], `Plasma.2012.12wk`=ddsFpkm[,"89PL"], `Plasma.2012.36wk`=ddsFpkm[,"91PL"], `Plasma.12wk`=rowMeans(ddsFpkm[,rn.12wk]), `Plasma.36wk`=rowMeans(ddsFpkm[,rn.36wk]), `Plasma.all`=rowMeans(ddsFpkm[,rn.case]), `Placenta`=rowMeans(ddsFpkm[,rn.control]))
	dt.meanFpkm.NR<-data.table(`ensembl_gene_id`=rownames(ddsFpkm.NR), `Plasma.2008.12wk`=ddsFpkm.NR[,"85PL"], `Plasma.2008.36wk`=ddsFpkm.NR[,"87PL"], `Plasma.2012.12wk`=ddsFpkm.NR[,"89PL"], `Plasma.2012.36wk`=ddsFpkm.NR[,"91PL"], `Plasma.12wk`=rowMeans(ddsFpkm.NR[,rn.12wk]), `Plasma.36wk`=rowMeans(ddsFpkm.NR[,rn.36wk]), `Plasma.all`=rowMeans(ddsFpkm.NR[,rn.case]), `Placenta`=rowMeans(ddsFpkm.NR[,rn.control]))
	#dt.meanFpm<-data.table(`ensembl_gene_id`=rownames(ddsFpm), `Plasma.2008.12wk`=ddsFpm[,"85PL"], `Plasma.2008.36wk`=ddsFpm[,"87PL"], `Plasma.2012.12wk`=ddsFpm[,"89PL"], `Plasma.2012.36wk`=ddsFpm[,"91PL"], `Plasma.12wk`=rowMeans(ddsFpm[,rn.12wk]), `Plasma.36wk`=rowMeans(ddsFpm[,rn.36wk]), `Plasma.all`=rowMeans(ddsFpm[,rn.case]), `Placenta`=rowMeans(ddsFpm[,rn.control]))
	#dt.meanCount<-data.table(`ensembl_gene_id`=rownames(counts(dds)), `Plasma.2008.12wk`=counts(dds)[,"85PL"], `Plasma.2008.36wk`=counts(dds)[,"87PL"], `Plasma.2012.12wk`=counts(dds)[,"89PL"], `Plasma.2012.36wk`=counts(dds)[,"91PL"], `Plasma.12wk`=rowMeans(counts(dds)[,rn.12wk]), `Plasma.36wk`=rowMeans(counts(dds)[,rn.36wk]), `Plasma.all`=rowMeans(counts(dds)[,rn.case]), `Placenta`=rowMeans(counts(dds)[,rn.control]))

	# set up data.table for FPKM by sample 
	#dt.Fpkm<-data.table(cbind(`ensembl_gene_id`=rownames(ddsFpkm), ddsFpkm))
	#dt.Fpkm[ensembl_gene_id %in% my.dt.meanFpkm[,ensembl_gene_id]]

	minFPKM=1
	###############################
	# 1. Plasma samples chrY only #
	###############################
	# chrY expressed either plasma or placenta (per group)
	my.dt.meanFpkm<-merge(dt.meanFpkm[ensembl_gene_id %in% my.anno[my.anno$chromosome_name=="Y","ensembl_gene_id"] & Plasma.all + Placenta >0],my.anno, by="ensembl_gene_id")
	my.dt.meanFpkm.NR<-merge(dt.meanFpkm.NR[ensembl_gene_id %in% my.anno[my.anno$chromosome_name=="Y","ensembl_gene_id"] & Plasma.all + Placenta >0],my.anno, by="ensembl_gene_id")
	my.dt.meanCount<-merge(dt.meanCount[ensembl_gene_id %in% my.anno[my.anno$chromosome_name=="Y","ensembl_gene_id"] & Plasma.all + Placenta >0],my.anno, by="ensembl_gene_id")

	write.csv(my.dt.meanFpkm[Plasma.all>minFPKM][order(-Plasma.all)], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.more.than.",minFPKM,".per.group.csv.gz"))))
	# chrY expressed either plasma or placenta (per sample)
	write.csv(ddsFpkm[my.dt.meanFpkm[Plasma.all + Placenta >0 ,ensembl_gene_id],], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.more.than.0.per.sample.csv.gz"))))
	# chrY expressed from plasma (per sample)
	write.csv(ddsFpkm[my.dt.meanFpkm[Plasma.all >0 ,ensembl_gene_id],], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".chrY.plasma.fpkm.more.than.0.per.sample.csv.gz"))))

	# VennDiagram between Plasma and Placent chrY 
	library(VennDiagram)
	venn.list=list()
	venn.list[["Plasma"]]=my.dt.meanFpkm[Plasma.all>minFPKM,ensembl_gene_id]
	venn.list[["Placenta"]]=my.dt.meanFpkm[Placenta>minFPKM,ensembl_gene_id]
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".venn.plasma.pt.fpkm.more.than.",minFPKM,".tiff"))
	venn.diagram(
		x=venn.list,
		filename = my.filename,
		col = "black",
		fill = c("dodgerblue", "goldenrod1"),
		alpha = 0.50,
		cat.col = c("dodgerblue", "goldenrod1"),
		cat.cex = 1.1,
		cat.fontface = "bold",
		margin = 0.05
	)

	# drop genes of all zero fpkm across samples
	myFpkm<-ddsFpkm[my.anno[my.anno$chromosome_name=="Y","ensembl_gene_id"],] # chrY only

	# heatmap of any chrY expressed
	#keep = rowSums(myFpkm) > 0  
	#myTopExp <- getTopExp(myFpkm[keep,], iqr=FALSE, nrow(myFpkm[keep,])) # isa 'list'

	# heatmap by samples
	# plasma chrY expressed sorted by rowMeans (N=47)
	myFpkm<-ddsFpkm[my.dt.meanFpkm[Plasma.all>0,ensembl_gene_id],]
	myTopExp <- getTopExp(myFpkm, iqr=FALSE, nrow(myFpkm)) # isa 'list'

	my.filename <- file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.plasma.heatmap"))
	tiff(filename=paste0(my.filename,".tiff"),width=10,height=12,units="in",res=300, compression = 'lzw')
	heatmap.2(as.matrix(log(myTopExp$fpkm+0.01)), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="log(FPKM)", main="Plasma ChrY Genes of FPKM>0", margins=c(4,8))
	dev.off()

	my.filename <- file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.plasma.heatmap.tr"))
	tiff(filename=paste0(my.filename,".tiff"),width=13,height=10,units="in",res=300, compression = 'lzw')
	heatmap.2(as.matrix(log(t(myTopExp$fpkm)+0.01)), Colv=FALSE, scale="none", trace="none", dendrogram="row", col = hmcol, keysize=1, key.xlab="log(FPKM)", main="Plasma ChrY Genes of FPKM>0", margins=c(8,4), srtCol=45)
	dev.off()

	# heatmap of plasma chrY expressed clustered by genes and samples (N=47)
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.plasma.row.dendro.heatmap"))
	tiff(filename=paste0(my.filename,".tiff"),width=10,height=12,units="in",res=300, compression = 'lzw')
	heatmap.2(as.matrix(log(myTopExp$fpkm+0.01)), Rowv=TRUE, scale="none", trace="none", dendrogram="both", col = hmcol, keysize=1, key.xlab="log(FPKM)", main="Plasma ChrY Genes of FPKM>0", margins=c(4,8))
	dev.off()

	# heatmap sorted by plasma 
	foo<-myTopExp$fpkm[names(sort(rowMeans(myTopExp$fpkm[,rn.case]),decreasing=T)),] # re-order by FPKM of Plasma
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.plasma.sorted.heatmap"))
	tiff(filename=paste0(my.filename,".tiff"),width=10,height=12,units="in",res=300, compression = 'lzw')
	heatmap.2(as.matrix(log(foo+0.01)), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="log(FPKM)", main="Plasma ChrY Genes of FPKM>0", margins=c(4,8))
	dev.off()

	# heatmap by group
	for(i in c(0, 0.1, 0.5, 1)){
		foo<-as.matrix(my.dt.meanFpkm[Plasma.all>i, .(Plasma.2008.12wk, Plasma.2008.36wk, Plasma.2012.12wk, Plasma.2012.36wk, Plasma.12wk, Plasma.36wk, Plasma.all, Placenta)])
		rownames(foo)<-my.dt.meanFpkm[Plasma.all>i, ensembl_gene_id]
		myTopExp <- getTopExp(foo, iqr=FALSE, nrow(foo)) # isa 'list'

		my.filename <- file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.more.than.",i,".plasma.heatmap"))
		tiff(filename=paste0(my.filename,".tiff"),width=10,height=12,units="in",res=300, compression = 'lzw')
		heatmap.2(as.matrix(log(myTopExp$fpkm+0.01)), Colv=FALSE, Rowv=FALSE, scale="none", trace="none", dendrogram="none", col = hmcol, keysize=1, key.xlab="log(FPKM)", main=paste("ChrY Genes of Avg. Plasma FPKM>",i), margins=c(8,10),  srtCol=45, cellnote=round(myTopExp$fpkm,4), notecol="black")
		dev.off()
	}



	# exploratory plots for plasma chrY genes FPKM>0
	ggplot(reshape2::melt(my.dt.meanFpkm[Plasma.all>0,.(ensembl_gene_id,Plasma.all, Placenta)], id="ensembl_gene_id", variable.name="Tissue", value.name="FPKM"), aes(Tissue, FPKM)) + geom_boxplot()
	ggplot(reshape2::melt(my.dt.meanFpkm[Plasma.all>0,.(ensembl_gene_id,Plasma.all, Placenta)], id="ensembl_gene_id", variable.name="Tissue", value.name="FPKM"), aes(FPKM)) + geom_density(aes(col=Tissue))

	################################
	## 2. Placenta Specific Genes ##
	################################
	# from Steve v1.
	my.ensg=c(`PSG1`="ENSG00000231924", `PSG2`="ENSG00000242221", `PSG3`="ENSG00000221826", `PSG4`="ENSG00000243137",
			`PSG5`="ENSG00000204941", `PSG6`="ENSG00000170848", `PSG7`="ENSG00000221878", `PSG8`="ENSG00000124467",
			`PSG9`="ENSG00000183668", `PSG10P`="ENSG00000248257", `PSG11`="ENSG00000243130", 
			`CGA`="ENSG00000135346", `CGB1`="ENSG00000267631", `CGB2`="ENSG00000104818", `CGB3`="ENSG00000104827",
			`CGB5`="ENSG00000189052", `CGB7`="ENSG00000196337", `CGB8`="ENSG00000213030", `CSH1`="ENSG00000136488",
			`CSH2`="ENSG00000213218", `PLAC1`="ENSG00000170965", `PLAC4`="ENSG00000280109", `PLAC8`="ENSG00000145287",
			`PAPPA`="ENSG00000182752", `PAPPA2`="ENSG00000116183"
			)
	# from Steve v2.
	my.ensg<-read.table("~/Pipelines/data/RNA-Seq/placenta.specific.gene.txt", header=F,col.names=c("ensembl_gene_id"),stringsAsFactors=F)$ensembl_gene_id
	# from Steve v3.
	my.ensg<-read.table("~/Pipelines/data/RNA-Seq/placenta.specific.coding.gene.txt", header=F,col.names=c("ensembl_gene_id"),stringsAsFactors=F)$ensembl_gene_id
	# ensg available from GRCh38
	dummy<-my.anno[my.anno$ensembl_gene_id %in% dt.meanFpkm[ensembl_gene_id %in% my.ensg & Plasma.all>0,ensembl_gene_id],c("ensembl_gene_id","gene_biotype")]
	# protein_coding gene only
	my.ensg<-dummy[dummy$gene_biotype=="protein_coding","ensembl_gene_id"]
	#heatmap by group
	for(i in c(0, 0.1, 0.5, 1)){
		foo<-as.matrix(dt.meanFpkm[ensembl_gene_id %in% my.ensg & Plasma.all>i, .(Plasma.2008.12wk, Plasma.2008.36wk, Plasma.2012.12wk, Plasma.2012.36wk, Plasma.12wk, Plasma.36wk, Plasma.all, Placenta)])
		rownames(foo)<-dt.meanFpkm[ensembl_gene_id %in% my.ensg & Plasma.all>i, ensembl_gene_id]
		myTopExp <- getTopExp(foo, iqr=FALSE, nrow(foo)) # isa 'list'

		my.filename <- file.path(cluster.dir, paste0(sampleType, ".placenta.specific.genes.fpkm.more.than.",i,".by.group.heatmap"))
		tiff(filename=paste0(my.filename,".tiff"),width=10,height=12,units="in",res=300, compression = 'lzw')
		heatmap.2(as.matrix(log(myTopExp$fpkm+0.01)), Colv=FALSE, Rowv=FALSE, scale="none", trace="none", dendrogram="none", col = hmcol, keysize=1, key.xlab="log(FPKM)", main=paste("Placenta Specific Genes of Avg Plasma FPKM>",i), margins=c(9,9),  srtCol=45, cellnote=round(myTopExp$fpkm,2), notecol="black")
		dev.off()
	}


	my.dt.meanFpkm<-merge(dt.meanFpkm[ensembl_gene_id %in% my.ensg & Plasma.all>0],my.anno, by="ensembl_gene_id")
	my.dt.meanCount<-merge(dt.meanCount[ensembl_gene_id %in% my.ensg & Plasma.all>0],my.anno, by="ensembl_gene_id")

	# FPKM per group 
	write.csv(my.dt.meanFpkm[Plasma.all>0][order(-Plasma.all)], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".placenta.specific.genes.fpkm.more.than.0.per.group.csv.gz"))))
	write.csv(my.dt.meanCount[Plasma.all>0][order(-Plasma.all)], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".placenta.specific.genes.read.count.more.than.0.per.group.csv.gz"))))
} # end of myProject=="Plasma"

if(FALSE){
	library(vsn) # for 'meanSdPlot'
	#resultsNames(dds) # inspect the result name
						# ‘resultsNames’ returns
						# the names of the estimated effects (coefficents) of the model

	#colnames(mcols(dds)) # display columns 
	#mcols(dds)[1,] # display the first gene
	#> table(mcols(dds)$allZero) # number of genes having zero count for all samples 
	#FALSE  TRUE 
	#42615 21062 
	# Scatter plot of sample 2 vs sample 1
	cat("Scatter plot ...","\n")
	par( mfrow = c( 1, 2 ) )
	plot( log2( 1+counts(dds, normalized=T)[, 1:2] ), col="#00000020", pch=20, cex=0.3 ) # log2 transform
	plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 ) # rlog transform

	# Per-gene standard deviation
	cat("Per-gene standard deviation...","\n")
	par(mfrow=c(1,3))
	notAllZero <- (rowSums(counts(dds))>0)
	meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1), ylim = c(0,2.5)) # requires library(vsn)
	meanSdPlot(assay(rld[notAllZero,]), ylim = c(0,2.5))
	meanSdPlot(assay(vsd[notAllZero,]), ylim = c(0,2.5))

	# without independant filtering
	if(grepl("^Boy.Girl",myProject)){
		resNoFilt <- results(ddsFlt, independentFiltering=FALSE) # this will disable independant filtering
	}else{
		resNoFilt <- results(dds, independentFiltering=FALSE) # this will disable independant filtering
	}
	table(filtering=(my.res$padj < .1), noFiltering=(resNoFilt$padj < .1))
	#         noFiltering
	#filtering FALSE  TRUE
	#    FALSE 19754     0
	#    TRUE     18    15


	# count genes satisfying cuffoff
	cat("counting genes satisfying cuffoff","\n")
	sum(my.res$pvalue < myPValue, na.rm=T)
	sum(my.res$padj < myFDR, na.rm=T)

	# plot : Mean counts as a filter statistic 
	plot(my.res$baseMean+1, -log10(my.res$pvalue),
		log="x", xlab="mean of normalized counts",
		ylab=expression(-log[10](pvalue)),
		ylim=c(0,30),
		cex=.4, col=rgb(0,0,0,.3))

	# browse the independent filtering
	attr(my.res,"filterThreshold")
	#68.92687% : % of genes having less than the number below
	#2.703189 : a cutoff meanBase
	use <- my.res$baseMean > attr(my.res,"filterThreshold")
	table(use)
	h1 <- hist(my.res$pvalue[!use], breaks=0:50/50, plot=FALSE)
	h2 <- hist(my.res$pvalue[use], breaks=0:50/50, plot=FALSE)
	colori <- c(`do not pass`="khaki", `pass`="powderblue")
	barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori, space = 0, main = "", ylab="frequency")
	text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
	legend("topright", fill=rev(colori), legend=rev(names(colori)))


	# Independent filtering. DESeq2 automatically determines a threshold, filtering on mean
	# normalized count, which maximizes the number of genes which will have an adjusted p value less than
	# a critical value.
	cat("plot baseMean","\n")
	plot(attr(my.res,"filterNumRej"),type="b", xlab="quantiles of 'baseMean'",ylab="number of rejections")

	# Ratio of small p values for groups of genes binned by mean normalized count
	# create bins using the quantile function
	qs <- c( 0, quantile( my.res$baseMean[my.res$baseMean > 0], 0:7/7 ) )
	# "cut" the genes into the bins
	bins <- cut( my.res$baseMean, qs )
	# rename the levels of the bins using the middle point
	levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
	# calculate the ratio of p-values less than .01 for each bin
	ratios <- tapply( my.res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
	# plot these ratios
	barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")

	# 4. More information on results columns
	plotDispEsts( dds, ylim = c(1e-6, 1e1) ) # Plot of dispersion estimates
	hist( my.res$pvalue, breaks=100, col="grey" ) # Histogram of the p values returned by the test for differential expression 
	cat("mcol(res)$description ","\n")
	mcols(my.res)$description
}#end of Plasma

# plot FC over read-count
cat("plot FC over read-count","\n")

###########################################
# conventional filtering (FDR<0.1 & FC>2)
###########################################
deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj) & (my.res$log2FoldChange > 1 | my.res$log2FoldChange < -1) & my.res$padj <= 0.1, ]
if(nrow(deseq.top.deg) > 1){
	plotMA(my.res, main=paste0("conventional filtering (FDR<=0.1 & |logFC|>1): ", nrow(deseq.top.deg), " genes (", table(my.res$log2FoldChange<0)["TRUE"], "<0, ", table(my.res$log2FoldChange<0)["FALSE"], ">=0)"), ylim=c(-2,2)) # default q-value: alpha = 0.1
}

#################
# for all genes
#################
cat("no filtering","\n")
deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange), ]
plotMA(my.res, main=paste0("No filtering: ", nrow(deseq.top.deg), " genes (", table(my.res$log2FoldChange<0)["TRUE"], "<0, ", table(my.res$log2FoldChange<0)["FALSE"], ">=0)"), ylim=c(-2,2)) # default q-value: alpha = 0.1
abline(v=c(0.01,1,minRead), h=c(-1,1), col="blue")

hist(my.res$pvalue, breaks=100, xlab='P-value', main=paste0('P value distribution of ', nrow(my.res[ !is.na(my.res$pvalue),]), ' genes'))
hist(my.res[my.res$baseMean>=minRead,]$pvalue, breaks=100, xlab='P-value', main=paste0('P value distribution of ', nrow(my.res[ my.res$baseMean>=minRead,]), ' genes  of minRead>=', minRead))

#######################################
# filter genes by mean count & p-value
#######################################
cat("filter by mean count & p-value\n")
deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$pvalue) & my.res$baseMean>=minRead & my.res$pvalue<=myPValue, ]
if(nrow(deseq.top.deg) > 1){
	hist(deseq.top.deg$padj, main=paste0("FDR distribution of ", nrow(deseq.top.deg) ," genes having P-value<=",myPValue))

	plotMA(deseq.top.deg, main=paste0(nrow(deseq.top.deg), " genes (",table(deseq.top.deg$log2FoldChange<0)["TRUE"] , "<0, ", table(deseq.top.deg$log2FoldChange<0)["FALSE"], ">=0) having mean>=",minRead, "& P-value<",myPValue), ylim=c(-2,2), alpha=0) # default q-value: alpha = 0.1
	abline(v=c(minRead), col="blue")
}

##############################################
# filter genes by mean count & p-value & FDR
##############################################
cat("filter by mean count & p-value & FDR\n")
deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$pvalue) & !is.na(my.res$padj) & my.res$baseMean>=minRead & my.res$pvalue<=myPValue & my.res$padj<=myFDR, ]
if(nrow(deseq.top.deg) > 1){
	plotMA(deseq.top.deg, main=paste0(nrow(deseq.top.deg), " genes (",table(deseq.top.deg$log2FoldChange<0)["TRUE"] , "<0, ", table(deseq.top.deg$log2FoldChange<0)["FALSE"], ">=0) having mean>=",minRead, "& P-value<=",myPValue," & FDR<=",myFDR), ylim=c(-2,2), alpha=0) # default q-value: alpha = 0.1
	abline(v=c(minRead), col="blue")
}

################
## Annotation ##
################
if(FALSE){
	cat("making deseq.top.deg.ens\n")
	fields <- c("ensembl_gene_id", my.id,"description", "gene_biotype")
	deseq.top.deg.ens=getBM(attributes = fields, filters = "ensembl_gene_id", values = rownames(deseq.top.deg), mart = myMart)

	cat("making deseq.top.deg.pheno\n")
	fields <- c("ensembl_gene_id", my.id,"description", "gene_biotype", "phenotype_description")
	deseq.top.deg.pheno=getBM(attributes = fields, filters = "ensembl_gene_id", values = rownames(deseq.top.deg), mart = myMart)

	cat("making deseq.top.deg.go\n")
	fields <- c("ensembl_gene_id", my.id, "go_id", "name_1006", "definition_1006", "namespace_1003")
	deseq.top.deg.go=getBM(attributes = fields, filters = "ensembl_gene_id", values = rownames(deseq.top.deg), mart = myMart)

	# 6. join two data.frame by a common column
	deseq.top.deg$my.filter<-rownames(deseq.top.deg)
	deseq.top.deg.anno=as.data.frame(merge(deseq.top.deg, deseq.top.deg.ens, by=my.filter))
	minDeseqlogFC=round(min(abs(deseq.top.deg.anno$log2FoldChange)),3)
	maxDeseqFDR=round(max(deseq.top.deg$padj),3)

	####################
	## Plot Gene Name ##
	####################
	# assign color
	# myCol defined within the config
	deseq.top.deg.anno<-assignFDRcol(deseq.top.deg.anno) # assignFDRcol defined within the config
	plot(y=deseq.top.deg.anno$log2FoldChange, x=log(deseq.top.deg.anno$baseMean), main=paste0("CPM>=1: ",nrow(deseq.top.deg)," top DEGs for ",sampleType," (logFC>=",minDeseqlogFC," FDR<=",maxDeseqFDR,")"), ylab='logFC', xlab='log10(baseMean)', cex=0.2)
	text(y=deseq.top.deg.anno$log2FoldChange, x=log(deseq.top.deg.anno$baseMean), labels=deseq.top.deg.anno[[my.id]], col=deseq.top.deg.anno$col, cex=0.75)
	abline(h=c(-minDeseqlogFC, minDeseqlogFC), col="blue")
	#legend("topright", legend=names(myCol), fill=myCol)

	# plot3D
	#library(rgl)
	#plot3d(y=deseq.top.deg.anno$padj, z=deseq.top.deg.anno$log2FoldChange, x=log(deseq.top.deg.anno$baseMean), ylab='FDR', zlab='logFC', xlab='log10(baseMean)')
	#text3d(y=deseq.top.deg.anno$padj, z=deseq.top.deg.anno$log2FoldChange, x=log(deseq.top.deg.anno$baseMean), deseq.top.deg.anno[[my.id]], col=deseq.top.deg.anno$col, cex=2)
}#end of FALSE

#########################################
## Top DEGs (meanBase>=10 & padj<=0.05) ##
#########################################
cat("Top DEGs (meanBase>=10) of padj<=0.05\n")
#deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj) & my.res$baseMean>=minRead & my.res$padj<=0.05, ] # at least 20 read
deseq.top.deg=my.res[!is.na(my.res$padj) & my.res$padj<=0.05, ]
if (nrow(deseq.top.deg)>0){getTopDeseq(deseq.top.deg,"padj.05")}
#########################################
## Top DEGs (meanBase>=10 & padj<=0.1) ##
#########################################
cat("Top DEGs (meanBase>=10) of padj<=0.1\n")
#deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj) & my.res$baseMean>=minRead & my.res$padj<=0.1, ] # at least 20 read
deseq.top.deg=my.res[!is.na(my.res$padj) & my.res$padj<=0.1, ]
if(nrow(deseq.top.deg)>0){ 
	getTopDeseq(deseq.top.deg,"padj.1")
	my.entry<-rownames(deseq.top.deg[order(deseq.top.deg$pvalue),])[1]
	plotExp(my.entry,"FPM",my.contrast,my.box=FALSE,my.save=T)
	plotExp(my.entry,"FPM",my.contrast,my.box=TRUE,my.save=T)
}
#########################################
## Top DEGs (meanBase>=20 & padj<=0.2) ##
#########################################
cat("Top DEGs (meanBase>=20) of padj<=0.2\n")
#deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj) & my.res$baseMean>=minRead & my.res$padj<=0.2, ] # at least 20 read
deseq.top.deg=my.res[!is.na(my.res$padj) & my.res$padj<=0.2, ]
if (nrow(deseq.top.deg)>0) getTopDeseq(deseq.top.deg,"padj.2")
#########################################
## Top DEGs (meanBase>=20 & padj<=0.3) ##
#########################################
cat("Top DEGs (meanBase>=20) of padj<=0.3\n")
#deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj) & my.res$baseMean>=minRead & my.res$padj<=0.3, ] # at least 20 read
deseq.top.deg=my.res[!is.na(my.res$padj) & my.res$padj<=0.3, ]
if (nrow(deseq.top.deg)>0) getTopDeseq(deseq.top.deg,"padj.3")
#########################################
## Top DEGs (meanBase>=20 & padj<=0.4) ##
#########################################
cat("Top DEGs (meanBase>=20) of padj<=0.4\n")
#deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj) & my.res$baseMean>=minRead & my.res$padj<=0.4, ] # at least 20 read
deseq.top.deg=my.res[!is.na(my.res$padj) & my.res$padj<=0.4, ]
if (nrow(deseq.top.deg)>0) getTopDeseq(deseq.top.deg,"padj.4")
##################################
## Top 100 genes (meanBase>=20) ##
##################################
cat("Top100 genes (meanBase>=20) by p-value\n")
#deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj) & my.res$baseMean>=minRead, ] # at least 20 read
deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj), ] # at least 20 read
deseq.top.deg<-head(deseq.top.deg[order(deseq.top.deg$pvalue),],n=100)
if (nrow(deseq.top.deg)>0) getTopDeseq(deseq.top.deg,"top100")
##################################
## Top 200 genes (meanBase>=20) ##
##################################
cat("Top200 genes (meanBase>=20) by p-value\n")
#deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj) & my.res$baseMean>=minRead, ] # at least 20 read
deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj), ] # at least 20 read
deseq.top.deg<-head(deseq.top.deg[order(deseq.top.deg$pvalue),],n=200)
if (nrow(deseq.top.deg)>0) getTopDeseq(deseq.top.deg,"top200")

#########################################
# save the result & FPKM into .csv file #
#########################################
cat("writing res and fpkm to file\n")
foo<-as.data.frame(my.res); foo[[my.filter]]<-rownames(my.res) # ensembl_gene_id|mirbase_id|pirbase_id  baseMean log2FoldChange      lfcSE       stat    pvalue      padj
bar<-data.frame(meanFpkm=rowMeans(ddsFpkm,na.rm=T), meanFpkm.case=rowMeans(ddsFpkm[,rn.case],na.rm=T), meanFpkm.control=rowMeans(ddsFpkm[,rn.control],na.rm=T)); bar[[my.filter]]<-rownames(bar)
write.csv(merge(foo, bar, all.x=TRUE), file=gzfile(paste0(deseq.dir,"/toptags_all_deseq.",sampleType,".csv.gz")))

cat("how the result looks like...\n")
head(my.res)

#################
# Samples Info ##
#################
write.csv(colData(dds), file=gzfile(file.path(deseq.dir,paste0(sampleType,".samples.txt.gz"))))


dev.off()
cat("All is done\n")
