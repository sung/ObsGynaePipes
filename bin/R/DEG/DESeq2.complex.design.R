#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# For detial, see http://www.nature.com/nprot/journal/v8/n9/full/nprot.2013.099.html
# complex design 14.D based on the publication above
# make sure to change my.local.db from config/Annotation.R

myCaller='DESeq2'
TR_PREFIX='GRCh38' # GRCh37 | GRCh38 
ENS_VER=82 # 75 (for igenome), 82 (for manual downlaod) 90 (for Cholestatis project)
source("~/Pipelines/config/DEG.R") # load config

cat(paste0("Running ", myProject, ":",sampleType,"...\n"))
# Set RData file
deseq.RData <- paste0(deseq.dir, "/deseq.",sampleType,'.RData')
# Load RData if exists
if(file.exists(deseq.RData)){
	cat("loading DESeq2 RData...\n")
	load(deseq.RData) # 
	cat("dds, res, rld?, loaded\n")
# create RData for the first time 
}else{
	####################
	# Read HTSeq count #
	####################
	cat("DESeqDataSetFromHTSeqCount()...\n")
	print(system.time(ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = samples, directory="", design= deseq.design ))) # isa 'DESeqDataSet'

	################## 
	## Filter genes ##
	################## 
	cat("filter out genes with no annotation...\n")
	#ddsKeep <- rownames(ddsHTSeq) %in% dt.ensg$ensembl_gene_id
	ddsKeep <- rownames(ddsHTSeq) %in% names(my.target.list) # remove entries of no GRange info
    if(!all(ddsKeep)){
	    newColData<- colData(ddsHTSeq)
        newCounts<- counts(ddsHTSeq)[ddsKeep,]
	    ddsHTSeq<- DESeqDataSetFromMatrix(newCounts, newColData, deseq.design) # isa 'DESeqDataSet'
    }
	##########################
	# Filter genes from chrY #
	##########################
	if(grepl("^Boy.Girl",myProject) | grepl("^GTEx",myProject)){
		if(grepl("^Boy.Girl",myProject) & sampleType!="BR"){
			cat("1. collapsing replicates...\n")
			# Error: sum(assay(object)) == sum(assay(collapsed)) is not TRUE (R version 3.1.1 (2014-07-10), DESeq2_1.6.3)
			# try with R 3.2.1 (DESeq2_1.10.1)
			# or collapseReplicatesBugFix
            if(packageVersion("DESeq2")<1.10){
                print(system.time(ddsCollapsed <- collapseReplicatesBugFix(ddsHTSeq, groupby=colData(ddsHTSeq)$CRN, run=colData(ddsHTSeq)$Source))) # isa 'DESeqDat# 'groupby' will be rownames of colData(dds)
																											# 'groupby' will be rownames of colData(dds)
																											# 'run' will be join by ','
            }else{
                print(system.time(ddsCollapsed <- collapseReplicates(ddsHTSeq, groupby=colData(ddsHTSeq)$CRN, run=colData(ddsHTSeq)$Source))) # isa 'DESeqDataSet'

            }
		}else{
            ddsCollapsed<-ddsHTSeq
		}
		# Boy.Girl.exon.GRCh38 (or GRCh37)
		if(grepl("exon",myProject)|grepl("iRNA",myProject)){ 
			ddsKeep <- rep(TRUE, nrow(ddsCollapsed))
		}else{
			cat("2. removing chrY...\n")
			#ddsKeep <- rownames(ddsCollapsed) %in% dt.ensg[chromosome_name!="Y",ensembl_gene_id]
			ddsKeep <- rep(TRUE, nrow(ddsCollapsed))
		}

		cat("DESeq from filtered count data...\n")
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
		cat("main DESeq2 pipeline: calculating dds...\n")
		print(system.time(dds<- DESeq(ddsMatrix, parallel=TRUE))) # isa 'DESeqDataSet'
    #################
    # RoadMap Exon ##
    #################
	}else if(grepl("^RoadMap.exon",myProject)){ 
		cat("1. collapsing replicates...\n")
        if(packageVersion("DESeq2")<1.10){
		    print(system.time(ddsCollapsed <- collapseReplicatesBugFix(ddsHTSeq, groupby=colData(ddsHTSeq)$Tissue))) # isa 'DESeqDat# 'groupby' will be rownames of colData(dds)
        }else{
            print(system.time(ddsCollapsed <- collapseReplicates(ddsHTSeq, groupby=colData(ddsHTSeq)$Tissue, run=colData(ddsHTSeq)$Source))) # isa 'DESeqDataSet'

        }
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
		cat("main DESeq2 pipeline: calculating dds...\n")
		print(system.time(dds<- DESeq(ddsMatrix))) # isa 'DESeqDataSet'
    ###########################
    ## Standard DEG analysis ##
    ###########################
	}else{
		#make sure the base level
		cat("main DESeq2 pipeline: calculating dds...\n")
		print(system.time(dds <- DESeq(ddsHTSeq, parallel=TRUE))) # isa 'DESeqDataSet'
	}
	cat("\ncalculating results...\n")
	print(system.time(res <- results(dds, parallel=TRUE))) # by default indepedant filtering at FDR (alpha) 0.1 (10%) 
											# independant filtering will set baseMean to filter genes automatically
											# by maximizing the number of genes satisfying < 0.1 FDR
											# isa 'DESeqResults'
	#print(system.time(res <- results(dds, parallel=TRUE, lfcThreshold=2))) # by default indepedant filtering at FDR (alpha) 0.1 (10%) 
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
			ddsKeep <- rownames(ddsHTSeq) %in% dt.ensg[chromosome_name!="Y",ensembl_gene_id] # remove ensg of chrY
		}else{
			ddsKeep = rowSums(ddsCpms >1) >= nrow(samples)/2
		}
		cat("\nDESeq from filtered count data ...\n")
		newCounts <- counts(ddsHTSeq)[ddsKeep,]
		ddsMatrix <- DESeqDataSetFromMatrix(newCounts, colData(ddsHTSeq), deseq.design) # isa 'DESeqDataSet'
		ddsFlt <- DESeq(ddsMatrix, parallel=TRUE) # isa 'DESeqDataSet'
		resFlt <- results(ddsFlt, parallel=TRUE) # by default indepedant filtering at FDR (alpha) 0.1 (10%) 
	}#end of FALSE

	####################################
	# set exon size per gene for FPKM  #
	####################################
	cat("\nsetting rowRanges(dds)...\n")
    dummy<-mcols(dds)
    if(packageVersion("DESeq2")<1.10){
        rowData(dds) <- my.target.list[rownames(dds)] # sorted by the original rownames
	    mcols(rowData(dds))<-dummy
    }else{
        rowRanges(dds) <- my.target.list[rownames(dds)] # sorted by the original rownames
	    rowData(dds)<-dummy
    }

	cat("\nsaving DESeq2 RData\n")
	save(dds,res, file=deseq.RData) #  save 'dds'

	#cat("calculating rld and vsd...\n")
	#print(system.time(rld <- rlog(dds)))# blind=TRUE (default)
										# it takes quite log time!	
	#rld.nb <- rlog(ddsHTSeq, blind=FALSE) # use thie design matrix 
	#vsd <- varianceStabilizingTransformation(ddsHTSeq)
	#save(dds,res,rld, file=deseq.RData) #  save 'dds'
}# end of loading RData

#if(packageVersion("DESeq2")>=1.16){} # CHANGES IN VERSION 1.15.9 adding prototype function lfcShrink().
                                        #‘>=1.16’, the default (betaPrior)is set to ‘FALSE’
if(priorInfo(res)$type=="none"){
    cat("Applying lfcShrink()...\n")
    print(system.time(res <- lfcShrink(dds,contrast=c(my.contrast,levels(colData(dds)[[my.contrast]])[2],levels(colData(dds)[[my.contrast]])[1]),res=res,parallel=TRUE)))
}
my.res <- res
dt.res<-data.table(data.frame(res)); dt.res[[my.filter]]<-rownames(res); 
rn=rownames(colData(dds)) #as.character(samples[["SampleName"]]) # samples$SampleName (or samples[,c("SampleName")])
rn.control=rownames(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1==0,]) 
			#as.character(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1==0,c("SampleName")]) #as.character(samples[as.numeric(samples[[my.contrast]])-1==0,c("SampleName")]) # control only (AGA or Male), which is background
rn.case=rownames(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1!=0,])
			#as.character(samples[as.numeric(samples[[my.contrast]])-1!=0,c("SampleName")]) # case only (SGA or female)
cat("setting FPKM...\n")
cat("fpkm(dds)...\n")
ddsFpkm <-fpkm(dds) # isa 'matrix'
cat("fpm(dds)...\n")
ddsFpm <-fpm(dds) # isa 'matrix'

######################
## Data Exploration ##
######################
if(FALSE){
	if(!exists("rld")){
		cat("rlog(dds)...\n")
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
	cat("making heatmap of all sampleDists...\n")
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".sample.clustering.heatmap"))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=9,units="in",res=300)
    }
	print(heatmap.2(mat, trace="none", col = hmcol, key.xlab="Distance", key.title="", main="Sample Distance"))
	dev.off()

	if(!grepl("^Retosiban",myProject)){
		cat("making sampleDists of case...\n")
		sampleDists<- dist(t(rlogMat[,rn.case])) # isa 'dist'
		mat <- as.matrix(sampleDists) # NxN distance matrix (row: samples, col: samples)
		cat("making heatmap of case sampleDists...\n")
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".sample.clustering.heatmap.case"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
        }else{
	        jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=9,units="in",res=300)
        }
		print(heatmap.2(mat, trace="none", col = hmcol, key.xlab="Distance", key.title="", main="Sample Distance of case"))
		dev.off()

		cat("making sampleDists of control...\n")
		sampleDists<- dist(t(rlogMat[,rn.control])) # isa 'dist'
		mat <- as.matrix(sampleDists) # NxN distance matrix (row: samples, col: samples)
		cat("making heatmap of control sampleDists...\n")
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".sample.clustering.heatmap.control"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
        }else{
	        jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=9,units="in",res=300)
        }
		print(heatmap.2(mat, trace="none", col = hmcol, key.xlab="Distance", key.title="", main="Sample Distance of control"))
		dev.off()

		# summary of boxplot of rlog diff by pair
		if(any(grepl("Pair", deseq.design))){
			diff.rlogMat<-rlogMat[,rn.case]-rlogMat[,rn.control]
			#summary(rlogMat[,rn.case]-rlogMat[,rn.control])
			my.filename <- file.path(cluster.dir, paste0(sampleType, ".boxplot.rlog.case-control"))
            if(capabilities()[["tiff"]]){
			    tiff(filename=paste0(my.filename,".tiff"),width=16,height=9,units="in",res=300, compression = 'lzw')
            }else{
			    jpeg(filename=paste0(my.filename,".jpeg"),width=16,height=9,units="in",res=300)
            }
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
    #select<-!rowData(dds)[["allZero"]] # genes having counts from at least one sample (not allZero)
	cat("making pca based on rlogMat...\n")
	pca<-prcomp(t(rlogMat[select,]), scale.=TRUE) # input (row:samples, col:genes)
	percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

	foo<-merge(data.frame(pca$x[,1:3], `SampleName`=rownames(pca$x)), samples[,c("SampleName","Library","Condition","Pair")])
	p<- ggplot(foo, aes(PC1,PC2,col=Condition)) + 
		geom_text(aes(label=SampleName), size=8) + 
		xlab(paste0("PC1 (",round(percentVar[1]*100,2)," % variance)")) + 
		ylab(paste0("PC2 (",round(percentVar[2]*100,2)," % variance)")) + 
		geom_hline(aes(yintercept=0), size=.2) + geom_vline(aes(xintercept=0), size=.2) +
		theme_Publication()
	cat("making pca plot...\n")
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".sample.pca"))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=9,units="in",res=300)
    }
	print(p)
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

	if(is.BM){
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
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=12,units="in",res=300)
    }
	print(heatmap.2(as.matrix(my.table), scale="row", trace="none", dendrogram="column", col = hmcol, keysize=1, key.title="", margin=c(5,7), offsetRow=0.1, main="Top Variable Genes"))
	dev.off()

	if(!grepl("^Retosiban",myProject)){
		cat("making heatmap of gene clustering from case...\n")
		select<- head( order(rowVars(rlogMat[,rn.case]), decreasing=TRUE ), n=20 ) # top 20 variable genes (requires 'genefilter' library)
		my.table<-rlogMat[select,rn.case]
		my.table<-as.data.frame(my.table[order(rownames(my.table)),]) # order by rownames (ensembl_gene_id or mirbase_id)

		if(is.BM){
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
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=12,units="in",res=300)
        }
		print(heatmap.2(as.matrix(my.table), scale="row", trace="none", dendrogram="column", col = hmcol, keysize=1, key.title="", margin=c(5,7), offsetRow=0.1, main="Top Variable Genes of case"))
		dev.off()

		cat("making heatmap of gene clustering from control...\n")
		#select<- head( order(rowVars(rlogMat[,rn.control]), decreasing=TRUE ), n=round(nrow(rlogMat[,rn.control])*0.1*0.01) ) # top 0.1% variable genes (requires 'genefilter' library)
		select<- head( order(rowVars(rlogMat[,rn.control]), decreasing=TRUE ), n=20 ) # top 20 variable genes (requires 'genefilter' library)
		my.table<-rlogMat[select,rn.control]
		my.table<-as.data.frame(my.table[order(rownames(my.table)),]) # order by rownames (ensembl_gene_id or mirbase_id)

		if(is.BM){
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
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=12,units="in",res=300)
        }
		print(heatmap.2(as.matrix(my.table), scale="row", trace="none", dendrogram="column", col = hmcol, keysize=1, key.title="", margin=c(5,7), offsetRow=0.1, main="Top Variable Genes of control"))
		dev.off()
	}

	###############################################################
	## Box plot of FPKM of genes having > 1 CPM for all samples  ##
	## All samples should have more than 1 CPM/FPM               ##
	###############################################################
	cat("making ddsFpkm from dds(ExonicGene)...\n")
	keep = rowSums(ddsFpm >1) >= ncol(ddsFpm)
	#cat("making heatmap of boxplot of FPKM in log2 scale...\n")
	#boxplot(log2(1+ddsFpkm[keep,]), xlab="Samples", ylab="log2(1+FPKM)", main=paste0("FPKM boxplots of ",nrow(ddsFpkm[keep,]), " genes having > 1 CPM across all samples"))

	cat("making heatmap of boxplot of FPKM...\n")
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".fpkm.boxplot"))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=9,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=9,units="in",res=300)
    }
	print(boxplot(ddsFpkm[keep,], xlab="Samples", ylab="FPKM", main=paste0("FPKM boxplots of ",nrow(ddsFpkm[keep,]), " genes having > 1 CPM across samples")))
	dev.off()
	if(!grepl("^Retosiban",myProject)){
		cat("making heatmap of boxplot of case FPKM...\n")
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".fpkm.boxplot.case"))
        if(capabilities()[["tiff"]]){
	    	tiff(filename=paste0(my.filename,".tiff"),width=10,height=9,units="in",res=300, compression = 'lzw')
        }else{
	    	jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=9,units="in",res=300)
        }
		print(boxplot(ddsFpkm[keep,rn.case], xlab="case Samples", ylab="FPKM", main=paste0("FPKM boxplots of ",nrow(ddsFpkm[keep,rn.case]), " genes having > 1 CPM across case samples")))
		dev.off()

		cat("making heatmap of boxplot of control FPKM...\n")
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".fpkm.boxplot.control"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=10,height=9,units="in",res=300, compression = 'lzw')
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=9,units="in",res=300)
        }
		print(boxplot(ddsFpkm[keep,rn.control], xlab="control Samples", ylab="FPKM", main=paste0("FPKM boxplots of ",nrow(ddsFpkm[keep,rn.control]), " genes having > 1 CPM across control samples")))
		dev.off()
	}

	###################
	## Top FPKM genes #
	###################
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
	if(!grepl("^Retosiban",myProject)){
		# Top 100 FPKM Genes
		cat("making heatmap of Top FPKM based on IQR mean...\n")
		myTopExp <- getTopExp(ddsFpkm[keep,], iqr=FALSE, 100) # isa 'list'
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.iqr.fpkm.gene.heatmap.all"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=12,units="in",res=300)
        }
		print(heatmap.2(head(as.matrix(log(myTopExp$fpkm+0.001)),n=round(nrow(ddsFpkm)*0.1*0.01)), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="FPKM", main="Top FPKM Genes"))
		dev.off()

		myTopExp$fpkm$mean<-rowMeans(myTopExp$fpkm,na.rm=T) # mean over all the samples
		myTopExp$fpkm$meanIQR<-apply(myTopExp$fpkm, 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]])) # mean over trimmed samples
		myTopExp$fpkm$SD<-rowSds(myTopExp$fpkm,na.rm=T) # mean over all the samples
		if(grepl("piRNA",myProject)){ # If piRNA
			my.fpkm<-round(myTopExp$fpkm,1)
		}else{
			my.fpkm<-merge(cbind(round(myTopExp$fpkm,1), hgnc_symbol=rownames(myTopExp$fpkm)), myTopExp$gene, all.x=TRUE)
		}
		write.csv(my.fpkm[order(my.fpkm$meanIQR,decreasing=T),], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".top.100.iqr.fpkm.csv.gz"))), row.names=F)

		###################
		# Top 100 case FPKM
		###################
		cat("making heatmap of Top FPKM based on IQR mean of case...\n")
		myTopExp.case <- getTopExp(ddsFpkm[keep,rn.case], iqr=TRUE, 100) # isa 'list'
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.iqr.fpkm.gene.heatmap.case"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=12,units="in",res=300)
        }
		print(heatmap.2(head(as.matrix(log(myTopExp.case$fpkm+0.001)),n=round(nrow(ddsFpkm)*0.1*0.01)), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="FPKM", main="Top IQR(FPKM) of case"))
		dev.off()

		myTopExp.case$fpkm$mean<-rowMeans(myTopExp.case$fpkm,na.rm=T) # mean over all the samples
		myTopExp.case$fpkm$meanIQR<-apply(myTopExp.case$fpkm, 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]])) # mean over trimmed samples
		myTopExp.case$fpkm$SD<-rowSds(myTopExp$fpkm,na.rm=T) # mean over all the samples
		if(grepl("piRNA",myProject)){ # If piRNA
			my.fpkm<-round(myTopExp.case$fpkm,1)
		}else{
			my.fpkm<-merge(cbind(round(myTopExp.case$fpkm,1), hgnc_symbol=rownames(myTopExp.case$fpkm)), myTopExp.case$gene, all.x=TRUE)
		}
		write.csv(my.fpkm[order(my.fpkm$meanIQR,decreasing=T),], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".top.100.iqr.fpkm.case.csv.gz"))))

		######################
		# Top 100 control FPKM
		######################
		cat("making heatmap of Top FPKM based on IQR mean of control...\n")
		myTopExp.ctl<-getTopExp(ddsFpkm[keep,rn.control], iqr=TRUE, 100) # isa 'list' 
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.iqr.fpkm.gene.heatmap.control"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=12,units="in",res=300)
        }
		print(heatmap.2(head(as.matrix(log(myTopExp.ctl$fpkm+0.001)),n=round(nrow(ddsFpkm)*0.1*0.01)), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="FPKM", main="Top IQR(FPKM) of control"))
		dev.off()

		myTopExp.ctl$fpkm$mean<-rowMeans(myTopExp.ctl$fpkm,na.rm=T) # mean over all the samples
		myTopExp.ctl$fpkm$meanIQR<-apply(myTopExp.ctl$fpkm, 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]])) # mean over trimmed samples
		myTopExp.ctl$fpkm$SD<-rowSds(myTopExp$fpkm,na.rm=T) # mean over all the samples
		if(grepl("piRNA",myProject)){ # If piRNA
			my.fpkm<-round(myTopExp.ctl$fpkm,1)
		}else{
			my.fpkm<-merge(cbind(round(myTopExp.ctl$fpkm,1), hgnc_symbol=rownames(myTopExp.ctl$fpkm)), myTopExp.ctl$gene, all.x=TRUE)
		}
		write.csv(my.fpkm[order(my.fpkm$meanIQR,decreasing=T),], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".top.100.iqr.fpkm.control.csv.gz"))))
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
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=9,height=10,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=10,units="in",res=300)
    }
	plot(log2(1+apply(ddsFpkm[keep,], 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]]))), log2(1+apply(ddsFpkm[keep,], 1, function(i) sd(i[i>=quantile(i)[2] & i<=quantile(i)[4]]))), xlab="log2(1+mean(IQR(FPKM)))", ylab="log2(1+SD(IQR(FPKM)))", main="log2(mean)-log2(SD) plot") # SD-Mean plot
	dev.off()

	if(!grepl("^Retosiban",myProject)){
		cat("making Sdplot based on ddsFpkm of case...\n")
		mean.iqr.fpkm<-apply(ddsFpkm[keep,rn.case], 1, function(i) ifelse(length(i)>=4,mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]]),mean(i))) # isa numeric vector
		sd.iqr.fpkm<-apply(ddsFpkm[keep,rn.case], 1, function(i) ifelse(length(i)>=4,sd(i[i>=quantile(i)[2] & i<=quantile(i)[4]]),sd(i))) # isa numeric vector
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".iqr.fpkm.mean.sd.case"))
        if(capabilities()[["tiff"]]){
	    	tiff(filename=paste0(my.filename,".tiff"),width=9,height=10,units="in",res=300, compression = 'lzw')
        }else{
	    	jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=10,units="in",res=300)
        }
		plot(mean.iqr.fpkm, sd.iqr.fpkm, xlab="mean(IQR(FPKM))", ylab="SD(IQR(FPKM))", main="Mean-SD plot of case") # SD-Mean plot
		text(cbind(mean.iqr.fpkm, sd.iqr.fpkm)[myTopExp.case$select[1:5],], labels=rownames(myTopExp.case$fpkm)[1:5], pos=2) # label to the left (pos=2)
		dev.off()

		cat("making Sdplot based on ddsFpkm of control...\n")
		mean.iqr.fpkm<-apply(ddsFpkm[keep,rn.control], 1, function(i) ifelse(length(i)>=4,mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]]),mean(i))) # isa numeric vector
		sd.iqr.fpkm<-apply(ddsFpkm[keep,rn.control], 1, function(i) ifelse(length(i)>=4,sd(i[i>=quantile(i)[2] & i<=quantile(i)[4]]),sd(i))) # isa numeric vector
		my.filename <- file.path(cluster.dir, paste0(sampleType, ".iqr.fpkm.mean.sd.control"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=9,height=10,units="in",res=300, compression = 'lzw')
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=10,units="in",res=300)
        }
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

# Boy.Girl chrX only (FG & JD)
if(FALSE){ 
	##############
	## Outliers ##
	##############
	df.cooks<-assays(dds)$cooks
	top.cooks.genes<-names(head(sort(rowSums(df.cooks), decreasing=T), n=10)) # top cooks genes
	ddsFpkm[top.cooks.genes,] # top cooks genes
	top.cooks.samples<-sort(colSums(ddsFpkm[top.cooks.genes,]), decreasing=T) # top cooks samples
	round(ddsFpkm[top.cooks.genes, head(names(top.cooks.samples),n=10)], 2)
	round(sort(ddsFpkm[top.cooks.genes[1],], decreasing=T))

	dt.ensg[ensembl_gene_id %in% top.cooks.genes]

	#IGFBP1: 84C, 88, 80C, 06C
	#SNORD3C: 02C, 02C
	dummy<-rbind(
		merge(cbind(reshape2::melt(ddsFpm[,rn.case]),`Sex`=c("F")), unique(dt.ensg[,.(ensembl_gene_id,chromosome_name)]), by.x="Var1", by.y="ensembl_gene_id"), # case only (SGA or feale)
		merge(cbind(reshape2::melt(ddsFpm[,rn.control]),`Sex`=c("M")), unique(dt.ensg[,.(ensembl_gene_id,chromosome_name)]), by.x="Var1", by.y="ensembl_gene_id") # control only (AGA or male)
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
	foo<-merge(cbind(as.data.frame(res), `ensembl_gene_id`=rownames(res)), unique(dt.ensg[,.(ensembl_gene_id,chromosome_name)]))
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
	dt.res[ensembl_gene_id==my.entry]
	ddsFpm[my.entry,rn.case]
	ddsFpm[my.entry,rn.control]
	wilcox.test(ddsFpm[my.entry,rn.case], ddsFpm[my.entry,rn.control])$p.value
	# up/down 1Mb
	my.upflank=10^6 # 1Mb up from target TSS
	my.downflank=10^6 # 1Mb down from target TSS 
	my.gr.ext<-promoters(gr.ensg[my.entry], upstream=my.upflank, downstream=end(gr.ensg[my.entry])-start(gr.ensg[my.entry])+1+my.downflank) # includes +- flanking from TSS and TES
	my.overlap=as.matrix( findOverlaps(gr.ensg, my.gr.ext, ignore.strand=TRUE, minoverlap=10L) ) #By default, any overlap is accepted
	my.target.ensg<-names(gr.ensg[my.overlap[,"queryHits"]])
	dt.ensg[ensembl_gene_id %in% my.target.ensg]

	ggplot(dt.csmd.fpkm, aes(tissue, FPKM, fill=Gender)) + geom_bar(stat="identity", position="dodge") + scale_fill_manual(name="Sex",values=my.col[["Gender"]]) + theme_Publication()


	##
	## chrX
	##
	dt.res<-data.table(mirbase_id=rownames(res),data.frame(res))
	dt.chr<-as.data.table(getBM(attributes =c("mirbase_id","ensembl_gene_id","chromosome_name"), filters = "mirbase_id", values = dt.res[,unique(mirbase_id)], mart = myMart)) # is a data.frame
	dt.res.chrX=merge(dt.res, dt.chr, by="mirbase_id", all.x=TRUE, allow.cartesian=TRUE)[chromosome_name=="X"]
	dt.res.chrX<-dt.res.chrX[,`new.padj` := p.adjust(pvalue, method="BH")]

} #myProject=="Boy.Girl"


# FG plasma samples 2016
if(myProject=="Plasma.2016"){
	library(data.table)

	# get chr and gene name
	cat("getting gene annotations...\n")
	# get 12wk and 36wk plasma samples
	rn.12wk=rownames(colData(dds)[colData(dds)$GA=="12",])
	rn.36wk=rownames(colData(dds)[colData(dds)$GA=="36",])

	# set up data.table for FPKM by group
	#dt.meanFpkm<-data.table(`ensembl_gene_id`=rownames(ddsFpkm), `Plasma.12wk`=rowMeans(ddsFpkm[,rn.12wk]), `Plasma.36wk`=rowMeans(ddsFpkm[,rn.36wk]), `Plasma.all`=rowMeans(ddsFpkm[,rn.case]), `Placenta`=rowMeans(ddsFpkm[,rn.control]))
	dt.meanFpkm<-data.table(`ensembl_gene_id`=rownames(ddsFpkm), `Plasma.2008.12wk`=ddsFpkm[,"85PL"], 
							`Plasma.2008.36wk`=ddsFpkm[,"87PL"], `Plasma.2012.12wk`=ddsFpkm[,"89PL"], 
							`Plasma.2012.36wk`=ddsFpkm[,"91PL"], `Plasma.12wk`=rowMeans(ddsFpkm[,rn.12wk]), 
							`Plasma.36wk`=rowMeans(ddsFpkm[,rn.36wk]), `Plasma.all`=rowMeans(ddsFpkm[,rn.case]), 
							`Placenta`=rowMeans(ddsFpkm[,rn.control]))

	#dt.meanFpkm.NR<-data.table(`ensembl_gene_id`=rownames(ddsFpkm.NR), `Plasma.2008.12wk`=ddsFpkm.NR[,"85PL"], `Plasma.2008.36wk`=ddsFpkm.NR[,"87PL"], `Plasma.2012.12wk`=ddsFpkm.NR[,"89PL"], `Plasma.2012.36wk`=ddsFpkm.NR[,"91PL"], `Plasma.12wk`=rowMeans(ddsFpkm.NR[,rn.12wk]), `Plasma.36wk`=rowMeans(ddsFpkm.NR[,rn.36wk]), `Plasma.all`=rowMeans(ddsFpkm.NR[,rn.case]), `Placenta`=rowMeans(ddsFpkm.NR[,rn.control]))
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
	my.dt.meanFpkm<-merge(dt.meanFpkm[ensembl_gene_id %in% dt.ensg[chromosome_name=="Y","ensembl_gene_id"] & Plasma.all + Placenta >0], dt.ensg, by="ensembl_gene_id")
	my.dt.meanFpkm.NR<-merge(dt.meanFpkm.NR[ensembl_gene_id %in% dt.ensg[chromosome_name=="Y","ensembl_gene_id"] & Plasma.all + Placenta >0], dt.ensg, by="ensembl_gene_id")
	my.dt.meanCount<-merge(dt.meanCount[ensembl_gene_id %in% dt.ensg[chromosome_name=="Y","ensembl_gene_id"] & Plasma.all + Placenta >0], dt.ensg, by="ensembl_gene_id")

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
    if(capabilities()[["tiff"]]){
	    my.filename <- file.path(cluster.dir, paste0(sampleType, ".venn.plasma.pt.fpkm.more.than.",minFPKM,".tiff"))
    }else{
	    my.filename <- file.path(cluster.dir, paste0(sampleType, ".venn.plasma.pt.fpkm.more.than.",minFPKM,".jpeg"))
    }
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
	myFpkm<-ddsFpkm[dt.ensg[chromosome_name=="Y","ensembl_gene_id"],] # chrY only

	# heatmap of any chrY expressed
	#keep = rowSums(myFpkm) > 0  
	#myTopExp <- getTopExp(myFpkm[keep,], iqr=FALSE, nrow(myFpkm[keep,])) # isa 'list'

	# heatmap by samples
	# plasma chrY expressed sorted by rowMeans (N=47)
	myFpkm<-ddsFpkm[my.dt.meanFpkm[Plasma.all>0,ensembl_gene_id],]
	myTopExp <- getTopExp(myFpkm, iqr=FALSE, nrow(myFpkm)) # isa 'list'

	my.filename <- file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.plasma.heatmap"))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=12,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=12,units="in",res=300)
    }
	heatmap.2(as.matrix(log(myTopExp$fpkm+0.01)), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="log(FPKM)", main="Plasma ChrY Genes of FPKM>0", margins=c(4,8))
	dev.off()

	my.filename <- file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.plasma.heatmap.tr"))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=13,height=10,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=13,height=10,units="in",res=300)
    }
	heatmap.2(as.matrix(log(t(myTopExp$fpkm)+0.01)), Colv=FALSE, scale="none", trace="none", dendrogram="row", col = hmcol, keysize=1, key.xlab="log(FPKM)", main="Plasma ChrY Genes of FPKM>0", margins=c(8,4), srtCol=45)
	dev.off()

	# heatmap of plasma chrY expressed clustered by genes and samples (N=47)
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.plasma.row.dendro.heatmap"))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=12,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=12,units="in",res=300)
    }
	heatmap.2(as.matrix(log(myTopExp$fpkm+0.01)), Rowv=TRUE, scale="none", trace="none", dendrogram="both", col = hmcol, keysize=1, key.xlab="log(FPKM)", main="Plasma ChrY Genes of FPKM>0", margins=c(4,8))
	dev.off()

	# heatmap sorted by plasma 
	foo<-myTopExp$fpkm[names(sort(rowMeans(myTopExp$fpkm[,rn.case]),decreasing=T)),] # re-order by FPKM of Plasma
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.plasma.sorted.heatmap"))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=12,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=12,units="in",res=300)
    }
	heatmap.2(as.matrix(log(foo+0.01)), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="log(FPKM)", main="Plasma ChrY Genes of FPKM>0", margins=c(4,8))
	dev.off()

	# heatmap by group
	for(i in c(0, 0.1, 0.5, 1)){
		foo<-as.matrix(my.dt.meanFpkm[Plasma.all>i, .(Plasma.2008.12wk, Plasma.2008.36wk, Plasma.2012.12wk, Plasma.2012.36wk, Plasma.12wk, Plasma.36wk, Plasma.all, Placenta)])
		rownames(foo)<-my.dt.meanFpkm[Plasma.all>i, ensembl_gene_id]
		myTopExp <- getTopExp(foo, iqr=FALSE, nrow(foo)) # isa 'list'

		my.filename <- file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.more.than.",i,".plasma.heatmap"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=10,height=12,units="in",res=300, compression = 'lzw')
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=12,units="in",res=300)
        }
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
	dummy<-dt.ensg[ensembl_gene_id %in% dt.meanFpkm[ensembl_gene_id %in% my.ensg & Plasma.all>0,ensembl_gene_id],.(ensembl_gene_id,gene_biotype)]
	# protein_coding gene only
	my.ensg<-dummy[dummy$gene_biotype=="protein_coding","ensembl_gene_id"]
	#heatmap by group
	for(i in c(0, 0.1, 0.5, 1)){
		foo<-as.matrix(dt.meanFpkm[ensembl_gene_id %in% my.ensg & Plasma.all>i, .(Plasma.2008.12wk, Plasma.2008.36wk, Plasma.2012.12wk, Plasma.2012.36wk, Plasma.12wk, Plasma.36wk, Plasma.all, Placenta)])
		rownames(foo)<-dt.meanFpkm[ensembl_gene_id %in% my.ensg & Plasma.all>i, ensembl_gene_id]
		myTopExp <- getTopExp(foo, iqr=FALSE, nrow(foo)) # isa 'list'

		my.filename <- file.path(cluster.dir, paste0(sampleType, ".placenta.specific.genes.fpkm.more.than.",i,".by.group.heatmap"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=10,height=12,units="in",res=300, compression = 'lzw')
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=12,units="in",res=300)
        }
		heatmap.2(as.matrix(log(myTopExp$fpkm+0.01)), Colv=FALSE, Rowv=FALSE, scale="none", trace="none", dendrogram="none", col = hmcol, keysize=1, key.xlab="log(FPKM)", main=paste("Placenta Specific Genes of Avg Plasma FPKM>",i), margins=c(9,9),  srtCol=45, cellnote=round(myTopExp$fpkm,2), notecol="black")
		dev.off()
	}

	my.dt.meanFpkm<-merge(dt.meanFpkm[ensembl_gene_id %in% my.ensg & Plasma.all>0],dt.ensg, by="ensembl_gene_id")
	my.dt.meanCount<-merge(dt.meanCount[ensembl_gene_id %in% my.ensg & Plasma.all>0],dt.ensg, by="ensembl_gene_id")

	# FPKM per group 
	write.csv(my.dt.meanFpkm[Plasma.all>0][order(-Plasma.all)], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".placenta.specific.genes.fpkm.more.than.0.per.group.csv.gz"))))
	write.csv(my.dt.meanCount[Plasma.all>0][order(-Plasma.all)], file=gzfile(file.path(cluster.dir, paste0(sampleType, ".placenta.specific.genes.read.count.more.than.0.per.group.csv.gz"))))

	#######################################
	## Top 5 highly expressed chrY genes ## 
	#######################################
	my.file.name<- paste0(cluster.dir,'/top.expressed.chrY.',sampleType)
	pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title=paste0(myProject,":",sampleType)) # A4 size
	# 1. top5 FPKM (protein-coding only)
	dt.foo<-my.dt.meanFpkm[gene_biotype=="protein_coding"][order(-Plasma.all)][1:5]
	dt.foo<-melt(dt.foo, id=c("ensembl_gene_id","hgnc_symbol"), measure=c("Plasma.12wk","Plasma.36wk"),  variable.name="GA", value.name="FPKM"); dt.foo$Gender="M"
	lapply(dt.foo[,unique(hgnc_symbol)][1], function(i) ggplot(dt.foo[hgnc_symbol==i], aes(GA, FPKM, group=Gender)) + geom_point(size=5) + geom_line(linetype=2) + ggtitle(i) + scale_x_discrete(limits=c("Plasma.12wk","Plasma.36wk"), labels=c("12WK", "36WK"))  + theme_Publication() )

	# 2. top5 FPKM
	dt.bar<-my.dt.meanFpkm[order(-Plasma.all)][1:5]
	dt.foo<-melt(dt.bar, id=c("ensembl_gene_id","hgnc_symbol"), measure=c("Plasma.12wk","Plasma.36wk"),  variable.name="GA", value.name="FPKM"); dt.foo$Gender="M"
	lapply(dt.foo[,unique(hgnc_symbol)][1], function(i) ggplot(dt.foo[hgnc_symbol==i], aes(GA, FPKM, group=Gender)) + geom_point(size=5) + geom_line(linetype=2) + ggtitle(i) + scale_x_discrete(limits=c("Plasma.12wk","Plasma.36wk"), labels=c("12WK", "36WK"))  + theme_Publication() )

	dev.off()	
} # end of myProject=="Plasma.2016"

# mapping and annotation based on GRCh38
if(myProject=="Plasma.2017"){
	library(data.table)

	#ddsFpkm <-fpkm(dds) # isa 'matrix'
	#ddsFpm <-fpm(dds) # isa 'matrix'

	cat("getting gene annotations...\n")
	#dt.ensg <- data.table(getBM(attributes = my.fields, filters = my.filter, values = rownames(dds), mart = myMart))
	#dt.ensg from config/Annotation.R

	# All plasma samples are breech (of JD)
	dt.samples<-data.table(as.data.frame(colData(dds)), SampleID=rownames(colData(dds)))

	# FPKM, FPM & Counts
	dt.Fpkm<-data.table(`ensembl_gene_id`=rownames(ddsFpkm), ddsFpkm)
	dt.Fpm<-data.table(`ensembl_gene_id`=rownames(ddsFpm), ddsFpm)
	dt.Count<-data.table(`ensembl_gene_id`=rownames(counts(dds)), counts(dds))

	##############################################################################
	## Top expressed chrY genes and some placenta-specific genes including PSG7 ##
	##############################################################################
	# placenta-specific genes
	my.ensg=c(`PSG1`="ENSG00000231924", `PSG2`="ENSG00000242221", `PSG3`="ENSG00000221826", `PSG4`="ENSG00000243137",
			`PSG5`="ENSG00000204941", `PSG6`="ENSG00000170848", `PSG7`="ENSG00000221878", `PSG8`="ENSG00000124467",
			`PSG9`="ENSG00000183668", `PSG10P`="ENSG00000248257", `PSG11`="ENSG00000243130", 
			`CGA`="ENSG00000135346", `CGB1`="ENSG00000267631", `CGB2`="ENSG00000104818", `CGB3`="ENSG00000104827",
			`CGB5`="ENSG00000189052", `CGB7`="ENSG00000196337", `CGB8`="ENSG00000213030", `CSH1`="ENSG00000136488",
			`CSH2`="ENSG00000213218", `PLAC1`="ENSG00000170965", `PLAC4`="ENSG00000280109", `PLAC8`="ENSG00000145287",
			`PAPPA`="ENSG00000182752", `PAPPA2`="ENSG00000116183",`SMS`="ENSG00000102172"
			)
	#my.ensg=c(`PSG7`="ENSG00000221878")
	# chrY and target genes only
	dt.target<-dt.ensg[chromosome_name=="Y" | hgnc_symbol %in% names(my.ensg),.(ensembl_gene_id,hgnc_symbol,chromosome_name)]

	dt.melt.Fpkm<-melt(dt.Fpkm[ensembl_gene_id %in% dt.target[,ensembl_gene_id]], id="ensembl_gene_id", variable.name="SampleID", value.name="FPKM")
	dt.melt.Fpm<-melt(dt.Fpm[ensembl_gene_id %in% dt.target[,ensembl_gene_id]], id="ensembl_gene_id", variable.name="SampleID", value.name="FPM")
	dt.melt.Count<-melt(dt.Count[ensembl_gene_id %in% dt.target[,ensembl_gene_id]], id="ensembl_gene_id", variable.name="SampleID", value.name="Count")
	#       ensembl_gene_id SampleID              FPKM
	#    1: ENSG00000012817      413                 0
	#	 2: ENSG00000067048      413                 0

	# FPKM of each samples
	dt.foo.fpkm<-merge(dt.melt.Fpkm, dt.samples) 
	dt.ga.sex.fpkm<-merge(dt.foo.fpkm[,.(`FPKM`=mean(FPKM)),"ensembl_gene_id,GA,SeqType,Sex"], dt.ensg) # mean FPKM by GA and Sex for a gene
	dt.sample.fpkm<-merge(dt.foo.fpkm, dt.ensg, by="ensembl_gene_id") # FPKM of each samples with gene names

	# FPM of each samples
	dt.foo.fpm<-merge(dt.melt.Fpm, dt.samples) 
	dt.ga.sex.fpm<-merge(dt.foo.fpm[,.(`FPM`=mean(FPM)),"ensembl_gene_id,GA,SeqType,Sex"], dt.ensg) # mean FPM by GA and Sex for a gene
	dt.sample.fpm<-merge(dt.foo.fpm, dt.ensg, by="ensembl_gene_id") # FPM of each samples with gene names

	# Count of each samples
	dt.foo.count<-merge(dt.melt.Count, dt.samples) 
	dt.ga.sex.count<-merge(dt.foo.count[,.(`Count`=mean(Count)),"ensembl_gene_id,GA,SeqType,Sex"], dt.ensg) # mean Count by "GA, SeqType and Sex"
	dt.sample.count<-merge(dt.foo.count, dt.ensg, by="ensembl_gene_id") # Count of each samples with gene names

	##
	## SMS
	##ggplot(merge(dt.foo.fpkm[BarCode!="NoIndex",.(ensembl_gene_id, Sex,GA,FPKM)], dt.ensg[hgnc_symbol=="SMS",.(ensembl_gene_id)], by="ensembl_gene_id"), aes(GA, FPKM)) + geom_boxplot(aes(fill=Sex),alpha=.5)
	##ggplot(merge(dt.foo.fpm[BarCode!="NoIndex",.(ensembl_gene_id, Sex,GA,FPM)], dt.ensg[hgnc_symbol=="SMS",.(ensembl_gene_id)], by="ensembl_gene_id"), aes(GA, FPM)) + geom_boxplot(aes(fill=Sex),alpha=.5)
	##ggplot(merge(dt.foo.count[BarCode!="NoIndex",.(ensembl_gene_id, Sex,GA,Count)], dt.ensg[hgnc_symbol=="SMS",.(ensembl_gene_id)], by="ensembl_gene_id"), aes(GA, Count)) + geom_boxplot(aes(fill=Sex),alpha=.5)

	dt.plasma.fpkm<-merge(
			# Highly expressed chrY genes by GA & Sex 
			data.table::dcast(dt.ga.sex.fpkm[chromosome_name=="Y" | hgnc_symbol %in% c("PSG7","SMS")], ensembl_gene_id+hgnc_symbol+gene_biotype~GA+SeqType+Sex, value.var="FPKM", sep="_"), # long (row-based) to wide (column-based)
			# Highly expressed chrY genes across all samples
			#dt.sample.fpkm[chromosome_name=="Y" | hgnc_symbol=="PSG7", .(all=mean(FPKM)),"ensembl_gene_id"][all>0][order(-all)],
			dt.sample.fpkm[chromosome_name=="Y" | hgnc_symbol %in% c("PSG7","SMS"), .(all=mean(FPKM)),"ensembl_gene_id"],
			by="ensembl_gene_id")[order(-all)]

	write.csv(dt.plasma.fpkm, file=file.path(cluster.dir, paste0(myProject, ".top.chrY.PSG7.fpkm.per.group.csv")))
	write.csv(dt.Fpkm[ensembl_gene_id %in% dt.plasma.fpkm$ensembl_gene_id], file=gzfile(file.path(cluster.dir, paste0(myProject, ".top.chrY.PSG7.fpkm.per.sample.csv.gz"))))

	write.csv(merge(dt.Count, dt.ensg[,.(ensembl_gene_id,chromosome_name)])[chromosome_name=="Y"], file=gzfile(file.path(cluster.dir, paste0(myProject, ".all.chrY.PSG7.count.per.sample.csv.gz"))))
	write.csv(merge(dt.Fpm, dt.ensg[,.(ensembl_gene_id,chromosome_name)])[chromosome_name=="Y"], file=gzfile(file.path(cluster.dir, paste0(myProject, ".all.chrY.PSG7.fpm.per.sample.csv.gz"))))

	dt.plasma.count<-merge(
			# Highly expressed chrY genes by GA & Sex 
			data.table::dcast(dt.ga.sex.count[chromosome_name=="Y" | hgnc_symbol=="PSG7"], ensembl_gene_id+hgnc_symbol+gene_biotype~GA+SeqType+Sex, value.var="Count", sep="_"), # long (row-based) to wide (column-based)
			# Highly expressed chrY genes across all samples
			#dt.sample.count[chromosome_name=="Y" | hgnc_symbol=="PSG7", .(all=mean(Count)),"ensembl_gene_id"][all>0][order(-all)],
			dt.sample.count[chromosome_name=="Y" | hgnc_symbol=="PSG7", .(all=mean(Count)),"ensembl_gene_id"][order(-all)],
			by="ensembl_gene_id")[order(-all)]

	write.csv(dt.plasma.count, file=file.path(cluster.dir, paste0(myProject, ".top.chrY.PSG7.per.group.count.csv")))
	write.csv(dt.Count[ensembl_gene_id %in% dt.plasma.count$ensembl_gene_id], file=gzfile(file.path(cluster.dir, paste0(myProject, ".top.chrY.PSG7.count.per.sample.csv.gz"))))

	# merge FPKM & Count for this gene
	#merge(dt.sample.count[hgnc_symbol=="TMSB4Y",.(SampleID,Count)], dt.sample.fpkm[hgnc_symbol=="TMSB4Y",.(SampleID,FPKM)])[order(FPKM)]

	#############
	## Figures ##
	#############
	my.file.name<- paste0(cluster.dir,'/top.expressed.chrY.',sampleType)
	pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title=paste0(myProject,":",sampleType)) # A4 size
	# genes of interests
	my.genes=c("PSG7", "RNA5-8SP6", "TMSB4Y")
	for(my.gene in my.genes){
		# Histogram of log2(FPKM)
		if(my.gene=="PSG7"){
			p0<-ggplot(dt.sample.fpkm[hgnc_symbol==my.gene,.(SampleID,Sex,GA,FPKM)], aes(x=log2(FPKM))) +  
				geom_histogram(colour="black", fill="white") + 
				ggtitle(my.gene) +
				theme_Publication()
			print(p0)

			p0.1<-ggplot(dt.sample.fpkm[hgnc_symbol==my.gene,.(SampleID,Sex,GA,FPKM)], aes(x=FPKM)) +  
				geom_histogram(colour="black", fill="white") + 
				ggtitle(my.gene) +
				theme_Publication()
			print(p0.1)

			p0.2<-ggplot(dt.sample.count[hgnc_symbol==my.gene,.(SampleID,Sex,GA,Count)], aes(x=Count)) +  
				geom_histogram(colour="black", fill="white") + 
				ggtitle(my.gene) +
				theme_Publication()
			print(p0.2)

			# boxplot of FPKM vs GA
			p2<-ggplot(dt.sample.fpkm[hgnc_symbol==my.gene], aes(GA, FPKM)) + 
				geom_boxplot(outlier.shape=1, outlier.size=5, size=1, width=.5) +
				ggtitle(my.gene) + 
				scale_colour_manual(values=my.col[["Sex"]]) +
				theme_Publication()
			print(p2)
		}else{
			# mean FPKM vs GA
			p1<-ggplot(dt.ga.sex.fpkm[hgnc_symbol==my.gene], aes(GA, FPKM, group=Sex)) + 
				geom_point(aes(col=Sex), size=6) + 
				geom_line(lty=2) + 
				ggtitle(my.gene) + 
				scale_colour_manual(values=my.col[["Sex"]]) +
				theme_Publication()
			print(p1)

			# mean FPKM vs GA
			p1.1<-ggplot(dt.ga.sex.fpkm[hgnc_symbol==my.gene], aes(GA, FPKM, group=Sex)) + 
				geom_point(aes(shape=Sex), size=8) + 
				geom_line() + 
				ggtitle(my.gene) + 
				scale_shape_discrete(solid=F) +
				theme_Publication()
			print(p1.1)

			# boxplot of FPKM vs GA
			p2<-ggplot(dt.sample.fpkm[hgnc_symbol==my.gene], aes(GA, FPKM)) + 
				geom_boxplot(aes(col=Sex), outlier.shape=1, outlier.size=5, size=1, width=.5) +
				ggtitle(my.gene) + 
				scale_colour_manual(values=my.col[["Sex"]]) +
				theme_Publication()
			print(p2)
		}
	}
	dev.off()

	#############################
	## Placenta breech samples ##
	## based on HT-Seq         ##
	## N=24 (n(F)=12, n(M)=12) ##
	## Plasma samples are BR   ##
	#############################
	load("~/results/RNA-Seq/Boy.Girl.FG.JD.GRCh38/DESeq2/BR/deseq.BR.RData")
	# samples
	dt.pt.samples<-data.table(as.data.frame(colData(dds)), SampleID=rownames(colData(dds)))
	pt.ddsFpkm <-fpkm(dds) # isa 'matrix'
	pt.ddsFpm <-fpm(dds) # isa 'matrix'

	# FPKM, FPM & Counts
	dt.pt.Fpkm<-data.table(`ensembl_gene_id`=rownames(pt.ddsFpkm), pt.ddsFpkm)
	dt.pt.Fpm<-data.table(`ensembl_gene_id`=rownames(pt.ddsFpm), pt.ddsFpm)
	dt.pt.Count<-data.table(`ensembl_gene_id`=rownames(counts(dds)), counts(dds))
	# target genes only (see dt.target above)
	dt.melt.pt.Fpkm<-melt(dt.pt.Fpkm[ensembl_gene_id %in% dt.target[,ensembl_gene_id]], id="ensembl_gene_id", variable.name="SampleID", value.name="FPKM")
	dt.melt.pt.Fpm<-melt(dt.pt.Fpm[ensembl_gene_id %in% dt.target[,ensembl_gene_id]], id="ensembl_gene_id", variable.name="SampleID", value.name="FPM")
	dt.melt.pt.Count<-melt(dt.pt.Count[ensembl_gene_id %in% dt.target[,ensembl_gene_id]], id="ensembl_gene_id", variable.name="SampleID", value.name="Count")
	
	#    1: ENSG00000012817      413                 0
	#	 2: ENSG00000067048      413                 0
	dt.foo.pt.fpkm<-merge(dt.melt.pt.Fpkm, dt.pt.samples) # FPKM of each samples
	dt.sample.pt.fpkm<-merge(dt.foo.pt.fpkm, dt.ensg, by="ensembl_gene_id") # FPKM of each samples with gene names

	dt.foo.pt.fpm<-merge(dt.melt.pt.Fpm, dt.pt.samples) # FPM of each samples
	dt.sample.pt.fpm<-merge(dt.foo.pt.fpm, dt.ensg, by="ensembl_gene_id") # FPM of each samples with gene names

	dt.foo.pt.count<-merge(dt.melt.pt.Count, dt.pt.samples) # Count of each samples
	dt.sample.pt.count<-merge(dt.foo.pt.count, dt.ensg, by="ensembl_gene_id") # Count of each samples with gene names

	################################################################################
	## Top expressed plasma chrY genes and the expression level from the placenta ##
	################################################################################
	dt.pt.avg.fpkm<-data.table::dcast(dt.sample.pt.fpkm[CRN %in% dt.pt.samples[dt.pt.samples$CRN %in% dt.samples$CRN]$CRN # placenta samples available with plasma only (n=22)
						& ensembl_gene_id %in% dt.plasma.fpkm[hgnc_symbol!='PSG7']$ensembl_gene_id,list(.N,FPKM=mean(FPKM)),"Sex,ensembl_gene_id"] # chrY and some targets only
						, ensembl_gene_id ~ Sex, value.var="FPKM")
	dt.pt.avg.count<-data.table::dcast(dt.sample.pt.count[CRN %in% dt.pt.samples[dt.pt.samples$CRN %in% dt.samples$CRN]$CRN # placenta samples available with plasma only
						& ensembl_gene_id %in% dt.plasma.fpkm[hgnc_symbol!='PSG7']$ensembl_gene_id,list(.N,Count=mean(Count)),"Sex,ensembl_gene_id"] # chrY and some targets only
						, ensembl_gene_id ~ Sex, value.var="Count")
	dt.pl.pt.fpkm<-merge(dt.plasma.fpkm, dt.pt.avg.fpkm, by="ensembl_gene_id")[order(-all)]
	dt.pl.pt.count<-merge(dt.plasma.count, dt.pt.avg.count, by="ensembl_gene_id")[order(-all)]

	# select genes of interest
	my.target.ensg<-dt.pt.avg.count[M>=10 & M/(F+1)>=10,ensembl_gene_id]


	my.mat<-list()
	my.mat[["fpkm"]]<-as.matrix(dt.pl.pt.fpkm[ensembl_gene_id %in% my.target.ensg,list(`12_PE75_M`, `20_PE75_M`, `28_PE75_M`, `36_PE75_M`, `36_PE150_M`, `12_PE75_F`, `20_PE75_F`, `28_PE75_F`, `36_PE75_F`, M, F)][order(-M)])	
	colnames(my.mat[["fpkm"]])<-c("12_PE75_M", "20_PE75_M", "28_PE75_M", "36_PE75_M", "36_PE150_M", "12_PE75_F", "20_PE75_F", "28_PE75_F", "36_PE75_F", "Placenta_M", "Placenta_F")
	rownames(my.mat[["fpkm"]])<-as.character(dt.pl.pt.fpkm[ensembl_gene_id %in% my.target.ensg][order(-M)]$hgnc_symbol)

	my.mat[["count"]]<-as.matrix(dt.pl.pt.count[ensembl_gene_id %in% my.target.ensg,list(`12_PE75_M`, `20_PE75_M`, `28_PE75_M`, `36_PE75_M`, `36_PE150_M`, `12_PE75_F`, `20_PE75_F`, `28_PE75_F`, `36_PE75_F`, M, F)][order(-M)])	
	colnames(my.mat[["count"]])<-c("12_PE75_M", "20_PE75_M", "28_PE75_M", "36_PE75_M", "36_PE150_M", "12_PE75_F", "20_PE75_F", "28_PE75_F", "36_PE75_F", "Placenta_M", "Placenta_F")
	rownames(my.mat[["count"]])<-as.character(dt.pl.pt.count[ensembl_gene_id %in% my.target.ensg][order(-M)]$hgnc_symbol)

	my.filename <- file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.plasma.heatmap.pt.10"))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=8.27,height=11.7,units="in",res=300, compression = 'lzw') # A4 size
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=8.27,height=11.7,units="in",res=300) # A4 size
    }
	heatmap.2(log(my.mat[["fpkm"]]+0.001), Rowv=FALSE, Colv=FALSE, scale="none", trace="none", dendrogram="none", col = hmcol, keysize=.7, key.xlab="log(FPKM)", density.info="none", main="Plasma ChrY Genes of FPKM>0", cellnote=round(my.mat[["fpkm"]],3), notecol="black", notecex=.7,  margins=c(6,7), srtCol=45, colCol=c(rep(my.col$Sex[["M"]], 5), rep(my.col$Sex[["F"]], 4), my.col$Sex))
	dev.off()

	my.filename <- file.path(cluster.dir, paste0(sampleType, ".chrY.read.count.plasma.heatmap.pt.10"))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=8.27,height=11.7,units="in",res=300, compression = 'lzw') # A4 size
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=8.27,height=11.7,units="in",res=300) # A4 size
    }
	heatmap.2(log(my.mat[["count"]]+0.001), Rowv=FALSE, Colv=FALSE, scale="none", trace="none", dendrogram="none", col = hmcol, keysize=.7, key.xlab="log(Count)", density.info="none", main="Plasma ChrY Genes of Read-Count>0", cellnote=round(my.mat[["count"]],1), notecol="black",  notecex=.7, margins=c(6,7), srtCol=45, colCol=c(rep(my.col$Sex[["M"]], 5), rep(my.col$Sex[["F"]], 4), my.col$Sex))
	dev.off()

	write.csv(my.mat[["fpkm"]], file=file.path(cluster.dir, paste0(sampleType, ".chrY.fpkm.plasma.pt.10.csv")))
	write.csv(my.mat[["count"]], file=file.path(cluster.dir, paste0(sampleType, ".chrY.read.count.plasma.pt.10.csv")))


	####################
	## Per-individual ##
	dt.all.fpkm<-merge(
					dt.sample.pt.fpkm[,.(CRN,Sex,ensembl_gene_id,hgnc_symbol,FPKM)], # placenta
					dcast(dt.sample.fpkm[,.(CRN,GA,SeqType,ensembl_gene_id,FPKM)], CRN+ensembl_gene_id~GA+SeqType, value.var="FPKM"), # plasma
					by=c("ensembl_gene_id","CRN")
					)
	dt.all.count<-merge(
					dt.sample.pt.count[,.(CRN,Sex,ensembl_gene_id,hgnc_symbol,Count)], # placenta
					dcast(dt.sample.count[,.(CRN,GA,SeqType,ensembl_gene_id,Count)], CRN+ensembl_gene_id~GA+SeqType, value.var="Count"), # plasma
					by=c("ensembl_gene_id","CRN")
					)



	#########################################################################
	## Sample-wise Comparsion between Plasma (GA=12,20,28,36) and Placenta ##
	########################################################################
	my.target=dt.target[hgnc_symbol=="PSG7",as.character(ensembl_gene_id)]
	dl.target<-list()
	for(my.target in dt.target[chromosome_name!="Y",ensembl_gene_id]){
		dt.pt.pl<-rbind(
			dt.sample.fpkm[ensembl_gene_id==my.target & SeqType=="PE75",.(CRN,Sex,Source,SampleID,SeqType,GA,FPKM)],
			dt.sample.fpkm[ensembl_gene_id==my.target & SeqType=="PE150",.(CRN,Sex,Source,SampleID,SeqType,GA,FPKM)],
			#dt.sample.pt.fpkm[ensembl_gene_id==my.target,.(CRN,Sex,Source='Placenta',SampleID,SeqType="SE125",GA='Placenta',FPKM)]
			dt.sample.pt.fpkm[ensembl_gene_id==my.target & CRN %in% dt.samples$CRN,.(CRN,Sex,Source='Placenta',SampleID,SeqType="SE125",GA=36,FPKM)]
			)
		dl.target[[my.target]]<-data.table::dcast(dt.pt.pl, CRN+Sex~Source+GA+SeqType,value.var="FPKM")[order(-Placenta_36_SE125)] # long (row-based) to wide (column-based)
	}
		#dt.target[[my.target]][,CRN:=NULL] # remove CRN
		#this.col=c("12","20","28","36","Placenta")
		#dt.dummy[[my.target]][, (this.col) := lapply(.SD,as.character), .SDcols=this.col] # numeric to character for the sake of NA to ""
		## http://stackoverflow.com/questions/20535505/replacing-all-missing-values-in-r-data-table-with-a-value
		#for (i in seq_along(dt.dummy[[my.target]])) set(dt.dummy[[my.target]], i=which(is.na(dt.dummy[[my.target]][[i]])), j=i, value="")
		#dt.dummy[[my.target]][,Sex:=ifelse(Sex=="F",1,0)] # 1 for F; 0 for M
		#write.csv(dt.dummy[[my.target]], file=file.path(cluster.dir,'Ulla', paste(myProject,my.target, "plasma.breech.csv",sep=".")), row.names=F)

	# PSG7
	my.file.name<- paste0(cluster.dir,'/PSG7.ROC.by.GA.',sampleType)
	pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=8.3, height=8.3, title=paste0(myProject,":",sampleType)) # A4 size

	dt.test<-dl.target[["ENSG00000221878"]][,.(`outcome`=factor(ifelse(Placenta_36_SE125>100,"High","Low")), `Placenta`=Placenta_36_SE125, `12WK`=Plasma_12_PE75, `20WK`=Plasma_20_PE75, `28WK`=Plasma_28_PE75, `36WK`=Plasma_36_PE75, `36KW.Deep`=Plasma_36_PE150)]
	my.roc.12<-plot.roc(dt.test$outcome, dt.test$`12WK`, main="12WK", legacy.axes=TRUE, print.auc=TRUE, ci=T)
	my.roc.20<-plot.roc(dt.test$outcome, dt.test$`20WK`, main="20WK",legacy.axes=TRUE, print.auc=TRUE, ci=T)
	my.roc.28<-plot.roc(dt.test$outcome, dt.test$`28WK`, main="28WK",legacy.axes=TRUE, print.auc=TRUE, ci=T)
	my.roc.36<-plot.roc(dt.test$outcome, dt.test$`36WK`, main="36WK",legacy.axes=TRUE, print.auc=TRUE, ci=T)

	my.roc.12<-plot.roc(dt.test$outcome, dt.test$`12WK`, main="AUC of PSG7 by GA ", legacy.axes=TRUE, ci=T)
	my.line.20<-lines.roc(dt.test$outcome, dt.test$`20WK`, col=cbPalette[2])
	my.line.28<-lines.roc(dt.test$outcome, dt.test$`28WK`, col=cbPalette[3])
	my.line.36<-lines.roc(dt.test$outcome, dt.test$`36WK`, col=cbPalette[4])
	text(.7, .85, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.roc.12$auc)), 2)))
	text(.7, .95, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.line.20$auc)), 2)), col=cbPalette[2])
	text(.65, .1, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.line.28$auc)), 2)), col=cbPalette[3])
	text(.65, .2, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.line.36$auc)), 2)), col=cbPalette[4])
	legend("bottomright", legend=c("12Wk","20Wk","28Wk","36Wk"), col=cbPalette2, lwd=2)

	dev.off()

	write.csv(dt.test, file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'FPKM','csv', sep ='.'), row.names=F) 


	if(FALSE){
		dt.foo<-rbind(dt.Fpkm[,list(ensembl_gene_id, SampleID="420","PE75"=`420`,"PE150"=`420a`)],dt.Fpkm[,list(ensembl_gene_id, SampleID="710","PE75"=`710`,"PE150"=`710a`)])
		dt.foo<-merge(dt.foo, dt.ensg)
		
		p.fpkm2<-ggplot(dt.foo[PE75>0 & PE150>0], aes(log(PE75),log(PE150))) + 
				geom_point(size=5,alpha=.1) + 
				labs(x="PE75",y="PE150") +
				ggtitle("log(FPKM)") + 
				facet_grid(~SampleID) + 
				theme_Publication()
		my.filename <- file.path(cluster.dir, paste0(sampleType, "PE75.PE150.fpkm"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=16,height=8,units="in",res=300, compression = 'lzw') # A4 size
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=16,height=8,units="in",res=300) # A4 size
        }
		print(p.fpkm2)
		dev.off()

		p.fpkm3<-ggplot(dt.foo[chromosome_name=="Y"], aes(log(PE75),log(PE150))) + 
				geom_point(size=5,alpha=.1) + 
				labs(x="PE75",y="PE150") +
				ggtitle("log(FPKM)") + 
				facet_grid(~SampleID) + 
				theme_Publication()

		p.fpkm4<-ggplot(dt.foo[chromosome_name=="Y"], aes(PE75,PE150,label=hgnc_symbol)) + 
				geom_point(size=5,alpha=.5) + 
				geom_text(check_overlap = TRUE, hjust =1.2, nudge_x=0.05, size=5) +
				ylim(c(0,13)) + xlim(c(0,13)) + 
				labs(x="PE75",y="PE150") +
				ggtitle("FPKM (chrY only)") + 
				facet_grid(~SampleID) + 
				theme_Publication()

		my.filename <- file.path(cluster.dir, paste0(sampleType, "PE75.PE150.fpkm.chrY"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=16,height=8,units="in",res=300, compression = 'lzw') # A4 size
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=16,height=8,units="in",res=300) # A4 size
        }
		print(p.fpkm4)
		dev.off()

		p.fpkm5<-ggplot(dt.foo[chromosome_name=="Y"], aes(PE75,PE150,label=hgnc_symbol)) + 
				geom_point(size=5,alpha=.5) + 
				geom_text(check_overlap = TRUE, hjust =1.2, nudge_x=0.05, size=5) +
				coord_cartesian(ylim=c(0,1), xlim=c(0,1)) + 
				labs(x="PE75",y="PE150") +
				ggtitle("FPKM (chrY & fpkm<1 only)") + 
				facet_grid(~SampleID) + 
				theme_Publication()

		my.filename <- file.path(cluster.dir, paste0(sampleType, "PE75.PE150.fpkm.chrY.fpkm.1"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=16,height=8,units="in",res=300, compression = 'lzw') # A4 size
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=16,height=8,units="in",res=300) # A4 size
        }
		print(p.fpkm5)
		dev.off()

		my.filename <- file.path(cluster.dir, paste0(sampleType, "PE75.PE150.fpkm.chrY.both"))
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=16,height=16,units="in",res=300, compression = 'lzw') # A4 size
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=16,height=16,units="in",res=300) # A4 size
        }
		multiplot(p.fpkm4, p.fpkm5)
		dev.off()
	}
} # end of Plasma.2017"

if(grepl("RPT",myProject)){ # re-constructed trascriptome based on our own RNA-Seq
	library(VennDiagram)
	dt.target<-as.data.table(as.data.frame(my.target))
	# re-constructed transcriptome
	dt.res.rpt<-merge(dt.res, unique(dt.target[type=="exon",.(gene_name,gene_id)]))

	# ref-based (based on featureCount)
    if(FALSE){
        load("~/results/RNA-Seq/PET.GRCh38.featureCount/DESeq2/ALL/deseq.ALL.RData")
        dt.res.fc<-data.table(ensembl_gene_id=rownames(res),data.frame(res))
        dt.res.fc<-merge(dt.res.fc, dt.ensg[,.(ensembl_gene_id,hgnc_symbol,gene_biotype)])
    }
	# ref-based (based on salmon)
	#load("~/results/RNA-Seq/PET.GRCh38.salmon/DESeq2/ALL/deseq.ALL.RData")
	load("~/results/RNA-Seq/SGA.AGA.SE125.GRCh38.salmon/DESeq2/ALL/deseq.ALL.RData")
	dt.res.salmon<-data.table(ensembl_gene_id=rownames(res),data.frame(res))
	dt.res.salmon<-merge(dt.res.salmon, dt.ensg[,.(ensembl_gene_id,hgnc_symbol,gene_biotype)])[order(pvalue)]

	# 1. RPT (novel via salmon) vs Ref (exon-union by featureCount)
	# merge placentome and the legacy (ref-based)
	# legacy either P>0.01 or NA
    if(FALSE){
        dt.foo<-merge(dt.res.fc, dt.res.rpt, by.x="hgnc_symbol", by.y="gene_name", all.y=T)
        dt.foo<-rbind(dt.foo[padj.x>0.01 & padj.y<0.01], dt.foo[padj.y<0.01 & is.na(padj.x)] )
        fields <- c("ensembl_gene_id","description")
        dt.foo<-merge(dt.foo, as.data.table(getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.foo$ensembl_gene_id, mart = myMart)))
        write.csv(dt.foo, file=file.path(deseq.dir, paste0(sampleType, ".novel10.not.in.reference.csv")))

        venn.list=list()
        venn.list[["Novel.10"]]=dt.res.rpt[padj<0.01,unique(gene_name)]
        venn.list[["Reference"]]=dt.res.fc[padj<0.01,unique(hgnc_symbol)]
        if(capabilities()[["tiff"]]){
            my.filename <- file.path(cluster.dir, paste0(sampleType, ".venn.novel10.vs.legacy.tiff"))
        }else{
            my.filename <- file.path(cluster.dir, paste0(sampleType, ".venn.novel10.vs.legacy.jpeg"))
        }
        venn.diagram(
            x=venn.list,
            filename = my.filename,
            col = "black",
            fill = c(cbPalette[1],cbPalette[2]),
            alpha = 0.5,
            cex = 1.4,
            fontfamily = "sans",
            cat.cex = 1.1,
            cat.fontface = "bold",
            margin = 0.05
        )
    }

	# 2. RPT (novel via salmon) vs Ref (transcript-aware by salmon))
    my.pval=0.05
	dt.foo<-merge(dt.res.salmon, dt.res.rpt, by.x="hgnc_symbol", by.y="gene_name", all.y=T)
	dt.foo<-rbind(dt.foo[padj.x>my.pval & padj.y<my.pval], dt.foo[padj.y<my.pval & is.na(padj.x)] )
	fields <- c("ensembl_gene_id","description")
	#dt.foo<-merge(dt.foo, as.data.table(getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.foo$ensembl_gene_id, mart = myMart)))
	write.csv(dt.foo, file=file.path(deseq.dir, paste0(sampleType, ".novel10.not.in.reference.salmon.csv")), row.names=FALSE, quote=FALSE)

	venn.list=list()
	venn.list[["Novel.10"]]=dt.res.rpt[padj<my.pval,unique(gene_name)]
	venn.list[["Reference\n(salmon)"]]=dt.res.salmon[padj<my.pval,unique(hgnc_symbol)]
    if(capabilities()[["tiff"]]){
	    my.filename <- file.path(cluster.dir, paste0(sampleType, ".venn.novel10.vs.legacy.salmon.tiff"))
        my.imagetype="tiff"
    }else{
	    my.filename <- file.path(cluster.dir, paste0(sampleType, ".venn.novel10.vs.legacy.salmon.png"))
        my.imagetype="png"
    }
	venn.diagram(
		x=venn.list,
		filename = my.filename,
        imagetype= my.imagetype,
		col = "black",
		fill = c(cbPalette[1],cbPalette[3]),
		alpha = 0.5,
		cex = 1.4,
		fontfamily = "sans",
		cat.cex = 1.1,
		cat.fontface = "bold",
		main.fontface = "bold",
        main.cex= 1.2,
		margin = 0.05
	)

    if(FALSE){
        venn.list=list()
        venn.list[["Novel.10"]]=dt.res.rpt[padj<0.01,unique(gene_name)]
        venn.list[["Reference"]]=dt.res.fc[padj<0.01,unique(hgnc_symbol)]
        venn.list[["Reference\n(salmon)"]]=dt.res.salmon[padj<0.01,unique(hgnc_symbol)]
        my.filename <- file.path(cluster.dir, paste0(sampleType, ".venn.novel10.vs.legacy.fc.vs.salmon.tiff"))
        venn.diagram(
            x=venn.list,
            filename = my.filename,
            col = "black",
            fill = c(cbPalette[1],cbPalette[2],cbPalette[3]),
            alpha = 0.5,
            cex = 1.4,
            fontfamily = "sans",
            cat.cex = 1.1,
            cat.fontface = "bold",
            margin = 0.05
        )
    }
}# end of RPT

# plot FC over read-count
cat("plot FC over read-count","\n")

#################
# Volcanio Plot #
#################
foo<-data.table(meanFpkm=rowMeans(ddsFpkm,na.rm=T), 
                meanFpkm.case=rowMeans(ddsFpkm[,rn.case],na.rm=T), 
                meanFpkm.control=rowMeans(ddsFpkm[,rn.control],na.rm=T))
foo[[my.filter]]<-rownames(ddsFpkm)
deseq.anno=merge(dt.res, foo, all.x=TRUE)
# miRNA, piRNA, exon, RPT
if(exists("my.target")){
    bar<-unique(with(as.data.frame(my.target), data.frame(seqnames,Name)))
    colnames(bar)<-c("chromosome_name",my.filter)
    deseq.anno=merge(deseq.anno, bar, all.x=TRUE)[order(pvalue)]
# ensembl_gene_id based analysis
}else{
    deseq.anno<-merge(deseq.anno, dt.ensg[,.(ensembl_gene_id,chromosome_name,hgnc_symbol,gene_biotype)], all.x=TRUE)[order(pvalue)]
}
# 1. all genes
p1 <- ggplot(deseq.anno, aes(log2FoldChange,-log10(padj))) + 
    geom_point(data=deseq.anno[padj>=0.05],alpha=0.5,size=1) + 
    geom_point(data=deseq.anno[padj<0.05 & log2FoldChange>=0],alpha=0.5,size=1,col="red") + 
    geom_point(data=deseq.anno[padj<0.05 & log2FoldChange<0],alpha=0.5,size=1,col="blue") + 
    geom_hline(aes(yintercept=-log10(0.05)),col="grey", linetype="dashed") +
    theme_Publication()
print(p1)
# 2. padj selected
my.pval=0.05
p2 <- ggplot(deseq.anno[padj<my.pval], aes(log2FoldChange,-log10(padj))) + 
    geom_point(alpha=0.8,size=1.5) + 
    geom_hline(yintercept=-log10(my.pval),col="grey", linetype="dashed") +
    geom_vline(xintercept=c(-log2(1.5),log2(1.5)),col="grey", linetype="dashed") +
	geom_text(data=deseq.anno[padj<my.pval & log2FoldChange>=0],aes_string(label=my.id),col="red",hjust=0,nudge_x=0.05, nudge_y=0.02,size=4.5) + 
	geom_text(data=deseq.anno[padj<my.pval & log2FoldChange<0],aes_string(label=my.id),col="blue",hjust=0,nudge_x=0.05, nudge_y=0.02,size=4.5) + 
    ggtitle(paste(nrow(deseq.anno[padj<my.pval]),"DEGs of padj<",my.pval)) +
    theme_Publication()
print(p2)
# 3. top 100 
p3 <- ggplot(deseq.anno[1:100], aes(log2FoldChange,-log10(padj))) + 
    geom_point(alpha=0.5,size=1.5) + 
    geom_hline(yintercept=-log10(0.05),col="grey", linetype="dashed") +
    geom_vline(xintercept=c(-log2(1.5),log2(1.5)),col="grey", linetype="dashed") +
	geom_text(data=deseq.anno[padj<0.05 & log2FoldChange>=0],aes_string(label=my.id),col="red",hjust=0,nudge_x=0.05, nudge_y=0.02,size=4.5) + 
	geom_text(data=deseq.anno[padj<0.05 & log2FoldChange<0],aes_string(label=my.id),col="blue",hjust=0,nudge_x=0.05, nudge_y=0.02,size=4.5) + 
    ggtitle("Top 100 DEGs") +
    theme_Publication()
print(p3)

#########################################
# save the result & FPKM into .csv file #
#########################################
cat("writing res and fpkm to file\n")
write.csv(deseq.anno, file=gzfile(paste0(deseq.dir,"/toptags_all_deseq.",sampleType,".csv.gz")), row.names=F, quote=F)

###########################################
# conventional filtering (FDR<0.1 & FC>2) #
###########################################
if(FALSE){
#    deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj) & (my.res$log2FoldChange > 1 | my.res$log2FoldChange < -1) & my.res$padj <= 0.1, ]
    deseq.top.deg=dt.res[ !is.na(log2FoldChange) & !is.na(padj) & (log2FoldChange > 1 | log2FoldChange < -1) & padj <= 0.1]
    if(nrow(deseq.top.deg) > 1){
        plotMA(res, main=paste0("conventional filtering (FDR<=0.1 & |logFC|>1): ", nrow(deseq.top.deg), " genes (", nrow(dt.res[log2FoldChange<0]), "<0, ", nrow(dt.res[log2FoldChange>=0]), ">=0)"), ylim=c(-2,2)) # default q-value: alpha = 0.1
    }
}

#################
# for all genes
#################
if(FALSE){
    cat("no filtering","\n")
    deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange), ]
    plotMA(my.res, main=paste0("No filtering: ", nrow(deseq.top.deg), " genes (", table(my.res$log2FoldChange<0)["TRUE"], "<0, ", table(my.res$log2FoldChange<0)["FALSE"], ">=0)"), ylim=c(-2,2)) # default q-value: alpha = 0.1
    abline(v=c(0.01,1,minRead), h=c(-1,1), col="blue")

    hist(my.res$pvalue, breaks=100, xlab='P-value', main=paste0('P value distribution of ', nrow(my.res[ !is.na(my.res$pvalue),]), ' genes'))
    hist(my.res[my.res$baseMean>=minRead,]$pvalue, breaks=100, xlab='P-value', main=paste0('P value distribution of ', nrow(my.res[ my.res$baseMean>=minRead,]), ' genes  of minRead>=', minRead))
}

#######################################
# filter genes by mean count & p-value
#######################################
if(FALSE){
    cat("filter by mean count & p-value\n")
    deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$pvalue) & my.res$baseMean>=minRead & my.res$pvalue<=myPValue, ]
    if(nrow(deseq.top.deg) > 1){
        hist(deseq.top.deg$padj, main=paste0("FDR distribution of ", nrow(deseq.top.deg) ," genes having P-value<=",myPValue))

        plotMA(deseq.top.deg, main=paste0(nrow(deseq.top.deg), " genes (",table(deseq.top.deg$log2FoldChange<0)["TRUE"] , "<0, ", table(deseq.top.deg$log2FoldChange<0)["FALSE"], ">=0) having mean>=",minRead, "& P-value<",myPValue), ylim=c(-2,2), alpha=0) # default q-value: alpha = 0.1
        abline(v=c(minRead), col="blue")
    }
}
##############################################
# filter genes by mean count & p-value & FDR
##############################################
if(FALSE){
    cat("filter by mean count & p-value & FDR\n")
    deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$pvalue) & !is.na(my.res$padj) & my.res$baseMean>=minRead & my.res$pvalue<=myPValue & my.res$padj<=myFDR, ]
    if(nrow(deseq.top.deg) > 1){
        plotMA(deseq.top.deg, main=paste0(nrow(deseq.top.deg), " genes (",table(deseq.top.deg$log2FoldChange<0)["TRUE"] , "<0, ", table(deseq.top.deg$log2FoldChange<0)["FALSE"], ">=0) having mean>=",minRead, "& P-value<=",myPValue," & FDR<=",myFDR), ylim=c(-2,2), alpha=0) # default q-value: alpha = 0.1
        abline(v=c(minRead), col="blue")
    }
}

################################################
## Plot expression level of a gene of interest #
################################################
my.entry<-dt.res[order(pvalue)][[my.filter]][1]
plotExp(my.entry,"FPM",my.contrast,my.box=FALSE,my.save=T)
plotExp(my.entry,"Count",my.contrast,my.box=FALSE,my.save=T)
#plotExp(my.entry,"FPM",my.contrast,my.box=TRUE,my.save=T)

##########################
## Top DEGs padj<=0.01  ##
##########################
cat("Top DEGs of padj<=0.01\n")
deseq.top.deg=my.res[!is.na(my.res$padj) & my.res$padj<=0.01, ]
if(nrow(deseq.top.deg)>0){getTopDeseq(deseq.top.deg,"padj.01")}
##########################
## Top DEGs padj<=0.05  ##
##########################
cat("Top DEGs of padj<=0.05\n")
deseq.top.deg=my.res[!is.na(my.res$padj) & my.res$padj<=0.05, ]
if(nrow(deseq.top.deg)>0){getTopDeseq(deseq.top.deg,"padj.05")}

########################
## Top DEGs padj<=0.2 ##
########################
if(FALSE){
    cat("Top DEGs of padj<=0.2\n")
    #deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj) & my.res$baseMean>=minRead & my.res$padj<=0.2, ] # at least 20 read
    deseq.top.deg=my.res[!is.na(my.res$padj) & my.res$padj<=0.2, ]
    if (nrow(deseq.top.deg)>0) getTopDeseq(deseq.top.deg,"padj.2")
}
###################
## Top 100 genes ##
###################
cat("Top100 genes by p-value\n")
#deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj) & my.res$baseMean>=minRead, ] # at least 20 read
deseq.top.deg=my.res[ !is.na(my.res$log2FoldChange) & !is.na(my.res$padj), ]
deseq.top.deg<-head(deseq.top.deg[order(deseq.top.deg$pvalue),],n=100)
if (nrow(deseq.top.deg)>0) getTopDeseq(deseq.top.deg,"top100")

#################
# Samples Info ##
#################
write.csv(colData(dds), file=gzfile(file.path(deseq.dir,paste0(sampleType,".samples.txt.gz"))), row.names=F, quote=F)
dev.off()
cat("All is done\n")

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
}#end of ANNOTATION 

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
		resNoFilt <- results(ddsFlt, parallel=TRUE, independentFiltering=FALSE) # this will disable independant filtering
	}else{
		resNoFilt <- results(dds, parallel=TRUE, independentFiltering=FALSE) # this will disable independant filtering
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
}# end of FALSE

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
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.exp.exon.pval.by.sex.dual.",TR_PREFIX))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=8,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=8,units="in",res=300)
    }
	ggplot_dual_axis(p.cnt,p.cnt.pval,"y")
	dev.off()

	# 2. log-count & p-val
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.log.cnt.pval.per.exon.by.sex.dual3.",TR_PREFIX))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=8,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=8,units="in",res=300)
    }
	ggplot_dual_axis(p.cnt3,p.cnt.meth.pval,"y") # p-val (exp) & p-val (meth)
	#ggplot_dual_axis(p.cnt3,p.cnt.pval,"y")
	#ggplot_dual_axis(p.cnt.pval,p.cnt3,"y")
	dev.off()
		
	# 3. meth & p-val
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.exon.pval.by.sex.dual.",TR_PREFIX))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=8,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=8,units="in",res=300)
    }
	ggplot_dual_axis(p.meth,p.cnt.pval,"y")
	dev.off()

	# 4. meth & p-val (meth)
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.exon.pval2.by.sex.dual.",TR_PREFIX,))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=8,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=8,units="in",res=300)
    }
	ggplot_dual_axis(p.meth,p.meth.pval,"y")
	dev.off()

	# 5-1. cnt & meth
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.exp.exon.by.sex.dual.",TR_PREFIX))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=7,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=7,units="in",res=300)
    }
	ggplot_dual_axis(p.cnt2,p.meth2,"y")
	dev.off()

	# 5-2. meth & cnt
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.exp2.exon.by.sex.dual.",TR_PREFIX))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=7,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=7,units="in",res=300)
q   }
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
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.diff.exp.by.exon.by.sex.",TR_PREFIX))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=8.5,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=8.5,units="in",res=300)
    }
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
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.exp.exon.pval.by.sex.",TR_PREFIX))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=4,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=4,units="in",res=300)
    }
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
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.exp.cnt.exon.pval.by.sex.",TR_PREFIX))
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=10,height=4,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=4,units="in",res=300)
    }
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

    if(capabilities()[["tiff"]]){
	    my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.domain.fpkm.pheatmap.by.sex.",TR_PREFIX,".tiff"))
    }else{
	    my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.domain.fpkm.pheatmap.by.sex.",TR_PREFIX,".jpeg"))
    }
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

    if(capabilities()[["tiff"]]){
	    my.filename<-"~/results/RNA-Seq/RoadMap.exon.GRCh37/Cluster/Consortium/RoadMap.POPS.CSMD1.domain.fpkm.pheatmap.by.tissue.GRCh37.tiff"
    }else{
	    my.filename<-"~/results/RNA-Seq/RoadMap.exon.GRCh37/Cluster/Consortium/RoadMap.POPS.CSMD1.domain.fpkm.pheatmap.by.tissue.GRCh37.jpeg"
    }
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

		my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.domain.fpkm.heatmap.by.sex.",TR_PREFIX,".v2"))
f
        if(capabilities()[["tiff"]]){
		    tiff(filename=paste0(my.filename,".tiff"),width=10,height=3,units="in",res=300, compression = 'lzw')
        }else{
		    jpeg(filename=paste0(my.filename,".jpeg"),width=10,height=3,units="in",res=300)
        }
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

	my.filename <- file.path("~/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/Cluster/AGA/CSMD1.PT.RoadMap.by.sex")
    if(capabilities()[["tiff"]]){
	    tiff(filename=paste0(my.filename,".tiff"),width=12,height=10,units="in",res=300, compression = 'lzw')
    }else{
	    jpeg(filename=paste0(my.filename,".jpeg"),width=12,height=10,units="in",res=300)
    }
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
} # end of exon

###################### 
## Elabela for Steve #
## PET datasets      #
## GRCh38.82         #
###################### 
if(FALSE){
    my.gene=c("ENSG00000248329","ENSG00000171388", "ENSG00000134817","ENSG00000102755","ENSG00000106991")
    dt.ensg[ensembl_gene_id %in% my.gene]

    write.csv(t(ddsFpm[my.gene,]), file=file.path(deseq.dir,"FPM.APEL.etc.csv"),quote=F)
    #APELA
    wilcox.test(ddsFpm["ENSG00000248329",rn.case], ddsFpm["ENSG00000248329",rn.control], paired=TRUE)$p.value
    #APLN
    wilcox.test(ddsFpm["ENSG00000248329",rn.case], ddsFpm["ENSG00000248329",rn.control], paired=TRUE)$p.value
    #APLNR
    wilcox.test(ddsFpm["ENSG00000134817",rn.case], ddsFpm["ENSG00000134817",rn.control], paired=TRUE)$p.value
    #FLT1
    wilcox.test(ddsFpm["ENSG00000102755",rn.case], ddsFpm["ENSG00000102755",rn.control], paired=TRUE)$p.value
    #ENG
    wilcox.test(ddsFpm["ENSG00000106991",rn.case], ddsFpm["ENSG00000106991",rn.control], paired=TRUE)$p.value
}
