#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
myCaller="Clustering"
source ("~/Pipelines/config/DEG.R") # load config

# pdf output filename 
pdf_file_name<- paste0(cluster.dir,'/cluster.',sampleType)
pdf(file=paste(pdf_file_name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

# DESeq2 RData file
deseq.RData <- paste0(deseq.dir, "/deseq.",sampleType,'.RData')
if(!file.exists(deseq.RData)){
	stop("RData file not exists\n")
}

#####################
## Sample Outliers ##
#####################
#cat("Sample Outlier\n")
#rlogMat[topVarGenes,c("88","75", "96")]
#assay(dds)[topVarGenes,c("88","75", "96")]
#samples.all[samples.all$SampleName %in% c("88","75", "96"),]

###################
## Loading RData ##
###################
cat("loading DESeq2 RData...\n")
load (deseq.RData) # dds,ddsFlt,res,resFlt,rld,vsd

#rld.nb <- rlog(dds,blind=FALSE)
#rlogMat.nb <- assay(rld.nb)
rlogMat <- assay(rld) #
#vstMat <- assay(vsd)

readCount<-counts(dds) # get the raw count
keep= rowSums(readCount >= 1) >= ncol(readCount)/2 # half of the sample should be at least >=1
cat("NO. of genes having at least>=1 read from half of the sample:\n")
print(table(keep))
readCount<-readCount[keep,]
rlogMat<-rlogMat[keep,] # apply the filter

#samples<-samples.all[!(samples.all$Pair %in% samples.all[samples.all$SampleName %in% c("88","75", "96"),]$Pair),] # remove outliers
rn=samples[,c("SampleName")]
rn.sga=samples[samples$Condition==1,c("SampleName")] # SGA only
rn.aga=samples[samples$Condition==0,c("SampleName")] # AGA only

# summary of boxplot
diff.rlogMat<-rlogMat[,rn.sga]-rlogMat[,rn.aga]
summary(rlogMat[,rn.sga]-rlogMat[,rn.aga])

# print Boxplot
# unaffected samples only
if(sampleType=="ALL.unaffected"){
	par(mfrow=c(3,1))
	boxplot(diff.rlogMat[,1:10], main=paste0("Boxplot for rlog(SGA)-rlog(AGA)"), xlab='SGA', ylab='rlog(SGA)-rlog(AGA)')
	boxplot(diff.rlogMat[,11:20], main=paste0("Boxplot for rlog(SGA)-rlog(AGA)"), xlab='SGA', ylab='rlog(SGA)-rlog(AGA)') 
	boxplot(diff.rlogMat[,21:30], main=paste0("Boxplot for rlog(SGA)-rlog(AGA)"), xlab='SGA', ylab='rlog(SGA)-rlog(AGA)')
	par(mfrow=c(1,1))
}

# outliers ( genes where abs(diff)>=xx )
outliers<-rowSums(abs(diff.rlogMat)>=9) >=1
outliers<-diff.rlogMat[outliers,] #rlogMat[outliers,]
outliers<-rownames(outliers[order(rowVars(outliers), decreasing=TRUE ),]) # sort gene by descending order of variation
# diff.rlogMat[outliers,]; rlogMat[outliers,]
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
fields <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description")
outlier.genes<-getBM(attributes = fields, filters = "ensembl_gene_id", values = outliers, mart = ensembl)
rownames(outlier.genes)<-outlier.genes$ensembl_gene_id
outlier.genes[outliers,]

# get DEGs
deseq.result = read.csv(paste0(deseq.dir,"/filtered_toptags_deseq.padj.1.",sampleType,".csv"), stringsAsFactors=FALSE)
deseq.result<-deseq.result[order(deseq.result$pvalue),]

#samplePhen<-samples.phen[samples.phen$SampleName %in% rn.sga,c(7:ncol(samples.phen))]
samplePhen<-samples.phen[samples.phen$SampleName %in% rn.sga,c('GVel','UAD','UBD','BW_PCT')]-samples.phen[samples.phen$SampleName %in% rn.aga,c('GVel','UAD','UBD','BW_PCT')]
samplePhen<-cbind(samplePhen, samples.phen[samples.phen$SampleName %in% rn.sga,c('FS','Smoke','PT','MD')])
for(i in colnames(samplePhen[c("MD","Smoke","FS","PT")]))(samplePhen[i]<-droplevels(samplePhen[i])) # adjust levels

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

makeCluster<-function(i){
	cat(paste0("PCA for Top ",i," %...\n"))
	#select<-deseq.result$ensg # Top DES names 
	#select <- order(rowMeans(counts(dds,normalized=T)),decreasing=TRUE)[1:30] # top 30 highly expressed genes
	#select<- head( order(rowVars(rlogMat[,rn.sga]-rlogMat[,rn.aga]), decreasing=TRUE ),  n=round(nrow(rlogMat)*i*0.01) ) # requires 'genefilter' library
	select<- head( order(rowVars(rlogMat), decreasing=TRUE ), n=round(nrow(rlogMat)*i*0.01) ) # requires 'genefilter' library

	#rlogMat[head(select,n=10),] # first 10 top variables

	#############################################
	# Heatmap of the sample-to-sample distances #
	#############################################
	cat("Heatmap of the sample-to-sample distances\n")
	hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
	mat<-rlogMat[select,rn.sga]-rlogMat[select,rn.aga] # isa matrix (row: ENSG, col: Pair(SGA-AGA))
	#colnames(mat)<-paste(colnames(mat), samples[samples$Condition==1,c("Code")], sep=":")
	#colnames(mat)<-samples[samples$Condition==1,c("PAPPA")]
	#colnames(mat)<-samples[samples$Condition==1,c("HT")]
	#colnames(mat)<-samples[samples$Condition==1,c("GVel")]

	#sampleDists <- as.matrix(dist(t( mat )))
	#heatmap.2(sampleDists, trace="none", col = hmcol, cexRow=0.7, cexCol=0.7)

	#Outliers (based on manual inspection of the Boxplot)
	if(FALSE){
		dummy<-c(
			names(head( mat[order(mat[,"59"], decreasing=TRUE),"59"],n=1 )),
			names(tail( mat[order(mat[,"63b"], decreasing=TRUE),"63b"],n=1 )),
			names(tail( mat[order(mat[,"71"], decreasing=TRUE),"71"],n=1 )),
			names(head( mat[order(mat[,"75"], decreasing=TRUE),"75"],n=7 )),
			names(tail( mat[order(mat[,"87"], decreasing=TRUE),"87"],n=8 )),
			names(head( mat[order(mat[,"95"], decreasing=TRUE),"95"],n=1 )),
			names(head( mat[order(mat[,"105"], decreasing=TRUE),"105"],n=1 )),
			names(head( mat[order(mat[,"53"], decreasing=TRUE),"53"],n=1 )),
			names(head( mat[order(mat[,"115"], decreasing=TRUE),"115"],n=1 )),
			names(tail( mat[order(mat[,"125"], decreasing=TRUE),"125"],n=1 )),
			names(head( mat[order(mat[,"129"], decreasing=TRUE),"129"],n=1 )),
			names(head( mat[order(mat[,"139"], decreasing=TRUE),"139"],n=1 ))
		)
		dummy<-unique(dummy)
	}
	if(FALSE){
	}

	mat<-mat[!(rownames(mat) %in% dummy),] # remove those outliers

	# all samples
	par(mfrow=c(2,1))
	boxplot(mat[,1:10], main=paste0("Top ",i,"%: Boxplot for rlog(SGA)-rlog(AGA)"), xlab='SGA', ylab='rlog(SGA)-rlog(AGA)')
	boxplot(mat[,11:20], main=paste0("Top ",i,"%: Boxplot for rlog(SGA)-rlog(AGA)"), xlab='SGA', ylab='rlog(SGA)-rlog(AGA)') 
	boxplot(mat[,21:30], main=paste0("Top ",i,"%: Boxplot for rlog(SGA)-rlog(AGA)"), xlab='SGA', ylab='rlog(SGA)-rlog(AGA)')
	boxplot(mat[,31:40], main=paste0("Top ",i,"%: Boxplot for rlog(SGA)-rlog(AGA)"), xlab='SGA', ylab='rlog(SGA)-rlog(AGA)')
	boxplot(mat[,41:50], main=paste0("Top ",i,"%: Boxplot for rlog(SGA)-rlog(AGA)"), xlab='SGA', ylab='rlog(SGA)-rlog(AGA)')
	boxplot(mat[,51:56], main=paste0("Top ",i,"%: Boxplot for rlog(SGA)-rlog(AGA)"), xlab='SGA', ylab='rlog(SGA)-rlog(AGA)')
	par(mfrow=c(1,1))

	#######
	# PCA #
	#######
	cat("PCA based on DEG profile\n")
	# 2. Paired Clustering (SGA-AGA)
	#pca<-prcomp(t(mat))
	pca<-prcomp(t(mat), scale.=TRUE)
	percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

	barplot(round(percentVar*100,2), main=paste0("% variance of Top ",i, "% (",length(select),") variable genes from ",sampleType), xlab="Principal Component", ylab="% variance")
	# print sample name
	plot(pca$x[,1:2], main=paste0("PCA of Top ",i, "% (",length(select),") variable genes from ",sampleType), cex=0.1, col=myDotCol, xlab=paste0("PC1: ",round(percentVar[1]*100,2)," % variance"), ylab=paste0("PC2: ",round(percentVar[2]*100,2)," % variance"))
	text(pca$x[,1:2], labels=rownames(pca$x))

	# print phen encoding
	#plot(pca$x[,1:2], main=paste0("PCA of Top ",i, "% (",length(select),") variable genes from ",sampleType), cex=0.1, col=myDotCol, xlab=paste0("PC1: ",round(percentVar[1]*100,2)," % variance"), ylab=paste0("PC2: ",round(percentVar[2]*100,2)," % variance"))
	#text(pca$x[,1:2], labels=samples[samples$Condition==1,c("Code")])

	if(i==100 | i==50 | i==20 | i==10){
		# print by UBD 
		plot(pca$x[,1:2], main=paste0("UBD(SGA-AGA) of Top ",i, "% (",length(select),") variable genes from ",sampleType), cex=0.1, col=myDotCol, xlab=paste0("PC1: ",round(percentVar[1]*100,2)," % variance"), ylab=paste0("PC2: ",round(percentVar[2]*100,2)," % variance"))
		text(pca$x[,1:2], labels=round(samplePhen[,c("UBD"),],2))
	
		#PC1 & UBD
		dummy<-cor(cbind(pca$x[,"PC1"], samplePhen[,'UBD']), use="complete.obs")
		dummy.pvalue<-round(cor.test(pca$x[,"PC1"], samplePhen[,'UBD'])$p.value,3)
		plot(x=pca$x[,"PC1"], y=samplePhen[,'UBD'], main=paste0("Top ", i, " r=",round(dummy[1,2],3), " (P-value=",dummy.pvalue,")"), xlab="PC1", ylab='UBD(SGA-AGA)')
		abline(lm(samplePhen[,'UBD']~pca$x[,"PC1"]), col="blue", lwd=2)
	}
	if(i==100){
		# print by UAD 
		plot(pca$x[,1:2], main=paste0("UAD(SGA-AGA) of Top ",i, "% (",length(select),") variable genes from ",sampleType), cex=0.1, col=myDotCol, xlab=paste0("PC1: ",round(percentVar[1]*100,2)," % variance"), ylab=paste0("PC2: ",round(percentVar[2]*100,2)," % variance"))
		text(pca$x[,1:2], labels=round(samplePhen[,c("UAD"),],2))

		#PC2 & UAD
		dummy<-cor(cbind(pca$x[,"PC2"], samplePhen[,'UAD']), use="complete.obs")
		dummy.pvalue<-round(cor.test(pca$x[,"PC2"], samplePhen[,'UAD'])$p.value,3)
		plot(x=pca$x[,"PC2"], y=samplePhen[,'UAD'], main=paste0("Top ", i, " r=",round(dummy[1,2],3), " (P-value=",dummy.pvalue,")"), xlab="PC2", ylab='UAD(SGA-AGA)')
		abline(lm(samplePhen[,'UAD']~pca$x[,"PC2"]), col="blue", lwd=2)

		# print by FS
		plot(pca$x[,1:2], main=paste0("By fetal sex of Top ",i, "% (",length(select),") variable genes from ",sampleType), cex=0.1, col=myDotCol, xlab=paste0("PC1: ",round(percentVar[1]*100,2)," % variance"), ylab=paste0("PC2: ",round(percentVar[2]*100,2)," % variance"))
		text(pca$x[,1:2], labels=samplePhen[,c("FS")])

		# print by Smoke
		plot(pca$x[,1:2], main=paste0("By 'smoking' of Top ",i, "% (",length(select),") variable genes from ",sampleType), cex=0.1, col=myDotCol, xlab=paste0("PC1: ",round(percentVar[1]*100,2)," % variance"), ylab=paste0("PC2: ",round(percentVar[2]*100,2)," % variance"))
		text(pca$x[,1:2], labels=samplePhen[,c("Smoke")])

		# print by PT 
		plot(pca$x[,1:2], main=paste0("By 'smoking' of Top ",i, "% (",length(select),") variable genes from ",sampleType), cex=0.1, col=myDotCol, xlab=paste0("PC1: ",round(percentVar[1]*100,2)," % variance"), ylab=paste0("PC2: ",round(percentVar[2]*100,2)," % variance"))
		text(pca$x[,1:2], labels=samplePhen[,c("PT")])
	}
	if(i==50){
		# print by UAD 
		plot(pca$x[,1:2], main=paste0("UAD(SGA-AGA) of Top ",i, "% (",length(select),") variable genes from ",sampleType), cex=0.1, col=myDotCol, xlab=paste0("PC1: ",round(percentVar[1]*100,2)," % variance"), ylab=paste0("PC2: ",round(percentVar[2]*100,2)," % variance"))
		text(pca$x[,1:2], labels=round(samplePhen[,c("UAD"),],2))

		#PC3 & UAD
		dummy<-cor(cbind(pca$x[,"PC3"], samplePhen[,'UAD']), use="complete.obs")
		dummy.pvalue<-round(cor.test(pca$x[,"PC3"], samplePhen[,'UAD'])$p.value,3)
		plot(x=pca$x[,"PC3"], y=samplePhen[,'UAD'], main=paste0("Top ", i, " r=",round(dummy[1,2],3), " (P-value=",dummy.pvalue,")"), xlab="PC3", ylab='UAD(SGA-AGA)')
		abline(lm(samplePhen[,'UAD']~pca$x[,"PC3"]), col="blue", lwd=2)
	}

	#print by MD
	plot(pca$x[,1:2], main=paste0("By mode of deliver from Top ",i, "% (",length(select),") variable genes from ",sampleType), pch=19, 
		xlab=paste0("PC1: ",round(percentVar[1]*100,2)," % variance"), ylab=paste0("PC2: ",round(percentVar[2]*100,2)," % variance"), col=samplePhen[,c("MD")])
	legend("topright", legend=levels(samplePhen[,c("MD")]), col=1:length(samplePhen[,c("MD")]), cex = 0.8, pch = 19)

	#sum((pca$sdev)^2) # should be ncol(x) if scale=T
	#sum((pca$rotation)^2) # should be ncol(x)

	# code from methylKit
	# https://github.com/al2na/methylKit/blob/master/R/batchControl.R
	# get association p-values using different tests
	phenAs=list()
	for(i in 1:ncol(samplePhen)){
		cat(paste0(names(samplePhen[i]),"\n"))
 		# for factors do kruskal.wallis or wilcox test
		if(is.factor(samplePhen[,i]) | is.character(samplePhen[,i]) | is.logical(samplePhen[,i])){
			annot=as.factor(samplePhen[,i])
			phenAs[[names(samplePhen)[i]]]=apply(pca$x,2,function(x){ # column-wise
				if(length(unique(annot))>2 ){
					kruskal.test(split(x,annot))$p.value
            	}else{
            		wilcox.test(split(x,annot)[[1]],split(x,annot)[[2]])$p.value 
            	}
			})
		}else{# for factors do cor.test
			annot=  samplePhen[,i]
			phenAs[[names(samplePhen)[i]]]=apply(pca$x,2,function(x){ # column-wise 
				cor.test(x,annot)$p.value
      		})
		}
	} # end of for i

	return(do.call("rbind",phenAs))
} # end of makeCluster

############################################
phenPvalue=list()
pct<-c(100, 50, 20, 10, 5, 1)
#sapply(pct, makeCluster)
for(i in pct){
	phenPvalue[[paste0("top",i)]]<-makeCluster(i)
}

lapply(phenPvalue, function(i) i[,1:10]) 
############################################

#for Ulla
if(FALSE){
	rlogMat<-as.data.frame(rlogMat)
	fg.rld.diff<-subset(rlogMat,select=samples[samples$Condition==1,c("SampleName")]) - subset(rlogMat,select=samples[samples$Condition==0,c("SampleName")]) # gene-wise
	fg.rld.diff<-t(fg.rld.diff) # sample-wise (isa 'matrix')
	fg.rld.diff<-cbind('sampleName'=rownames(fg.rld.diff), fg.rld.diff)
	write.csv(fg.rld.diff, "~/scratch/results/RNA-Seq/SGA.AGA/Ulla/fg.rld.diff.csv", row.names=FALSE, quote=FALSE) # gene-wise

	read.counts<-as.data.frame(counts(dds))
	fg.count.diff<-subset(read.counts,select=samples[samples$Condition==1,c("SampleName")])-subset(read.counts,select=samples[samples$Condition==0,c("SampleName")])
	fg.count.diff<-t(fg.count.diff) # sample-wise
	fg.count.diff<-cbind('sampleName'=rownames(fg.count.diff), fg.count.diff)
	write.csv(fg.count.diff, "~/scratch/results/RNA-Seq/SGA.AGA/Ulla/fg.count.diff.csv", row.names=FALSE, quote=FALSE) # gene-wise

	pct<-c(1, 5, 10, 20, 50) 
	makeFile<-function(i){
		select<- head( order(rowVars(rlogMat), decreasing=TRUE ), n=round(nrow(rlogMat)*i*0.01) ) # requires 'genefilter' library

		#select<- head( order( rowVars(fg.count.diff), decreasing=TRUE ), round(nrow(fg.count.diff)*i*0.01)) # requires 'genefilter' library
		dummy<-fg.count.diff[,select]
		dummy<-cbind('sampleName'=rownames(dummy), dummy)
		write.csv(dummy, paste0("~/scratch/results/RNA-Seq/SGA.AGA/Ulla/fg.count.diff.top.",i,".csv"), row.names=FALSE, quote=FALSE) # sample-wise

		#select<- head( order( rowVars(fg.rld.diff), decreasing=TRUE ), round(nrow(fg.rld.diff)*i*0.01)) # requires 'genefilter' library
		dummy<-fg.rld.diff[,select]
		dummy<-cbind('sampleName'=rownames(dummy), dummy)
		write.csv(dummy, paste0("~/scratch/results/RNA-Seq/SGA.AGA/Ulla/fg.rld.diff.top.",i,".csv"), row.names=FALSE, quote=FALSE) # sample-wise
	}
	sapply(pct, makeFile)
}

dev.off()
cat("All is done\n")
