#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

myCaller="Merger"
source ("~/Pipelines/config/DEG.R") # load config

sampleType1<-sampleType
#sampleType2<-"HT"

if(!file.exists(deg.dir)){dir.create(deg.dir)}

mergeDEG("all","all")
mergeDEG("FDR",4)
mergeDEG("FDR",3)
mergeDEG("FDR",2)
mergeDEG("FDR",1)

mergeDEG("top",200)
mergeDEG("top",100)

############
# All Degs #
############
if(FALSE){
edgeR.result1 = read.csv(paste0(edgeR.dir,"/toptags_all_edgeR.",sampleType1,".csv"), stringsAsFactors=FALSE)
deseq.result1 = read.csv(paste0(deseq.dir,"/toptags_all_deseq.",sampleType1,".csv"), stringsAsFactors=FALSE)

names(edgeR.result1)[names(edgeR.result1) == 'X'] <- 'ensg'
names(deseq.result1)[names(deseq.result1) == 'X'] <- 'ensg'

union.top.deg1=merge(edgeR.result1, deseq.result1, by="ensg", all=T) # isa 'data.frame'   # union between two
common.top.deg1=union.top.deg1[!is.na(union.top.deg1$logFC) & !is.na(union.top.deg1$log2FoldChange), ] # common DEG (intersection)

write.csv(union.top.deg1[order(union.top.deg1$pvalue, union.top.deg1$PValue),], file=paste0(deg.dir,"/unified.all.deg.",sampleType1,".csv"))
write.csv(common.top.deg1[order(common.top.deg1$pvalue, common.top.deg1$PValue),], file=paste0(deg.dir,"/common.all.deg.",sampleType1,".csv"))

##############
# Top200/100 #
##############
edgeR.result1 = read.csv(paste0(edgeR.dir,"/filtered_toptags_edgeR.top200.",sampleType1,".csv"), stringsAsFactors=FALSE)
deseq.result1 = read.csv(paste0(deseq.dir,"/filtered_toptags_deseq.top200.",sampleType1,".csv"), stringsAsFactors=FALSE)
union.top.deg1=merge(edgeR.result1, deseq.result1, by="ensg", all=T) # isa 'data.frame'   # union between two
common.top.deg1=union.top.deg1[!is.na(union.top.deg1$X.x) & !is.na(union.top.deg1$X.y), ] # common DEG (intersection)
write.csv(union.top.deg1, file=paste0(deg.dir,"/unified.top.deg.top200.",sampleType1,".csv"))
write.csv(common.top.deg1, file=paste0(deg.dir,"/common.top.deg.top200.",sampleType1,".csv"))

edgeR.result1 = read.csv(paste0(edgeR.dir,"/filtered_toptags_edgeR.top100.",sampleType1,".csv"), stringsAsFactors=FALSE)
deseq.result1 = read.csv(paste0(deseq.dir,"/filtered_toptags_deseq.top100.",sampleType1,".csv"), stringsAsFactors=FALSE)

union.top.deg1=merge(edgeR.result1, deseq.result1, by="ensg", all=T) # isa 'data.frame'   # union between two
common.top.deg1=union.top.deg1[!is.na(union.top.deg1$X.x) & !is.na(union.top.deg1$X.y), ] # common DEG (intersection)

write.csv(union.top.deg1, file=paste0(deg.dir,"/unified.top.deg.top100.",sampleType1,".csv"))
write.csv(common.top.deg1, file=paste0(deg.dir,"/common.top.deg.top100.",sampleType1,".csv"))

###################
# DEGs by q-value #
###################
edgeR.result1 = read.csv(paste0(edgeR.dir,"/filtered_toptags_edgeR.",sampleType1,".csv"), stringsAsFactors=FALSE)
deseq.result1 = read.csv(paste0(deseq.dir,"/filtered_toptags_deseq.",sampleType1,".csv"), stringsAsFactors=FALSE)

union.top.deg1=merge(edgeR.result1, deseq.result1, by="ensg", all=T) # isa 'data.frame'   # union between two
common.top.deg1=union.top.deg1[!is.na(union.top.deg1$X.x) & !is.na(union.top.deg1$X.y), ] # common DEG (intersection)

write.csv(union.top.deg1, file=paste0(deg.dir,"/unified.top.deg.",sampleType1,".csv"))
write.csv(common.top.deg1, file=paste0(deg.dir,"/common.top.deg.",sampleType1,".csv"))
}


if(exists("sampleType2")){
	edgeR.dir2=paste0(edgeRHome,'/',sampleType2)
	deseq.dir2=paste0(deseqHome,'/',sampleType2)
	deg.dir2=paste0(DegHome,'/',sampleType2)

	if(!file.exists(deg.dir2)){dir.create(deg.dir2)}
	####################
	## edgeR & DESeq2 ##
	####################
	edgeR.result2 = read.csv(paste0(edgeR.dir2,"/filtered_toptags_edgeR.",sampleType2,".csv"), stringsAsFactors=FALSE)
	deseq.result2 = read.csv(paste0(deseq.dir2,"/filtered_toptags_deseq.",sampleType2,".csv"), stringsAsFactors=FALSE)

	union.top.deg2=merge(edgeR.result2, deseq.result2, by="ensg", all=T) # isa 'data.frame'
	common.top.deg2=union.top.deg2[!is.na(union.top.deg2$X.x) & !is.na(union.top.deg2$X.y), ] # common DEG (intersection)

	write.csv(union.top.deg2, file=paste0(deg.dir2,"/unified.top.deg.",sampleType2,".csv"))
	write.csv(common.top.deg2, file=paste0(deg.dir2,"/common.top.deg.",sampleType2,".csv"))

	##########################################
	## edgeR from sampleType1 & sampleType2 ##
	##########################################
	deg.dir.sub=paste0(DegHome,'/',sampleType1,'/',sampleType2)
	if(!file.exists(deg.dir.sub)){dir.create(deg.dir.sub)}

	union.top.edgeR=merge(edgeR.result1, edgeR.result2, by="ensg", all=T) # isa 'data.frame'
	common.top.edgeR=union.top.edgeR[!is.na(union.top.edgeR$X.x) & !is.na(union.top.edgeR$X.y), ] # common DEG (intersection)

	write.csv(union.top.edgeR, file=paste0(deg.dir.sub,"/unified.top.deg.edgeR.",sampleType1,".",sampleType2,".csv"))
	write.csv(common.top.edgeR, file=paste0(deg.dir.sub,"/common.top.deg.edgeR.",sampleType1,".",sampleType2,".csv"))

	cat("edgeR of ",sampleType1,":\n")
	print(nrow(edgeR.result1))
	cat("edgeR of ",sampleType2,":\n")
	print(nrow(edgeR.result2))
	cat("Union of edgeR between ",sampleType1," and ", sampleType2,":\n")
	print(nrow(union.top.edgeR))
	cat("Intersection of edgeR between ",sampleType1," and ", sampleType2,":\n")
	print(nrow(common.top.edgeR))

	# FDR<=0.4
	print(nrow(edgeR.result1[edgeR.result1$FDR<=myFDR,]))
	print(nrow(edgeR.result2[edgeR.result2$FDR<=myFDR,]))
	print(nrow(union.top.edgeR[(!is.na(union.top.edgeR$FDR.x) & union.top.edgeR$FDR.x<=myFDR) | (!is.na(union.top.edgeR$FDR.y) & union.top.edgeR$FDR.y<=myFDR),]))
	print(nrow(common.top.edgeR[common.top.edgeR$FDR.x<=myFDR & common.top.edgeR$FDR.y<=myFDR,]))
	###########################################
	## DESeq2 from sampleType1 & sampleType2 ##
	###########################################
	union.top.deseq=merge(deseq.result1, deseq.result2, by="ensg", all=T) # isa 'data.frame'
	common.top.deseq=union.top.deseq[!is.na(union.top.deseq$X.x) & !is.na(union.top.deseq$X.y), ] # common DEG (intersection)

	write.csv(union.top.deseq, file=paste0(deg.dir.sub,"/unified.top.deg.DESeq2.",sampleType1,".",sampleType2,".csv"))
	write.csv(common.top.deseq, file=paste0(deg.dir.sub,"/common.top.deg.DESeq2.",sampleType1,".",sampleType2,".csv"))

	cat("DESeq2 of ",sampleType1,":\n")
	print(nrow(deseq.result1))
	cat("DESeq2 of ",sampleType2,":\n")
	print(nrow(deseq.result2))
	cat("Union of DESeq2 between ",sampleType1," and ", sampleType2,":\n")
	print(nrow(union.top.deseq))
	cat("Intersection of DESeq2 between ",sampleType1," and ", sampleType2,":\n")
	print(nrow(common.top.deseq))

	# FDR<=0.4
	print(nrow(deseq.result1[deseq.result1$padj<=myFDR,]))
	print(nrow(deseq.result2[deseq.result2$padj<=myFDR,]))
	print(nrow(union.top.deseq[(!is.na(union.top.deseq$padj.x) & union.top.deseq$padj.x<=myFDR) | (!is.na(union.top.deseq$padj.y) & union.top.deseq$padj.y<=myFDR),]))
	print(nrow(common.top.deseq[common.top.deseq$padj.x<=myFDR & common.top.deseq$padj.y<=myFDR,]))
} # sampleType2

if(FALSE){
	# Plot common DEG (intersection): coloured coded by the first
	dummy.top.deg<-common.top.deg
	# edgeR
	minlogFC=min(abs(union.top.deg[!is.na(union.top.deg$logFC.x),]$logFC.x))
	plot(y=dummy.top.deg$logFC.x, x=dummy.top.deg$logCPM.x, cex=0.5, ylab='logFC', xlab='logCPM', cex.lab=1.8)
	text(y=dummy.top.deg$logFC.x, x=dummy.top.deg$logCPM.x, cex=2, labels=dummy.top.deg$hgnc_symbol.x, col=dummy.top.deg$col.x)
	abline(h=c(-minlogFC, minlogFC), col="blue")
	legend("topright", legend=names(myCol), fill=myCol, cex=2 )
	# DESeq2
	minlogFC=min(abs(union.top.deg[!is.na(union.top.deg$log2FoldChange.x),]$log2FoldChange.x))
	plot(y=dummy.top.deg$log2FoldChange.x, x=log(dummy.top.deg$baseMean.x), cex=0.5, ylab='logFC', xlab='log10(baseMean)', cex.lab=1.7)
	text(y=dummy.top.deg$log2FoldChange.x, x=log(dummy.top.deg$baseMean.x), cex=2, labels=dummy.top.deg$hgnc_symbol.x, col=dummy.top.deg$col.x)
	abline(h=c(-minlogFC, minlogFC), col="blue")
	legend("topright", legend=names(myCol), fill=myCol, cex=2 ) 

	# Plot common DEG (intersection): coloured coded by the second
	dummy.top.deg<-common.top.deg
	# edgeR
	minlogFC=min(abs(union.top.deg[!is.na(union.top.deg$logFC.y),]$logFC.y))
	plot(y=dummy.top.deg$logFC.y, x=dummy.top.deg$logCPM.y, cex=0.5, ylab='logFC', xlab='logCPM', cex.lab=1.8)
	text(y=dummy.top.deg$logFC.y, x=dummy.top.deg$logCPM.y, cex=2, labels=dummy.top.deg$hgnc_symbol.y, col=as.character(dummy.top.deg$col.y))
	abline(h=c(-minlogFC, minlogFC), col="blue")
	legend("topright", legend=names(myCol), fill=myCol, cex=2 )
	# DESeq2
	minlogFC=min(abs(union.top.deg[!is.na(union.top.deg$log2FoldChange.y),]$log2FoldChange.y))
	plot(y=dummy.top.deg$log2FoldChange.y, x=log(dummy.top.deg$baseMean.y), cex=0.5, ylab='logFC', xlab='log10(baseMean)', cex.lab=1.7)
	text(y=dummy.top.deg$log2FoldChange.y, x=log(dummy.top.deg$baseMean.y), cex=2, labels=dummy.top.deg$hgnc_symbol.y, col=dummy.top.deg$col.y)
	abline(h=c(-minlogFC, minlogFC), col="blue")
	legend("topright", legend=names(myCol), fill=myCol, cex=2 ) 

	#################################
	# DEGs called only by the first #
	#################################
	dummy.top.deg<-union.top.deg[!is.na(union.top.deg$X.x) & is.na(union.top.deg$X.y), ]
	# edgeR
	minlogFC=min(abs(union.top.deg[!is.na(union.top.deg$logFC.x),]$logFC.x))
	plot(y=dummy.top.deg$logFC.x, x=dummy.top.deg$logCPM.x, cex=0.5, ylab='logFC', xlab='logCPM', cex.lab=1.8)
	text(y=dummy.top.deg$logFC.x, x=dummy.top.deg$logCPM.x, cex=2, labels=dummy.top.deg$hgnc_symbol.x, col=dummy.top.deg$col.x)
	abline(h=c(-minlogFC, minlogFC), col="blue")
	legend("topright", legend=names(myCol), fill=myCol, cex=2 )
	# DESeq2
	minlogFC=min(abs(union.top.deg[!is.na(union.top.deg$log2FoldChange.x),]$log2FoldChange.x))
	plot(y=dummy.top.deg$log2FoldChange.x, x=log(dummy.top.deg$baseMean.x), cex=0.5, main=paste0("DESeq2: ",nrow(dummy.top.deg)," DEGs from NHT only (FDR<=0.4)"), ylab='logFC', xlab='log10(baseMean)', cex.lab=1.7)
	text(y=dummy.top.deg$log2FoldChange.x, x=log(dummy.top.deg$baseMean.x), cex=2, labels=dummy.top.deg$hgnc_symbol.x, col=as.character(dummy.top.deg$col.x))
	abline(h=c(-minlogFC, minlogFC), col="blue")
	legend("topright", legend=names(myCol), fill=myCol, cex=2 ) 

	##################################
	# DEGs called only by the second #
	##################################
	dummy.top.deg<-union.top.deg[is.na(union.top.deg$X.x) & !is.na(union.top.deg$X.y), ]
	# edgeR
	minlogFC=min(abs(union.top.deg[!is.na(union.top.deg$logFC.y),]$logFC.y))
	plot(y=dummy.top.deg$logFC.y, x=dummy.top.deg$logCPM.y, cex=0.5, ylab='logFC', xlab='logCPM', cex.lab=1.8)
	text(y=dummy.top.deg$logFC.y, x=dummy.top.deg$logCPM.y, cex=2, labels=dummy.top.deg$hgnc_symbol.y, col=as.character(dummy.top.deg$col.y))
	abline(h=c(-minlogFC, minlogFC), col="blue")
	legend("topright", legend=names(myCol), fill=myCol, cex=2 )
	# DESeq2
	minlogFC=min(abs(union.top.deg[!is.na(union.top.deg$log2FoldChange.y),]$log2FoldChange.y))
	plot(y=dummy.top.deg$log2FoldChange.y, x=log(dummy.top.deg$baseMean.y), cex=0.5, ylab='logFC', xlab='log10(baseMean)', cex.lab=1.7)
	text(y=dummy.top.deg$log2FoldChange.y, x=log(dummy.top.deg$baseMean.y), cex=2, labels=dummy.top.deg$hgnc_symbol.y, col=as.character(dummy.top.deg$col.y))
	abline(h=c(-minlogFC, minlogFC), col="blue")
	legend("topright", legend=names(myCol), fill=myCol, cex=2 ) 
}

cat("All is done", "\n")
