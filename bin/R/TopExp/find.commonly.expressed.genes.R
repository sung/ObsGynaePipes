#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

myCaller="commonExp"
source ("~/Pipelines/config/DEG.R") # load config
options(scipen=999) # diable scientific notation

# RData file
deseq.RData <- paste0(RDataDir, "/deseq.",sampleType,'.RData')

###################################
## collapse Technical Replicates ##
###################################
if(file.exists(deseq.RData)){
	cat("loading DESeq2 RData...\n")
	load (deseq.RData) # 
	cat("dds and ddsExonicGene loaded\n")
}else{
	if(myProject=="PDN"){collapsed=TRUE}else{collapsed=FALSE}
	initDDS(collapsed)  # this function defined within 'DEG.R'
						# this saves RData
	cat("loading DESeq2 RData...\n")
	load (deseq.RData) # 
	cat("dds and ddsExonicGene loaded\n")
}

if(!exists("rld")){
	rld <- rlog(dds) # isa 'SummarizedExperiment'
}

cat("making ddsFpkm from ddsExonicGene...\n")
ddsFpkm <-fpkm(ddsExonicGene) # isa 'matrix'
ddsFpm <-fpm(ddsExonicGene) # isa 'matrix'
#############################################
# Heatmap of the sample-to-sample distances #
# Sample Clustering Based on rlog           #
#############################################
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
rlogMat <- assay(rld) # isa 'matirx' (row: genes, col: samples)
cat("making sampleDists...\n")
sampleDists<- dist(t(rlogMat)) # isa 'dist'
mat <- as.matrix(sampleDists) # NxN distance matrix (row: samples, col: samples)

my.filename <- file.path(cluster.dir, paste0(sampleType, ".sample.clustering.heatmap"))
tiff(filename=paste0(my.filename,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
cat("making heatmap of sampleDists...\n")
print(heatmap.2(mat, trace="none", col = hmcol, key.xlab="Distance", key.title=""))
dev.off()

#########################################
## PCA Based on Sample-Sample Distance ##
#########################################
select<-!mcols(dds)[["allZero"]] # genes having counts from at least one sample (not allZero)
cat("making pca based on rlogMat...\n")
pca<-prcomp(t(rlogMat[select,]), scale.=TRUE) # input (row:samples, col:genes)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
if(FALSE){
	plot(pca)
	plot(pca$x, main="PCA based on Regularised Log Count (rlog)", cex=0.1,  xlab=paste0("PC1: ",round(percentVar[1]*100,2)," % variance"), ylab=paste0("PC2: ",round(percentVar[2]*100,2)," % variance"))
	text(pca$x[,1],pca$x[,2], labels=rownames(pca$x))
}
#ggplot way
plot <- ggplot(as.data.frame(pca$x), aes(PC1,PC2)) + geom_text(label=rownames(pca$x)) 
plot <- plot + xlab(paste0("PC1 (",round(percentVar[1]*100,2)," % variance)")) + ylab(paste0("PC2 (",round(percentVar[2]*100,2)," % variance)"))
plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)

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
select<- head( order(rowVars(rlogMat), decreasing=TRUE ), n=round(nrow(rlogMat)*0.1*0.01) ) # top 0.1% variable genes (requires 'genefilter' library)
fields <- c("ensembl_gene_id", "hgnc_symbol", "description") # biomart
#fields <- c("ensembl_gene_id", "hgnc_symbol", "description", "superfamily")
#fields <- c("ensembl_gene_id", "hgnc_symbol", "description", "family", "family_description")
anno<-getBM(attributes = fields, filters = "ensembl_gene_id", values = rownames(rlogMat[select,]), mart = grch37)
anno<-anno[order(anno$ensembl_gene_id),] # order by ensembl_gene_id
dummy<-rlogMat[select,]
dummy<-dummy[order(rownames(dummy)),] # order by rownames (ensembl_gene_id)
rownames(dummy) <- ifelse(anno$hgnc_symbol=="" | is.na(anno$hgnc_symbol), anno$ensembl_gene_id, anno$hgnc_symbol) # replace with ensembl id if no gene name

cat("making heatmap of gene clustering...\n")
my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.var.gene.clustering.heatmap"))
tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
print(heatmap.2(dummy, scale="row", trace="none", dendrogram="column", col = hmcol, keysize=1, key.title="", margin=c(5,7), offsetRow=0.1))
dev.off()

write.csv(anno[order(anno$hgnc_symbol),], file=file.path(cluster.dir, paste0(sampleType, ".top.var.gene.csv")), row.names = FALSE)

##########
## FPKM ##
##########
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)

dummy<-getTopExp(ddsFpkm, annotation=TRUE, iqr=FALSE, round(nrow(ddsFpkm)*0.1*0.01)) # isa 'list'

cat("making heatmap of Top FPKM...\n")
my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.fpkm.gene.heatmap"))
tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
print(heatmap.2(as.matrix(dummy$fpkm), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="FPKM"))
dev.off()

cat("making heatmap of Top FPKM in log2 scale...\n")
my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.log2.fpkm.gene.heatmap"))
tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
print(heatmap.2(as.matrix(log2(1+dummy$fpkm)), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="log2(1+FPKM)"))
dev.off()

###############################################################
## Box plot of FPKM of genes having > 1 CPM for all samples  ##
###############################################################
keep = rowSums(ddsFpm >1) >= ncol(ddsExonicGene)

cat("making heatmap of boxplot of FPKM in log2 scale...\n")
my.filename <- file.path(cluster.dir, paste0(sampleType, ".log2.fpkm.boxplot"))
tiff(filename=paste0(my.filename,".tiff"),width=10,height=9,units="in",res=300, compression = 'lzw')
print(boxplot(log2(1+ddsFpkm[keep,]), xlab="Samples", ylab="log2(1+FPKM)", main=paste0("FPKM boxplots of ",nrow(ddsFpkm[keep,]), " genes having > 1 CPM across all samples")))
dev.off()

cat("making heatmap of boxplot of FPKM...\n")
my.filename <- file.path(cluster.dir, paste0(sampleType, ".fpkm.boxplot"))
tiff(filename=paste0(my.filename,".tiff"),width=10,height=9,units="in",res=300, compression = 'lzw')
print(boxplot(ddsFpkm[keep,], xlab="Samples", ylab="FPKM", main=paste0("FPKM boxplots of ",nrow(ddsFpkm[keep,]), " genes having > 1 CPM across all samples")))
dev.off()

##################
## FPKM  of IQR ##
##################
# Top 200 FPKM
dummy<-getTopExp(ddsFpkm[keep,], annotation=TRUE, iqr=TRUE, 200) # top 200

cat("making heatmap of Top FPKM based on IQR mean...\n")
my.filename <- file.path(cluster.dir, paste0(sampleType, ".top.iqr.fpkm.gene.heatmap"))
tiff(filename=paste0(my.filename,".tiff"),width=9,height=12,units="in",res=300, compression = 'lzw')
print(heatmap.2(head(as.matrix(dummy$fpkm),n=round(nrow(ddsFpkm)*0.1*0.01)), Rowv=FALSE, scale="none", trace="none", dendrogram="column", col = hmcol, keysize=1, key.xlab="FPKM"))
dev.off()

write.csv(anno[order(anno$hgnc_symbol),], file=file.path(cluster.dir, paste0(sampleType, ".top.200.iqr.fpkm.gene.csv")), row.names = FALSE)
dummy$fpkm$mean<-rowMeans(dummy$fpkm) # mean over all the samples
dummy$fpkm$meanIQR<-apply(dummy$fpkm, 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]])) # mean over trimmed samples
write.csv(dummy$fpkm, file=file.path(cluster.dir, paste0(sampleType, ".top.200.iqr.fpkm.csv")))

###############################
## SD-Mean plot of IQR(FPKM) ##
###############################
cat("making Sdplot based on ddsFpkm...\n")
mean.iqr.fpkm<-apply(ddsFpkm[keep,], 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]])) # isa numeric vector
sd.iqr.fpkm<-apply(ddsFpkm[keep,], 1, function(i) sd(i[i>=quantile(i)[2] & i<=quantile(i)[4]])) # isa numeric vector

my.filename <- file.path(cluster.dir, paste0(sampleType, ".iqr.fpkm.boxplot"))
tiff(filename=paste0(my.filename,".tiff"),width=9,height=10,units="in",res=300, compression = 'lzw')
boxplot(mean.iqr.fpkm, ylab="mean(IRQ(FPKM))") # box plot of FPKM IQR
text(cbind(1, dummy$fpkm$meanIQR[1:5]), labels=rownames(dummy$fpkm)[1:5], pos=4) # label to the right (pos=4)
dev.off()

my.filename <- file.path(cluster.dir, paste0(sampleType, ".iqr.fpkm.mean.sd"))
tiff(filename=paste0(my.filename,".tiff"),width=9,height=10,units="in",res=300, compression = 'lzw')
plot(mean.iqr.fpkm, sd.iqr.fpkm, xlab="mean(IQR(FPKM))", ylab="SD(IQR(FPKM))" ) # SD-Mean plot
#[todo]
text(cbind(mean.iqr.fpkm, sd.iqr.fpkm)[dummy$select[1:5],], labels=rownames(dummy$fpkm)[1:5], pos=2) # label to the left (pos=2)
dev.off()

my.filename <- file.path(cluster.dir, paste0(sampleType, ".iqr.fpkm.log2mean.log2sd"))
tiff(filename=paste0(my.filename,".tiff"),width=9,height=10,units="in",res=300, compression = 'lzw')
plot(log2(apply(ddsFpkm[keep,], 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]]))), log2(apply(ddsFpkm[keep,], 1, function(i) sd(i[i>=quantile(i)[2] & i<=quantile(i)[4]]))), xlab="log2(mean(IQR(FPKM)))", ylab="log2(SD(IQR(FPKM)))" ) # SD-Mean plot
dev.off()

my.filename <- file.path(cluster.dir, paste0(sampleType, ".iqr.fpkm.log2.one.plus.mean.log2sd"))
tiff(filename=paste0(my.filename,".tiff"),width=9,height=10,units="in",res=300, compression = 'lzw')
plot(log2(1+apply(ddsFpkm[keep,], 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]]))), log2(1+apply(ddsFpkm[keep,], 1, function(i) sd(i[i>=quantile(i)[2] & i<=quantile(i)[4]]))), xlab="log2(1+mean(IQR(FPKM)))", ylab="log2(1+SD(IQR(FPKM)))" ) # SD-Mean plot
dev.off()

##################################
## FPKM from Protein Atlas Data ##
##################################
if(FALSE){
	p.atlas<-read.table("~/data/Annotation/protein.atlas.fpkm.new.header.txt", header=TRUE)
	rownames(p.atlas)<-p.atlas$ensg
	p.atlas$ensg<-NULL

	# split the table by the unique tissue types
	tissue.list<-unique(lapply(strsplit(colnames(p.atlas),"_"), function(i) i[1]))
	p.atlas.fpkm<-lapply(tissue.list, function(i) p.atlas[,grepl(i, colnames(p.atlas))])
	names(p.atlas.fpkm)<-unlist(tissue.list)

	# Top FPKM from placenta 
	#p.atlas.fpkm$placenta
	p.atlas.placenta<-getTopExp(p.atlas.fpkm$placenta, annotation=TRUE, iqr=TRUE, 60)
}
