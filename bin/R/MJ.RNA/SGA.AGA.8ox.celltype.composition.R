#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

# make sure 
# myProject="SGA.AGA"
# sampleType='OXBS' # MG 4 low PAPP-A SGA vs. 4 matched AGA 
myCaller="cellType"
source ("~/Pipelines/config/DEG.R") # load config

library(reshape) # 'melt'

options(scipen=999) # diable scientific notation

# pdf output filename 
pdf_file_name<-"~/scratch/results/RNA-Seq/SGA.AGA/MetExp/plot"
pdf(file=paste(pdf_file_name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

##########################
## Read MJ RNA-Seq data ##
##########################
#DESeq2 RData for this sampleType
deseq.RData <- paste0(deseq.dir, "/deseq.",sampleType,'.RData')
load(deseq.RData) #dds, ddsFlt, res, resFlt, rld, rld.nb, vsd

############################# 
## Cell-type Heterogeneity ##
############################# 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
cell.types<-read.csv("~/Pipelines/data/cell.type.expression.csv") #read 21 genes of 4 cell types
my.metrics<-c(`RLOG`="Regulaised Log", `NCPM`="Normalised Count Per Million", `CPM`="Count Per Million", `FPKM`="Fragments Per Kilobase per Million mapped fragments")
for(i in names(my.metrics)){
	#CPM-based
	if(i=="CPM"){
		ddsCpms<-as.data.frame(cpm(counts(dds))) 
		dummy<-as.data.frame(ddsCpms[rownames(ddsCpms) %in% cell.types$ensg,]) # isa 'matrix'
		#boxplot(scale(t(dummy))) # scale(t(dummy)) same with t((dummy-rowMeans(dummy))/rowSds(dummy))
	}else if(i=="NCPM"){
		ddsCpms<-as.data.frame(cpm(counts(dds))) 
		dummy<-ddsCpms[rownames(ddsCpms) %in% cell.types$ensg,] # isa 'matrix'
		dummy<-as.data.frame(t(scale(t(dummy))))  # scale CPM by each gene
	}else if(i=="RLOG"){
		rlogMat <- as.data.frame(assay(rld)) 
		dummy<-rlogMat[rownames(rlogMat) %in% cell.types$ensg,] # isa 'matrix'
	}else{
	#FPKM-based
		#annotation db
		my.local.db<-"~/data/Annotation/hg19ensGene.sqlite"
		# load annotation from local file
		hg19ensGene<- loadDb(my.local.db)
		exons.list.per.gene <- exonsBy(hg19ensGene,"gene") # isa "GRangeList"
		keep <- rownames(dds) %in% names(exons.list.per.gene) # remove ensg of no exons info
		newCounts <- counts(dds)[keep,]
		ddsMatrix <- DESeqDataSetFromMatrix(newCounts, colData(ddsCollapsed), deseq.design) # isa 'DESeqDataSet'
		dds.exonic.gene <- DESeq(ddsMatrix) # isa 'DESeqDataSet'
		rowData(dds.exonic.gene) <- exons.list.per.gene # set rowData for the sake of fpkm
		ddsFpkm <-fpkm(dds.exonic.gene) # isa 'matrix'

		dummy<-as.data.frame(ddsFpkm[rownames(ddsFpkm) %in% cell.types$ensg,]) # isa 'matrix'
	}

	# replace ENSG by GeneName
	rownames(dummy)<-sapply(rownames(dummy), function(i) cell.types[cell.types$ensg==i,c("gene")])
	dummy$cellType<-sapply(rownames(dummy), function(i) cell.types[cell.types$gene==i,c("cellType")])
	dummy$gene<-rownames(dummy)
	dummy<-dummy[order(dummy$cellType),]

	# ENSG00000261371(PECMA1) and ENSG00000267631(CGB1) not expressed 
	select<-!(dummy$gene %in% c("PECAM1", "CGB1"))
	dummy<-dummy[select,]

	# melt down columns
	melt.dummy<-melt(dummy,variable_name="sampleId")  # isa 'data.frame', melt from reshape library
	melt.dummy<-merge(melt.dummy, samples[,c("SampleName","Condition")], by.x="sampleId", by.y="SampleName") # add 'Conditioni'
	list.melt.dummy<-split(melt.dummy, melt.dummy$Condition)
	type=c(`0`="AGA",`1`="SGA")

	for(j in names(list.melt.dummy)){
		my.melt.dummy<-list.melt.dummy[[j]]	
		# boxplot vs. gene coloured by cellType
		my.filename=paste0(pdf_file_name,i,".",type[j],".by.gene")
		tiff(filename=paste0(my.filename,".tiff"),width=16,height=9,units="in",res=300, compression = 'lzw')
		print(ggplot(my.melt.dummy, aes(gene,value)) + geom_boxplot(aes(fill=cellType)) + geom_point() + scale_x_discrete(limits=dummy[order(dummy$cellType),c("gene")]) + ylab(my.metrics[i]) + ggtitle(type[j]))
		dev.off()
			
		# boxplot vs. cellType
		my.filename=paste0(pdf_file_name,i,".",type[j],".by.celltype")
		tiff(filename=paste0(my.filename,".tiff"),width=16,height=9,units="in",res=300, compression = 'lzw')
		print(ggplot(my.melt.dummy, aes(cellType,value)) + geom_boxplot(aes(fill=cellType)) + geom_point() + ylab(my.metrics[i]) + ggtitle(type[j]))
		dev.off()

		# boxplot vs. sampleId coloured by cellType
		my.filename=paste0(pdf_file_name,i,".",type[j],".by.sample")
		tiff(filename=paste0(my.filename,".tiff"),width=16,height=9,units="in",res=300, compression = 'lzw')
		print(ggplot(my.melt.dummy, aes(sampleId,value)) + geom_boxplot(aes(fill=factor(cellType))) + ylab(my.metrics[i]) + ggtitle(type[j]) )
		dev.off()
		#ggplot(melt.dummy, aes(cellType,value)) + geom_boxplot(aes(fill=cellType)) + facet_grid(sampleId ~ .)
	}
}

##############################
## Heatmap (gene vs sample) ##
##############################
my.filename=paste0(pdf_file_name,i,".heatmap.gene.sample")
#tiff(filename=paste0(my.filename,".tiff"),width=11,height=9,units="in",res=300, compression = 'lzw')
heatmap.2(as.matrix(dummy[,1:8]), col = hmcol, scale="none", trace="none", margin=c(10, 6))
#dev.off()

########################################
## Heatmap (sample-sample clustering) ##
########################################
sampleDists<- dist(t(dummy[,1:8]))
mat <- as.matrix(sampleDists) # row: samples, col: samples
#colnames(mat) <- rownames(mat)<- paste(rownames(mat), paste(samples$Code,samples$Pair,samples$Sex,sep="."), sep=":")
my.filename=paste0(pdf_file_name,i,".heatmap.sample.sample")
#tiff(filename=paste0(my.filename,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
heatmap.2(mat, trace="none", col = hmcol, margin=c(13, 13))
#dev.off()

#########
## PCA ##
#########
#pca<-prcomp(mat, scale=TRUE)
#pca<-prcomp(dummy[,1:8], scale=TRUE)
mat <- as.matrix(t(dummy[,1:8])) # row: samples, col: samples
#colnames(mat) <- rownames(mat)<- paste(rownames(mat), paste(samples$Code,samples$Pair,samples$Sex,sep="."), sep=":")
pca<-prcomp(t(dummy[,1:8]))
my.filename=paste0(pdf_file_name,i,".pca")
#tiff(filename=paste0(my.filename,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
plot(pca$x, main=paste0("PCA of ",sampleType), cex=0.1)
text(pca$x[,1], pca$x[,2], labels=rownames(pca$x))
#dev.off()

cat("All is done", "\n")
