#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
myCaller="VennDiagram"
source ("~/Pipelines/config/DEG.R") # load config
library(VennDiagram)

#my.type='uniq.SGA'
#my.type='exclusive.SGA'
my.type='abnormal.SGA'
#my.type='normal.SGA'
if(my.type=='abnormal.SGA'){
	sampleSubTypes=c("LOW.PAPPA", "LOW.GVEL", "HIGH.UAD", "HIGH.UBD", "HT")
}else if(my.type=='uniq.SGA'){
	sampleSubTypes=c("LOW.PAPPA.UNIQ", "LOW.GVEL.UNIQ", "HIGH.UAD.UNIQ", "HIGH.UBD.UNIQ", "HT.UNIQ")
}else if(my.type=='exclusive.SGA'){
	sampleSubTypes=c("LOW.PAPPA.EXC", "LOW.GVEL.EXC", "HIGH.UAD.EXC", "HIGH.UBD.EXC", "HT.EXC")
}else if(my.type=='normal.SGA'){
	sampleSubTypes=c("NORMAL.PAPPA", "NORMAL.GVEL", "NORMAL.UAD", "NORMAL.UBD", "NHT")
}else{
	stop(paste0(my.type," Not supported"))
}
deseq.result<-list()

############################
# padj (adjusted p-valued) #
############################
deseq.result<-lapply(sampleSubTypes, function(i) read.csv(paste0(deseqHome,"/",i,"/filtered_toptags_deseq.padj.4.",i,".csv"), stringsAsFactors=FALSE)) # read csv
names(deseq.result)<-sampleSubTypes # name each slice
deseq.result.padj4<-lapply(deseq.result, function(i) i<-unique(i[i$padj<=0.4,c(my.filter)])) # get ensg filter by FDR 
deseq.result.padj3<-lapply(deseq.result, function(i) i<-unique(i[i$padj<=0.3,c(my.filter)])) # get ensg filter by FDR 
deseq.result.padj2<-lapply(deseq.result, function(i) i<-unique(i[i$padj<=0.2,c(my.filter)])) # get ensg filter by FDR 
deseq.result.padj1<-lapply(deseq.result, function(i) i<-unique(i[i$padj<=0.1,c(my.filter)])) # get ensg filter by FDR 
deseq.result.padj05<-lapply(deseq.result, function(i) i<-unique(i[i$padj<=0.05,c(my.filter)])) # get ensg filter by FDR 

###########
# Top 100 #
###########
deseq.result<-lapply(sampleSubTypes, function(i) read.csv(paste0(deseqHome,"/",i,"/filtered_toptags_deseq.top100.",i,".csv"), stringsAsFactors=FALSE)) # read csv
names(deseq.result)<-sampleSubTypes # name each slice
deseq.result.top100<-lapply(deseq.result, function(i) i<-unique(i[,c(my.filter)])) # get only ensg

###########
# Top 200 #
###########
deseq.result<-lapply(sampleSubTypes, function(i) read.csv(paste0(deseqHome,"/",i,"/filtered_toptags_deseq.top200.",i,".csv"), stringsAsFactors=FALSE)) # read csv
names(deseq.result)<-sampleSubTypes # name each slice
deseq.result.top200<-lapply(deseq.result, function(i) i<-unique(i[,c(my.filter)])) # get only ensg

##################
## Venn Diagram ##
##################
cat("making Venn Diagram...\n")
###############
# DEG Overlap #
###############
draw.venn(deseq.result.padj4,paste0(deseqHome,"/Venn.",my.type,".padj4.tiff"))
draw.venn(deseq.result.padj3,paste0(deseqHome,"/Venn.",my.type,".padj3.tiff"))
draw.venn(deseq.result.padj2,paste0(deseqHome,"/Venn.",my.type,".padj2.tiff"))
draw.venn(deseq.result.padj1,paste0(deseqHome,"/Venn.",my.type,".padj1.tiff"))
draw.venn(deseq.result.padj05,paste0(deseqHome,"/Venn.",my.type,".padj05.tiff"))
draw.venn(deseq.result.top100,paste0(deseqHome,"/Venn.",my.type,".top100.tiff"))
draw.venn(deseq.result.top200,paste0(deseqHome,"/Venn.",my.type,".top200.tiff"))

if(FALSE){
	##################
	# Sample Overlap #
	##################
	if(my.type!='uniq.SGA'){
		subtype.list<-lapply(sampleSubTypes, selectSamples, samples.all) # selectSamples defined 
		names(subtype.list)<-sampleSubTypes # name each slice
		subtype.list<-lapply(subtype.list, function(i) i<-i[i$Condition==1,c("SampleName")]) # get SGA sample only
		draw.venn(subtype.list,paste0(deg.metaDir,"/Venn.",my.type,".samples.tiff"))
	}
}


##
## ALL vs. HIGH.UBD
##
sampleSubTypes=c("ALL", "HT")
deseq.result<-lapply(sampleSubTypes, function(i) read.csv(paste0(deseqHome,"/",i,"/filtered_toptags_deseq.padj.4.",i,".csv"), stringsAsFactors=FALSE)) # read csv
names(deseq.result)<-sampleSubTypes # name each slice
deseq.result.padj05<-lapply(deseq.result, function(i) i<-unique(i[i$padj<=0.05,c(my.filter)])) # get ensg filter by FDR 

venn.diagram(
	x=deseq.result.padj05,
	filename = paste0(deseqHome,"/Venn.",sampleSubTypes[1],"vs.",sampleSubTypes[2],".samples.tiff"),
	col = "black",
	fill = c("dodgerblue", "goldenrod1"),
	alpha = 0.50,
	cat.col = c("dodgerblue", "goldenrod1"),
	cat.cex = 1.1,
	cat.fontface = "bold",
	margin = 0.05
)

cat("All is done\n")
