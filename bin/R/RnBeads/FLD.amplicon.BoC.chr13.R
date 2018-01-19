#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

library(RnBeads)
library(methylKit) # for 'read.feature.flank'
library(reshape2) # for 'melt' and 'dcast'
library(GenomicRanges)
library(gplots) # heatmap.2
library(rtracklayer)
library(ggbio)

############
## Config ##
############
logger.start("Logger started", fname = NA )
source("~/Pipelines/config/RnBeads.sbs.R")
source("~/Pipelines/config/Annotation.R")

my.upflank=2000 # up from TSS
my.downflank=200 # down form TES

my.case="LPA" #"SGA"
my.control="CPA" #"AGA"

#my.fld.ver="FLD.v2"; my.assay<-"BS"
#my.fld.ver="FLD.v3"; my.assay<-"oxBS"
my.fld.ver="FLD.v4"; my.assay<-"oxBS"
my.slx<-c(`FLD.v2`="SLX-8772",`FLD.v3`="SLX-8773", `FLD.v4`="SLX-10295")

#########################
## RnBeads Sample Info ##
#########################
if(!file.exists(my.meta.file)){ 
	stop("my.meta.file not found. Stopped!")
}else{
	my.meta.sample<-read.csv(my.meta.file, stringsAsFactors=FALSE)
	group1<-my.meta.sample[my.meta.sample$Condition==my.control,c("SampleName")] # AGA
	group2<-my.meta.sample[my.meta.sample$Condition==my.case,c("SampleName")] # SGA
	names(group1)<-sapply(my.meta.sample[my.meta.sample$Condition==my.control,c("BedFile")], function(i) unlist(strsplit(unlist(strsplit(i, "/"))[8], "[.]"))[2])
	names(group2)<-sapply(my.meta.sample[my.meta.sample$Condition==my.case,c("BedFile")], function(i) unlist(strsplit(unlist(strsplit(i, "/"))[8], "[.]"))[2])
	my.samples=c(group1,group2)
}

###############################
## Fludigim Meta Sample Info ##
###############################
fld.file<-paste0("~/results/RnBeads/Meta/",my.fld.ver,"/",my.fld.ver,".sample.info.csv") # 141 FLD v2 samples 
if(!file.exists(fld.file)){
	stop("fld.file not found. Stopped!")
}else{
	fld.meta.sample<-read.csv(fld.file)
}
my.fld.meta.sample<-fld.meta.sample[fld.meta.sample$Barcode %in% names(my.samples),]


###########################################
### Fluidigm Methylation from BED Files ###
###########################################
if(!file.exists(preprocessed.rnb.set.dir)){ 
	gr.fld<-lapply(as.character(my.fld.meta.sample$BedFile), rtracklayer::import.bed) # isa 'GRangesList'
	names(gr.fld)<-my.fld.meta.sample$SampleName
	gr.fld.new<-lapply(names(gr.fld), function(i) 
			GRanges(
				seqnames=seqnames(gr.fld[[i]]), 
				ranges<-IRanges(start=start(gr.fld[[i]]),end=end(gr.fld[[i]])),
				strand=strand(gr.fld[[i]]), 
				mcols<-DataFrame(with(mcols(gr.fld[[i]]), data.frame(SampleName=i, Condition=ifelse(grepl(my.control,i), my.control, my.case), Sex=ifelse(grepl("M",i), "Boy", "Girl"), score=round(as.numeric(name)/100,2), coverage=score)))
				)
			)
	gr.fld.merged<-do.call(c, gr.fld.new) # isa 'GRanges'
	gr.fld.merged<- gr.fld.merged[mcols(gr.fld.merged)$coverage>=min.doc,] # min 10 depth  
	#fld.df<-dcast(as.data.frame(gr.fld.merged), start ~ SampleName, value.var="score"); rownames(fld.df)<-fld.df$start; fld.df$start<-NULL; 
	#gr.fld.merged<-gr.fld.merged[start(gr.fld.merged) %in% rownames(fld.df[rowSums(!is.na(fld.df))==ncol(fld.df),]),]
	fld.df<-dcast(with(as.data.frame(gr.fld.merged), data.frame(id=paste(seqnames,start,sep="."),SampleName, score)), id ~ SampleName, value.var="score"); rownames(fld.df)<-fld.df$id; fld.df$id<-NULL;
	gr.fld.merged<-gr.fld.merged[paste(seqnames(gr.fld.merged), start(gr.fld.merged), sep=".") %in% rownames(fld.df[rowSums(!is.na(fld.df))>=ncol(fld.df)-4,]),]
}else{
#########################################
### Fluidigm Methylation from RnBeads ###
#########################################
	# 1. load RnBSet object
	cat(paste0("Importing rnb.set (preprocessed RnBeads data) from ",preprocessed.rnb.set.dir,"...\n"))
	rnb.set<-load.rnb.set(path=preprocessed.rnb.set.dir) #isa 'RnBSet'
	gr.fld<-rnb.RnBSet.to.GRangesList(rnb.set) # isa 'GRangesList'

	# Adjust start and end 
	# RnBeads write CG sites (whereas 'C' only required)
	# Add 'Condition' and 'SampleName' 
	gr.fld.new<-lapply(names(gr.fld), function(i) 
			GRanges(
				seqnames=seqnames(gr.fld[[i]]), 
				ranges<-IRanges(
					start=ifelse(strand(gr.fld[[i]])=="+", start(gr.fld[[i]]), end(gr.fld[[i]])),
					end=ifelse(strand(gr.fld[[i]])=="+", start(gr.fld[[i]]), end(gr.fld[[i]]))
				), 
				strand=strand(gr.fld[[i]]), 
				mcols<-DataFrame(with(mcols(gr.fld[[i]]), data.frame(SampleName=i, Condition=ifelse(grepl(my.control,i), my.control, my.case), Sex=ifelse(grepl("M",i), "Boy", "Girl"), score, coverage)))
			)
		)
	gr.fld.merged<-do.call(c, gr.fld.new) # isa 'GRanges'
}

#######################################
### WGBS Methylation from BED Files ###
#######################################
wgbs.file<-c(`MAB21L1`="~/Pipelines/data/MAB21L1/MAB21L1.top2.grch37.WGBS.bed.txt", 
			`chr7.DMR`="~/Pipelines/data/Chr7.DMR/chr7.top2.grch37.WGBS.bed.txt")
dummy<-lapply(wgbs.file, read.csv, stringsAsFactors=FALSE) # isa 'list'

gr.chr13<-lapply(dummy[["MAB21L1"]]$BedFile, rtracklayer::import.bed) 
names(gr.chr13)<-dummy[["MAB21L1"]]$SampleName
gr.chr13.new<-lapply(names(gr.chr13), function(i) 
		GRanges(
			seqnames=seqnames(gr.chr13[[i]]), 
			ranges<-IRanges(start=start(gr.chr13[[i]]),end=end(gr.chr13[[i]])),
			strand=strand(gr.chr13[[i]]), 
			mcols<-DataFrame(with(mcols(gr.chr13[[i]]), data.frame(SampleName=i, Condition=ifelse(grepl("AGA",i), "AGA", "SGA"), Sex=ifelse(grepl("M",i), "Boy", "Girl"), score=round(as.numeric(name)/100,2), coverage=score)))
			)
		)
gr.chr13.merged<-do.call(c, gr.chr13.new) # isa 'GRanges'
gr.chr13.merged<- gr.chr13.merged[mcols(gr.chr13.merged)$coverage>=min.doc,] # min 10 depth  
#CpG called by all samples (dcast from 'reshape2')
wgbs.chr13<-dcast(as.data.frame(gr.chr13.merged), start ~ SampleName, value.var="score"); rownames(wgbs.chr13)<-wgbs.chr13$start; wgbs.chr13$start<-NULL; 
gr.chr13.merged<-gr.chr13.merged[start(gr.chr13.merged) %in% rownames(wgbs.chr13[rowSums(!is.na(wgbs.chr13))==ncol(wgbs.chr13),]),]

gr.chr7<-lapply(dummy[["chr7.DMR"]]$BedFile, rtracklayer::import.bed) 
names(gr.chr7)<-dummy[["chr7.DMR"]]$SampleName
gr.chr7.new<-lapply(names(gr.chr7), function(i) 
		GRanges(
			seqnames=seqnames(gr.chr7[[i]]), 
			ranges<-IRanges(start=start(gr.chr7[[i]]),end=end(gr.chr7[[i]])),
			strand=strand(gr.chr7[[i]]), 
			mcols<-DataFrame(with(mcols(gr.chr7[[i]]), data.frame(SampleName=i, Condition=ifelse(grepl("AGA",i), "AGA", "SGA"), Sex=ifelse(grepl("M",i), "Boy", "Girl"), score=round(as.numeric(name)/100,2), coverage=score)))
			)
		)
gr.chr7.merged<-do.call(c, gr.chr7.new) # isa 'GRanges'
gr.chr7.merged<- gr.chr7.merged[mcols(gr.chr7.merged)$coverage>=min.doc,] # min 10 depth  
#CpG called by all samples (dcast from 'reshape2')
wgbs.chr7<-dcast(as.data.frame(gr.chr7.merged), start ~ SampleName, value.var="score"); rownames(wgbs.chr7)<-wgbs.chr7$start; wgbs.chr7$start<-NULL; 
gr.chr7.merged<-gr.chr7.merged[start(gr.chr7.merged) %in% rownames(wgbs.chr7[rowSums(!is.na(wgbs.chr7))==ncol(wgbs.chr7),]),]

##############################
## Amplicon Coveage Heatmap ##
##############################
# 1. Breadth of Coverage Per Amplicon
#SLX=SLX-10295; VERSION=Homo_sapiens.v2
#echo "Barcode Amplicon Base Length Cov" > ~/results/$SLX.$VERSION/Coverage/FLD.BoC.per.amplicon.txt
#coverageFiles=$(ls ~/results/$SLX.$VERSION/Coverage/FLD*/$SLX.FLD*.r_1_val_1.fq.gz_bismark_bt2_pe.q1.sorted.bam.CoverageHist.per.amplicon)
#amplicons=$(awk '{print $4}' ~/Pipelines/data/MJ.FLD/$SLX/mj.amplicon.coordinates.sorted.bed)
#for i in $coverageFiles;do barcode=`echo $i| cut -d'/' -f7`;for k in $amplicons; do printf "$barcode $k "; awk "BEGIN{sum=0}\$4==\"$k\"{len=\$9}\$4==\"$k\" && \$7>=10{sum+=\$8}END{print sum,len,sum/len}" $i; done; done >> ~/results/$SLX.$VERSION/Coverage/FLD.BoC.per.amplicon.txt

# 2. Depth of Coverage Per Amplicon
#echo "Barcode Amplicon Base Length Cov" > ~/results/$SLX.$VERSION/Coverage/FLD.DoC.per.amplicon.txt
#for i in $coverageFiles;do barcode=`echo $i| cut -d'/' -f7`;for k in $amplicons; do printf "$barcode $k "; awk "BEGIN{sum=0}\$4==\"$k\"{len=\$9; sum+=\$7*\$8}END{print sum,len,sum/len}" $i; done; done >> ~/results/$SLX.$VERSION/Coverage/FLD.DoC.per.amplicon.txt

########################
## Breadh of Coverage ##
########################
my.boc.file=file.path("~/results",paste0(my.slx[my.fld.ver], ".Homo_sapiens.v2/Coverage/FLD.BoC.per.amplicon.txt"))
if(file.exists(my.boc.file)){
	fld.boc<-read.delim(my.boc.file, sep=" ")
}else{
	stop(paste0(fld.boc.file, " not found\n"))
}
BoC<-dcast(fld.boc, Barcode ~ Amplicon, value.var="Cov") # from 'reshape2'
rownames(BoC)<-BoC$Barcode
BoC$Barcode<-NULL # remove 'Barcode' column

# BoC Heatmap (for all samples)
#filter<-rownames(BoC) %in% fld.boc[fld.boc$assay==my.assay,c("Barcode")]
filter<-TRUE
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".amplicon.boc.heatmap")
tiff(filename=paste0(out.file.name,".tiff"),width=11.7,height=8.27,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("Avg. Breadth of Coverage Per Amplicons (>=10x ",my.assay,")")
heatmap.2(as.matrix(BoC[filter,])*100, col=hmcol, trace="none", dendrogram="none", cexCol=0.7, cexRow=0.5, key.xlab="Breadth of Coverage (%)", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, Rowv=F, Colv=F, main=my.main) # cexCol to adjust column labels
dev.off()

# for ggbio later
avg.boc<-rowMeans(t(BoC)) # mean coverage (across samples)
avg.boc.q1<-rowMeans(t(BoC[rownames(BoC) %in% names(my.samples),]))

######################
## Depth of Coveage ##
######################
my.doc.file=file.path("~/results",paste0(my.slx[my.fld.ver], ".Homo_sapiens.v2/Coverage/FLD.DoC.per.amplicon.txt"))
if(file.exists(my.doc.file)){
	fld.doc<-read.delim(my.doc.file, sep=" ")
}else{
	stop(paste0(fld.doc.file, " not found\n"))
}
DoC<-dcast(fld.doc, Barcode ~ Amplicon, value.var="Cov") # from 'reshape2'
rownames(DoC)<-DoC$Barcode
DoC$Barcode<-NULL # remove Barcode

# DoC Heatmap (for all samples)
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".amplicon.doc.heatmap")
tiff(filename=paste0(out.file.name,".tiff"),width=11.7,height=8.27,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("Avg. Depth of Coverage Per Amplicons (",my.assay,")")
heatmap.2(as.matrix(DoC[filter,]), col=hmcol, trace="none", dendrogram="none", cexCol=0.7, cexRow=0.5, key.xlab="Depth of Coverage", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, Rowv=F, Colv=F, main=my.main) # cexCol to adjust column labels
dev.off()

# log10(DoC) Heatmap (for all samples)
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".amplicon.log10.doc.heatmap")
tiff(filename=paste0(out.file.name,".tiff"),width=11.7,height=8.27,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("Avg. Log10(Depth of Coverage) (",my.assay,")")
heatmap.2(as.matrix(log10(DoC[filter,]+1)), col=hmcol, trace="none", dendrogram="none", cexCol=0.7, cexRow=0.5, key.xlab="log10(Doc)", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, Rowv=F, Colv=F, main=my.main) # cexCol to adjust column labels
dev.off()

# for ggbio later
avg.doc<-rowMeans(t(DoC)) # mean coverage (across samples)
avg.doc.q1<-rowMeans(t(DoC[rownames(DoC) %in% names(my.samples),]))


#######################
### Count TS-primer ###
#######################
#SLX=SLX-10295; VERSION=Homo_sapiens.v2
#BARCODES=$(awk '$6=="q3"{print $2}' $HOME/Pipelines/data/MJ.FLD/$SLX/RnBeads.sample.sheet.txt)
#for i in $BARCODES;do myFile=~/results/$SLX.$VERSION/Trim/$i/$SLX.$i.000000000-AH4LJ.s_1.trimming_report.txt; awk 'BEGIN{RS="=== First read: Adapter"}NR>1{print $1,$11}' $myFile > dummy1.txt; awk 'BEGIN{RS="=== Second read: Adapter"}NR>1{print $11}' $myFile > dummy2.txt; paste dummy1.txt dummy2.txt > ~/results/$SLX.$VERSION/Trim/$i/$SLX.$i.trim.cnt.txt; rm dummy1.txt dummy2.txt; done

#echo "Barcode Amplicon Read1 Read2" > ~/results/$SLX.$VERSION/Trim/FLD.amplicon.trim.cnt.txt
#trimCountFiles=$(ls ~/results/$SLX.$VERSION/Trim/FLD*/$SLX.FLD*.trim.cnt.txt)
#amplicons=$(awk '{print $4}' ~/Pipelines/data/MJ.FLD/$SLX/mj.amplicon.coordinates.sorted.bed)
#for i in $trimCountFiles;do barcode=`echo $i| cut -d'/' -f7`;for k in $amplicons; do printf "$barcode $k "; awk "/$k/{print \$2,\$3}" $i; done; done >> ~/results/$SLX.$VERSION/Trim/FLD.amplicon.trim.cnt.txt
myFile=file.path("~/results",paste0(my.slx[my.fld.ver], ".Homo_sapiens.v2/Trim/FLD.amplicon.trim.cnt.txt"))
if(file.exists(myFile)){
	fld.trim<-read.delim(myFile, sep=" ")
}else{
	stop(paste0(myFile, " not found\n"))
}
trimRead1<-dcast(fld.trim, Barcode ~ Amplicon, value.var="Read1") # from 'reshape2'
rownames(trimRead1)<-trimRead1$Barcode
trimRead1$Barcode<-NULL # remove Barcode

trimRead2<-dcast(fld.trim, Barcode ~ Amplicon, value.var="Read2") # from 'reshape2'
rownames(trimRead2)<-trimRead2$Barcode
trimRead2$Barcode<-NULL # remove Barcode

out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".read.cnt.by.amplicon.primer1.heatmap")
tiff(filename=paste0(out.file.name,".tiff"),width=11.7,height=8.27,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("NO. of Read with the Target-Specific Primer1 (",my.assay,")")
heatmap.2(as.matrix(trimRead1), col=hmcol, trace="none", dendrogram="none", cexCol=0.7, cexRow=0.5, key.xlab="Read Count", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, Rowv=F, Colv=F, main=my.main) # cexCol to adjust column labels
dev.off()

out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".read.cnt.by.amplicon.primer2.heatmap")
tiff(filename=paste0(out.file.name,".tiff"),width=11.7,height=8.27,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("NO. of Read with the Target-Specific Primer2 (",my.assay,")")
heatmap.2(as.matrix(trimRead2), col=hmcol, trace="none", dendrogram="none", cexCol=0.7, cexRow=0.5, key.xlab="Read Count", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, Rowv=F, Colv=F, main=my.main) # cexCol to adjust column labels
dev.off()

#############################
## PCA based on BoC or DoC ##
#############################
cat("making pca based on coverage...\n")
## Based on Breadth of Coverage
select<-!rownames(BoC) %in% c("FLD0008", "FLD0035") # these failed on chr7
pca<-prcomp(BoC[select,], scale.=F) # input (row:samples, col:amplicons)
#pca<-prcomp(BoC, scale.=T) # it failes
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

dummy<-cbind(pca$x[,1:3], fld.meta.sample[fld.meta.sample$Barcode %in% rownames(BoC[select,]),],`breadth`=rowMeans(BoC[select,]))
pca.plot <- ggplot(dummy, aes(PC1,PC2)) +
		geom_text(aes(label=rownames(pca$x),size=breadth)) + scale_size(range=c(3,10)) +
		#geom_point(aes(color=PCR.batch),size=7,alpha=0.5) + 
		xlab(paste0("PC1 (",round(percentVar[1]*100,2)," % variance)")) + ylab(paste0("PC2 (",round(percentVar[2]*100,2)," % variance)")) + 
		geom_hline(aes(0), size=.2) + 
		geom_vline(aes(0), size=.2)

cat("making pca plot...\n")
out.file.name<-"~/results/RnBeads/FLD.v2.pca.by.breadh.coverage"
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".pca.by.breadth.coverage")
tiff(filename=paste0(out.file.name,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
print(pca.plot)
dev.off()

## Based on Depth of Coverage
select<-!rownames(DoC) %in% c("FLD0008", "FLD0035") # these failed on chr7
pca<-prcomp(DoC[select,], scale.=F) # input (row:samples, col:amplicons)
#pca<-prcomp(DoC[select,], scale.=T) # input (row:samples, col:amplicons)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

dummy<-cbind(pca$x[,1:3], fld.meta.sample[fld.meta.sample$Barcode %in% rownames(DoC[select,]),],`depth`=rowMeans(DoC[select,]))
pca.plot <- ggplot(dummy, aes(PC1,PC2)) +
		geom_point(aes(color=coll.batch,size=depth),alpha=0.5) + scale_size(range=c(3,10)) +
		#geom_text(aes(label=rownames(pca$x),color=factor(PCR.batch),size=depth)) + scale_size(range=c(3,10)) +
		xlab(paste0("PC1 (",round(percentVar[1]*100,2)," % variance)")) + ylab(paste0("PC2 (",round(percentVar[2]*100,2)," % variance)")) + 
		geom_hline(aes(0), size=.2) + 
		geom_vline(aes(0), size=.2)

cat("making pca plot...\n")
#out.file.name<-"~/results/RnBeads/FLD.v2.pca.by.depth.coverage.no.scale"
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".pca.by.depth.coverage.with.scale.dots")
tiff(filename=paste0(out.file.name,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
print(pca.plot)
dev.off()

#############################
## PCA based on Methylation #
#############################
cat("making pca based on % methylation...\n")

############################
## Amplicon Coveage GGBIO ##
############################
# genomic region of interesg (target)
my.target<-'MAB21L1' # 'chr7.DMR' #'MAB21L1' #NKX1-2
my.chr<-'chr13' #'chr7' #'chr13'
if(my.target=='chr7.DMR'){
	my.contrast<-"Sex"
	my.col<-my.col.Sex
	my.gr.ensg <- rtracklayer::import.bed("~/Pipelines/data/Chr7.DMR/chr7.top2.bed")
	my.gr.ensg <- reduce(my.gr.ensg)
	names(my.gr.ensg)<-my.target
}else{
	my.contrast<-"Condition"
	my.col=c("#00BFC4", "#F8766D"); names(my.col)<-c(my.control, my.case)
	my.ensg=getBM(attributes = c("ensembl_gene_id"), filters = "hgnc_symbol", values=my.target, mart = grch37) # is a data.frame
	my.gr.ensg<-gr.ensg[mcols(gr.ensg)$gene_id %in% my.ensg] # GRanges of my.target gene 
}
my.gr.ext.ensg<-promoters(my.gr.ensg, upstream=my.upflank, downstream=end(my.gr.ensg)-start(my.gr.ensg)+1+my.downflank) # includes 2000 +TSS and 100 +TES 

# overlap between cpgi and my.target 
cpg.obj=read.feature.flank(location=my.cgi.file, feature.flank.name=c("CpGi","shores")) # a ‘GenomicRangesList’ contatining one GRanges object for flanks and one for GRanges object for the main feature.
my.overlap=as.matrix( findOverlaps(cpg.obj[["CpGi"]], my.gr.ext.ensg) )
my.cpgi<-cpg.obj$CpGi[my.overlap[,"queryHits"]]

# Amplicons 
my.AmpBedFile=file.path("~/Pipelines/data/MJ.FLD",my.fld.ver,"mj.amplicon.coordinates.sorted.bed") # FLD amplicon coordinate
if(!file.exists(my.AmpBedFile)){
	stop(paste0(my.AmpBedFile, " not found\n"))
}
#foo<-read.delim(my.AmpBedFile, header=F)
#colnames(foo)<-c('chr','start','end','id','score','strand')
#gr.target <- with(foo, GRanges(chr, IRanges(start, end), strand, score, id=id))
gr.target.all <- rtracklayer::import.bed(my.AmpBedFile)
gr.target <- gr.target.all[as.character(seqnames(gr.target.all))==my.chr,] # my.chr only
names(gr.target)<-mcols(gr.target)$name # assign amplicon name

mcols(gr.target)<-DataFrame(with(mcols(gr.target), 
					data.frame(
						name=name, 
						boc=avg.boc[mcols(gr.target)$name]*100, 
						boc.q1=avg.boc.q1[mcols(gr.target)$name]*100, 
						doc=avg.doc[mcols(gr.target)$name], 
						doc.q1=avg.doc.q1[mcols(gr.target)$name], 
						order=seq(5,200,by=3)[1:nrow(mcols(gr.target))] # amplicon left-most order 
					)
				)) 

# gr for forward/reverse strand separately
gr.target.fwd <- gr.target.all[as.character(seqnames(gr.target.all))==my.chr & strand(gr.target.all)=="+",]
mcols(gr.target.fwd)<-DataFrame(with(mcols(gr.target.fwd), data.frame(name=name, order=seq(5,200,by=3)[1:nrow(mcols(gr.target.fwd))]))) 
names(gr.target.fwd)<-mcols(gr.target.fwd)$name # assign amplicon name

gr.target.rev<- gr.target.all[as.character(seqnames(gr.target.all))==my.chr & strand(gr.target.all)=="-",]
mcols(gr.target.rev)<-DataFrame(with(mcols(gr.target.rev), data.frame(name=name, order=seq(5,200,by=3)[1:nrow(mcols(gr.target.rev))]))) 
names(gr.target.rev)<-mcols(gr.target.rev)$name# assign amplicon name

# extends +-100 from start/end of the target region
gr.target.ext<-reduce(gr.target + 100) 
gr.target.ext<-gr.target.ext[strand(gr.target.ext)=="+",]

##############
# chromosome #
##############
p0 <- Ideogram(genome = "hg19", xlabel=TRUE) + xlim(ranges(gr.target.ext))
####################
# transcript track # 
####################
p1<-autoplot(hg19ensGene, which=my.gr.ext.ensg, fill = "orange", color = "grey") 
############################
# cpgi track for my.target #
############################
p2<-autoplot(my.cpgi, fill = "darkgreen", color = "grey") 
###########################################
# amplicon breadth of coveage (geom_rect) #
###########################################
p3<-ggplot(gr.target) +  
	geom_rect(stat="identity", rect.height = 1, aes(fill=boc, y=order), col="black", ylab=NULL) +  
	scale_fill_gradient2(name="coveage(%)",low=hmcol[1], mid=hmcol[round(length(hmcol)/2)], high=hmcol[length(hmcol)], midpoint=50) + 
	theme(axis.ticks.y = element_blank(), axis.text.y = element_blank() ) +
	geom_text(aes(x=end(gr.target)+10, y=order), label=names(gr.target), size=2.6)
#########################################
# amplicon depth of coveage (geom_rect) #
#########################################
p8<-ggplot(gr.target) +  
	geom_rect(stat="identity", rect.height = 1, aes(fill=doc, y=order), col="black", ylab=NULL) +  
	scale_fill_gradient2(name="depth",low=hmcol.GnBu[1], mid=hmcol.GnBu[round(length(hmcol)/2)], high=hmcol.GnBu[length(hmcol.GnBu)], midpoint=median(avg.doc)/2) + 
	theme(axis.ticks.y = element_blank(), axis.text.y = element_blank() ) +
	geom_text(aes(x=end(gr.target)+10, y=order), label=names(gr.target), size=2.6)
# amplicon relative coverage (based on amplion target bed file) 
p4<-ggplot(gr.target)+stat_coverage()
# amplicon loci and ther names (geom_arrow)
p5<-ggplot(gr.target) + 
	geom_arrow(stat="identity", aes(y=order), col="black", ylab=NULL) +  
	theme(axis.ticks.y = element_blank(), axis.text.y = element_blank() ) +
	geom_text(aes(x=end(gr.target)+10, y=order), label=names(gr.target), size=3)
# amplicon loci and ther names colored by strand 
p7<-ggplot(gr.target) + 
	geom_rect(stat="identity", rect.height = 1, aes(fill=strand, y=order), col="black", ylab=NULL) +  
	theme(axis.ticks.y = element_blank(), axis.text.y = element_blank() ) +
	geom_text(aes(x=end(gr.target)+10, y=order), label=names(gr.target), size=2.6)
p7.fwd<-ggplot(gr.target.fwd) + 
	geom_rect(stat="identity", rect.height = 1, aes(fill=strand, y=order), col="black", ylab=NULL) +  
	theme(axis.ticks.y = element_blank(), axis.text.y = element_blank() ) +
	geom_text(aes(x=end(gr.target.fwd)+10, y=order), label=names(gr.target.fwd), size=2.6)
p7.rev<-ggplot(gr.target.rev) + 
	geom_rect(stat="identity", rect.height = 1, aes(fill=strand, y=order), col="black", ylab=NULL) +  
	theme(axis.ticks.y = element_blank(), axis.text.y = element_blank() ) +
	geom_text(aes(x=end(gr.target.rev)+10, y=order), label=names(gr.target.rev), size=2.6)
##########################
# chr13 WGBS methylation #
##########################
p9<-ggplot(gr.chr13.merged, aes(x=start(gr.chr13.merged),y=score,col=Condition)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col.Condition)
p9.fwd<-ggplot(gr.chr13.merged[strand(gr.chr13.merged)=="+",], aes(x=start(gr.chr13.merged[strand(gr.chr13.merged)=="+",]),y=score,col=Condition)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col.Condition)
p9.rev<-ggplot(gr.chr13.merged[strand(gr.chr13.merged)=="-",], aes(x=start(gr.chr13.merged[strand(gr.chr13.merged)=="-",]),y=score,col=Condition)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col.Condition)

#########################
# chr7 WGBS methylation #
#########################
p10<-ggplot(gr.chr7.merged, aes(x=start(gr.chr7.merged),y=score,col=Sex)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)
p10.fwd<-ggplot(gr.chr7.merged[strand(gr.chr7.merged)=="+",], aes(x=start(gr.chr7.merged[strand(gr.chr7.merged)=="+",]),y=score,col=Sex)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)
p10.rev<-ggplot(gr.chr7.merged[strand(gr.chr7.merged)=="-",], aes(x=start(gr.chr7.merged[strand(gr.chr7.merged)=="-",]),y=score,col=Sex)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)

select<-mcols(gr.chr7.merged)$Condition=="AGA"; my.gr.chr7<-as.data.frame(gr.chr7.merged[select,])
p10.aga<-ggplot(my.gr.chr7, aes(x=start,y=score,col=Sex)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)
select<-mcols(gr.chr7.merged)$Condition=="SGA"; my.gr.chr7<-as.data.frame(gr.chr7.merged[select,])
p10.sga<-ggplot(my.gr.chr7, aes(x=start,y=score,col=Sex)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)

select<-mcols(gr.chr7.merged)$Condition=="AGA"  & strand(gr.chr7.merged)=="+"; my.gr.chr7<-as.data.frame(gr.chr7.merged[select,])
p10.aga.fwd<-ggplot(my.gr.chr7, aes(x=start,y=score,col=Sex)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)
select<-mcols(gr.chr7.merged)$Condition=="AGA"  & strand(gr.chr7.merged)=="-"; my.gr.chr7<-as.data.frame(gr.chr7.merged[select,])
p10.aga.rev<-ggplot(my.gr.chr7, aes(x=start,y=score,col=Sex)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)

select<-mcols(gr.chr7.merged)$Condition=="SGA"  & strand(gr.chr7.merged)=="+"; my.gr.chr7<-as.data.frame(gr.chr7.merged[select,])
p10.sga.fwd<-ggplot(my.gr.chr7, aes(x=start,y=score,col=Sex)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)
select<-mcols(gr.chr7.merged)$Condition=="SGA"  & strand(gr.chr7.merged)=="-"; my.gr.chr7<-as.data.frame(gr.chr7.merged[select,])
p10.sga.rev<-ggplot(my.gr.chr7, aes(x=start,y=score,col=Sex)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)

########################
# Fluidigm methylation #
########################
select<-seqnames(gr.fld.merged)==my.chr; my.gr.fld<-as.data.frame(gr.fld.merged[select,])
p11<-ggplot(my.gr.fld, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)

select<-seqnames(gr.fld.merged)==my.chr & strand(gr.fld.merged)=="+"; my.gr.fld<-as.data.frame(gr.fld.merged[select,])
p11.fwd<-ggplot(my.gr.fld, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)

select<-seqnames(gr.fld.merged)==my.chr & strand(gr.fld.merged)=="-"; my.gr.fld<-as.data.frame(gr.fld.merged[select,])
p11.rev<-ggplot(my.gr.fld, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)

select<-seqnames(gr.fld.merged)==my.chr & mcols(gr.fld.merged)$Condition==my.control; my.gr.fld<-as.data.frame(gr.fld.merged[select,])
p11.ctrl<-ggplot(my.gr.fld, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)

select<-seqnames(gr.fld.merged)==my.chr & mcols(gr.fld.merged)$Condition==my.case; my.gr.fld<-as.data.frame(gr.fld.merged[select,])
p11.case<-ggplot(my.gr.fld, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)

select<-seqnames(gr.fld.merged)==my.chr & mcols(gr.fld.merged)$Condition==my.control & strand(gr.fld.merged)=="+"; my.gr.fld<-as.data.frame(gr.fld.merged[select,])
p11.ctrl.fwd<-ggplot(my.gr.fld, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)
select<-seqnames(gr.fld.merged)==my.chr & mcols(gr.fld.merged)$Condition==my.control & strand(gr.fld.merged)=="-"; my.gr.fld<-as.data.frame(gr.fld.merged[select,])
p11.ctrl.rev<-ggplot(my.gr.fld, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)

select<-seqnames(gr.fld.merged)==my.chr & mcols(gr.fld.merged)$Condition==my.case & strand(gr.fld.merged)=="+"; my.gr.fld<-as.data.frame(gr.fld.merged[select,])
p11.case.fwd<-ggplot(my.gr.fld, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)
select<-seqnames(gr.fld.merged)==my.chr & mcols(gr.fld.merged)$Condition==my.case & strand(gr.fld.merged)=="-"; my.gr.fld<-as.data.frame(gr.fld.merged[select,])
p11.case.rev<-ggplot(my.gr.fld, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 1) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)

#####################################################################################################################################################################################
my.gr.fld<-as.data.frame(gr.fld.merged[seqnames(gr.fld.merged)=="chr7",])
# colour methylation by Sample 
out.file.name<-file.path("~/results/RnBeads",my.project,my.cpg.context,paste("fld.high.met",my.target,my.assay,sep="."))
tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
ggplot(my.gr.fld[my.gr.fld$Sex=="Girl" & my.gr.fld$Condition=="SGA",], aes(x=start,y=score,col=SampleName)) + geom_point(alpha=0.5, size=3) + ylim(0, 1) + geom_smooth(se=FALSE,)
dev.off()
# colour methylation by technical and biological factors for chr7
ggplot(merge(my.gr.fld, my.fld.meta.sample), aes(x=start,y=score, col=coll.batch)) + geom_point(alpha=0.5,size=3) + ylim(0, 1) + geom_smooth(se=FALSE,)
ggplot(merge(my.gr.fld, my.fld.meta.sample), aes(x=start,y=score, col=factor(tissue.colleciton.time))) + geom_point(alpha=0.5,size=3) + ylim(0, 1) + geom_smooth(se=FALSE,)
#####################################################################################################################################################################################

# 1. Amplicon Breadth and Depth of Coverage
out.file.name<-file.path("~/results/RnBeads",my.project,my.cpg.context,paste("amplicon.boc.doc",my.target,my.assay,"ggbio",sep="."))
tiff(filename=paste0(out.file.name,".tiff"),width=8.27, height=11.7 ,units="in",res=300, compression = 'lzw') #A4 size 
if(my.target=='chr7.DMR'){
	tracks(CPGi= p2, `Avg. Breadth of Coverage`=p3, `Avg. Depth of Coverage`=p8, Amplicons=p7, heights = c(1,5,5,5)) + xlim(ranges(gr.target.ext))  + theme_alignment()
}else{
	tracks(Transcripts = p1, CPGi= p2, `Avg. Breadth of Coverage`=p3, `Avg. Depth of Coverage`=p8, `Amplicons`=p7, heights = c(1,0.5,5,5,5)) + xlim(ranges(gr.target.ext))  + theme_alignment() # grid=FALSE, border = T) # plot ggbio
	#tracks(Transcripts = p1, CPGi= p2, `Avg. Breadth of Coverage`=p3, `Avg. Depth of Coverage`=p8, `Amplicons`=p7, heights = c(1,0.5,5,5,5))  # this will break coordinate - always use xlim(range(gr))
} 
dev.off()

# 4. WGBS & Fludigm methylation profile over the amplicon regions
if(my.target=='chr7.DMR'){
	###################
	# Control Samples #
	###################
	#tracks(Transcripts = p1, p2, `WGBS (AGA)`=p10.aga, `Fluidigm (AGA)`=p11.ctrl, Amplicon=p7, heights = c(1.3,0.5,3,3,3)) + xlim(ranges(my.gr.ext.ensg))
	# zoom-in amplicon regions
	out.file.name<-file.path("~/results/RnBeads",my.project,my.cpg.context,paste("wgbs.fld.by.amplicon.AGA",my.target,my.assay,"zoom.ggbio",sep="."))
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
	tracks(p2, `WGBS (AGA)`=p10.aga, `Fluidigm (AGA)`=p11.ctrl, Amplicon=p7, heights = c(0.5,3,3,3)) + xlim(ranges(gr.target.ext))
	dev.off()

	################
	# Case samples #
	################
	#tracks(Transcripts = p1, p2, `WGBS (SGA)`=p10.sga, `Fluidigm (SGA)`=p11.case, Amplicon=p7, heights = c(1.3,0.5,3,3,3)) + xlim(ranges(my.gr.ext.ensg))
	# zoom-in amplicon regions
	out.file.name<-file.path("~/results/RnBeads",my.project,my.cpg.context,paste("wgbs.fld.by.amplicon.SGA",my.target,my.assay,"zoom.ggbio",sep="."))
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
	tracks(p2, `WGBS (SGA)`=p10.sga, `Fluidigm (SGA)`=p11.case, `Avg. Breadth of Coverage`=p3, heights = c(0.5,3,3,3)) + xlim(ranges(gr.target.ext))
	dev.off()

	#################
	# Fluidigm Set3 #
	#################
	out.file.name<-file.path("~/results/RnBeads",my.project,my.cpg.context,paste("wgbs.fld.by.amplicon",my.target,my.assay,"zoom.ggbio",sep="."))
	tiff(filename=paste0(out.file.name,".tiff"),width=8.27, height=11.7,units="in",res=300, compression = 'lzw') #A4 size 
	tracks(p2, `WGBS (AGA)`=p10.aga, `Fluidigm (Set3)`=p11, `Fluidigm (Set3-controls)`=p11.ctrl, `Fluidigm (Set3-case)`=p11.case, Amplicon=p7, heights = c(0.5,3,3,3,3,3)) + xlim(ranges(gr.target.ext))
	dev.off()
}else{
	# extended transcript region
	out.file.name<-file.path("~/results/RnBeads",my.project,my.cpg.context,paste("wgbs.fld.by.amplicon",my.target,my.assay,"ggbio",sep="."))
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
	tracks(Transcripts = p1, p2, WGBS=p9, Fluidigm=p11, Amplicon=p7, heights = c(1.3,0.5,3,3,3)) + xlim(ranges(gr.target.ext))
	dev.off()

	# zoom-in amplicon regions
	out.file.name<-file.path("~/results/RnBeads",my.project,my.cpg.context,paste("wgbs.fld.by.amplicon",my.target,my.assay,"zoom.ggbio",sep="."))
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
	tracks(Transcripts = p1, p2, WGBS=p9, Fluidigm=p11, Amplicon=p7, heights = c(1.3,0.5,3,3,3)) + xlim(ranges(gr.target.ext))
	dev.off()

	# zoom-in fwd amplicon regions
	out.file.name<-file.path("~/results/RnBeads",my.project,my.cpg.context,paste("wgbs.fld.by.amplicon",my.target,"fwd",my.assay,"zoom.ggbio",sep="."))
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
	tracks(Transcripts = p1, p2, WGBS=p9.fwd, Fluidigm=p11.fwd, Amplicon=p7.fwd, heights = c(1.3,0.5,3,3,2)) + xlim(ranges(gr.target.ext))
	dev.off()

	# zoom-in rev amplicon regions
	out.file.name<-file.path("~/results/RnBeads",my.project,my.cpg.context,paste("wgbs.fld.by.amplicon",my.target,"rev",my.assay,"zoom.ggbio",sep="."))
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
	tracks(Transcripts = p1, p2, WGBS=p9.rev, Fluidigm=p11.rev, Amplicon=p7.rev, heights = c(1.3,0.5,3,3,2)) + xlim(ranges(gr.target.ext))
	dev.off()
}

################################################################
# Methylation profile re-formatted (row: CpG, column: Samples) #
################################################################
# dcast from 'reshape2'
wgbs.chr7[rowSums(!is.na(wgbs.chr7))==ncol(wgbs.chr7),]
wgbs.chr13[rowSums(!is.na(wgbs.chr13))==ncol(wgbs.chr13),]

fld.chr7<-dcast(as.data.frame(gr.fld.merged[seqnames(gr.fld.merged)=="chr7",]), start ~ SampleName, value.var="score"); rownames(fld.chr7)<-fld.chr7$start; fld.chr7$start<-NULL;
fld.chr13<-dcast(as.data.frame(gr.fld.merged[seqnames(gr.fld.merged)=="chr13",]), start ~ SampleName, value.var="score"); rownames(fld.chr13)<-fld.chr13$start; fld.chr13$start<-NULL; 

############################################
## Concordance between Bismark and BisSNP ##
############################################
my.cmd<-"bedtools intersect -a ~/results/SLX-8772.Homo_sapiens.v1/Bismark/MethylCall/FLD0193/SLX-8772.FLD0193.r_1_val_1.fq.gz_bismark_bt2_pe.q1.ontarget.sorted.bed -b ~/results/SLX-8772.Homo_sapiens.v1/BisSNP/FLD0193/SLX-8772.FLD0193.cpg.filtered.CG.bed -wa -wb -s  | awk '{print $4,$5,$10,$11}'"
agreement<-system(my.cmd, intern=T)
bismark.bissnp<-do.call("rbind", lapply(strsplit(agreement, " "), function(i) as.numeric(cbind(i[1],i[3])))) # bismark bissnp

ggplot(as.data.frame(bismark.bissnp), aes(V1,V2)) + geom_point(size=3) + geom_smooth(method=lm,se=FALSE) + labs(x="Bismark",y="BisSNP") + ggtitle("FLD0193")
smoothScatter(bismark.bissnp)

##########################################
## Concordance between WGBS and Fludigm ## 
##########################################
# see README how to make this file
my.cmd<-paste0("ls ~/results/",my.slx[my.fld.ver],".Homo_sapiens.v2/BisSNP/FLD*/",my.slx[my.fld.ver],".FLD*.cpg.concordance.bed")

dummy.bed<-system(my.cmd, intern=T)
barcode.lst<-lapply(dummy.bed, function(i) unlist(strsplit(i, "/"))[7])
foo<-lapply(dummy.bed, read.delim, header=F)
names(foo)<- do.call(c, barcode.lst)
wgbs.fld<-lapply(foo, function(i) i[i$V5>=10 & i$V16>=10,c("V1","V2","V3","V4","V5","V15","V16","V6")]) # V4:Fludigm %met, V15:WGBS %met

# Colour by Strand
out.file.name<-file.path("~/results/RnBeads",my.project,my.cpg.context,"FLD.WGBS.concordance.by.strand")
pdf(file=paste(out.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title=my.project)
lapply(names(wgbs.fld), 
	function(i) ggplot(wgbs.fld[[i]], aes(x=V4, y=V15, color=V6)) + 
				geom_point(size=3,alpha=0.7) + ylim(0, 100) + xlim(0,100) +
				geom_smooth(method=lm,se=FALSE) + 
				labs(x=paste0("Fludigm (NO. of CpG = ",nrow(wgbs.fld[[i]]),")"),y="WGBS") + 
				ggtitle(paste0(my.samples[i]," (R^2=",round(summary(lm(wgbs.fld[[i]][,c("V4","V15")]))$r.squared,2),")")) +
				scale_color_discrete(name="Strand")
) 
dev.off()

# Colour by Chromosome 
out.file.name<-file.path("~/results/RnBeads",my.project,my.cpg.context,"FLD.WGBS.concordance.by.chr")
pdf(file=paste(out.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'),  width=11.7, height=8.3, title=my.project)
lapply(names(wgbs.fld), 
	function(i) ggplot(wgbs.fld[[i]], aes(x=V4, y=V15, color=V1)) + 
				geom_point(size=3,alpha=0.7) + ylim(0, 100) + xlim(0,100) + 
				geom_smooth(method=lm,se=FALSE) + 
				labs(x=paste0("Fludigm (NO. of CpG = ",nrow(wgbs.fld[[i]]),")"),y="WGBS") + 
				ggtitle(paste0(my.samples[i]," (R^2=",round(summary(lm(wgbs.fld[[i]][,c("V4","V15")]))$r.squared,2),")")) +
				scale_color_discrete(name="Chromosome")
)
dev.off()

#ggplot(wgbs.fld[["FLD0194"]], aes(x=V4, y=V15)) + geom_point(aes(color=name),size=3) + geom_smooth(method=lm,se=FALSE) + labs(x="Fludigm",y="WGBS") + ggtitle(my.samples["FLD0194"])
#ggplot(wgbs.fld[["FLD0194"]], aes(x=V4, y=V15,color=V6)) + geom_point(size=3,alpha=0.8) + geom_smooth(method=lm,se=FALSE,fullrange=TRUE) + labs(x="Fludigm",y="WGBS") + ggtitle(my.samples["FLD0194"]) # by strand
#ggplot(wgbs.fld[["FLD0194"]], aes(x=V4, y=V15,color=V1)) + geom_point(size=3,alpha=0.8) + geom_smooth(method=lm,se=FALSE,fullrange=TRUE) + labs(x="Fludigm",y="WGBS") + ggtitle(my.samples["FLD0194"]) # by chr 
#smoothScatter(wgbs.fld[[1]][,c("V4","V15")], cex=5, xlab="Fludigm", ylab="WGBS")
#cor(wgbs.fld[[1]][,c("V4","V15")]) isa 'matrix'
#lm(wgbs.fld[[1]][,c("V4","V15")]) # isa 'lm'

##################################################
## NO of Mapped Read by Chromosome (Fludigm v2) ##
##################################################
my.cmd<-"samtools view ~/results/SLX-8772.Homo_sapiens.v1/Bismark/Alignment/FLD0194/SLX-8772.FLD0194.000000000-AF7HP.s_1.r_1_val_1.fq.gz_bismark_bt2_pe.q1.bam | awk '{print $3}' | sort | uniq -c | awk '{print substr($2,4),$1/2}' | sort -k1,1n"
dummy<-system(my.cmd, intern=T)
strsplit(dummy, " ")
simplify2array(strsplit(dummy, " ")) # equivalent: sapply(strsplit(dummy, " "), function(i) cbind(i[1],i[2]))
do.call("rbind", strsplit(dummy, " "))

foo<-lapply(names(my.samples), 
				function(i){
					my.cmd<-paste0("samtools view ~/results/SLX-8772.Homo_sapiens.v1/Bismark/Alignment/",i,"/SLX-8772.",i,".000000000-AF7HP.s_1.r_1_val_1.fq.gz_bismark_bt2_pe.q1.bam | awk '{print $3}' | sort | uniq -c | awk '{print substr($2,4),$1/2}' | sort -k1,1n"); 
					cbind(`Barcode`=i,do.call("rbind",strsplit(system(my.cmd, intern=T)," "))) 
				}
	)
read.by.chr<-do.call("rbind",foo)
read.by.chr<-as.data.frame(read.by.chr)
colnames(read.by.chr)<-c("Barcode","Chr","NO")
read.by.chr$NO<-as.numeric(as.character(read.by.chr$NO)) # factor to numeric must be proceeded by factor to character, then numeric
#> head(read.by.chr)
#  Barcode Chr   NO
#1 FLD0193   X   90
#2 FLD0193   1 2590
#3 FLD0193   2  146
#4 FLD0193   3  126
#5 FLD0193   4  102
#6 FLD0193   5  136
ggplot(read.by.chr, aes(x=Chr,y=NO, color=Barcode)) + geom_point()
ggplot(read.by.chr, aes(x=Chr,y=log(NO), color=Barcode)) + geom_point()

bar<-dcast(read.by.chr, Barcode ~ Chr)
#> dcast(read.by.chr, Barcode ~ Chr)
#  Barcode    1    10  11  12 13  14  15  16  17  18  19   2 20 21 22   3   4   5   6   7   8   9   X  Y
#1 FLD0193 2590 32437  95 114 56  75  67  90 101  50 110 146 61 21 37 126 102 136 120 105 150  98  90 NA
#2 FLD0194 2088 30196 103 130 82  70  71 106 139  67 121 209 62 28 64 152 124 136 120 136 180 103  96 NA
rownames(bar)<-bar$Barcode
bar$Barcode<-NULL
heatmap.2(as.matrix(bar), col=hmcol, trace="none", dendrogram="none", key.xlab="NO. of Mapped Read", keysize=1, cexCol=0.8, cexRow=0.8, srtCol=0, Rowv=F, Colv=F, cellnote=as.matrix(bar), notecol="black", notecex=0.7)


logger.completed()
logger.close()
cat("All is done", "\n")
