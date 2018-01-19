#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

##############################
## Amplicon Coveage Heatmap ##
##############################
library(reshape2) # for 'melt' and 'dcast'
library(gplots) # heatmap.2
library(RColorBrewer) # for brewer.pal
#hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(10)
hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)

#my.fld.ver="FLD.v2"; my.assay<-"BS"
#my.fld.ver="FLD.v3"; my.assay<-"oxBS"
my.fld.ver="FLD.v4"; my.assay<-"oxBS"
my.slx<-c(`FLD.v2`="SLX-8772",`FLD.v3`="SLX-8773", `FLD.v4`="SLX-10295")

## bin/R/RnBeads/make.FLD.meta.sample.R
fld.info.file<-file.path("~/Pipelines/data/MJ.FLD",my.fld.ver,"mj.FLD.RnBeads.sample.sheet.txt")
fld.info<-read.table(fld.info.file, header=T)

#fld.bed.file<-file.path(paste0("~/Pipelines/data/MJ.FLD/",my.fld.ver,"/",my.fld.ver,".chr10.BedFile.txt"))
fld.bed.file<-file.path(paste0("~/Pipelines/data/MJ.FLD/",my.fld.ver,"/",my.fld.ver,".BedFile.txt"))
fld.bed<-read.table(fld.bed.file,header=T)

fld.file<-paste0("~/results/RnBeads/Meta/",my.fld.ver,"/",my.fld.ver,".sample.info.csv") # 141 FLD v2 samples 
if(file.exists(fld.file)){
	fld.meta.sample<-read.csv(fld.file)
	# from bin/R/RnBeads/make.FLD.meta.sample.R
}else{
	stop(paste0(fld.file, " not found\n"))
}

if(my.fld.ver=="FLD.v2"){
	# Technical Sample Info
	batch.file<-"~/Pipelines/data/MJ.FLD/FLD.v2/mj.FLD.technical.batches.txt"
	batch.info<-read.table(batch.file, header=T)
	rownames(batch.info)<-batch.info$Barcode
	select<-!is.na(fld.meta.sample$PairGroup) & fld.meta.sample$PairGroup=='q1' 
}else{
	# Analysis Design Sample Info
	select<-!is.na(fld.meta.sample$PairGroup) & fld.meta.sample$PairGroup=='q1' & fld.meta.sample$assay==my.assay
}
group1<-as.character(fld.meta.sample[select & fld.meta.sample$Condition=="AGA",c("SampleName")])
names(group1)<-as.character(fld.meta.sample[select & fld.meta.sample$Condition=="AGA",c("Barcode")])

group2<-as.character(fld.meta.sample[select & fld.meta.sample$Condition=="SGA",c("SampleName")])
names(group2)<-as.character(fld.meta.sample[select & fld.meta.sample$Condition=="SGA",c("Barcode")])

my.samples=c(group1,group2)

# 1. Breadth of Coverage Per Amplicon
#SLX=SLX-8772
#echo "Barcode Amplicon Base Length Cov" > ~/results/$SLX.Homo_sapiens.v1/Coverage/FLD.BoC.per.amplicon.txt
#coverageFiles=$(ls ~/results/$SLX.Homo_sapiens.v1/Coverage/FLD*/$SLX.FLD*.r_1_val_1.fq.gz_bismark_bt2_pe.q1.sorted.bam.CoverageHist.per.amplicon)
#amplicons=$(awk '{print $4}' ~/Pipelines/data/MJ.FLD/FLD.v2/mj.amplicon.coordinates.sorted.bed)
#for i in $coverageFiles;do barcode=`echo $i| cut -d'/' -f7`;for k in $amplicons; do printf "$barcode $k "; awk "BEGIN{sum=0}\$4==\"$k\"{len=\$9}\$4==\"$k\" && \$7>=10{sum+=\$8}END{print sum,len,sum/len}" $i; done; done >> ~/results/$SLX.Homo_sapiens.v1/Coverage/FLD.BoC.per.amplicon.txt

# 2. Depth of Coverage Per Amplicon
#echo "Barcode Amplicon Base Length Cov" > ~/results/$SLX.Homo_sapiens.v1/Coverage/FLD.DoC.per.amplicon.txt
#coverageFiles=$(ls ~/results/$SLX.Homo_sapiens.v1/Coverage/FLD*/$SLX.FLD*.r_1_val_1.fq.gz_bismark_bt2_pe.q1.sorted.bam.CoverageHist.per.amplicon)
#amplicons=$(awk '{print $4}' ~/Pipelines/data/MJ.FLD/FLD.v2/mj.amplicon.coordinates.sorted.bed)
#for i in $coverageFiles;do barcode=`echo $i| cut -d'/' -f7`;for k in $amplicons; do printf "$barcode $k "; awk "BEGIN{sum=0}\$4==\"$k\"{len=\$9; sum+=\$7*\$8}END{print sum,len,sum/len}" $i; done; done >> ~/results/$SLX.Homo_sapiens.v1/Coverage/FLD.DoC.per.amplicon.txt

########################
## Breadh of Coverage ##
########################
my.boc.file=file.path("~/results",paste0(my.slx[my.fld.ver], ".Homo_sapiens.v1/Coverage/FLD.BoC.per.amplicon.txt"))
if(file.exists(my.boc.file)){
	fld.boc<-read.delim(my.boc.file, sep=" ")
}else{
	stop(paste0(fld.boc.file, " not found\n"))
}
BoC<-dcast(fld.boc, Barcode ~ Amplicon, value.var="Cov")
rownames(BoC)<-BoC$Barcode
BoC$Barcode<-NULL # remove Barcode
BoC<-BoC[grepl("NKX",names(BoC))] # NKX amplicons only
amplicons<-colnames(BoC)
colnames(BoC)<-as.numeric(unlist(lapply(strsplit(colnames(BoC),"[.]"), function(i) i[2]))) # remove "NKX1-2." from "NKX1-2.2"
BoC<-BoC[,order(as.numeric(colnames(BoC)))]
colnames(BoC)<-paste0('NKX1-2.',colnames(BoC))

# BoC Heatmap (for all samples)
filter<-rownames(BoC) %in% fld.boc[fld.boc$assay==my.assay,c("Barcode")]
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".amplicon.boc.heatmap")
tiff(filename=paste0(out.file.name,".tiff"),width=9,height=11.7,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("Avg. Breadth of Coverage (>=10x ",my.assay,")")
heatmap.2(as.matrix(BoC[filter,])*100, col=hmcol, trace="none", dendrogram="row", cexCol=0.7, cexRow=0.5, key.xlab="Breadth of Coverage (%)", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, Colv=F, main=my.main) # cexCol to adjust column labels
dev.off()

# BoC Heatmap (for set 1 samples only)
filter<-names(my.samples)
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".set1.amplicon.boc.heatmap")
tiff(filename=paste0(out.file.name,".tiff"),width=11.7,height=8.27,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("Avg. Breadth of Coverage (>=10x ",my.assay,")")
heatmap.2(as.matrix(BoC[filter,])*100, col=hmcol, trace="none", dendrogram="none", cexCol=0.7, cexRow=0.8, key.xlab="Breadth of Coverage (%)", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, Rowv=F, Colv=F, main=my.main)
dev.off()

avg.boc<-rowMeans(t(BoC)) # mean coverage (across samples)
avg.boc<-avg.boc[grepl("NKX",names(avg.boc))] # NKX amplicons only

avg.boc.q1<-rowMeans(t(BoC[rownames(BoC) %in% names(my.samples),]))
avg.boc.q1<-avg.boc.q1[grepl("NKX",names(avg.boc.q1))] # NKX amplicons only

######################
## Depth of Coveage ##
######################
my.doc.file=file.path("~/results",paste0(my.slx[my.fld.ver], ".Homo_sapiens.v1/Coverage/FLD.DoC.per.amplicon.txt"))
if(file.exists(my.doc.file)){
	fld.doc<-read.delim(my.doc.file, sep=" ")
}else{
	stop(paste0(fld.doc.file, " not found\n"))
}
DoC<-dcast(fld.doc, Barcode ~ Amplicon, value.var="Cov")
rownames(DoC)<-DoC$Barcode
DoC$Barcode<-NULL # remove Barcode
DoC<-DoC[grepl("NKX",names(DoC))] # NKX amplicons only
amplicons<-colnames(DoC)
colnames(DoC)<-as.numeric(unlist(lapply(strsplit(colnames(DoC),"[.]"), function(i) i[2]))) # remove "NKX1-2." from "NKX1-2.2"
DoC<-DoC[,order(as.numeric(colnames(DoC)))]
colnames(DoC)<-paste0('NKX1-2.',colnames(DoC))

# DoC Heatmap (for all samples)
filter<-rownames(DoC) %in% fld.doc[fld.doc$assay==my.assay,c("Barcode")]
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".amplicon.doc.heatmap")
tiff(filename=paste0(out.file.name,".tiff"),width=9,height=11.7,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("Avg. Depth of Coverage (",my.assay,")")
heatmap.2(as.matrix(DoC[filter,]), col=hmcol, trace="none", dendrogram="row", cexCol=0.7, cexRow=0.5, key.xlab="Depth of Coverage", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, Colv=F, main=my.main) # cexCol to adjust column labels
dev.off()

# DoC Heatmap (for set 1 samples)
filter<-names(my.samples)
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".set1.amplicon.doc.heatmap")
tiff(filename=paste0(out.file.name,".tiff"),width=11.7,height=9,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("Avg. Depth of Coverage (",my.assay,")")
heatmap.2(as.matrix(DoC[filter,]), col=hmcol, trace="none", dendrogram="none", cexCol=0.7, cexRow=0.8, key.xlab="Depth of Coverage", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, Rowv=F,Colv=F, main=my.main) # cexCol to adjust column labels
dev.off()

# log10(DoC) Heatmap (for all samples)
filter<-rownames(DoC) %in% fld.doc[fld.doc$assay==my.assay,c("Barcode")]
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".amplicon.log10.doc.heatmap")
tiff(filename=paste0(out.file.name,".tiff"),width=9,height=11.7,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("Avg. Log10(Depth of Coverage) (",my.assay,")")
heatmap.2(as.matrix(log10(DoC[filter,]+1)), col=hmcol, trace="none", dendrogram="row", cexCol=0.7, cexRow=0.5, key.xlab="log10(Doc)", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, Colv=F, main=my.main) # cexCol to adjust column labels
dev.off()

# log10(DoC) Heatmap (for set 1 samples)
filter<-names(my.samples)
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,".set1.amplicon.log10.doc.heatmap")
tiff(filename=paste0(out.file.name,".tiff"),width=11.7,height=8.27,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("Avg. Log10(Depth of Coverage) (",my.assay,")")
heatmap.2(as.matrix(log10(DoC[filter,]+1)), col=hmcol, trace="none", dendrogram="none", cexCol=0.7, cexRow=0.5, key.xlab="log10(Doc)", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, Rowv=F,Colv=F, main=my.main) # cexCol to adjust column labels
dev.off()


# for ggbio later
avg.doc<-rowMeans(t(DoC)) # mean coverage (across samples)
avg.doc<-avg.doc[grepl("NKX",names(avg.doc))] # NKX amplicons only

avg.doc.q1<-rowMeans(t(DoC[rownames(DoC) %in% names(my.samples),]))
avg.doc.q1<-avg.doc.q1[grepl("NKX",names(avg.doc.q1))] # NKX amplicons only

#############################
## PCA based on BoC or DoC ##
#############################
cat("making pca based on coverage...\n")
pca<-prcomp(BoC, scale.=F) # input (row:samples, col:amplicons)
#pca<-prcomp(BoC, scale.=T) # it failes
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

pca.plot <- ggplot(cbind(pca$x[,1:3], batch.info,`breadth`=rowMeans(BoC)), aes(PC1,PC2)) +
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

pca<-prcomp(DoC, scale.=T) # input (row:samples, col:amplicons)
#pca<-prcomp(DoC, scale.=F) # input (row:samples, col:amplicons)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

pca.plot <- ggplot(cbind(pca$x[,1:3], batch.info,`depth`=rowMeans(DoC)), aes(PC1,PC2)) +
		geom_point(aes(color=PCR.batch,size=depth),alpha=0.5) + scale_size(range=c(3,10)) +
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

############################
## Amplicon Coveage GGBIO ##
############################
library(methylKit) # for 'read.feature.flank'
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(biomaRt)
library(ggbio)

my.local.db<-"~/data/Annotation/hg19ensGene.sqlite" #annotation db
my.cgi.file='~/data/Annotation/cpgi.hg19.bed'
my.upflank=2000 # up from TSS
my.downflank=200 # down form TES

hg19ensGene<- loadDb(my.local.db) # from GenomicFeatures
gr.ensg<-genes(hg19ensGene)

cpg.obj=read.feature.flank(location=my.cgi.file, feature.flank.name=c("CpGi","shores")) # a ‘GenomicRangesList’ contatining one GRanges object for flanks and one for GRanges object for the main feature.

my.target<-'MAB21L1' #NKX1-2
my.chr<-'chr13'
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
my.ensg=getBM(attributes = c("ensembl_gene_id"), filters = "hgnc_symbol", values=my.target, mart = grch37) # is a data.frame
my.gr.ensg<-gr.ensg[mcols(gr.ensg)$gene_id %in% my.ensg] # GRanges of my.target gene 
my.gr.ext.ensg<-promoters(my.gr.ensg, upstream=my.upflank, downstream=end(my.gr.ensg)-start(my.gr.ensg)+1+my.downflank) # includes 2000 +TSS and 100 +TES 

# overlap between cpgi and my.target 
my.overlap=as.matrix( findOverlaps(cpg.obj[["CpGi"]], my.gr.ext.ensg) )
my.cpgi<-cpg.obj$CpGi[my.overlap[,"queryHits"]]

my.AmpBedFile=file.path("~/Pipelines/data/MJ.FLD",my.fld.ver,"mj.amplicon.coordinates.sorted.bed") # FLD amplicon coordinate
if(!file.exists(my.AmpBedFile)){
	stop(paste0(my.AmpBedFile, " not found\n"))
}

#foo<-read.delim(my.AmpBedFile, header=F)
#colnames(foo)<-c('chr','start','end','id','score','strand')
#gr.target <- with(foo, GRanges(chr, IRanges(start, end), strand, score, id=id))
gr.target <- rtracklayer::import.bed(my.AmpBedFile)
gr.target <- gr.target[as.character(seqnames(gr.target))==my.chr,] # my.chr only
names(gr.target)<-mcols(gr.target)$name # assign amplicon name

if(exists("avg.boc") & exists("avg.doc")){
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
}else{
	mcols(gr.target)<-DataFrame(with(mcols(gr.target), 
						data.frame(
							name=name, 
							order=seq(5,200,by=3)[1:nrow(mcols(gr.target))] # amplicon left-most order 
						)
					)) 
}

if(my.fld.ver!="FLD.v4"){
	## Top 10 CpGs from WGBS
	my.top.cpg.bed<-"~/Pipelines/data/NKX1-2/MJ.WGBS.SGA.AGA.top.chr10.cpg.named.bed" # top ranked CpG from MJ WGBS
	gr.top.cpg<- rtracklayer::import.bed(my.top.cpg.bed)
	names(gr.top.cpg)<-mcols(gr.top.cpg)$name # assign amplicon name
	gr.top10.cpg<-gr.top.cpg[mcols(gr.top.cpg)$score<=10,] # this score is just index for the region
}

# chromosome 
p0 <- Ideogram(genome = "hg19") + xlim(gr.target)
# transcript track
p1<-autoplot(hg19ensGene, which=my.gr.ext.ensg, fill = "orange", color = "grey") 
# cpgi track for my.target
p2<-autoplot(my.cpgi, fill = "darkgreen", color = "grey") 
# amplicon breadth of coveage (geom_rect)
if(exists("avg.boc")){
p3<-ggplot(gr.target) +  
	geom_rect(stat="identity", rect.height = 1, aes(fill=boc, y=order), col="black", ylab=NULL) +  
	scale_fill_gradient2(name="coveage(%)",low=hmcol[1], mid=hmcol[round(length(hmcol)/2)], high=hmcol[length(hmcol)], midpoint=50) + 
	theme(axis.ticks.y = element_blank(), axis.text.y = element_blank() ) +
	geom_text(x=start(gr.target)+250, y=seq(5,200,by=3)[1:nrow(mcols(gr.target))], label=names(gr.target), size=2.6)
}
# amplicon relative coverage (based on amplion target bed file) 
p4<-ggplot(gr.target)+stat_coverage()
# amplicon loci and ther names (geom_arrow)
p5<-ggplot(gr.target) + 
	geom_arrow(stat="identity", aes(y=order), col="black", ylab=NULL) +  
	theme(axis.ticks.y = element_blank(), axis.text.y = element_blank() ) +
	geom_text(x=start(gr.target)+100, y=seq(5,200,by=3)[1:nrow(mcols(gr.target))]+1, label=names(gr.target), size=3)
# top 10 CpGs
if(exists("gr.top10.cpg")){
p6<-ggplot(gr.top10.cpg) + geom_rect()
}
# amplicon loci and ther names colored by strand 
p7<-ggplot(gr.target) + 
	geom_rect(stat="identity", rect.height = 1, aes(fill=strand, y=order), col="black", ylab=NULL) +  
	theme(axis.ticks.y = element_blank(), axis.text.y = element_blank() ) +
	geom_text(x=start(gr.target)+250, y=seq(5,200,by=3)[1:nrow(mcols(gr.target))], label=names(gr.target), size=2.6)
# amplicon depth of coveage (geom_rect)
if(exists("avg.doc")){
p8<-ggplot(gr.target) +  
	geom_rect(stat="identity", rect.height = 1, aes(fill=doc, y=order), col="black", ylab=NULL) +  
	scale_fill_gradient2(name="depth",low=hmcol[1], mid=hmcol[round(length(hmcol)/2)], high=hmcol[length(hmcol)], midpoint=max(avg.doc)/2) + 
	theme(axis.ticks.y = element_blank(), axis.text.y = element_blank() ) +
	geom_text(x=start(gr.target)+250, y=seq(5,200,by=3)[1:nrow(mcols(gr.target))], label=names(gr.target), size=2.6)
}
# 1. Amplicon Breadth of Coverage
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,"amplicon.breadh.coverage.ggbio")
tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
#tracks(Transcripts = p1, CPGi= p2, `Amplicon Relative Coveage`=p4, `Avg. Breadh of Coverage (per amplicon)`=p3, heights = c(1.3,0.5,0.5,8)) + xlim(ranges(my.gr.ext.ensg))  + theme_alignment() # grid=FALSE, border = T) # plot ggbio
tracks(Transcripts = p1, CPGi= p2, `Top CpGs`=p6, `Avg. Breadh of Coverage (per amplicon)`=p3, heights = c(1.3,0.5,0.5,8)) + xlim(ranges(my.gr.ext.ensg))  + theme_alignment() # grid=FALSE, border = T) # plot ggbio
dev.off()

# 2. Amplicon Depth of Coverage
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,"amplicon.depth.coverage.ggbio")
tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
tracks(Transcripts = p1, CPGi= p2, `Top CpGs`=p6, `Avg. Depth of Coverage (per amplicon)`=p8, heights = c(1.3,0.5,0.5,8)) + xlim(ranges(my.gr.ext.ensg))  + theme_alignment() # grid=FALSE, border = T) # plot ggbio
dev.off()

# 3. Amplicon by Strand
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,"amplicon.by.strand.ggbio")
tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw')
tracks(Transcripts = p1, CPGi= p2, `Top CpGs`=p6, `Amplicon`=p7, heights = c(1.3,0.5,0.5,8)) + xlim(ranges(my.gr.ext.ensg))  + theme_alignment() # grid=FALSE, border = T) # plot ggbio
dev.off()

# 4. Top 10 CpGs loci
out.file.name<-paste0("~/results/RnBeads/",my.fld.ver,".",my.assay,"top10.cpg")
tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size
tracks(Transcripts = p1, `Top 10 CpGs`=p6,`Amplicons`=p5, heights=c(1,1,2)) + xlim(ranges(gr.target[mcols(gr.target)$name %in% c("NKX1-2.14", "NKX1-2.16","NKX1-2.17","NKX1-2.18","NKX1-2.19","NKX1-2.20","NKX1-2.21","NKX1-2.22"),])) 
dev.off()

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
my.cmd<-"ls ~/results/SLX-8773.Homo_sapiens.v2/BisSNP/FLD*/SLX-8773.FLD*.cpg.concordance.bed"# TS-primer remained
#my.cmd<-"ls ~/results/SLX-8772.Homo_sapiens.v2/BisSNP/FLD*/SLX-8772.FLD*.cpg.concordance.bed" # TS-primer removed
dummy.bed<-system(my.cmd, intern=T)
barcode.lst<-lapply(dummy.bed, function(i) unlist(strsplit(i, "/"))[7])
foo<-lapply(dummy.bed, read.delim, header=F)
names(foo)<- do.call(c, barcode.lst)
wgbs.fld<-lapply(foo, function(i) i[i$V5>=10 & i$V16>=10,c("V1","V2","V3","V4","V5","V15","V16","V6")]) # V4:Fludigm %met, V15:WGBS %met
wgbs.fld<-lapply(wgbs.fld, function(i) merge(i, as.data.frame(gr.top10.cpg), by.x=c("V1","V3","V6"), by.y=c("seqnames","end","strand"), all.x=TRUE)[,c("V1","V2","V3","V4","V5","V15","V16","V6","name","score")]) # join with top10 CpG

# Highlight Top 10 CpG
#out.file.name<-"~/results/RnBeads/FLD.v2.the.eight/CpG/all/FLD.the.eight.condordance.with.wgbs"
#out.file.name<-"~/results/RnBeads/FLD.v2.the.eight/CpG/all/FLD.the.eight.condordance.with.wgbs.ts.primer.removed" # SLX-8772.Homo_sapiens.v2
out.file.name<-"~/results/RnBeads/FLD.v3.the.eight/CpG/all/FLD.the.eight.condordance.with.wgbs"
pdf(file=paste(out.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))
lapply(names(wgbs.fld), 
	function(i) ggplot(wgbs.fld[[i]], aes(x=V4, y=V15)) + 
				geom_point(aes(color=name),size=3,alpha=0.7) + 
				geom_smooth(method=lm,se=FALSE) + 
				labs(x=paste0("Fludigm (NO. of CpG = ",nrow(wgbs.fld[[i]]),")"),y="WGBS") + 
				ggtitle(paste0(i," (R^2=",round(summary(lm(wgbs.fld[[i]][,c("V4","V15")]))$r.squared,2),")"))
) 
dev.off()

# Colour by Strand
out.file.name<-"~/results/RnBeads/FLD.v2.the.eight/CpG/all/FLD.the.eight.condordance.with.wgbs.by.strand"
#out.file.name<-"~/results/RnBeads/FLD.v2.the.eight/CpG/all/FLD.the.eight.condordance.with.wgbs.by.strand.ts.primer.removed" # SLX-8772.Homo_sapiens.v2
pdf(file=paste(out.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))
lapply(names(wgbs.fld), 
	function(i) ggplot(wgbs.fld[[i]], aes(x=V4, y=V15, color=V6)) + 
				geom_point(size=3,alpha=0.7) + 
				geom_smooth(method=lm,se=FALSE) + 
				labs(x=paste0("Fludigm (NO. of CpG = ",nrow(wgbs.fld[[i]]),")"),y="WGBS") + 
				ggtitle(paste0(my.samples[i]," (R^2=",round(summary(lm(wgbs.fld[[i]][,c("V4","V15")]))$r.squared,2),")"))
) 
dev.off()

# Colour by Chromosome 
out.file.name<-"~/results/RnBeads/FLD.v2.the.eight/CpG/all/FLD.the.eight.condordance.with.wgbs.by.chr" # SLX-8772.Homo_sapiens.v2
pdf(file=paste(out.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))
lapply(names(wgbs.fld), 
	function(i) ggplot(wgbs.fld[[i]], aes(x=V4, y=V15, color=V1)) + 
				geom_point(size=3,alpha=0.7) + 
				geom_smooth(method=lm,se=FALSE) + 
				labs(x=paste0("Fludigm (NO. of CpG = ",nrow(wgbs.fld[[i]]),")"),y="WGBS") + 
				ggtitle(paste0(my.samples[i]," (R^2=",round(summary(lm(wgbs.fld[[i]][,c("V4","V15")]))$r.squared,2),")"))
) 
dev.off()

#ggplot(wgbs.fld[["FLD0194"]], aes(x=V4, y=V15)) + geom_point(aes(color=name),size=3) + geom_smooth(method=lm,se=FALSE) + labs(x="Fludigm",y="WGBS") + ggtitle(my.samples["FLD0194"])
#ggplot(wgbs.fld[["FLD0194"]], aes(x=V4, y=V15,color=V6)) + geom_point(size=3,alpha=0.8) + geom_smooth(method=lm,se=FALSE,fullrange=TRUE) + labs(x="Fludigm",y="WGBS") + ggtitle(my.samples["FLD0194"]) # by strand
#ggplot(wgbs.fld[["FLD0194"]], aes(x=V4, y=V15,color=V1)) + geom_point(size=3,alpha=0.8) + geom_smooth(method=lm,se=FALSE,fullrange=TRUE) + labs(x="Fludigm",y="WGBS") + ggtitle(my.samples["FLD0194"]) # by chr 
#smoothScatter(wgbs.fld[[1]][,c("V4","V15")], cex=5, xlab="Fludigm", ylab="WGBS")
#cor(wgbs.fld[[1]][,c("V4","V15")]) isa 'matrix'
#lm(wgbs.fld[[1]][,c("V4","V15")]) # isa 'lm'

#############################################
## Bismark Top 10 Methylation HeatMap from ##
## Fludigm Only                            ##
#############################################
my.cmd<-paste0("ls ~/results/",my.slx[my.fld.ver],".Homo_sapiens.v1/Bismark/MethylCall/*/",my.slx[my.fld.ver],".*.r_1_val_1.fq.gz_bismark_bt2_pe.q1.ontarget.sorted.top10.bed")
bismark.bed<-system(my.cmd, intern=T)
barcode.lst<-lapply(bismark.bed, function(i) unlist(strsplit(i, "/"))[8])
foo<-lapply(bismark.bed, function(i){j<-read.delim(i, header=F)[,c("V1","V2","V3","V4","V5","V6","V10")]; colnames(j)<-c("chr","start","end","met","depth","str","rank");return(j)})
names(foo)<- do.call(c, barcode.lst)
foo.doc<-lapply(names(foo), function(i){j<-foo[[i]][,c("chr","start","end","depth","str","rank")]; colnames(j)<-c("chr","start","end",i,"str","rank");return(j)})
foo.met<-lapply(names(foo), function(i){j<-foo[[i]][,c("chr","start","end","met","str","rank")]; colnames(j)<-c("chr","start","end",i,"str","rank");return(j)})

bismark.cpg<-Reduce(function(x,y) merge(x,y,by=c("chr","start","end","str","rank")), foo.met) # Merge a list to data.frame
bismark.cpg<-bismark.cpg[order(bismark.cpg$start),]
rownames(bismark.cpg)<-bismark.cpg$rank

bismark.doc<-Reduce(function(x,y) merge(x,y,by=c("chr","start","end","str","rank")), foo.doc)
bismark.doc<-bismark.doc[order(bismark.doc$start),]
rownames(bismark.doc)<-bismark.doc$rank

filter<-names(my.samples[names(my.samples) %in% as.character(barcode.lst)])
out.file.name=file.path("~/results/RnBeads",paste0(my.fld.ver,".the.eight/CpG/all/",my.fld.ver,".",my.assay,".the.eight.top10.cpg.methylation.bismark"))
tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("% Methylation by Bismark (",my.assay,")")
heatmap.2(t(as.matrix(bismark.cpg[,filter])), col=hmcol, trace="none", dendrogram="none", key.xlab="% methylation", keysize=1, cexCol=0.8, cexRow=0.8, srtCol=30, Rowv=F, Colv=F, 
	cellnote=round(t(as.matrix(bismark.cpg[,filter])),2), notecol="black", notecex=0.7, main=my.main)
dev.off()

out.file.name=file.path("~/results/RnBeads",paste0(my.fld.ver,".the.eight/CpG/all/",my.fld.ver,".",my.assay,".the.eight.top10.cpg.depth.bismark"))
tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size
my.main=paste0("Depth of Coverage by Bismark (",my.assay,")")
heatmap.2(t(as.matrix(bismark.doc[,filter])), col=hmcol, trace="none", dendrogram="none", key.xlab="depth of coverage", keysize=1, cexCol=0.8, cexRow=0.8, srtCol=30, Rowv=F, Colv=F, 
	cellnote=round(t(as.matrix(bismark.doc[,filter])),2), notecol="black", notecex=0.7, main=my.main)
dev.off()

##################################################
## NO of Mapped Read by Chromosome (Fludigm v2) ##
##################################################
my.cmd<-"samtools view ~/results/SLX-8772.Homo_sapiens.v1/Bismark/Alignment/FLD0194/SLX-8772.FLD0194.000000000-AF7HP.s_1.r_1_val_1.fq.gz_bismark_bt2_pe.q1.bam | awk '{print $3}' | sort | uniq -c | awk '{print substr($2,4),$1/2}' | sort -k1,1n"
dummy<-system(my.cmd, intern=T)
strsplit(dummy, " ")
simplify2array(strsplit(dummy, " ")) # equivalent: sapply(strsplit(dummy, " "), function(i) cbind(i[1],i[2]))
do.call("rbind", strsplit(dummy, " "))

foo<-lapply(names(my.samples), function(i){my.cmd<-paste0("samtools view ~/results/SLX-8772.Homo_sapiens.v1/Bismark/Alignment/",i,"/SLX-8772.",i,".000000000-AF7HP.s_1.r_1_val_1.fq.gz_bismark_bt2_pe.q1.bam | awk '{print $3}' | sort | uniq -c | awk '{print substr($2,4),$1/2}' | sort -k1,1n"); cbind(`Barcode`=i,do.call("rbind",strsplit(system(my.cmd, intern=T)," "))) })
read.by.chr<-do.call("rbind",foo)
read.by.chr<-as.data.frame(read.by.chr)
colnames(read.by.chr)<-c("Barcode","Chr","NO")
read.by.chr$NO<-as.numeric(as.character(read.by.chr$NO))
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


cat("All is done", "\n")
