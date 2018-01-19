#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

options(scipen=999) # diable scientific notation
# make sure 
myCaller="Expression"
TR_PREFIX='GRCh37' # GRCh37|GRCh38
source ("~/Pipelines/config/DEG.R") # load config

# Define up/down flanking size
my.upflank=10^6 # 1Mb up from target TSS
my.downflank=10^6 # 1Mb down from target TSS 

#######################
## Read RNA-Seq data ##
#######################
# load DESeq2 RData for this sampleType
deseq.RData <- paste0(deseq.dir, "/deseq.",sampleType,'.RData')
if(file.exists(deseq.RData)){
	cat("loading DESeq2 RData...\n")
	load(deseq.RData) #dds, ddsFlt, res, resFlt, rld, rld.nb, vsd
	cat("dds, res, rld, and ddsExonicGene loaded\n")
}else{
	stop(paste(deseq.RData, "does not exist\n"))
}

cat("making fpm...\n")
ddsFpm <-fpm(dds)
cat("making fpkm...\n")
ddsFpkm <-fpkm(ddsExonicGene) # isa 'matrix'

# set the output dir
ExpHome <- paste0(RNAresultDir,'/TopDMR'); if(!file.exists(ExpHome)){dir.create(ExpHome)}
expDir=paste0(ExpHome,'/',sampleType); if(!file.exists(expDir)){dir.create(expDir)}

############################################
# Target region (i.e. DMR) to interrorage ##
############################################
cat("reading target file...\n")
#my.target.file="~/Pipelines/data/POPS.top.DMR.txt" # BED format
my.target.file="~/Pipelines/data/CSMD1/CSMD1.boy.girl.dmr.bed"
if(FALSE){
	# below reads .bed file which is 0-based for Start, 1-based for End
	# However, GRange Start, End are both 1-based
	my.targets<-read.delim(my.target.file, header=FALSE)
	colnames(my.targets)<-c('chr','start','end','id','score','strand'); my.gr=makeGRangesFromDataFrame(my.targets,keep.extra.columns=TRUE); genome(my.gr)<-"hg19"
}
my.gr <- rtracklayer::import.bed(my.target.file, genome="hg19")
#my.gr <- with(my.gr, GRanges(seqnames,IRanges(start(my.gr),end(my.gr)),strand="*",name, score)) #ignore strand info
# extend the target regions using the promoters function which is within the GenomicRanges
cat("extending +- 1MB...\n")
my.gr.ext<-promoters(my.gr, upstream=my.upflank, downstream=my.downflank) # includes +- flanking from TSS
#lapply(my.gr, function(i) promoters(i, upstream=my.upflank, downstream=end(i)-start(i)+1+my.downflank)) # includes +- my.up(down)flank from the target

my.overlap=as.matrix( findOverlaps(gr.ensg, my.gr.ext, ignore.strand=TRUE, minoverlap=10L) ) #By default, any overlap is accepted
my.target.ensg<-gr.ensg[my.overlap[,"queryHits"]]

#####################################
## Get Annotations of Target Genes ##
#####################################
deseq.res<-as.data.frame(res[names(my.target.ensg),])
deseq.res$ensg<-rownames(deseq.res)

#grch37
fields <- c("ensembl_gene_id","hgnc_symbol","description","gene_biotype")
anno.grch37=getBM(attributes = fields, filters = "ensembl_gene_id", values = mcols(my.target.ensg)[,"gene_id"], mart = grch37)
my.target.anno.grch37<-merge(as.data.frame(my.target.ensg), anno.grch37, by.x="gene_id",by.y="ensembl_gene_id", all.x=TRUE)
my.target.anno.grch37<-merge(my.target.anno.grch37, deseq.res, by.x="gene_id", by.y="ensg", all.x=TRUE)

select<-my.target.anno.grch37$gene_id %in% rownames(ddsFpkm)
dummy<-as.data.frame(cbind(fpkm.mean=rowMeans(ddsFpkm[my.target.anno.grch37[select,]$gene_id,]), fpkm.var=rowVars(ddsFpkm[my.target.anno.grch37[select,]$gene_id,])))
dummy$gene_id<-rownames(ddsFpkm[my.target.anno.grch37[select,]$gene_id,])
my.target.anno.grch37<-merge(my.target.anno.grch37, dummy, by="gene_id", all.x=TRUE)

cat("writing output...\n")
my.file.grch37<-paste0(expDir,"/grch37.top.DMR.",sampleType,".csv")
write.csv(my.target.anno.grch37[order(my.target.anno.grch37$seqnames, my.target.anno.grch37$start),], file=my.file.grch37)

#grch38
anno.grch38=getBM(attributes = fields, filters = "ensembl_gene_id", values = mcols(my.target.ensg)[,"gene_id"], mart = grch38)
my.target.anno.grch38<-merge(as.data.frame(my.target.ensg), anno.grch38, by.x="gene_id",by.y="ensembl_gene_id", all.x=TRUE)
my.target.anno.grch38<-merge(my.target.anno.grch38, deseq.res, by.x="gene_id", by.y="ensg", all.x=TRUE)

select<-my.target.anno.grch38$gene_id %in% rownames(ddsFpkm)
dummy<-as.data.frame(cbind(fpkm.mean=rowMeans(ddsFpkm[my.target.anno.grch38[select,]$gene_id,]), fpkm.var=rowVars(ddsFpkm[my.target.anno.grch38[select,]$gene_id,])))
dummy$gene_id<-rownames(ddsFpkm[my.target.anno.grch38[select,]$gene_id,])
my.target.anno.grch38<-merge(my.target.anno.grch38, dummy, by="gene_id", all.x=TRUE)

cat("writing output...\n")
my.file.grch38<-paste0(expDir,"/grch38.top.DMR.",sampleType,".csv")
write.csv(my.target.anno.grch38[order(my.target.anno.grch38$seqnames, my.target.anno.grch38$start),], file=my.file.grch38)

####################################################
## Get Expression Level (ie. FPM) of Target Genes ##
####################################################
#barplot(counts(dds)["ENSG00000164669",], main="INTS4L1", xlab="sample", ylab="counts", col=c("green","green","orange","orange"))
#plotExpMulti("ENSG00000146757") 
#plotExpMulti("ENSG00000146757","FPM") # make sure 'ddsFpm' defined?
if(FALSE){
	my.file.name<-paste0(expDir,"/FPM.top.DMR.2Mbp.",sampleType)
	pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=22, height=9, title=paste0("FPM.top.DMR.2Mbp.",sampleType))
	sapply(my.target.anno.grch37[order(my.target.anno.grch37$seqnames, my.target.anno.grch37$start),c("gene_id")], plotExpMulti,"FPM") 
	dev.off()

	my.file.name<-paste0(expDir,"/Count.top.DMR.2Mbp.",sampleType)
	pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=22, height=9, title=paste0("Count.top.DMR.2Mbp.",sampleType))
	sapply(my.target.anno.grch37[order(my.target.anno.grch37$seqnames, my.target.anno.grch37$start),c("gene_id")], plotExpMulti) 
	dev.off()
}

cat("All is done", "\n")
