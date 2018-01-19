#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# Last modified 16/Mar/2016

TR_PREFIX='GRCh37'

library(ggbio)
############
## Config ##
############
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R") # load colour config
source("~/Pipelines/bin/R/RoadMap/local.R")
cat("All source loaded\n")

my.upflank=1500 # up from TSS
my.downflank=300 # down form TES

# ignore strand info
subject.keys=c("seqnames","start","end")
query.key=c("V1","V2","End")

#########################
## Input Target Region ##
#########################
#import target regions of interest
my.target="CSMD1" # KDM5C, KDM6A, CSMD1, NKX1-2, XIST, AVPR2, Chr7.DMR, PNPLA4
if(my.target=="NKX1-2"){
	my.gr.ensg <- rtracklayer::import.bed("~/Pipelines/data/NKX1-2/NKX1-2.up.1500.down.500.bed", genome="hg19")
}else if(my.target=="CSMD1"){
	my.gr.ensg <- rtracklayer::import.bed("~/Pipelines/data/CSMD1/CSMD1.boy.girl.dmr.bed", genome="hg19")
}else if(my.target=="Chr7.DMR"){
	my.gr.ensg <- GenomicRanges::reduce(rtracklayer::import.bed("~/Pipelines/data/Chr7.DMR/chr7.top2.bed"))
	my.gr.ensg<-GenomicRanges::promoters(my.gr.ensg, upstream=500, downstream=end(my.gr.ensg)-start(my.gr.ensg)+1+500) # extends 500 up and down 
}else{
	gr.subject<-get.region("Genes")
	#select<-!is.na(gr.subject$symbol) & grepl(my.target, gr.subject$symbol)
	select<-!is.na(gr.subject$symbol) & gr.subject$symbol==my.target # exact name match
	my.gr.ensg<-gr.subject[select,]
	my.gr.ensg<-GenomicRanges::promoters(my.gr.ensg, upstream=my.upflank, downstream=end(my.gr.ensg)-start(my.gr.ensg)+1+my.downflank) # extends 1500 +TSS and 300 +TES 
}

subject<-as.data.frame(my.gr.ensg)
cat("\tRemoving leading 'chr'...\n")
#subject$seqnames<-simplify2array(strsplit(as.character(subject$seqnames), "chr"))[2,] # chrX = > X
subject$seqnames<-substr(subject$seqnames, 4, 5) #chrX=>X
cat("\tconvert to data.table...\n")
dt.subject<-as.data.table(subject) # isa data.table
setkeyv(dt.subject,subject.keys) # regions of interests 
rm(subject)

cat("Subject loaded\n")
####################
# transcript track # 
####################
if(my.target!="Chr7.DMR"){
	p.tr<-ggbio::autoplot(hg19ensGene, which=my.gr.ensg, fill = "orange", color = "grey") # hg19ensGene defined in 'Annotation.R'
}

##############
# chromosome #
##############
p.chr <- Ideogram(genome = "hg19", xlabel=TRUE) + xlim(ranges(my.gr.ensg))

############################
# cpgi track for my.target #
############################
#cpg.obj=methylKit::read.feature.flank(location=my.cgi.file, feature.flank.name=c("CpGi","shores")) # a ‘GenomicRangesList’ contatining one GRanges object for flanks and one for GRanges object for the main feature.
cpg.obj<-Reduce(c,RnBeads::rnb.get.annotation(type="cpgislands", assembly = "hg19")) # isa 'GRanges'
my.overlap=as.matrix( findOverlaps(cpg.obj, my.gr.ensg) )
my.cpgi <-cpg.obj[my.overlap[,"queryHits"]]
p.cpgi <-ggbio::autoplot(my.cpgi, fill = "darkgreen", color = "grey") 

############
## Tissue ##
############
my.cpg.type="CG" # CG, CHH, CHG
dt.tissue.locus<-list()
for(my.tissue in c('PT',avail.tissues[!avail.tissues %in% "PA"])){
	## Load this RoadMap tissue
	dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)

	cat("\tFinding CpGs overlap...\n")
	system.time(dt.overlap<-foverlaps(dt.query, dt.subject, by.x=query.key, type="any", nomatch=0L))
	#> dt.overlap
	#   V1    Start      End Strand       V2 V3  V4 V7 V5.x V6.x V5.y V6.y    i.End
	#1:  X  2746554  2748235      *  2746590  - CGA  1    6   17    2   10  2746590
	#2:  X  2835950  2836236      *  2836235  + CGG  1   14   14    9   10  2836235
	if(nrow(dt.overlap)){
		dt.meth<-rbind(dt.overlap[,.(V1,V2,V5.x/V6.x*100,"Female",my.tissue)], dt.overlap[,.(V1,V2,V5.y/V6.y*100,"Male",my.tissue)]) # isa 'data.table'
		setnames(dt.meth,c("chr","start","meth","Gender","Tissue"))
		#dt.tissue.locus[[my.tissue]]<-ggplot(dt.meth, aes(start, meth, col=Gender)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col[["Gender"]])
		dt.tissue.locus[[my.tissue]]<-dt.meth
	}
}# end of my.tissue 

# Save the image 
my.file.name<-file.path("~/results/RoadMap/BS-Seq",paste(my.target,"meth.profile.roadmap.tiff",sep="."))
tiff(filename=my.file.name,width=8.27, height=11.7 ,units="in",res=300, compression = 'lzw') #A4 size 
if(my.target=="NKX1-2"){
	p.tissue<-ggplot(rbindlist(dt.tissue.locus), aes(start, meth)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + ggtitle(my.target) + theme_bw() + facet_grid(Tissue ~ .)
	ggbio::tracks(Transcripts = p.tr, CPGi=p.cpgi, Tissues=p.tissue, heights = c(1,0.5,9)) + xlim(ranges(my.gr.ensg))
}else{
	# Gender-specific
	p.tissue<-ggplot(rbindlist(dt.tissue.locus), aes(start, meth, col=Gender)) + geom_point(size=0.8,alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col[["Gender"]]) + ggtitle(my.target) + theme_bw() + facet_grid(Tissue ~ .)
	#print(p.tissue)
	if(length(my.cpgi)){
		if(my.target=="Chr7.DMR"){
			print(tracks(CPGi=p.cpgi, Tissues=p.tissue, heights = c(0.5,9)) + xlim(ranges(my.gr.ensg)))
		}else{
			print(tracks(Transcripts = p.tr, CPGi=p.cpgi, Tissues=p.tissue, heights = c(1,0.5,9)) + xlim(ranges(my.gr.ensg)))
		}
	}else{
		print(tracks(Transcripts = p.tr, Tissues=p.tissue, heights = c(1,9)) + xlim(ranges(my.gr.ensg))) 
	}
}
dev.off()

if(my.target=="NKX1-2"){
	my.file.name<-file.path("~/results/RoadMap/BS-Seq",paste(my.target,"meth.profile.roadmap.hz.tiff",sep="."))
	tiff(filename=my.file.name,width=11.7, height=8.27,units="in",res=300, compression = 'lzw') #A4 size 
	p.tissue<-ggplot(rbindlist(dt.tissue.locus), aes(start, meth, col=Tissue)) + geom_point(size=1.5,alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette) + theme_bw() + ggtitle(my.target)
	ggbio::tracks(Transcripts = p.tr, CPGi=p.cpgi, Tissues=p.tissue, heights = c(1,0.5,7)) + xlim(ranges(my.gr.ensg))
	dev.off()
}
