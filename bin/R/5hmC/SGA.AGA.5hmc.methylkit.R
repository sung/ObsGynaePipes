#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# data prepared from 'bin/R/5hmC/SGA.AGA.5hmc.R'
# basic 5hmC diff met analysis using mehtylkit
# also plot 5hmC level by BS, oxBS, 5hmC tracks using ggbio
# ~/results/RnBeads/5hmC.SGA.AGA/SGA.AGA.bs.oxbs.*

library(methylKit)
library(GenomicRanges)
library(biomaRt)
library(GenomicFeatures)
library(ggbio)

top.dir<-"/home/ssg29/results/SLX-8074.SLX-8080.v1"

pdf.file.name<- paste0(top.dir,'/MethylKit/5hmC.plot')
pdf(file=paste(pdf.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

my.depth<-10 # at least 10 depth
my.diff<-10 # 10% difference
my.pvalue<-0.05
# cores to use for the diffmet
my.np<-4
# promoter definitions
my.promoter.up=1500
my.promoter.down=500
# annotation files 
my.gene.file='~/data/Annotation/ENSG.hg19.bed'
my.cgi.file='~/data/Annotation/cpgi.hg19.bed'
#annotation db
my.local.db<-"~/data/Annotation/hg19ensGene.sqlite"
#up and down flank size
my.upflank=2000 # 2000bp up from TSS
my.downflank=100 # 100bp down form TES

# load annotation from local file
hg19ensGene<- loadDb(my.local.db)
gr.ensg<-genes(hg19ensGene)
gr.enst<-transcripts(hg19ensGene)
gr.promoters<-promoters(hg19ensGene, upstream=my.promoter.up, downstream=my.promoter.down)

# gene object
gene.obj=read.transcript.features(location=my.gene.file, up.flank=my.promoter.up, down.flank=my.promoter.down) # a ‘GRangesList’ containing locations of exons/introns/promoters/TSSes
#gene.obj=read.transcript.features(location=my.gene.file, up.flank=my.promoter.up, down.flank=my.promoter.down, unique.prom=FALSE) # a ‘GRangesList’ containing locations of exons/introns/promoters/TSSes
cpg.obj=read.feature.flank(location=my.cgi.file, feature.flank.name=c("CpGi","shores")) # a ‘GenomicRangesList’ contatining one GRanges object for flanks and one for GRanges object for the main feature.

#https://support.bioconductor.org/p/62347/
#ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

####################################
## 5hmC (BS-oxBS) for AGA2F,SGA2F ##
####################################
#my.filtered.5hmC.data<-file.path(top.dir,"MethylKit/AGA2F.SGA2F.non.neg.5hmC.RData") # problematic 5hmC data 
my.filtered.5hmC.data<-file.path(top.dir,"MethylKit/AGA2F.SGA2F.non.neg.5hmC.new.RData") # new 5hmC data (depth problem solved)
if(file.exists(my.filtered.5hmC.data)){
	cat("Loading my.filtered.5hmc.data...\n")
	load(my.filtered.5hmC.data)
	cat("mk.5hmC, mk.filtered.5hmC, mk.5hmC.meth, mk.filtered.5hmC.meth loaded\n")
}else{
	#my.5hmc.data<-file.path(top.dir,"MethylKit/AGA2F.SGA2F.5hmC.RData") # there could be CpGs where BS(>=10) whereas oxBS(<=10)
	my.5hmc.data<-file.path(top.dir,"MethylKit/AGA2F.SGA2F.5hmC.new.RData")
	if(file.exists(my.5hmc.data)){
		cat("Loading my.5hmc.data...\n")
		load(my.5hmc.data)
		cat("AGA2F.5hmC.strand, SGA2F.5hmC.strand, AGA2F.5hmC.destrand and SGA2F.5hmC.destrand loaded\n") 

		mk.5hmC =new("methylRawList",list(SGA2F.5hmC.strand,AGA2F.5hmC.strand),treatment=c(1,0)) # sample.ids=c("SGA2F.bs", "AGA2F.bs")

		# apply filters (+ve numC & depth>=10)
		dummy.list=lapply(mk.5hmC, function(i){dummy<-getData(i); select=dummy$numCs>=0 & dummy$coverage>=my.depth; dummy<-dummy[select,]; new("methylRaw",dummy,sample.id=i@sample.id,assembly=i@assembly,context=i@context,resolution=i@resolution) } ) 
		mk.filtered.5hmC =new("methylRawList",dummy.list,treatment=c(1,0))
	}else{
		file.list=list(
			file.path(top.dir,"BisSNP/SGA2F.5hmC.strand.methylkit.txt"), # SGA2F 5hmC methylKit file
			file.path(top.dir,"BisSNP/AGA2F.5hmC.strand.methylkit.txt")  # AGA2F 5hmC methylKit file
		)
		mk.5hmC=read(file.list, sample.id=list("SGA", "AGA"), assembly="hg19",treatment=c(1,0),context="CpG") # is a 'methylRawList'
		# Filter by depth
		mk.filtered.5hmC<-filterByCoverage(mk.5hmC,lo.count=my.depth, lo.perc=NULL, hi.count=NULL,hi.perc=NULL) # is a 'methylRawList'
	}
	# unite
	mk.5hmC.meth<-unite(mk.5hmC,destrand=FALSE) # a methylBase object
	mk.filtered.5hmC.meth<-unite(mk.filtered.5hmC,destrand=FALSE) # a methylBase object
	# Base-level Difference 
	#mk.5hmC.diff=calculateDiffMeth(mk.5hmC.meth, num.cores=my.np) # Fisher's exact test because 1 sample per group
																# slim=TRUE,weighted.mean=TRUE by default
																# weighted.mean only for multi-sample(>=2) within a group
																# diff for SGA-AGA because SGA listed firstly followed by AGA
	#mk.filtered.5hmC.diff=calculateDiffMeth(mk.filtered.5hmC.meth, num.cores=my.np) # Fisher's exact test because 1 sample per group
	#save(mk.filtered.5hmC, mk.filtered.5hmC.meth, mk.5hmC.diff, file=my.filtered.5hmC.data)
	save(mk.5hmC, mk.filtered.5hmC, mk.5hmC.meth, mk.filtered.5hmC.meth, file=my.filtered.5hmC.data)
}

################################
## 5mC (oxBS) for AGA2F,SGA2F ##
################################
my.filtered.5mC.data<-file.path(top.dir,"MethylKit/AGA2F.SGA2F.5mC.RData") 
if(file.exists(my.filtered.5mC.data)){
	cat("Loading my.filtered.5mc.data...\n")
	load(my.filtered.5mC.data)
	cat("mk.filtered.5mC, mk.5mC.meth, mk.5mC.diff loaded\n")
}else{
	oxBS.2F.file.list=list(
		file.path(top.dir,"BisSNP/A003/SLX-8074.SLX-8080.A003.cpg.filtered.homo.CG.methylkit.txt"), #oxBS SGA2F
		file.path(top.dir,"BisSNP/A001/SLX-8074.SLX-8080.A001.cpg.filtered.homo.CG.methylkit.txt")  #oxBS AGA2F
	)
	mk.5mC=read(oxBS.2F.file.list, sample.id=list("SGA2F.ox", "AGA2F.ox"), assembly="hg19",treatment=c(1,0),context="CpG") # is a 'methylRawList'
	# Filter by depth
	mk.filtered.5mC<-filterByCoverage(mk.5mC,lo.count=my.depth, lo.perc=NULL, hi.count=NULL,hi.perc=NULL) # is a 'methylRawList'
	# unite
	mk.5mC.meth<-unite(mk.filtered.5mC,destrand=FALSE) # a methylBase object
	# Base-level Difference 
	mk.5mC.diff=calculateDiffMeth(mk.5mC.meth, num.cores=my.np) # Fisher's exact test because 1 sample per group
	save(mk.filtered.5mC, mk.5mC.meth, mk.5mC.diff, file=my.filtered.5mC.data)
}

###################################
## 5mC+5hmC (BS) for AGA2F,SGA2F ##
###################################
my.filtered.bs.data<-file.path(top.dir,"MethylKit/AGA2F.SGA2F.BS.RData") 
if(file.exists(my.filtered.bs.data)){
	cat("Loading my.filtered.bs.data...\n")
	load(my.filtered.bs.data)
	cat("mk.filtered.bs, mk.bs.meth, mk.bs.diff loaded\n")
}else{
	BS.2F.file.list=list(
		file.path(top.dir,"BisSNP/A012/SLX-8074.SLX-8080.A012.cpg.filtered.homo.CG.methylkit.txt"), #BS SGA2F
		file.path(top.dir,"BisSNP/A010/SLX-8074.SLX-8080.A010.cpg.filtered.homo.CG.methylkit.txt")  #BS AGA2F
	)
	mk.bs=read(BS.2F.file.list, sample.id=list("SGA2F.bs", "AGA2F.bs"), assembly="hg19",treatment=c(1,0),context="CpG") # is a 'methylRawList'
	# Filter by depth
	mk.filtered.bs<-filterByCoverage(mk.bs,lo.count=my.depth, lo.perc=NULL, hi.count=NULL,hi.perc=NULL) # is a 'methylRawList'
	# unite
	mk.bs.meth<-unite(mk.filtered.bs,destrand=FALSE) # a methylBase object
	# Base-level Difference 
	mk.bs.diff=calculateDiffMeth(mk.bs.meth, num.cores=my.np) # Fisher's exact test because 1 sample per group
	save(mk.filtered.bs, mk.bs.meth, mk.bs.diff, file=my.filtered.bs.data)
}

################################
## Basic Exploratory Analysis ##
################################
# 1. Plot
#getMethylationStats(mk.filtered.5hmC[[1]],plot=T,both.strands=F)
#getMethylationStats(mk.filtered.5hmC[[2]],plot=T,both.strands=F)
getCorrelation(mk.filtered.5hmC.meth, plot=T)

# 2. get top diff
mk.top.dmc=get.methylDiff(mk.5hmC.diff,difference=my.diff,qvalue=my.pvalue) # isa methylDiff (resulution: base)

# 3. Annotate differentially methylated CpG bases 
# over exons/introns/promoters/TSSes
# Resolution: base
# promoter/exon/intron/TSS
top.met.anno=annotate.WithGenicParts(mk.top.dmc, gene.obj) # isa 'annotationByGenicParts'
if(FALSE){
	top.met.anno@annotation        
	top.met.anno@dist.to.TSS       
	top.met.anno@members           
	top.met.anno@no.of.OlapFeat    
	top.met.anno@num.annotation    
	top.met.anno@num.precedence    
	top.met.anno@perc.of.OlapFeat  
	top.met.anno@precedence
}
getTargetAnnotationStats(top.met.anno, percentage=TRUE, precedence=TRUE)
getFeatsWithTargetsStats(top.met.anno, percentage=TRUE)
plotTargetAnnotation(top.met.anno, precedence=TRUE, main="differential methylation annotation")

## 4. get annotations of top.dmc
top.dmc<-getData(mk.top.dmc) # isa data.frame
top.dmc[,c("enst","enst.strand","dist.to.feature")]<-getAssociationWithTSS(top.met.anno)[,c("feature.name","feature.strand","dist.to.feature")] # add a new column (assuming the order or record is the same between dmc & dmc.enst)
fields <- c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol","transcript_biotype") # get the gene names
dummy.enst=getBM(attributes = fields, filters = "ensembl_transcript_id", values = top.dmc$enst, mart = grch37) # is a data.frame
top.dmc=merge(top.dmc, dummy.enst, by.x="enst", by.y="ensembl_transcript_id") # get annotations 
top.dmc=merge(top.dmc, getData(mk.filtered.5hmC.meth), by=c("chr","start","end","strand")) # get number of 5hmC and their coverage
top.dmc<-top.dmc[order(top.dmc$qvalue),] # order by qvalue 
write.csv(head(top.dmc,n=100), file=file.path(top.dir,"top.100.5hmc.site.csv"))

#################################################
# methylKit to GRanges (for the sake of ggbio) ##
#################################################
gr.5hmC.meth<-as(mk.filtered.5hmC.meth,"GRanges")
gr.5mC.meth<-as(mk.5mC.meth,"GRanges")
gr.bs.meth<-as(mk.bs.meth,"GRanges")
#gr.top.dmc<-as(mk.top.dmc,"GRanges")
#gr.5hmC.diff<-as(mk.5hmC.diff,"GRanges")

# 1-a. set GRange of target transcripts (use my.gr.ensg)
#my.enst<-c("ENST00000451024","ENST00000440536","ENST00000604581") # first two from "NKX1-2" (ENSG00000229544) and the last from "RP13-238F13.3" (ENSG00000271434)
#my.enst<-c("ENST00000441415") # RP11-460N20.5 ENSG00000226002 (chr7 locus)
#my.gr.enst<-gr.enst[mcols(gr.enst)$tx_name %in% my.enst] # NKX1-2

# 1-b. set GRange of target genes
my.target<-'NKX1-2'
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
my.ensg=getBM(attributes = c("ensembl_gene_id"), filters = "hgnc_symbol", values=my.target, mart = grch37) # is a data.frame
my.gr.ensg<-gr.ensg[mcols(gr.ensg)$gene_id %in% my.ensg] # GRanges of my.target gene 
my.gr.ext.ensg<-promoters(my.gr.ensg, upstream=my.upflank, downstream=end(my.gr.ensg)-start(my.gr.ensg)+1+my.downflank) # includes 2000 +TSS and 100 +TES 
#my.gr.promoter<-promoters(my.gr.ensg, upstream=my.upflank, downstream=my.downflank) # includes 2000 +TSS and 100 +TSS

# 1-c. set GRange of target region (if no ensg/enst avaible) 
my.target<-'chr7.DMR'
my.chr='chr7'
my.start=64541001 - 100
my.end=  64542000 + 2000
my.strand='*'
my.range=IRanges( start=my.start, end=my.end, names=my.target)
my.gr.ext.ensg <- GRanges(seqnames<-Rle(my.chr) , ranges<-my.range, strand<-Rle(my.strand), mcols<-DataFrame(data.frame(name=my.target)))

# overlap between united 5hmC SGA/AGA methylation data (methylBase) and my.target 
my.overlap=as.matrix( findOverlaps(gr.5hmC.meth, my.gr.ext.ensg) )
gr.5hmC.sga<-gr.5hmC.meth[my.overlap[,"queryHits"]] # for SGA
gr.5hmC.aga<-gr.5hmC.meth[my.overlap[,"queryHits"]] # for AGA
mcols(gr.5hmC.sga)<-DataFrame(with(mcols(gr.5hmC.sga), data.frame(Condition="SGA", seqType="5hmC (BS-oxBS)", methylation=numCs1/coverage1))) # SGA
mcols(gr.5hmC.aga)<-DataFrame(with(mcols(gr.5hmC.aga), data.frame(Condition="AGA", seqType="5hmC (BS-oxBS)", methylation=numCs2/coverage2))) # AGA

# overlap between united 5mC SGA/AGA methylation data (methylBase) and my.target 
my.overlap=as.matrix( findOverlaps(gr.5mC.meth, my.gr.ext.ensg) )
gr.5mC.sga<-gr.5mC.meth[my.overlap[,"queryHits"]] # for SGA
gr.5mC.aga<-gr.5mC.meth[my.overlap[,"queryHits"]] # for AGA
mcols(gr.5mC.sga)<-DataFrame(with(mcols(gr.5mC.sga), data.frame(Condition="SGA", seqType="5mC (oxBS)", methylation=numCs2/coverage2))) # SGA
mcols(gr.5mC.aga)<-DataFrame(with(mcols(gr.5mC.aga), data.frame(Condition="AGA", seqType="5mC (oxBS)", methylation=numCs1/coverage1))) # AGA

# overlap between united BS SGA/AGA methylation data (methylBase) and my.target 
my.overlap=as.matrix( findOverlaps(gr.bs.meth, my.gr.ext.ensg) )
gr.bs.sga<-gr.bs.meth[my.overlap[,"queryHits"]] # for SGA
gr.bs.aga<-gr.bs.meth[my.overlap[,"queryHits"]] # for AGA
mcols(gr.bs.sga)<-DataFrame(with(mcols(gr.bs.sga), data.frame(Condition="SGA", seqType="5mC+5hmC (BS)", methylation=numCs1/coverage1))) # SGA
mcols(gr.bs.aga)<-DataFrame(with(mcols(gr.bs.aga), data.frame(Condition="AGA", seqType="5mC+5hmC (BS)", methylation=numCs2/coverage2))) # AGA

# combine data
my.5hmC.meth=c(gr.5hmC.sga, gr.5hmC.aga) # combine 5hmC SGA and AGA
my.5mC.meth=c(gr.5mC.sga, gr.5mC.aga) # combine oxBS SGA and AGA
my.bs.meth=c(gr.bs.sga, gr.bs.aga) # combine BS SGA and AGA
my.sga.meth=c(gr.5hmC.sga, gr.5mC.sga, gr.bs.sga)
my.aga.meth=c(gr.5hmC.aga, gr.5mC.aga, gr.bs.aga)

# overlap between cpgi and my.target 
my.overlap=as.matrix( findOverlaps(cpg.obj[["CpGi"]], my.gr.ext.ensg) )
my.cpgi<-cpg.obj$CpGi[my.overlap[,"queryHits"]]

#p1<-autoplot(hg19ensGene, which=my.gr.enst, fill = "orange", color = "grey") # transcript track for my.target
p1<-autoplot(hg19ensGene, which=my.gr.ext.ensg, fill = "orange", color = "grey") # transcript track for my.target
p2<-autoplot(my.cpgi, fill = "darkgreen", color = "grey") # cpgi track for my.target

p3<-ggplot(my.5hmC.meth) + geom_point(aes(x=start(my.5hmC.meth), y=methylation, col=Condition)) + geom_smooth(aes(x=start(my.5hmC.meth), y=methylation, col=Condition)) # BS-oxBS (5hmC) % methylation
p4<-ggplot(my.5mC.meth) + geom_point(aes(x=start(my.5mC.meth), y=methylation, col=Condition)) + geom_smooth(aes(x=start(my.5mC.meth), y=methylation, col=Condition))    # oxBS (5mC)  % methylation
p5<-ggplot(my.bs.meth) + geom_point(aes(x=start(my.bs.meth), y=methylation, col=Condition)) + geom_smooth(aes(x=start(my.bs.meth), y=methylation, col=Condition))       # BS (5mC+5hmC) % methylation

p6<-ggplot(my.sga.meth) + geom_point(aes(x=start(my.sga.meth), y=methylation, col=seqType)) + geom_smooth(aes(x=start(my.sga.meth), y=methylation, col=seqType))       # BS (5mC+5hmC) % methylation
p7<-ggplot(my.aga.meth) + geom_point(aes(x=start(my.aga.meth), y=methylation, col=seqType)) + geom_smooth(aes(x=start(my.aga.meth), y=methylation, col=seqType))       # BS (5mC+5hmC) % methylation

p8<-ggplot(gr.bs.sga) + geom_point(aes(x=start(gr.bs.sga), y=methylation, col=Condition)) + geom_smooth(aes(x=start(gr.bs.sga), y=methylation, col=Condition))       # BS (5mC+5hmC) % methylation
p9<-ggplot(gr.5mC.sga) + geom_point(aes(x=start(gr.5mC.sga), y=methylation, col=Condition)) + geom_smooth(aes(x=start(gr.5mC.sga), y=methylation, col=Condition))       # BS (5mC+5hmC) % methylation
p10<-ggplot(gr.5hmC.sga) + geom_point(aes(x=start(gr.5hmC.sga), y=methylation, col=Condition)) + geom_smooth(aes(x=start(gr.5hmC.sga), y=methylation, col=Condition))       # BS (5mC+5hmC) % methylation

p11<-ggplot(gr.bs.aga) + geom_point(aes(x=start(gr.bs.aga), y=methylation, col=Condition)) + geom_smooth(aes(x=start(gr.bs.aga), y=methylation, col=Condition))       # BS (5mC+5hmC) % methylation
p12<-ggplot(gr.5mC.aga) + geom_point(aes(x=start(gr.5mC.aga), y=methylation, col=Condition)) + geom_smooth(aes(x=start(gr.5mC.aga), y=methylation, col=Condition))       # BS (5mC+5hmC) % methylation
p13<-ggplot(gr.5hmC.aga) + geom_point(aes(x=start(gr.5hmC.aga), y=methylation, col=Condition)) + geom_smooth(aes(x=start(gr.5hmC.aga), y=methylation, col=Condition))       # BS (5mC+5hmC) % methylation

my.filename=paste0("~/results/RnBeads/5hmC.SGA.AGA/SGA.AGA.bs.oxbs.",my.target,".profile.by.condition.new")
tiff(filename=paste0(my.filename,".tiff"),width=16,height=9,units="in",res=300, compression = 'lzw')
tracks(transcripts = p1, cpgi= p2, `BS (5mc+5hmC)`=p5, `oxBS (5mC)`=p4, `BS - oxBS (5hmC)`=p3, heights = c(3,1,4,4,4)) # + theme_alignment(grid=FALSE, border = FALSE) # plot ggbio
dev.off()

my.filename=paste0("~/results/RnBeads/5hmC.SGA.AGA/SGA.AGA.bs.oxbs.",my.target,".profile.by.assay.new")
tiff(filename=paste0(my.filename,".tiff"),width=16,height=9,units="in",res=300, compression = 'lzw')
tracks(transcripts = p1, cpgi= p2, `SGA`=p6, `AGA`=p7, heights = c(3,1,4,4)) # + theme_alignment(grid=FALSE, border = FALSE) # plot ggbio
dev.off()

my.filename=paste0("~/results/RnBeads/5hmC.SGA.AGA/SGA.AGA.bs.oxbs.",my.target,".profile.SGA.only.new")
tiff(filename=paste0(my.filename,".tiff"),width=16,height=9,units="in",res=300, compression = 'lzw')
tracks(transcripts = p1, cpgi= p2, `BS (5mc+5hmC)`=p8, `oxBS (5mC)`=p9, `BS - oxBS (5hmC)`=p10, heights = c(3,1,4,4,4)) # + theme_alignment(grid=FALSE, border = FALSE) # plot ggbio
dev.off()

#######################################
# Regional Analysis                   #
# promoter as a unit of analysis      # 
# rather than individual methylated-C #
# Resolution: region (promoter)       #
#######################################
mk.promoters=regionCounts(mk.filtered.5hmC, gene.obj$promoters, cov.bases=1, strand.aware=TRUE) # is a methylRawList object with length(filtered.myobj) methylRaw objects
#mk.promoters=regionCounts(mk.filtered.5hmC, gr.promoters, cov.bases=1, strand.aware=TRUE) # is a methylRawList object with length(filtered.myobj) methylRaw objects
mk.promoters.meth=unite(mk.promoters,destrand=FALSE) # is a 'methylBase'
# Difference 
mk.promoters.diff=calculateDiffMeth(mk.promoters.meth, num.cores=my.np) # Fisher's exact test because 1 sample per group
# top diff
mk.promoters.top.dmr=get.methylDiff(mk.promoters.diff,difference=my.diff,qvalue=my.pvalue) # isa methylDiff (resolution: region)
# methylKit to GRanges
gr.promoters.meth<-as(mk.promoters.meth,"GRanges")
gr.promoters.diff<-as(mk.promoters.diff,"GRanges")
gr.promoters.top.dmr<-as(mk.promoters.top.dmr,"GRanges")

#promoters.top.dmr=getData(mk.promoters.top.dmr) 
# visualize the distribution of hypo/hyper-methylated bases/regions per chromosome
#diffMethPerChr(promoters.top.dmr,plot=F, qvalue.cutoff=myQvalue, meth.cutoff=myDiffvalue)
#diffMethPerChr(promoters.top.dmr,plot=T, qvalue.cutoff=myQvalue, meth.cutoff=myDiffvalue)

# overlap between top diff promoters and the promoter gr 
my.overlap=as.matrix( findOverlaps(gr.promoters.top.dmr, gr.promoters, minoverlap=1990L) )
gr.promoters.top.dmr<-gr.promoters.top.dmr[my.overlap[,"queryHits"]]
mcols(gr.promoters.top.dmr)<-DataFrame(mcols(gr.promoters.top.dmr), mcols(gr.promoters[my.overlap[,"subjectHits"]]))

## get gene name of gr.promoters.top.dmr 
top.prompter.dmr<-as.data.frame(gr.promoters.top.dmr) # isa data.frame
fields <- c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol","transcript_biotype") # get the gene names
dummy.enst=getBM(attributes = fields, filters = "ensembl_transcript_id", values = top.promoter.dmr$tx_name, mart = grch37) # is a data.frame
top.promoter.dmr=merge(top.promoter.dmr, dummy.enst, by.x="tx_name", by.y="ensembl_transcript_id") # get gene names 
unique.top.promoter.dmr<-unique(top.promoter.dmr[,c("seqnames","start","end","strand","qvalue","meth.diff","hgnc_symbol","transcript_biotype")])
unique.top.promoter.dmr<-unique.top.promoter.dmr[order(unique.top.promoter.dmr$qvalue),] # order by qvalue 
write.csv(head(unique.top.promoter.dmr,n=100), file=file.path(top.dir,"top.100.5hmc.promoters.csv"))


if(FALSE){
	################################################
	# Annotate differentially methylated CpG bases #
	# over CPGi
	# Resolution: base
	################################################
	cpg.obj=read.feature.flank(location=my.cgi.file, feature.flank.name=c("CpGi","shores")) # a ‘GenomicRangesList’ contatining one GRanges object for flanks and one for GRanges object for the main feature.
	met.cgi.anno=annotate.WithFeature.Flank(mk.top.dmc, cpg.obj$CpGi, cpg.obj$shores, feature.name="CpGi",flank.name="shores")
	met.cgi.anno

	getTargetAnnotationStats(met.cgi.anno, percentage=TRUE, precedence=TRUE)
	getFeatsWithTargetsStats(met.cgi.anno, percentage=TRUE)
	plotTargetAnnotation(met.cgi.anno,col=c("green","gray","white"), main="differential methylation annotation")

	# USEr-defined annotations here
	#diffmyAnn=annotate.WithFeature(mk.top.dmc, my.anno.obj, feature.name="my_feature")


	# CPGi
	promoters.cgi.anno=annotate.WithFeature.Flank(promoters.top.dmr, cpg.obj$CpGi, cpg.obj$shores, feature.name="CpGi",flank.name="shores")
}

dev.off()
cat("All done\n")
