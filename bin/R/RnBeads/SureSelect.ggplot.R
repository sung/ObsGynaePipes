#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

library(methylKit) # for 'read.feature.flank'
library(reshape2) # for 'melt' and 'dcast'
library(GenomicRanges)
library(gplots) # heatmap.2
library(rtracklayer)
library(ggbio)

############
## Config ##
############
TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R")
source("~/Pipelines/config/graphic.R") # cbPalette
source("~/Pipelines/bin/R/RnBeads/local.R") # defines 'prep.sureselect' functiton

my.project="SureSelect"
my.slx<-"SLX-10409.v3"

my.upflank=2000 # up from TSS
my.downflank=200 # down form TES

min.doc <- 10 # depth of coverage
min.missing.sample.ratio <- 0   # similar with filtering.missing.value.quantile of RnBeads
							# 0  : does not allow missing sample 
							# 0.2: allow 20% missing samples
							# 1  : allow any missing  
#######################################
### WGBS Methylation from BED Files ###
#######################################
			#`NKX1-2`="~/Pipelines/data/NKX1-2/NKX1-2.up.1500.down.500.grch37.WGBS.bed.txt",
			#`MAB21L1`="~/Pipelines/data/MAB21L1/MAB21L1.top2.grch37.WGBS.bed.txt", 
			#`chr7.DMR`="~/Pipelines/data/Chr7.DMR/chr7.top2.grch37.WGBS.bed.txt",
			#`VIPR2`="~/Pipelines/data/VIPR2/VIPR2.5kb.grch37.WGBS.bed.txt",
			#`C1orf216`="~/Pipelines/data/C1orf216/C1orf216.up.1500.down.500.grch37.WGBS.bed.txt",
			#`TICAM2`="~/Pipelines/data/TICAM2/TICAM2.grch37.bed.txt"
dummy<-c(
			`CSMD1`="~/Pipelines/data/CSMD1/CSMD1.225kb.grch37.WGBS.bed.txt",
			`CSMD1.promo`="~/Pipelines/data/CSMD1/CSMD1.putative.promoter.4K.grch37.WGBS.bed.txt"
		)
dummy<-lapply(dummy, read.csv, stringsAsFactors=FALSE) # isa 'list'

gl.wgbs.locus<-list() # a list of GRanges
df.wgbs.locus<-list() # a list of data.frame
for(i in names(dummy)){
	gl.dummy<-lapply(dummy[[i]]$BedFile, rtracklayer::import.bed) 
	names(gl.dummy)<-dummy[[i]]$SampleName # e.g. AGA2F SGA2F
	gl.dummy.new<-lapply(names(gl.dummy), function(i) 
			GRanges(
				seqnames=seqnames(gl.dummy[[i]]), 
				ranges<-IRanges(start=start(gl.dummy[[i]]),end=end(gl.dummy[[i]])),
				strand=strand(gl.dummy[[i]]), 
				mcols<-DataFrame(with(mcols(gl.dummy[[i]]), data.frame(SampleName=i, Condition=ifelse(grepl("AGA",i), "AGA", "SGA"), Gender=ifelse(grepl("M",i), "Male", "Female"), score=as.numeric(name), coverage=score))) # score in percentage methylation
				)
			) # isa GRangeList
	gr.dummy <- do.call(c, gl.dummy.new) # isa 'GRanges'
	gr.dummy <- gr.dummy[mcols(gr.dummy)$coverage>=min.doc,] # min 10 depth  
	# filter CpG called by all samples (dcast from 'reshape2')
	df.dummy <- dcast(as.data.frame(gr.dummy), start ~ SampleName, value.var="score"); rownames(df.dummy) <- df.dummy$start; df.dummy$start<-NULL; 
	# CpG called across all samples
	gl.wgbs.locus[[i]] <- gr.dummy[start(gr.dummy) %in% rownames(df.dummy[rowSums(!is.na(df.dummy))==ncol(df.dummy),]),]
	df.wgbs.locus[[i]] <- df.dummy
}

#############################################
### SureSelect Methylation from BED Files ###
### Read the following first              ###
### ~/Pipelines/data/MJ.SureSelect/README ###
### bin/sureselect.on.target.sh           ###
#############################################
##########################
## Regions of Interests ##
##########################
# genomic region of interests (target)
#my.target='NKX1-2'
#my.assay <- 'Set1' # either 'Set1': original samples or 'Set2': validation samples
#my.targets<-c("NKX1-2","MAB21L1","chr7.DMR","CSMD1","VIPR2","C1orf216")
if(FALSE){
	my.targets<-c("CSMD1")
	my.assays<-c("Set1","Set2") # Set1: technical replicates (same samples); Set2: biological replicates (independent cohorts)
	for(my.target in my.targets){ 
		for(my.assay in my.assays){
			if(my.assay=="Set2" & my.target=="chr7.DMR"){
				min.missing.sample.ratio <- 0.5
			}
			#list of (`gl`=gl.sureselect.locus, `df`=df.sureselect.locus, `gender`=boy.girl.cpg)
			sureselect.v1<-prep.sureselect("SureSelect.v1",my.assay)
			sureselect.v2<-prep.sureselect("SureSelect.v2",my.assay)
			sureselect.v3<-prep.sureselect("SureSelect.v3",my.assay)

			#my.chr<-'chr10' #'chr7' #'chr13'
			if(my.target=="NKX1-2" || my.target=="MAB21L1"){ # Targets from SGA/AGA study
				my.contrast<-"Condition"
				my.col=c(`AGA`="#00BFC4", `SGA`="#F8766D")
				my.ensg=getBM(attributes = c("ensembl_gene_id"), filters = "hgnc_symbol", values=my.target, mart = myMart) # is a data.frame
				my.gr.ensg<-gr.ensg[mcols(gr.ensg)$gene_id %in% my.ensg] # GRanges of my.target gene (gr.ensg defined in 'Annotation.R') 
			}else{
				my.contrast<-"Gender"
				my.col<-c(`Female`="#00BFC4", `Male`="#F8766D")
				if(my.target=='chr7.DMR'){
					my.gr.ensg <- rtracklayer::import.bed("~/Pipelines/data/Chr7.DMR/chr7.top2.bed")
				}else if(my.target=='CSMD1'){
					my.gr.ensg <- rtracklayer::import.bed("~/Pipelines/data/CSMD1/CSMD1.boy.girl.dmr.bed")
				}else if(my.target=='CSMD1.promo'){
					my.gr.ensg <- rtracklayer::import.bed("~/Pipelines/data/CSMD1/CSMD1.placenta.putative.promoter.GRCh37.bed")
				}else if(my.target=='VIPR2'){
					my.gr.ensg <- rtracklayer::import.bed("~/Pipelines/data/VIPR2/VIPR2.top5k.grch37.bed")
				}else{
					my.ensg=getBM(attributes = c("ensembl_gene_id"), filters = "hgnc_symbol", values=my.target, mart = myMart) # is a data.frame
					my.gr.ensg<-gr.ensg[mcols(gr.ensg)$gene_id %in% my.ensg] # GRanges of my.target gene (gr.ensg defined in 'Annotation.R') 
				}
				my.gr.ensg <- reduce(my.gr.ensg)
				names(my.gr.ensg)<-my.target
			}
			my.gr.ext.ensg<-promoters(my.gr.ensg, upstream=my.upflank, downstream=end(my.gr.ensg)-start(my.gr.ensg)+1+my.downflank) # includes 2000 +TSS and 100 +TES 

			# overlap between cpgi and my.target 
			cpg.obj=read.feature.flank(location=my.cgi.file, feature.flank.name=c("CpGi","shores")) # a ‘GenomicRangesList’ contatining one GRanges object for flanks and one for GRanges object for the main feature.
			my.overlap=as.matrix( findOverlaps(cpg.obj[["CpGi"]], my.gr.ext.ensg) )
			my.cpgi <-cpg.obj$CpGi[my.overlap[,"queryHits"]]

			##############
			# chromosome #
			##############
			p.chr <- Ideogram(genome = "hg19", xlabel=TRUE) + xlim(ranges(my.gr.ext.ensg))

			####################
			# transcript track # 
			####################
			p.tr<-autoplot(hg19ensGene, which=my.gr.ext.ensg, fill = "orange", color = "grey") # hg19ensGene defined in 'Annotation.R'

			############################
			# cpgi track for my.target #
			############################
			p.cpgi <-autoplot(my.cpgi, fill = "darkgreen", color = "grey") 
			#p.cpgi <-ggplot(as.data.frame(my.cpgi)) + geom <- rect(mapping=aes(xmin=start,xmax=end,ymin=1,ymax=2), fill="darkgreen") + theme(axis.text.y=element <- blank(), axis.ticks.y=element <- blank())

			##########################
			# WGBS methylation GGBIO #
			##########################
			if(my.target=="NKX1-2" || my.target=="MAB21L1"){ # Targets from SGA/AGA study
				p.wgbs.locus <-ggplot(gl.wgbs.locus[[my.target]], aes_string(x="start",y="score",col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)
			}else{
				select<-mcols(gl.wgbs.locus[[my.target]])$Condition=="AGA"; my.gr.wgbs<-as.data.frame(gl.wgbs.locus[[my.target]][select,])
				p.wgbs.locus.aga<-ggplot(my.gr.wgbs, aes_string(x="start",y="score",col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)

				select<-mcols(gl.wgbs.locus[[my.target]])$Condition=="SGA"; my.gr.wgbs<-as.data.frame(gl.wgbs.locus[[my.target]][select,])
				p.wgbs.locus.sga<-ggplot(my.gr.wgbs, aes_string(x="start",y="score",col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)

				## NO. of CpG & P-value
				if(FALSE){
					dummy<-with(my.gr.wgbs, data.frame(start, SampleName, Gender, Cs=round(coverage*score/100), Ts=round(coverage*(1-score/100)), coverage))
					## NO. of CpG
					length(unique(dummy$start))
					## P-value
					my.mat<-matrix(c(sum(dummy[dummy$Gender=="Female", "Cs"]), sum(dummy[dummy$Gender=="Female", "Ts"]), sum(dummy[dummy$Gender=="Male", "Cs"]), sum(dummy[dummy$Gender=="Male", "Ts"])),nrow=2)
					#> my.mat
					#      [,1]  [,2]
					#[1,] 21306 38682
					#[2,] 24185 15558
					## P-value
					fisher.test(my.mat)$p.value
				}
			}
			################################
			# SureSelect methylation GGBIO #
			################################
			if(my.target=="NKX1-2" || my.target=="MAB21L1"){ # Targets from SGA/AGA study
				#################
				# sureselect.v1 #
				#################
				p.sureselect.v1.locus<-ggplot(sureselect.v1[["gl"]][[my.target]], aes(x=start(sureselect.v1[["gl"]][[my.target]]),y=score,col=Condition)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)
				p.sureselect.v1.locus.by.sample<-ggplot(sureselect.v1[["gl"]][[my.target]], aes(x=start(sureselect.v1[["gl"]][[my.target]]),y=score,col=SampleName)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				select.v1.aga<-mcols(sureselect.v1[["gl"]][[my.target]])$Condition=="AGA"
				p.sureselect.v1.locus.aga.by.sample<-ggplot(sureselect.v1[["gl"]][[my.target]][select.v1.aga,], aes(x=start(sureselect.v1[["gl"]][[my.target]][select.v1.aga,]),y=score,col=SampleName)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				select.v1.sga<-mcols(sureselect.v1[["gl"]][[my.target]])$Condition=="SGA"
				p.sureselect.v1.locus.sga.by.sample<-ggplot(sureselect.v1[["gl"]][[my.target]][select.v1.sga,], aes(x=start(sureselect.v1[["gl"]][[my.target]][select.v1.sga,]),y=score,col=SampleName)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				#################
				# sureselect.v2 #
				#################
				p.sureselect.v2.locus<-ggplot(sureselect.v2[["gl"]][[my.target]], aes(x=start(sureselect.v2[["gl"]][[my.target]]),y=score,col=Condition)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)
				p.sureselect.v2.locus.by.sample<-ggplot(sureselect.v2[["gl"]][[my.target]], aes(x=start(sureselect.v2[["gl"]][[my.target]]),y=score,col=SampleName)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				select.v2.aga<-mcols(sureselect.v2[["gl"]][[my.target]])$Condition=="AGA"
				p.sureselect.v2.locus.aga.by.sample<-ggplot(sureselect.v2[["gl"]][[my.target]][select.v2.aga,], aes(x=start(sureselect.v2[["gl"]][[my.target]][select.v2.aga,]),y=score,col=SampleName)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				select.v2.sga<-mcols(sureselect.v2[["gl"]][[my.target]])$Condition=="SGA"
				p.sureselect.v2.locus.sga.by.sample<-ggplot(sureselect.v2[["gl"]][[my.target]][select.v2.sga,], aes(x=start(sureselect.v2[["gl"]][[my.target]][select.v2.sga,]),y=score,col=SampleName)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)
				#################
				# sureselect.v3 #
				#################
				p.sureselect.v3.locus<-ggplot(sureselect.v3[["gl"]][[my.target]], aes(x=start(sureselect.v3[["gl"]][[my.target]]),y=score,col=Condition)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)
				p.sureselect.v3.locus.by.sample<-ggplot(sureselect.v3[["gl"]][[my.target]], aes(x=start(sureselect.v3[["gl"]][[my.target]]),y=score,col=SampleName)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				select.v3.aga<-mcols(sureselect.v3[["gl"]][[my.target]])$Condition=="AGA"
				p.sureselect.v3.locus.aga.by.sample<-ggplot(sureselect.v3[["gl"]][[my.target]][select.v3.aga,], aes(x=start(sureselect.v3[["gl"]][[my.target]][select.v3.aga,]),y=score,col=SampleName)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				select.v3.sga<-mcols(sureselect.v3[["gl"]][[my.target]])$Condition=="SGA"
				p.sureselect.v3.locus.sga.by.sample<-ggplot(sureselect.v3[["gl"]][[my.target]][select.v3.sga,], aes(x=start(sureselect.v3[["gl"]][[my.target]][select.v3.sga,]),y=score,col=SampleName)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				#################
				# sureselect.v3 #
				# plot pair-wise#
				#################
				this.pair<-sureselect.v3[["gl"]][[my.target]]$SampleName %in% colnames(sureselect.v3[["df"]][[my.target]])[1:2]
				p.14M=ggplot(sureselect.v3[["gl"]][[my.target]][this.pair,], aes(x=start(sureselect.v3[["gl"]][[my.target]][this.pair,]),y=score,col=Condition)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)

				this.pair<-sureselect.v3[["gl"]][[my.target]]$SampleName %in% colnames(sureselect.v3[["df"]][[my.target]])[3:4]
				p.16F=ggplot(sureselect.v3[["gl"]][[my.target]][this.pair,], aes(x=start(sureselect.v3[["gl"]][[my.target]][this.pair,]),y=score,col=Condition)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)

				this.pair<-sureselect.v3[["gl"]][[my.target]]$SampleName %in% colnames(sureselect.v3[["df"]][[my.target]])[5:6]
				p.19M=ggplot(sureselect.v3[["gl"]][[my.target]][this.pair,], aes(x=start(sureselect.v3[["gl"]][[my.target]][this.pair,]),y=score,col=Condition)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)

				this.pair<-sureselect.v3[["gl"]][[my.target]]$SampleName %in% colnames(sureselect.v3[["df"]][[my.target]])[7:8]
				p.20M=ggplot(sureselect.v3[["gl"]][[my.target]][this.pair,], aes(x=start(sureselect.v3[["gl"]][[my.target]][this.pair,]),y=score,col=Condition)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)

				this.pair<-sureselect.v3[["gl"]][[my.target]]$SampleName %in% colnames(sureselect.v3[["df"]][[my.target]])[9:10]
				p.21F=ggplot(sureselect.v3[["gl"]][[my.target]][this.pair,], aes(x=start(sureselect.v3[["gl"]][[my.target]][this.pair,]),y=score,col=Condition)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)

				this.pair<-sureselect.v3[["gl"]][[my.target]]$SampleName %in% colnames(sureselect.v3[["df"]][[my.target]])[11:12]
				p.22F=ggplot(sureselect.v3[["gl"]][[my.target]][this.pair,], aes(x=start(sureselect.v3[["gl"]][[my.target]][this.pair,]),y=score,col=Condition)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=my.col)
			}else{
				#################
				# sureselect.v1 #
				#################
				# sureselect.v1 AGA
				select.v1.aga<-start(sureselect.v1[["gl"]][[my.target]]) %in% sureselect.v1[["gender"]][[my.target]][["AGA"]] & mcols(sureselect.v1[["gl"]][[my.target]])$Condition=="AGA"
				my.gr.sureselect<-as.data.frame(sureselect.v1[["gl"]][[my.target]][select.v1.aga,])
				p.sureselect.v1.locus.aga<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)
				p.sureselect.v1.locus.aga.by.sample<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col="SampleName")) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				# sureselect.v1 SGA
				select.v1.sga<-start(sureselect.v1[["gl"]][[my.target]]) %in% sureselect.v1[["gender"]][[my.target]][["SGA"]] & mcols(sureselect.v1[["gl"]][[my.target]])$Condition=="SGA"
				my.gr.sureselect<-as.data.frame(sureselect.v1[["gl"]][[my.target]][select.v1.sga,])
				p.sureselect.v1.locus.sga<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)
				p.sureselect.v1.locus.sga.by.sample<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col="SampleName")) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				#################
				# sureselect.v2 #
				#################
				# sureselect.v2 AGA
				select.v2.aga<-start(sureselect.v2[["gl"]][[my.target]]) %in% sureselect.v2[["gender"]][[my.target]][["AGA"]] & mcols(sureselect.v2[["gl"]][[my.target]])$Condition=="AGA"
				my.gr.sureselect<-as.data.frame(sureselect.v2[["gl"]][[my.target]][select.v2.aga,])
				p.sureselect.v2.locus.aga<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)
				p.sureselect.v2.locus.aga.by.sample<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col="SampleName")) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				# sureselect.v2 SGA
				select.v2.sga<-start(sureselect.v2[["gl"]][[my.target]]) %in% sureselect.v2[["gender"]][[my.target]][["SGA"]] & mcols(sureselect.v2[["gl"]][[my.target]])$Condition=="SGA"
				my.gr.sureselect<-as.data.frame(sureselect.v2[["gl"]][[my.target]][select.v2.sga,])
				p.sureselect.v2.locus.sga<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)
				p.sureselect.v2.locus.sga.by.sample<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col="SampleName")) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				#################
				# sureselect.v3 #
				#################
				# sureselect.v3 AGA
				select.v3.aga<-start(sureselect.v3[["gl"]][[my.target]]) %in% sureselect.v3[["gender"]][[my.target]][["AGA"]] & mcols(sureselect.v3[["gl"]][[my.target]])$Condition=="AGA"
				my.gr.sureselect<-as.data.frame(sureselect.v3[["gl"]][[my.target]][select.v3.aga,])
				p.sureselect.v3.locus.aga<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)
				p.sureselect.v3.locus.aga.by.sample<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col="SampleName")) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				# sureselect.v3 SGA
				select.v3.sga<-start(sureselect.v3[["gl"]][[my.target]]) %in% sureselect.v3[["gender"]][[my.target]][["SGA"]] & mcols(sureselect.v3[["gl"]][[my.target]])$Condition=="SGA"
				my.gr.sureselect<-as.data.frame(sureselect.v3[["gl"]][[my.target]][select.v3.sga,])
				p.sureselect.v3.locus.sga<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col=my.contrast)) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE,) + scale_colour_manual(values=my.col)
				p.sureselect.v3.locus.sga.by.sample<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col="SampleName")) + geom_point(alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette)

				## NO. of CpG & P-value
				if(FALSE){
					dummy<-with(as.data.frame(sureselect.v1[["gl"]][[my.target]][select.v1.aga,]), data.frame(start, SampleName, Gender, Cs=round(coverage*score/100), Ts=round(coverage*(1-score/100)), coverage))
					dummy<-with(as.data.frame(sureselect.v2[["gl"]][[my.target]][select.v2.aga,]), data.frame(start, SampleName, Gender, Cs=round(coverage*score/100), Ts=round(coverage*(1-score/100)), coverage))
					dummy<-with(as.data.frame(sureselect.v3[["gl"]][[my.target]][select.v3.aga,]), data.frame(start, SampleName, Gender, Cs=round(coverage*score/100), Ts=round(coverage*(1-score/100)), coverage))
					## NO. of CpG
					length(unique(dummy$start))
					my.mat<-matrix(c(sum(dummy[dummy$Gender=="Female", "Cs"]), sum(dummy[dummy$Gender=="Female", "Ts"]), sum(dummy[dummy$Gender=="Male", "Cs"]), sum(dummy[dummy$Gender=="Male", "Ts"])),nrow=2)
					#Set1
					#> my.mat
					#     [,1] [,2]
					#[1,] 5452 8693
					#[2,] 3625 4392
					#Set2
					#> my.mat
					#      [,1]  [,2]
					#[1,] 11310 11823
					#[2,]  8766  4890
					## P-value
					fisher.test(my.mat)$p.value
				}
			}

			##############################################
			# PLOT WGBS & SureSelect methylation profile #
			##############################################
			# 1. Growth Restriction (aga/sga) analysis
			if(my.target=="NKX1-2" || my.target=="MAB21L1"){ # Targets from SGA/AGA study
				# extended transcript region
				out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect",my.target,my.assay,"ggbio",sep="."))
				tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
				print(tracks(Transcripts = p.tr, CPGi=p.cpgi, WGBS=p.wgbs.locus, SureSelect.v1=p.sureselect.v1.locus, SureSelect.v2=p.sureselect.v2.locus, SureSelect.v3=p.sureselect.v3.locus, heights = c(1.3,0.5,3,3,3,3)) + xlim(ranges(my.gr.ext.ensg)))
				dev.off()

				# 1:Tr, 2:CPGi, 3:WGBS, 4:SureSelect, 5:SureSelect.By.Sample, 5:SureSelect.SGA.By.Sample 
				out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect.by.sample",my.target,my.assay,"ggbio",sep="."))
				tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
				print(tracks(Transcripts = p.tr, CPGi=p.cpgi, WGBS=p.wgbs.locus, `SureSelect.v1`=p.sureselect.v1.locus.by.sample, `SureSelect.v2`=p.sureselect.v2.locus.by.sample, `SureSelect.v3`=p.sureselect.v3.locus.by.sample, heights = c(1.3,0.5,3,3,3,3)) + xlim(ranges(my.gr.ext.ensg)))
				dev.off()

				# per pair-wise
				out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect.by.pair",my.target,my.assay,"ggbio",sep="."))
				tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
				if(my.assay=="Set2"){
					print(tracks(Transcripts = p.tr, CPGi=p.cpgi, WGBS=p.wgbs.locus, `14M`=p.14M, `16F`=p.16F, `19M`=p.19M, `20M`=p.20M, `21F`=p.21F, `22F`=p.22F,  heights = c(1.5,0.7,3,3,3,3,3,3,3)) + xlim(ranges(my.gr.ext.ensg)))
				}
				dev.off()

				# zoom-in MAB21L1
				if(my.target=="MAB21L1"){
					# extended transcript region
					out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect",my.target,my.assay,"ggbio.zoom",sep="."))
					tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
					print(tracks(Transcripts = p.tr, CPGi=p.cpgi, WGBS=p.wgbs.locus, SureSelect.v1=p.sureselect.v1.locus, SureSelect.v2=p.sureselect.v2.locus, SureSelect.v3=p.sureselect.v3.locus, heights = c(1.3,0.5,3,3,3,3)) + xlim(reduce(ranges(my.cpgi[1,]) + 300)))
					dev.off()

					# by sample
					out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect.by.sample",my.target,my.assay,"ggbio.zoom",sep="."))
					tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
					print(tracks(Transcripts = p.tr, CPGi=p.cpgi, WGBS=p.wgbs.locus, `SureSelect.v1`=p.sureselect.v1.locus.by.sample, `SureSelect.v2`=p.sureselect.v2.locus.by.sample, `SureSelect.v3`=p.sureselect.v3.locus.by.sample, heights = c(1.3,0.5,3,3,3,3)) + xlim(reduce(ranges(my.cpgi[1,]) + 300)))
					dev.off()

					# per pair-wise
					out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect.by.pair",my.target,my.assay,"ggbio.zoom",sep="."))
					tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
					if(my.assay=="Set2"){
						print(tracks(Transcripts = p.tr, CPGi=p.cpgi, WGBS=p.wgbs.locus, `14M`=p.14M, `16F`=p.16F, `19M`=p.19M, `20M`=p.20M, `21F`=p.21F, `22F`=p.22F,  heights = c(1.5,0.7,3,3,3,3,3,3,3)) + xlim(reduce(ranges(my.cpgi[1,]) + 300))) 
					}
					dev.off()
				}
			# 2. Gender analysis
			}else{
				###############
				# AGA Samples #
				###############
				out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect.AGA",my.target,my.assay,"ggbio",sep="."))
				tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
				if(my.target=="chr7.DMR"){
					print(tracks(`CPGi`=p.cpgi, `WGBS (AGA)`=p.wgbs.locus.aga, `SureSelect.v1 (AGA)`=p.sureselect.v1.locus.aga, `SureSelect.v2 (AGA)`=p.sureselect.v2.locus.aga, `SureSelect.v3 (AGA)`=p.sureselect.v3.locus.aga, heights = c(0.5,3,3,3,3), xlim=my.gr.ensg))
				}else{
					print(tracks(Transcripts=p.tr, `CPGi`=p.cpgi, `WGBS (AGA)`=p.wgbs.locus.aga, `SureSelect.v1 (AGA)`=p.sureselect.v1.locus.aga, `SureSelect.v2 (AGA)`=p.sureselect.v2.locus.aga, `SureSelect.v3 (AGA)`=p.sureselect.v3.locus.aga, heights = c(2,0.5,3,3,3,3), xlim=my.gr.ext.ensg))
				}
				dev.off()

				# colour by sample
				out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect.AGA.by.sample",my.target,my.assay,"ggbio",sep="."))
				tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
				print(tracks(CPGi=p.cpgi, `WGBS (AGA)`=p.wgbs.locus.aga, `SureSelect.v1 (AGA)`=p.sureselect.v1.locus.aga.by.sample, `SureSelect.v2 (AGA)`=p.sureselect.v2.locus.aga.by.sample, `SureSelect.v3 (AGA)`=p.sureselect.v3.locus.aga.by.sample, heights = c(0.5,3,3,3,3)) + xlim(ranges(my.gr.ext.ensg)))
				dev.off()

				###############
				# SGA samples #
				###############
				out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect.SGA",my.target,my.assay,"ggbio",sep="."))
				tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
				if(my.target=="chr7.DMR"){
					print(tracks(`CPGi`=p.cpgi, `WGBS (SGA)`=p.wgbs.locus.sga, `SureSelect.v1 (SGA)`=p.sureselect.v1.locus.sga, `SureSelect.v2 (SGA)`=p.sureselect.v2.locus.sga, `SureSelect.v3 (SGA)`=p.sureselect.v3.locus.sga, heights = c(0.5,3,3,3,3), xlim=my.gr.ensg))
				}else{
					print(tracks(Transcripts=p.tr, `CPGi`=p.cpgi, `WGBS (SGA)`=p.wgbs.locus.sga, `SureSelect.v1 (SGA)`=p.sureselect.v1.locus.sga, `SureSelect.v2 (SGA)`=p.sureselect.v2.locus.sga, `SureSelect.v3 (SGA)`=p.sureselect.v3.locus.sga, heights = c(2,0.5,3,3,3,3), xlim=my.gr.ext.ensg))
				}
				dev.off()

				# colour by sample
				out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect.SGA.by.sample",my.target,my.assay,"ggbio",sep="."))
				tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
				print(tracks(CPGi=p.cpgi, `WGBS (SGA)`=p.wgbs.locus.sga, `SureSelect.v1 (SGA)`=p.sureselect.v1.locus.sga.by.sample, `SureSelect.v2 (SGA)`=p.sureselect.v2.locus.sga.by.sample, `SureSelect.v3 (SGA)`=p.sureselect.v3.locus.sga.by.sample, heights = c(0.5,3,3,3,3)) + xlim(ranges(my.gr.ext.ensg)))
				dev.off()
			}
		}#end of my.assays
	}# end of my.targets
}

#################################
## Publication Quality Figures ##
#################################
# R/3.1.1 (this ggbio conflicts with ggplot which may be >2.0)
# use R/3.2.1 (this version of ggbio works, but color does does not work)
if(TRUE){
	my.target<-"CSMD1.promo" # CSMD1
	sureselect.v3.set1<-prep.sureselect("SureSelect.v3","Set1")
	sureselect.v3.set2<-prep.sureselect("SureSelect.v3","Set2")

	my.contrast<-"Gender"
	mycol<-my.col[[my.contrast]]
	#my.gr.ensg <- rtracklayer::import.bed("~/Pipelines/data/CSMD1/CSMD1.boy.girl.dmr.bed")
	my.gr.ensg <- rtracklayer::import.bed("~/Pipelines/data/CSMD1/CSMD1.placenta.putative.promoter.4K.GRCh37.bed")
	my.gr.ensg <- reduce(my.gr.ensg)
	names(my.gr.ensg)<-my.target
	my.gr.ext.ensg<-promoters(my.gr.ensg, upstream=my.upflank, downstream=end(my.gr.ensg)-start(my.gr.ensg)+1+my.downflank) # includes 2000 +TSS and 100 +TES 

	# overlap between cpgi and my.target 
	cpg.obj=read.feature.flank(location=my.cgi.file, feature.flank.name=c("CpGi","shores")) # a ‘GenomicRangesList’ contatining one GRanges object for flanks and one for GRanges object for the main feature.
	my.overlap=as.matrix( findOverlaps(cpg.obj[["CpGi"]], my.gr.ext.ensg) )
	my.cpgi <-cpg.obj$CpGi[my.overlap[,"queryHits"]]

	##############
	# chromosome #
	##############
	p.chr <- Ideogram(genome = "hg19", xlabel=TRUE) + xlim(ranges(my.gr.ext.ensg)) + theme_Publication() 

	####################
	# transcript track # 
	####################
	p.tr<-autoplot(hg19ensGene, which=my.gr.ext.ensg, fill = "orange", color = "grey") + theme_Publication() # hg19ensGene defined in 'Annotation.R'

	############################
	# cpgi track for my.target #
	############################
	p.cpgi <-autoplot(my.cpgi, fill = "darkgreen", color = "grey") + theme_Publication() 

	##########################
	# WGBS methylation GGBIO #
	##########################
	select<-mcols(gl.wgbs.locus[[my.target]])$Condition=="AGA"; my.gr.wgbs<-as.data.frame(gl.wgbs.locus[[my.target]][select,])
	p.wgbs.locus.aga<-ggplot(my.gr.wgbs, aes_string(x="start",y="score",col=my.contrast)) + geom_point(size=0.8,alpha=0.5) + ylim(0, 100) + 
						geom_smooth(se=FALSE) + 
						scale_colour_manual(name="Sex", values=mycol) + 
						theme_Publication() + 
						theme(legend.position=c(0.91,0.7)) 

	#################
	# sureselect.v3 #
	#################
	# sureselect.v3.set1 
	select.v3.aga<-start(sureselect.v3.set1[["gl"]][[my.target]]) %in% sureselect.v3.set1[["gender"]][[my.target]][["AGA"]] & mcols(sureselect.v3.set1[["gl"]][[my.target]])$Condition=="AGA"
	my.gr.sureselect<-as.data.frame(sureselect.v3.set1[["gl"]][[my.target]][select.v3.aga,])
	p.sureselect.v3.set1.locus.aga<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col=my.contrast)) + 
									geom_point(size=0.8,alpha=0.5) + ylim(0, 100) + 
									geom_smooth(se=FALSE,) + 
									scale_colour_manual(values=mycol) + 
									theme_Publication() +
									theme(legend.position="none") 
	#p.sureselect.v3.set1.locus.aga.by.sample<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col="SampleName")) + geom_point(size=0.8,alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette) + theme_Publication()

	# sureselect.v3.set2
	select.v3.aga<-start(sureselect.v3.set2[["gl"]][[my.target]]) %in% sureselect.v3.set2[["gender"]][[my.target]][["AGA"]] & mcols(sureselect.v3.set2[["gl"]][[my.target]])$Condition=="AGA"
	my.gr.sureselect<-as.data.frame(sureselect.v3.set2[["gl"]][[my.target]][select.v3.aga,])
	p.sureselect.v3.set2.locus.aga<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col=my.contrast)) + 
									geom_point(size=0.8,alpha=0.5) + ylim(0, 100) + 
									geom_smooth(se=FALSE,) + 
									scale_colour_manual(values=mycol) + 
									theme_Publication() +
									theme(legend.position="none") 
	#p.sureselect.v3.set2.locus.aga.by.sample<-ggplot(my.gr.sureselect, aes_string(x="start",y="score", col="SampleName")) + geom_point(size=0.8,alpha=0.5) + ylim(0, 100) + geom_smooth(se=FALSE) + scale_colour_manual(values=cbPalette) + theme_Publication()

	##############################################
	# PLOT WGBS & SureSelect methylation profile #
	##############################################
	p.placenta<-tracks(Transcripts=p.tr, `WGBS`=p.wgbs.locus.aga, `Technical`=p.sureselect.v3.set1.locus.aga, `Biological`=p.sureselect.v3.set2.locus.aga, heights = c(3,3.5,3.5,3.5), xlim=my.gr.ext.ensg)

	out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect.AGA",my.target,"publication.ggbio",sep="."))
	tiff(filename=paste0(out.file.name,".tiff"),width=15, height=8.5 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(p.placenta)
	dev.off()

	#########################################################
	# PLOT WGBS, SureSelect and RoadMap methylation profile #
	#########################################################
	source("~/Pipelines/bin/R/RoadMap/local.R")
	subject.keys=c("seqnames","start","end")
	query.key=c("V1","V2","End")
	subject<-as.data.frame(my.gr.ensg)
	cat("\tRemoving leading 'chr'...\n")
	#subject$seqnames<-simplify2array(strsplit(as.character(subject$seqnames), "chr"))[2,] # chrX = > X
	subject$seqnames<-substr(subject$seqnames, 4, 5) #chrX=>X
	cat("\tconvert to data.table...\n")
	dt.subject<-as.data.table(subject) # isa data.table
	setkeyv(dt.subject,subject.keys) # regions of interests 
	rm(subject)
	cat("Subject loaded\n")

	## Load RoadMap Tissue
	my.cpg.type="CG" # CG, CHH, CHG
	dt.tissue.locus<-list()
	for(my.tissue in c(avail.tissues[!avail.tissues %in% "PA"])){ # except 'PT': placenta
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
			setnames(dt.meth,c("chr","start","methylation_level","Gender","Tissue"))
			dt.tissue.locus[[my.tissue]]<-dt.meth
		}
	}# end of my.tissue 
	# Gender-specific
	p.roadmap<-ggplot(rbindlist(dt.tissue.locus), aes(start, methylation_level, col=Gender)) + geom_point(size=0.8,alpha=0.5) + ylim(0, 100) + 
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				ggtitle(my.target) + 
				theme_Publication() + 
				theme(legend.position="none") + 
				facet_grid(Tissue ~ .)

	# wgbs, tech, bio, roadmap
	out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect.AGA",my.target,"with.roadmap.publication.ggbio",sep="."))
	tiff(filename=paste0(out.file.name,".tiff"),width=10, height=15 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(tracks(Transcripts=p.tr, `WGBS`=p.wgbs.locus.aga, `technical`=p.sureselect.v3.set1.locus.aga, `biological`=p.sureselect.v3.set2.locus.aga, `RoadMap`=p.roadmap, heights = c(3,3.5,3.5,3.5,28), xlim=my.gr.ext.ensg))
	dev.off()

	## multiple ##
	#p1
	out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect.AGA",my.target,"with.roadmap.publication.ggbio.p1",sep="."))
	tiff(filename=paste0(out.file.name,".tiff"),width=10, height=15 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(tracks(Transcripts=p.tr, `WGBS`=p.wgbs.locus.aga, `technical`=p.sureselect.v3.set1.locus.aga, `biological`=p.sureselect.v3.set2.locus.aga, heights = c(3,3.5,3.5,3.5), xlim=my.gr.ext.ensg))
	dev.off()
	#p2
	out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,paste("wgbs.sureselect.AGA",my.target,"with.roadmap.publication.ggbio.p2",sep="."))
	tiff(filename=paste0(out.file.name,".tiff"),width=9, height=15 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(tracks(Transcripts=p.tr, `RoadMap`=p.roadmap, heights = c(3,28), xlim=my.gr.ext.ensg))
	dev.off()
}

############################################
## Concordance between WGBS and SureSelect##
############################################
# see README how to make this file
if(FALSE){
	if(my.slx=="SLX-10409.v1"){
		my.cmd<-paste0("ls ~/results/SLX-10409.Homo_sapiens.v1/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.top.concordance.with.wgbs.bed")
	}else if(my.slx=="SLX-10409.v2"){
		my.cmd<-paste0("ls ~/results/SLX-10409.Homo_sapiens.v2/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.top.concordance.with.wgbs.bed")
	}else if(my.slx=="SLX-10409.v3"){
		my.cmd<-paste0("ls ~/results/SLX-10409.Homo_sapiens.v3/BisSNP/HAL*/SLX-10409.HAL*.cpg.filtered.top.concordance.with.wgbs.bed")
	}else{
		stop()
	}

	dummy.bed<-system(my.cmd, intern=T)
	barcode.lst<-lapply(dummy.bed, function(i) unlist(strsplit(i, "/"))[7])
	foo<-lapply(dummy.bed, read.delim, header=F)
	names(foo)<- do.call(c, barcode.lst)
	wgbs.sureselect<-lapply(foo, function(i) i[i$V5>=min.doc & i$V16>=min.doc ,c("V1","V2","V3","V4","V5","V15","V16","V6")]) # V4:SureSelect %met, V15:WGBS %met

	# Colour by Chromosome 
	out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,"SureSelect.WGBS.concordance.by.chr")
	pdf(file=paste(out.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'),  width=11.7, height=8.3, title=my.project)
	lapply(names(wgbs.sureselect), 
		function(i) ggplot(wgbs.sureselect[[i]], aes(x=V4, y=V15, color=V1)) + 
					geom_point(size=3,alpha=0.7) + ylim(0, 100) + xlim(0,100) + 
					geom_smooth(method=lm,se=FALSE) + 
					labs(x=paste0("SureSelect (NO. of CpG = ",nrow(wgbs.sureselect[[i]]),")"),y="WGBS") + 
					ggtitle(paste0(i," (R^2=",round(summary(lm(wgbs.sureselect[[i]][,c("V4","V15")]))$r.squared,2),")")) +
					scale_color_discrete(name="Chromosome")
	)
	dev.off()
}

# Colour by Strand
if(FALSE){
out.file.name<-file.path("~/results/RnBeads",my.project,my.slx,"SureSelect.WGBS.concordance.by.strand")
pdf(file=paste(out.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title=my.project)
lapply(names(wgbs.sureselect), 
	function(i) ggplot(wgbs.sureselect[[i]], aes(x=V4, y=V15, color=V6)) + 
				geom_point(size=3,alpha=0.7) + ylim(0, 100) + xlim(0,100) +
				geom_smooth(method=lm,se=FALSE) + 
				labs(x=paste0("SureSelect (NO. of CpG = ",nrow(wgbs.sureselect[[i]]),")"),y="WGBS") + 
				ggtitle(paste0(i," (R^2=",round(summary(lm(wgbs.sureselect[[i]][,c("V4","V15")]))$r.squared,2),")")) +
				scale_color_discrete(name="Strand")
) 
dev.off()
}

cat("All is done", "\n")
