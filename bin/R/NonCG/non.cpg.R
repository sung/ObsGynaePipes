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

meta.sample<-read.csv("/home/ssg29/Pipelines/data/WGBS.NON.CG/wgbs.non.cg.A001.csv")
#################################################
### WGBS NON-CG PER_CHR Methylation BED Files ###
#################################################
gr.fld<-lapply(as.character(meta.sample$BedFile), rtracklayer::import.bed) # isa 'GRangesList'
names(gr.fld)<-meta.sample$SampleName
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

