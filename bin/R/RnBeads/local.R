
library(reshape2)

# UCSC-style (e.g. chr8)
target.region<-list(
	`CSMD1.DMR`="~/Pipelines/data/CSMD1/CSMD1.boy.girl.dmr.bed", # Original CSMD1 DMR
	`CSMD1.DMR2`="~/Pipelines/data/CSMD1/CSMD1.boy.girl.dmr+2exons.bed", # above + 5' 2 exons 
	`CSMD1.legacy.promo`="~/Pipelines/data/CSMD1/CSMD1.legacy.promo.4K.GRCh37.bed",
	`CSMD1.putative.promo`="~/Pipelines/data/CSMD1/CSMD1.placenta.putative.promoter.4K.GRCh37.bed"
)	

prep.wgbs<-function(){
	gl.wgbs.locus<-list() # a list of GRanges
	df.wgbs.locus<-list() # a list of data.frame

	# WGBS Methylation from BED Files ###
	dummy<-c(
				`CSMD1.DMR`="~/Pipelines/data/CSMD1/CSMD1.225kb.grch37.WGBS.bed.txt",
				`CSMD1.DMR2`="~/Pipelines/data/CSMD1/CSMD1.247kb.grch37.WGBS.bed.txt",
				`CSMD1.legacy.promo`="~/Pipelines/data/CSMD1/CSMD1.legacy.promoter.4K.grch37.WGBS.bed.txt",
				`CSMD1.putative.promo`="~/Pipelines/data/CSMD1/CSMD1.putative.promoter.4K.grch37.WGBS.bed.txt"
			)
	dummy<-lapply(dummy, read.csv, stringsAsFactors=FALSE) # isa 'list'

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
	return(list(`gl`=gl.wgbs.locus, `df`=df.wgbs.locus))
}

prep.sureselect<-function(my.ver,my.assay){
	gl.sureselect.locus<-list() # a list of GRanges
	df.sureselect.locus<-list() # a list of data.frame
	boy.girl.cpg<-list()        # a list of list of numeric vector to store CpG position by SGA/AGA separately

	dummy <-c(
				#`NKX1-2`=paste0("~/Pipelines/data/MJ.SureSelect/",my.ver,"/SLX-10409.NKX1-2.sample.annotation.csv"),
				#`MAB21L1`=paste0("~/Pipelines/data/MJ.SureSelect/",my.ver,"/SLX-10409.MAB21L1.top2.sample.annotation.csv"),
				#`chr7.DMR`=paste0("~/Pipelines/data/MJ.SureSelect/",my.ver,"/SLX-10409.chr7.DMR.top2.sample.annotation.csv"),
				#`VIPR2`=paste0("~/Pipelines/data/MJ.SureSelect/",my.ver,"/SLX-10409.VIPR2.5kb.sample.annotation.csv"),
				#`C1orf216`=paste0("~/Pipelines/data/MJ.SureSelect/",my.ver,"/SLX-10409.C1orf216.sample.annotation.csv"),
				#`TICAM2`=paste0("~/Pipelines/data/MJ.SureSelect/",my.ver,"/SLX-10409.TICAM2.sample.annotation.csv")
				`CSMD1.DMR`=paste0("~/Pipelines/data/MJ.SureSelect/",my.ver,"/SLX-10409.CSMD1.225kb.sample.annotation.csv"),
				`CSMD1.DMR2`=paste0("~/Pipelines/data/MJ.SureSelect/",my.ver,"/SLX-10409.CSMD1.247kb.sample.annotation.csv"),
				`CSMD1.legacy.promo`=paste0("~/Pipelines/data/MJ.SureSelect/",my.ver,"/SLX-10409.CSMD1.legacy.promoter.sample.annotation.csv"),
				`CSMD1.putative.promo`=paste0("~/Pipelines/data/MJ.SureSelect/",my.ver,"/SLX-10409.CSMD1.placenta.putative.promoter.sample.annotation.csv")
			)
	dummy <-lapply(dummy, read.csv, stringsAsFactors=FALSE) # isa 'list'
	dummy <-lapply(dummy, function(i) i[i$SetNo==my.assay,]) # Set1 only (WGBS samples)

	# for each target
	for(i in names(dummy)){
		gl.dummy<-lapply(dummy[[i]]$BedFile, rtracklayer::import.bed) 
		names(gl.dummy)<-dummy[[i]]$SampleName
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

		# filter CpG called by all/some samples (dcast from 'reshape2')
		df.dummy <- dcast(as.data.frame(gr.dummy), start ~ SampleName, value.var="score"); rownames(df.dummy) <- df.dummy$start; df.dummy$start<-NULL; 
		if(i=="NKX1-2" || i=="MAB21L1"){ # Targets from SGA/AGA study
			gl.sureselect.locus[[i]] <- gr.dummy[start(gr.dummy) %in% rownames(df.dummy[rowSums(!is.na(df.dummy)) >= ncol(df.dummy)*(1-min.missing.sample.ratio),]),]
		}else{  # Targets from the boy/girl study
			boy.girl.cpg[[i]][["AGA"]]<-rownames(df.dummy[rowSums(!is.na(df.dummy[,grepl("AGA",colnames(df.dummy))])) >= ncol(df.dummy[,grepl("AGA",colnames(df.dummy))])*(1-min.missing.sample.ratio),grepl("AGA",colnames(df.dummy))]) # CpGs in AGA
			boy.girl.cpg[[i]][["SGA"]]<-rownames(df.dummy[rowSums(!is.na(df.dummy[,grepl("SGA",colnames(df.dummy))])) >= ncol(df.dummy[,grepl("SGA",colnames(df.dummy))])*(1-min.missing.sample.ratio),grepl("SGA",colnames(df.dummy))]) # CpGs in SGA
			gl.sureselect.locus[[i]] <- gr.dummy[start(gr.dummy) %in% union(boy.girl.cpg[[i]][["AGA"]],boy.girl.cpg[[i]][["SGA"]]),]
		}
		df.sureselect.locus[[i]] <- df.dummy
	}
	return(list(`gl`=gl.sureselect.locus, `df`=df.sureselect.locus, `gender`=boy.girl.cpg))
} # end of function 'prep.sureselect'

