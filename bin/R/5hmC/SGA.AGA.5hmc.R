#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# This R script generate two 5hmC files from two individuals each of which have one BS-Seq and one osBS-Seq
# AGA2F: AGA2F.bs + AGA2F.ox
# SGA2F: SGA2F.bs + SGA2F.ox
# min depth: 10
# output:
# /home/ssg29/results/SLX-8074.SLX-8080.v1/BisSNP/*.bed
# /home/ssg29/results/SLX-8074.SLX-8080.v1/BisSNP/*.methylkit

library(methylKit)
library(GenomicRanges)
top.dir<-"/home/ssg29/results/SLX-8074.SLX-8080.v1"

pdf.file.name<- paste0(top.dir,'/MethylKit/5hmC.plot')
pdf(file=paste(pdf.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

my.depth<-10 # at least 10 depth

#BS-Seq: 5mC+5hmC
#oxBS-Seq: 5mC

if(FALSE){
	my.unite.data<-file.path(top.dir,"MethylKit/AGA2F.SGA2F.unite.RData")
	if(file.exists(my.save.data)){
		load(my.unite.data)
		cat("AGA2F, SGA2F, AGA2F.unite.strand, SGA2F.unite.strand, AGA2F.unite.destrand, and SGA2F.unite.destrand loaded \n")
	}else{
		AGA2F.file.list=list(
			file.path(top.dir,"BisSNP/A010/SLX-8074.SLX-8080.A010.cpg.filtered.homo.CG.methylkit.txt"), #BS AGA2F
			file.path(top.dir,"BisSNP/A001/SLX-8074.SLX-8080.A001.cpg.filtered.homo.CG.methylkit.txt")  #oxBS AGA2F
		)
		SGA2F.file.list=list(
			file.path(top.dir,"BisSNP/A012/SLX-8074.SLX-8080.A012.cpg.filtered.homo.CG.methylkit.txt"), #BS SGA2F
			file.path(top.dir,"BisSNP/A003/SLX-8074.SLX-8080.A003.cpg.filtered.homo.CG.methylkit.txt")  #oxBS SGA2F
		)

		#methylRawList
		AGA2F=read(AGA2F.file.list, sample.id=list("AGA2F.bs", "AGA2F.ox"), assembly="hg19",treatment=c(0,0),context="CpG") # is a 'methylRawList'
		SGA2F=read(SGA2F.file.list, sample.id=list("SGA2F.bs", "SGA2F.ox"), assembly="hg19",treatment=c(1,1),context="CpG") # is a 'methylRawList'
		#filtered.AGA2F=filterByCoverage(AGA2F,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9) # is a 'methylRawList'  # depth filter will be applied later

		# unite BS & oxBS
		AGA2F.unite.strand<-unite(AGA2F,destrand=FALSE)
		SGA2F.unite.strand<-unite(SGA2F,destrand=FALSE)

		AGA2F.unite.destrand<-unite(AGA2F,destrand=TRUE)
		SGA2F.unite.destrand<-unite(SGA2F,destrand=TRUE)

		save(AGA2F, SGA2F, AGA2F.unite.strand, SGA2F.unite.strand, AGA2F.unite.destrand, SGA2F.unite.destrand, file=my.unite.data) #
	}
}

# function based on methylKit
get.5hmC<-function(united){
	merged<-getData(united)
	# %5hmC = (% of 5mC and 5mhC from BS) - (% of 5mC from oxBS)
	# = (merged$numCs1/merged$coverage1) - (merged$numCs2/merged$coverage2)
	# num(5hmC) = (merged$numCs1) - (merged$coverage1)*((merged$numCs2/merged$coverage2))

	# possible filter?
	# 1. coverage1 >=10 & coverage2 >=10
	# this should be done at this (or earlier via filterByCoverage) stage
	# otherwise, the depth info from 5mC (oxBS) will be gone
	select<-merged$coverage1>=my.depth & merged$coverage2>=my.depth 
	merged<-merged[select,]

	diff=(merged$numCs1)-round(merged$coverage1*(merged$numCs2/merged$coverage2))
	merged$numCs1=diff

	# possible filter?
	# 2. remove -ve value or set to 0
	# this is done later when saving the data in .bed or methylKit format 
	#diff[diff<0]=0

	merged$numTs1=merged$coverage1-merged$numCs1
	colnames(merged)[5:7]=c("coverage","numCs","numTs")
	new("methylRaw",merged[,1:7],sample.id=united@sample.ids[1], assembly=united@assembly, context =united@context, resolution=united@resolution)

}

#> head(getData(SGA2F.unite))
#   chr start   end strand coverage1 numCs1 numTs1 coverage2 numCs2 numTs2
#1 chr1 10470 10470      -         4      0      4         1      0      1
get.5hmC.new<-function(united){
	merged<-getData(united) # isa data.frame

	cat("calculating 5hmC...\n")
	merged[,"num5hmC"]<-apply(merged[,c("numCs1","coverage1","numCs2","coverage2")], 1, function(x){x[1]-round(x[2]*(x[3]/x[4]))}) # numCs1-round(coverage1*(numCs2/coverage2))

	cat("calculating pvalue...\n")
	merged[,"pvalue"]<-apply(merged[,c("numCs1","numTs1","numCs2","numTs2")],1, function(x){round(fisher.test(matrix(x, nrow=2), alternative="greater")$p.value,5)})
	#merged[,"pvalue"]<-with(merged, fisher.test(matrix(c(numCs1,numTs1,numCs2,numTs2), nrow=2), alternative="greater")$p.value) # segfault

	cat("calculating adjusted pvalue...\n")
	merged[,"pvalue.adj"]<-p.adjust(merged$pvalue,method="BH")

	#return(merged)
	return(
		GRanges(
				seqnames<-Rle(merged$chr), # isa factor Rle object
				ranges<-IRanges(start=merged$start, end=merged$end, names=paste0(merged$chr,'.',merged$start)), # isa IRanges
				strand<-Rle(merged$strand), # isa factor Rle object
				mcols<-DataFrame(with(merged,data.frame(bs.coverage=coverage1,bs.numCs=numCs1,bs.numTs=numTs1,ox.coverage=coverage2,ox.numCs=numCs2,ox.numTs=numTs2,num5hmC=num5hmC,pvalue=pvalue,pvalue.adj=pvalue.adj)))
		)
	)
}

# get 5hmC
#my.5hmc.data<-file.path(top.dir,"MethylKit/AGA2F.SGA2F.5hmC.RData") # there could be CpGs where BS(>=10) whereas oxBS(<=10)
my.5hmc.data<-file.path(top.dir,"MethylKit/AGA2F.SGA2F.5hmC.new.RData")
if(file.exists(my.5hmc.data)){
	cat("Loading my.5hmc.data...\n")
	load(my.5hmc.data)
	cat("AGA2F.5hmC.strand, SGA2F.5hmC.strand, AGA2F.5hmC.destrand and SGA2F.5hmC.destrand loaded\n")
}else{
	AGA2F.5hmC.strand<-get.5hmC(AGA2F.unite.strand)
	SGA2F.5hmC.strand<-get.5hmC(SGA2F.unite.strand)
	AGA2F.5hmC.destrand<-get.5hmC(AGA2F.unite.destrand)
	SGA2F.5hmC.destrand<-get.5hmC(SGA2F.unite.destrand)
	save(AGA2F.5hmC.strand, SGA2F.5hmC.strand, AGA2F.5hmC.destrand, SGA2F.5hmC.destrand, file=my.5hmc.data) # 8oxBS.meth.RData
}

plot.5hmc<-function(strand,destrand,my.type){
	# diagnostic
	dummy<-list()
	dummy[["strand"]]<-getData(strand)
	dummy[["destrand"]]<-getData(destrand)
	my.col<-c(`strand`="red",`destrand`="blue")

	# plot the proportion of negative 5hmC at depth 1,5,10,15,20..100
	# by stranded and unstranded CpG read depth
	#neg.5hmc<-data.frame(type=character(),depth=numeric(),ratio=numeric(),stringsAsFactors=FALSE)
	index<-0
	for(depth in c(1,seq(5,100,5))){
		for(i in names(dummy)){
			index=index+1

			filter<-dummy[[i]]$coverage>=depth # depth filter
			select<-dummy[[i]][filter,]$numCs<0 # negagive 5hmC at depth X

			depth.ratio<-round( table(filter)["TRUE"]/length(filter), 3) # ratio of CpG more than depth x out of total CpG
			neg.ratio<-round( table(select)["TRUE"]/length(select), 3) # ratio of negagive 5hmC out of total 5hmC at depth x

			cat(paste("type=",i,"depth=",depth,"%5hmC<0:",neg.ratio, "\n"))
			#neg.5hmc[index,]<-c(i,depth,neg.ratio)

			#plot the distribution via histogram
			#dummy2<-dummy[[i]][filter,]$numCs/dummy[[i]][filter,]$coverage
			#hist(dummy2, probability=TRUE, main=paste0(i,"ed CpG at depth=",depth), xlab="% 5hmC methylation", right=FALSE) # [start,end)

			if(depth==1 && i=="strand"){
				plot(depth, neg.ratio, xlim=c(0,100), ylim=c(0,1), main=my.type, xlab="CpG Depth", ylab="Ratio", col=my.col[i])
				par(new=TRUE); plot(depth, depth.ratio, xlim=c(0,100), ylim=c(0,1), xlab="", ylab="", axes=F, pch=4, col=my.col[i])
			}else{
				par(new=TRUE)
				# http://www.sixhat.net/plotting-multiple-data-series-in-r.html
				#NB: 1. xlim & ylim should be same as above
				#    2. xlab, ylab and axes
				plot(depth, neg.ratio, xlim=c(0,100), ylim=c(0,1), xlab="", ylab="", axes=F, col=my.col[i])
				par(new=TRUE); plot(depth, depth.ratio, xlim=c(0,100), ylim=c(0,1), xlab="", ylab="", axes=F, pch=4, col=my.col[i])
			}
		}
	}
	#legend("topright", legend=names(my.col), col=my.col, pch=c(1,1,))
	legend("topright", legend=c(names(my.col), NA, "CpG depth", "Negative 5hmC"), col=c(my.col,NA,"black","black"), text.col=c(my.col,NA,"black","black"), pch=c(NA,NA,NA,4,1) )
	par(new=FALSE)
}

#plot.5hmc(AGA2F.5hmC.strand, AGA2F.5hmC.destrand, "AGA")
#plot.5hmc(SGA2F.5hmC.strand, SGA2F.5hmC.destrand, "SGA")


##################
## write bed file
##################
#track name=SLX-8074.SLX-8080.A012.cpg.filtered.CG.bed type=bedDetail description='CpG methylation level' visibility=3
#chr1    10468   10469   0       1       +       10468   10469   0,240,0 0       4
#chr1    10469   10470   0       4       -       10469   10470   0,240,0 0       4
if(FALSE){
	options(scipen=999) # diable scientific notation
	# AGA
	dummy<-getData(AGA2F.5hmC.strand)
	select<-dummy$numCs>=0 # remove negagive 5hmC at depth X
	dummy<-with(dummy[select,], data.frame(chr=chr,start=start-1,end=end,met=round(numCs/coverage*100,3),coverage=coverage,strand=strand,start2=start-1,end2=end,rbg="0,240,0",score1=0,score2=4))

	write.table("track name=AGA.5hmC.strand type=bedDetail description='5hmC methylation level' visibility=3", file=file.path(top.dir,"BisSNP/AGA2F.5hmC.strand.bed"), quote=FALSE, row.names=F, col.names=F)
	write.table(dummy, file.path(top.dir,"BisSNP/AGA2F.5hmC.strand.bed"), append=TRUE, quote=FALSE, sep = "\t", row.names=F, col.names=F)

	# SGA
	dummy<-getData(SGA2F.5hmC.strand)
	select<-dummy$numCs>=0 # remove negagive 5hmC at depth X
	dummy<-with(dummy[select,], data.frame(chr=chr,start=start-1,end=end,met=round(numCs/coverage*100,3),coverage=coverage,strand=strand,start2=start-1,end2=end,rbg="0,240,0",score1=0,score2=4))
	write.table("track name=SGA.5hmC.strand type=bedDetail description='5hmC methylation level' visibility=3", file=file.path(top.dir,"BisSNP/SGA2F.5hmC.strand.bed"), quote=FALSE, row.names=F, col.names=F)
	write.table(dummy, file=file.path(top.dir,"BisSNP/SGA2F.5hmC.strand.bed"), append=TRUE, quote=FALSE, sep = "\t", row.names=F, col.names=F)


#########################
## write MethylKit file
#########################
#ChrBase chr base strand coverage freqC freqT
#chr1.10469 chr1 10469 F 6 0 100
#chr1.10470 chr1 10470 R 7 14.2857 85.7143
#chr1.10471 chr1 10471 F 12 8.33333 91.6667
	#AGA
	dummy<-getData(AGA2F.5hmC.strand)
	select<-dummy$numCs>=0 # remove negagive 5hmC at depth X
	dummy<-with(dummy[select,], data.frame(ChrBase=paste(chr,start,sep='.'), chr=chr, base=start, strand=ifelse(strand=='+',"F","R"), coverage=coverage, freqC=round(numCs/coverage*100,3), freqT=round(numTs/coverage*100,3)))
	write.table(dummy, file.path(top.dir,"BisSNP/AGA2F.5hmC.strand.methylkit.txt"), quote=FALSE, sep = " ", row.names=F, col.names=T)
	#SGA
	dummy<-getData(SGA2F.5hmC.strand)
	select<-dummy$numCs>=0 # remove negagive 5hmC at depth X
	dummy<-with(dummy[select,], data.frame(ChrBase=paste(chr,start,sep='.'), chr=chr, base=start, strand=ifelse(strand=='+',"F","R"), coverage=coverage, freqC=round(numCs/coverage*100,3), freqT=round(numTs/coverage*100,3)))
	write.table(dummy, file.path(top.dir,"BisSNP/SGA2F.5hmC.strand.methylkit.txt"), quote=FALSE, sep = " ", row.names=F, col.names=T)
}

#################
## get 5hmC.gr ##
#################
if(FALSE){
	my.gr.data<-file.path(top.dir,"MethylKit/AGA2F.SGA2F.gr.RData")
	if(file.exists(my.gr.data)){
		cat("Loading my.5hmc.gr.data...\n")
		load(my.gr.data)
		cat("AGA2F.gr.strand, SGA2F.gr.strand loaded\n")
	}else{
		AGA2F.gr.strand<-get.5hmC.new(AGA2F.unite.strand)
		SGA2F.gr.strand<-get.5hmC.new(SGA2F.unite.strand)
		save(AGA2F.gr.strand, SGA2F.gr.strand, file=my.gr.data) 

		#AGA2F.gr.destrand<-get.5hmC.new(AGA2F.unite.destrand)
		#SGA2F.gr.destrand<-get.5hmC.new(SGA2F.unite.destrand)
		#save(AGA2F.gr.strand, SGA2F.gr.strand, AGA2F.gr.destrand, SGA2F.gr.destrand, file=my.gr.data) 
	}
}

# How to adjust negagive 5hmC?

dev.off()
cat("All done\n")
