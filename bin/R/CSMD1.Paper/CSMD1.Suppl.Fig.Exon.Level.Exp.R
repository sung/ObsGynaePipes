	# part of ~/Pipelines/bin/R/CSMD1.Paper/CSMD1.paper.figures.R
	library(DESeq2)
	library(data.table)
	source("~/lib/ggplot_dual_axis.R") # only R>=3.2


	###########################
	## Figure xx.            ##
	## Exon-Level Expression ##
	## CSMD1
	###########################
	#exon definition
	my.target<-rtracklayer::import("/home/ssg29/data/genome/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/Homo_sapiens.GRCh37.75.CSMD1.exon.union.sorted.gtf")
	my.target.list<-split(my.target, mcols(my.target)[["exon_id"]]) # isa 'GRangeList'

	my.exon.RData<-paste0("~/results/CSMD1.2016.paper/CSMD1.",TR_PREFIX,".boy.girl.exon.RData")
	if(file.exists(my.exon.RData)){
		load(my.exon.RData)
	}else{
		## at this point, R version should be the same when deseq.RData was compiled.
		## probably 3.1
		deseq.RData <- file.path("~/results/RNA-Seq",paste0("Boy.Girl.exon.",TR_PREFIX,"/DESeq2/AGA/deseq.AGA.RData"))
		load (deseq.RData) # 

		my.contrast="Sex"

		dt.res<-as.data.frame(res); dt.res$ensg<-rownames(dt.res); dt.res<-as.data.table(dt.res)
		rn=rownames(colData(dds)) #as.character(samples[["SampleName"]]) # samples$SampleName (or samples[,c("SampleName")])
		rn.control=rownames(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1==0,]) 
		rn.case=rownames(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1!=0,])
		ddsFpkm <-fpkm(dds) # isa 'matrix'
		ddsFpm <-fpm(dds) # isa 'matrix'
		save(dt.res, rn, rn.control, rn.case, ddsFpkm, ddsFpm, file=my.exon.RData) 
	}

	##############################
	# import Ensembl exon region #
	##############################
	#my.target; my.target.list # defined within DEG.R
	subject.exon.keys=c("seqnames","strand","start","end")
	subject.exon<-as.data.frame(my.target)
	cat("\tconvert to data.table...\n")
	dt.ens.exon<-as.data.table(subject.exon) # isa data.table
	setkeyv(dt.ens.exon,subject.exon.keys) # regions of interests 
	cat("Ensembl Exon loaded\n")

	######################################
	# import reduced Ensembl exon region #
	######################################
	subject.rdx.exon<-as.data.frame(GenomicRanges::reduce(my.target))
	dt.rdx.exon<-as.data.table(subject.rdx.exon) # isa data.table
	dt.rdx.exon[,feature_id:=paste0("exon.",.N-.I+1)]
	dt.rdx.exon[,feature_type:="exon"]
	cat("Ensembl reduced Exon loaded\n")
	# exon.45 corresponds to the first expressed exon (ENSE00002127078)

	########################
	# import intron region #
	########################
	subject.intron<-as.data.frame(gaps(reduce(my.target))[2:length(reduce(my.target))])
	dt.intron<-as.data.table(subject.intron) # isa data.table
	dt.intron[,feature_id:=paste0("intron.",.N-.I+1)]
	dt.intron[,feature_type:="intron"]
	cat("Ensembl reduced Intron loaded\n")
	#rm(subject.exon, subject.rdx.exon, subject.intron)

	#################################
	# merge reduced exon and intron #
	#################################
	dt.rdx.feature=rbind(dt.rdx.exon, dt.intron)[order(start)]
	setkeyv(dt.rdx.feature,subject.exon.keys) # regions of interests 

	#####################################
	# methylation level at Ensembl exon #
	#####################################
	dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type) # ida 'data.table' for PT

	query.key=c("V1","V3","V2","End") # V1: chr, V3: strand, V2: start
	system.time(dt.overlap<-foverlaps(dt.query, dt.ens.exon, by.x=query.key, type="any", nomatch=0L))
	system.time(dt.dmr.exon<-dt.overlap[,list("num.sites"=.N,c.f=sum(V5.x), t.f=sum(V6.x)-sum(V5.x), c.m=sum(V5.y), t.m=sum(V6.y)-sum(V5.y),
										met=round(sum(V5.x,V5.y)/sum(V6.x,V6.y)*100,2), 
										met.f=round(sum(V5.x)/sum(V6.x)*100,2), 
										met.m=round(sum(V5.y)/sum(V6.y)*100,2), 
										p.value=fisher.test(matrix(c(sum(V5.x),sum(V6.x)-sum(V5.x),sum(V5.y),sum(V6.y)-sum(V5.y)),nrow=2))$p.value),
										by="exon_id"])
	########################################################
	# methylation level at redued features (exon & intron) #
	########################################################
	system.time(dt.overlap<-foverlaps(dt.query, dt.rdx.feature, by.x=query.key, type="any", nomatch=0L))
	system.time(dt.dmr.feature<-dt.overlap[,list("num.sites"=.N,c.f=sum(V5.x), t.f=sum(V6.x)-sum(V5.x), c.m=sum(V5.y), t.m=sum(V6.y)-sum(V5.y),
										met=round(sum(V5.x,V5.y)/sum(V6.x,V6.y)*100,2), 
										met.f=round(sum(V5.x)/sum(V6.x)*100,2), 
										met.m=round(sum(V5.y)/sum(V6.y)*100,2), 
										p.value=fisher.test(matrix(c(sum(V5.x),sum(V6.x)-sum(V5.x),sum(V5.y),sum(V6.y)-sum(V5.y)),nrow=2))$p.value),
										by="feature_id"])
	dt.dmr.feature<-merge(dt.dmr.feature, dt.rdx.feature, by="feature_id")[order(start)]
	dt.dmr.feature<-merge(dt.rdx.feature, dt.dmr.feature, by="feature_id", all.x=TRUE)[order(start)]

	# methylation level across reduced exons and intron
	p.feature.meth<-ggplot(dt.dmr.feature, aes(feature_id, met, group=1)) + 
				geom_point(aes(col=feature_type), size=rel(1)) +  
				geom_smooth(aes(group=feature_type,col=feature_type),se=FALSE) +
				scale_x_discrete(limits=dt.dmr.feature[,feature_id]) + 
				geom_line(linetype="dotted",na.rm=T) + 
				labs(x=" ", y="% Methylation") +
				theme_bw() +
				theme_Publication() +
				theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.title = element_text(face="italic",size=rel(.9)))

	# methylation level boxplot at exon and intron
	p.feature.meth.box<-ggplot(dt.dmr.feature[num.sites>5]) + 
						geom_boxplot(aes(feature_type, met)) + 
						theme_Publication()

	###########################
	# process expression data #
	###########################
	my.pval=0.05
	exon.pval<-sapply(dt.res[,ensg], function(i) wilcox.test(ddsFpm[i,rn.case], ddsFpm[i,rn.control])$p.value)
	dt.fpkm.exon<-data.table(`exon_id`=rownames(ddsFpkm), 
							 `FPKM.F`=rowMeans(ddsFpkm[,rn.case]), `FPKM.M`=rowMeans(ddsFpkm[,rn.control]), 
							 `FPM.F`=rowMeans(ddsFpm[,rn.case]), `FPM.M`=rowMeans(ddsFpm[,rn.control]),
							 `CNT.F`=rowMeans(counts(dds,normalized=T)[,rn.case]), `CNT.M`=rowMeans(counts(dds,normalized=T)[,rn.control]),
							 `pvalue`=exon.pval)
	# the index (row number) of my first exon
	dt.fpkm.exon[,.I[exon_id==my.first.exon_id]]


	##########################################
	# merge expression & methylation by exon #
	##########################################
	dt.exp.meth.exon<-merge(dt.fpkm.exon, dt.dmr.exon, by="exon_id", all.x=TRUE)[exon_id %in% dt.fpkm.exon[,exon_id]]
	dt.cnt.melt<-reshape2::melt(dt.exp.meth.exon, id=c("exon_id","pvalue"), measure=c("CNT.F","CNT.M"), variable.name="Gender", value.name="Count"); dt.cnt.melt[,Gender:=ifelse(Gender=="CNT.F","Female","Male")] # melt: col to row
	dt.cnt.melt[,significant:=ifelse(pvalue<my.pval,"Yes","No")]
	dt.meth.melt<-reshape2::melt(dt.exp.meth.exon, id=c("exon_id","num.sites","p.value"), measure=c("met.f","met.m"), variable.name="Gender", value.name="Meth"); dt.meth.melt[,Gender:=ifelse(Gender=="met.f","Female","Male")] # melt: col to row
	#merge(dt.cnt.melt, dt.meth.melt, by=c("exon_id","Gender"))

	dt.cnt.meth.diff<-rbind(dt.exp.meth.exon[,.(exon_id, `measure`="Expression", `diff`=CNT.F-CNT.M, `size`=(CNT.F+CNT.M)/2, `pval`=pvalue)],
							dt.exp.meth.exon[,.(exon_id, `measure`="Methylation",`diff`=met.f-met.m, `size`=num.sites, `pval`=p.value)])

	#####################
	## P-value plot    ##
	## with Read Count ##
	## with % Meth     ##
	#####################
	p.cnt.pval<-ggplot(dt.exp.meth.exon, aes(exon_id, pvalue, group=1)) + 
		geom_point(na.rm=T,size=rel(.5)) + 
		geom_line(aes(linetype="Read-Count Difference"),na.rm=T) + 
		geom_hline(aes(yintercept=my.pval)) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_y_continuous(breaks=c(0, my.pval, 0.25, .5, .75, 1)) +
		scale_linetype_manual(name="P-value", values=c(3)) + # dotted
		labs(x="Exon ID", y="P-value") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.position="top" ,legend.title = element_text(face="italic",size=rel(.9)))

	p.meth.pval<-ggplot(dt.exp.meth.exon, aes(exon_id, p.value, group=1)) + 
		geom_point(na.rm=T,size=rel(.5)) + 
		geom_line(linetype="dotted",na.rm=T) + 
		geom_hline(aes(yintercept=my.pval)) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_y_continuous(breaks=c(0, my.pval, 0.25, .5, .75, 1)) +
		labs(x="Exon ID", y="P-value") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.title = element_text(face="italic",size=rel(.9)))

	# p-value both from cnt & meth
	p.cnt.meth.pval<-ggplot(dt.cnt.meth.diff, aes(exon_id, pval, group=measure)) + 
		geom_point(aes(shape=measure), na.rm=T) + 
		geom_line(aes(linetype=measure),na.rm=T) +
		geom_hline(aes(yintercept=my.pval)) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_y_continuous(breaks=c(0, my.pval, 0.25, .5, .75, 1)) +
		scale_linetype_manual(name="Source of P-value",values=c(3,1)) +
		labs(x="Exon ID", y="P-value") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.title = element_text(face="italic",size=rel(.9)))

	# p-value both from cnt & meth with diff and size
	p.cnt.meth.pval2<-ggplot(dt.cnt.meth.diff, aes(exon_id, pval, group=measure)) + 
		geom_point(aes(size=size,col=diff,shape=measure), na.rm=T) + 
		geom_line(aes(linetype=measure),col='gray70',na.rm=T) +
		geom_hline(aes(yintercept=my.pval)) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_y_continuous(breaks=c(0, my.pval, 0.25, .5, .75, 1)) +
		scale_linetype_manual(name="Source of P-value",values=c(1,2)) +
		scale_colour_gradient2(name="Difference\n(Female-Male)",mid="gray48", high="red", low="blue", midpoint=0) +
		labs(x="Exon ID", y="P-value") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.title = element_text(face="italic",size=rel(.9)))

	# methylation level without no line and with No. CpG by shape
	p.meth<-ggplot(dt.meth.melt, aes(exon_id, Meth, group=Gender, colour=Gender)) +
		geom_point(aes(size=num.sites),shape=8,na.rm=T, alpha=0.9) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_size_continuous(name="No. CpG") +
		scale_colour_manual(name="% Methylation\n(by sex)",values=my.col[["Gender"]]) +
		labs(x="Exon ID", y="% Methylation") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.position="top", legend.title = element_text(face="italic",size=rel(.7)))
	# dotted-line with No. CpG by shape
	p.meth2<-ggplot(dt.meth.melt, aes(exon_id, Meth, group=Gender, colour=Gender)) +
		geom_point(aes(size=num.sites),shape=8,na.rm=T, alpha=0.9) + 
		geom_line(linetype="dotted",na.rm=T) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_size_continuous(name="No. CpG") +
		scale_colour_manual(name="% Methylation\n(by sex)",values=my.col[["Gender"]]) +
		labs(x="Exon ID", y="% Methylation") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.justification="top", legend.position=c(0.8,1), legend.title = element_text(face="italic",size=rel(.7)))

	# read-count without line with read-count by shape
	p.cnt<-ggplot(dt.cnt.melt, aes(exon_id, Count, colour=Gender)) +
		geom_point(aes(size=Count),na.rm=T, alpha=0.9) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_size_continuous(name="Read Count") +
		scale_colour_manual(name="Read Count\n(by sex)",values=my.col[["Gender"]]) +
		labs(x="Exon ID", y="Read Count") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.position="top", legend.title = element_text(face="italic",size=rel(.7)))

	# read-count solid-line with different shape by p-value threshold
	p.cnt2<-ggplot(dt.cnt.melt, aes(exon_id, Count, colour=Gender)) +
		geom_point(aes(shape=factor(significant)),size=3,na.rm=T, alpha=0.9) + 
		geom_line(aes(group=Gender),na.rm=T) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_colour_manual(name="Read Count\n(by sex)",values=my.col[["Gender"]]) +
		scale_shape(name=paste("P-value<",my.pval)) +
		labs(x="Exon ID", y="Read Count") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.justification="top", legend.position=c(0.8,1), legend.title = element_text(face="italic",size=rel(.7)))

	# log read-count without line with read-count by shape
	p.cnt3<-ggplot(dt.cnt.melt, aes(exon_id, log(Count+0.1), colour=Gender)) +
		geom_point(aes(size=Count),na.rm=T, alpha=0.9) + 
		scale_x_discrete(limits=dt.fpkm.exon[,exon_id]) + 
		scale_size_continuous(name="Read Count") +
		scale_colour_manual(name="Sex",values=my.col[["Gender"]]) +
		labs(x="Exon ID", y="Log(Read-Count + 0.1)") +
		theme_bw() +
		theme(axis.text.x = element_text(angle =-90, hjust = 1, size=rel(0.7)), axis.title=element_text(face="bold",size=rel(1.2)), legend.position="top", legend.title = element_text(face="italic",size=rel(.7)))

	# 1. count & p-val
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.exp.exon.pval.by.sex.dual.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=8,units="in",res=300, compression = 'lzw')
	ggplot_dual_axis(p.cnt,p.cnt.pval,"y")
	dev.off()

	# 2. log-count & p-val
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.log.cnt.pval.per.exon.by.sex.dual3.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=8,units="in",res=300, compression = 'lzw')
	ggplot_dual_axis(p.cnt3,p.cnt.meth.pval,"y") # p-val (exp) & p-val (meth)
	#ggplot_dual_axis(p.cnt3,p.cnt.pval,"y")
	#ggplot_dual_axis(p.cnt.pval,p.cnt3,"y")
	dev.off()
		
	# 3. meth & p-val
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.exon.pval.by.sex.dual.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=8,units="in",res=300, compression = 'lzw')
	ggplot_dual_axis(p.meth,p.cnt.pval,"y")
	dev.off()

	# 4. meth & p-val (meth)
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.exon.pval2.by.sex.dual.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=8,units="in",res=300, compression = 'lzw')
	ggplot_dual_axis(p.meth,p.meth.pval,"y")
	dev.off()

	# 5-1. cnt & meth
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.exp.exon.by.sex.dual.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=7,units="in",res=300, compression = 'lzw')
	ggplot_dual_axis(p.cnt2,p.meth2,"y")
	dev.off()

	# 5-2. meth & cnt
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.exp2.exon.by.sex.dual.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=7,units="in",res=300, compression = 'lzw')
	ggplot_dual_axis(p.meth2,p.cnt2,"y")
	dev.off()

	# 6. diff meth & diff exp
	hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
	#p.diff<-ggplot(dt.exp.meth.exon, aes(CNT.F-CNT.M, met.f-met.m)) + 
	p.diff<-ggplot(dt.exp.meth.exon[(CNT.F+CNT.M)/2>1 & num.sites>1], aes(CNT.F-CNT.M, met.f-met.m)) + 
		geom_point(aes(size=(CNT.F+CNT.M)/2, col=met)) +
		scale_size_continuous(name="Average\nRead Count") +
		scale_colour_gradientn(name="Average\n% Methylation",colors=hmcol) +
		labs(x="Read Count Difference (Female-Male)", y="% Methylation Difference (Female-Male)") +
		geom_hline(aes(yintercept=0)) +
		geom_vline(aes(xintercept=0)) +
		theme_Publication()
	my.filename <- file.path(cluster.dir, paste0(sampleType, ".CSMD1.diff.meth.diff.exp.by.exon.by.sex.",TR_PREFIX,".tiff"))
	tiff(filename=my.filename,width=10,height=8.5,units="in",res=300, compression = 'lzw')
	print(p.diff)
	dev.off()

