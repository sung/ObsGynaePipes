#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 30/Mar/2016
# Last modified 27/Jul/2016

TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/RoadMap/local.R")
library(gplots)

################
# get XCI list #
################
# this is name-based annotation (carrel.xci$gene)
my.RData="~/results/RoadMap/X-inactivation/RData/carrel.xci.genes.grch37.RData"
if(file.exists(my.RData)){
	load(my.RData)
	cat("xci.anno and dt.xci.anno loaded\n")
}else{
	fields <- c("chromosome_name", "ensembl_gene_id", "hgnc_symbol", "gene_biotype")
	# method1: by gene name
	carrel.xci=read.csv("~/data/X-inactivation/Carrel.science.2005.csv")
	xci<-getBM(attributes = fields, filters = "hgnc_symbol", values = carrel.xci$gene, mart = myMart) # is a data.frame
	xci.anno<-merge(carrel.xci, xci, by.x="gene", by.y="hgnc_symbol") # isa 'data.frame'
	xci.anno<-xci.anno[xci.anno$chromosome_name=="X",]
	dt.xci.anno<-as.data.table(xci.anno)
	save(xci.anno, dt.xci.anno, file="~/results/RoadMap/X-inactivation/RData/carrel.xci.genes.grch37.RData")
}

###################################################
# get methylation of Carrel for all tissue types ##
###################################################
#see bin/R/RoadMap/plot.roadmap.gender.meth.diff.by.XCI.Carrel.science.2005.R
#defined only following regions: cpg.contexts=c("Genes","Promo_15.05", "Promo_10.03","CPGi","CPGi-shores")
#for "CPGi" and "CPGi-shores", the overlapped genes are from Genes_20.20 (see the code)
# ensembl_gene_id from RnBeads (Ensembl 73 or gencode 18)
load("~/results/RoadMap/X-inactivation/RData/dt.meth.chrX.per.tissue.CG.Genes_20.20.RData")
cat("dt.meth.chrX.per.tissue loaded\n")
dt.meth.all.chrX<-data.table::rbindlist(dt.meth.chrX.per.tissue) 
## Note that ensembl_gene_id can be null <NA> because of left-join (all.x=TRUE)
dt.meth.xci<-merge(dt.meth.all.chrX, xci.anno, by="ensembl_gene_id")

################################
## get the RoadMap expression ##
################################
my.exp.RData="~/results/RoadMap/X-inactivation/RData/dt.exp.roadmap.RData"
if(file.exists(my.exp.RData)){
	cat(paste0("loading ",my.exp.RData,"...\n"))
	load (my.exp.RData) 
	cat("dt.tissue.exp.list loaded\n")
}else{
	roadmap.exp<-read.delim("~/data/RoadMap/Schultz2015/Supplementary_Information/RNA_matrix_FPKM.tsv", stringsAsFactors=FALSE)
	roadmap.exp.gene=simplify2array(strsplit(roadmap.exp$gene_id,"[.]"))[1,]
	#dt.roadmap.exp<-data.table::fread("~/data/RoadMap/Schultz2015/Supplementary_Information/RNA_matrix_FPKM.tsv")
	#roadmap.exp.gene=simplify2array(strsplit(dt.roadmap.exp[,gene_id],"[.]"))[1,]

	dt.tissue.exp.list=list()
	for(my.tissue in avail.tissues[!avail.tissues %in% "PA"]){
		this.tissue.female<-paste0(my.tissue, "_2_FPKM") # female
		this.tissue.male<-paste0(my.tissue, "_3_FPKM") # male

		dt.tissue.exp.list[[my.tissue]]=data.table(
				`ensembl_gene_id`=roadmap.exp.gene,
				`tissue`=my.tissue,
				`female.exp`=roadmap.exp[,this.tissue.female],
				`male.exp`=roadmap.exp[,this.tissue.male]
			)
	}
	#######################
	## get PT expression ##  
	#######################
	# Boy.Girl.FG.JD.GRCh37
	myCaller='DESeq2'
	source ("~/Pipelines/config/DEG.R") # load config
	deseq.RData <- paste0(deseq.dir, "/deseq.",sampleType,'.RData')
	load (deseq.RData) 
	cat("dds, res, rld, and (ddsExonicGene?) loaded\n")

	my.res<-res
	rn=as.character(samples[["SampleName"]]) # samples$SampleName (or samples[,c("SampleName")])
	rn.case=as.character(samples[as.numeric(samples[[my.contrast]])-1==1,c("SampleName")]) # case only (SGA or Male)
	rn.control=as.character(samples[as.numeric(samples[[my.contrast]])-1==0,c("SampleName")]) # control only (AGA or Female), which is background
	cat("setting FPKM...\n")
	if(grepl("piRNA$",myProject)){ # If piRNA from small RNA-Seq
		ddsFpkm <-fpkm(dds) # isa 'matrix'
		ddsFpm <-fpm(dds) # isa 'matrix'
	}else{ 
		ddsFpkm <-fpkm(ddsExonicGene) # isa 'matrix'
		ddsFpm <-fpm(ddsExonicGene) # isa 'matrix'
	}
	
	my.res<-res
	rn=rownames(colData(dds)) #as.character(samples[["SampleName"]]) # samples$SampleName (or samples[,c("SampleName")])
	rn.control=rownames(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1==0,])
	            #as.character(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1==0,c("SampleName")]) #as.character(samples[as.numeric(samples[[my.contrast]])-1==0,c("SampleName")]) # control only (AGA or Male), which is background
	rn.case=rownames(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1!=0,])
	            #as.character(samples[as.numeric(samples[[my.contrast]])-1!=0,c("SampleName")]) # case only (SGA or female)
	cat("setting FPKM...\n")
	keep <- rownames(dds) %in% names(my.target.list) # remove entries of no GRange info
	if(table(keep)["TRUE"]==length(keep)){
		    ddsFpkm <-fpkm(dds) # isa 'matrix'
	    ddsFpm <-fpm(dds) # isa 'matrix'
	}else{
		    ddsFpkm <-fpkm(ddsExonicGene) # isa 'matrix'
	    ddsFpm <-fpm(ddsExonicGene) # isa 'matrix'
	}


	#save this expression to list
	this.source="PT" # FG.AGA RNA-Seq (2/Aug/2016)
	dt.tissue.exp.list[[this.source]]<-data.table(
				`ensembl_gene_id`=rownames(ddsFpkm),
				`tissue`=this.source,
				`female.exp`=rowMeans(ddsFpkm[,rn.case]), # mean FPKM of female
				`male.exp`=rowMeans(ddsFpkm[,rn.control])    # mean FPKM of male
			)

	save(dt.tissue.exp.list, file=my.exp.RData) 
}

if(FALSE){
rbindlist(dt.tissue.exp.list)[ensembl_gene_id=="ENSG00000183117",.(tissue,female.exp,male.exp, ratio=female.exp/male.exp)] # CSMD1
rbindlist(dt.tissue.exp.list)[ensembl_gene_id=="ENSG00000176595",.(tissue,female.exp,male.exp, ratio=female.exp/male.exp)] # KBTBD11
rbindlist(dt.tissue.exp.list)[ensembl_gene_id=="ENSG00000022267",.(tissue,female.exp,male.exp, ratio=female.exp/male.exp)] # FHL1
rbindlist(dt.tissue.exp.list)[ensembl_gene_id=="ENSG00000047634",.(tissue,female.exp,male.exp, ratio=female.exp/male.exp)] # SCML1
rbindlist(dt.tissue.exp.list)[ensembl_gene_id=="ENSG00000067992",.(tissue,female.exp,male.exp, ratio=female.exp/male.exp)] # PDK3
rbindlist(dt.tissue.exp.list)[ensembl_gene_id=="ENSG00000229807",.(tissue,female.exp,male.exp, ratio=female.exp/male.exp)] # XIST
}

##############################################
## CPG-density of XCI Genes (Carrel et al.) ##
##############################################
if(FALSE){
	file.name<-file.path("~/results/RoadMap/BS-Seq/cpg.density.chrX.pdf")
	pdf(file=file.name, width=11.7, height=8.3, title="CpG Density @ chrX")
	## Note that ensembl_gene_id can be null <NA> because of left-join (all.x=TRUE)

	# 1. barplot
	ggplot(dt.meth.all.chrX[,list(density=sum(num.sites)/.N),by="region,tissue"], aes(x=region,y=density)) + 
	geom_bar(aes(fill=tissue), stat="identity", position="dodge") + 
	theme_bw() +
	scale_fill_manual(values=cbPalette, name="Tissues") + ggtitle("CpG-density @ chrX")

	# cpg density per 1K bases
	ggplot(dt.meth.all.chrX[,list(cpg.num=sum(num.sites), size=sum(end-start+1), num.cpg.per.1k=sum(num.sites)/sum(end-start+1)*1000, num.region=.N),by="region,tissue"], aes(x=region, y=num.cpg.per.1k)) +
	geom_bar(aes(fill=tissue), stat="identity", position="dodge") +
	theme_bw() +
	scale_fill_manual(values=cbPalette, name="Tissues") + ggtitle("CpG-density @ chrX")

	# 2. boxplot
	ggplot(dt.meth.all.chrX[,list(region, tissue, num.sites)], aes(x=region, y=num.sites)) + 
	geom_boxplot(aes(fill=tissue)) + 
	scale_fill_manual(values=cbPalette, name="Tissues") + 
	ggtitle("CpG-density @ chrX") + 
	theme_bw() +
	scale_y_log10()

	dev.off()
}
  
######################################
## 1. XIE Expression difference     ##
## log2(fc) of FPKM across tissue   ##
## log2(fpkm.female/fepml.male)     ##
######################################
if(FALSE){
	file.name<-file.path("~/results/RoadMap/X-inactivation/Heatmap",paste("XCI.expression.by.tissue",time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=11.7, height=8.3, title="XCI Expression by Tissue")

	dt.exp.xci<-merge(data.table::rbindlist(dt.tissue.exp.list), xci.anno, by="ensembl_gene_id")
	## dt.meth.all.chrX: Note that ensembl_gene_id can be null <NA> because of left-join (all.x=TRUE)
	dt.exp.chrX<-merge(data.table::rbindlist(dt.tissue.exp.list), dt.meth.all.chrX[!is.na(ensembl_gene_id),.(type='chrX'),ensembl_gene_id], by="ensembl_gene_id")

	# concat dt.exp.xci and dt.exp.chrX
	#rbind(dt.exp.xci[,.(ensembl_gene_id,tissue,female.exp,male.exp,type)], dt.exp.chrX)

	##############################################
	# boxplot of expression difference by gender #
	##############################################
	my.lim=c(-5,5)
	p<-ggplot(rbind(dt.exp.chrX, dt.exp.xci[,.(ensembl_gene_id,tissue,female.exp,male.exp,type)]), aes(x=type, y=log2(female.exp/male.exp))) + 
		geom_boxplot(aes(fill=tissue), outlier.shape=NA) + 
		scale_fill_manual(values=cbPalette, name="Tissues") +
		geom_hline(yintercept=0) +
		coord_cartesian(ylim=c(-5, 5)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
		ggtitle("XIE") +
		theme_bw() +
		scale_y_continuous(breaks=seq(from=my.lim[1],to=my.lim[2], by=1) )
		#expand_limits(y=c(-5,5))
		#ylim(-5,5)
	print(p)

	#####################################
	## female-male expression dot plot ##
	#####################################
	#ggplot(dt.exp.xci[type=="XIE"], aes(x=female.exp,y=male.exp, col=tissue)) + geom_point(size=3,alpha=.7) + scale_colour_manual(values=cbPalette, name="Tissues") + geom_abline(slope=1)
	#ggplot(dt.exp.xci[type=="XIE"], aes(x=female.exp,y=male.exp, col=tissue)) + geom_point(size=3,alpha=.7) + scale_colour_manual(values=cbPalette, name="Tissues") + scale_y_log10() + scale_x_log10() + geom_abline(slope=1)
	#ggplot(dt.exp.xci[type=="XIE"], aes(x=log2(female.exp),y=log2(male.exp), col=gene)) + geom_point(size=3,alpha=.7) + geom_abline(slope=1)
	p<-ggplot(dt.exp.xci[type=="XIE"], aes(x=log2(female.exp),y=log2(male.exp), col=tissue)) + 
		geom_point(size=3,alpha=.7) + 
		scale_colour_manual(values=cbPalette, name="Tissues") + 
		geom_abline(slope=1) + 
		ggtitle("XIE") +
		theme_bw()
	print(p)

	# male-female expression dot plot for all tissues
	#ggplot(dt.exp.xci[type=="XIE"], aes(x=log2(female.exp),y=log2(male.exp))) + 
	#geom_point(aes(shape=tissue, col=gene),size=3,alpha=.7) + 
	#geom_abline(slope=1)

	# male-female expression dot plot for selected genes across all tissues
	p<-ggplot(dt.exp.xci[type=="XIE" & (gene=="XIST"|gene=="MAOA"|gene=="AVPR2"|gene=="ARSD")], aes(x=log2(female.exp),y=log2(male.exp))) + 
		geom_point(aes(colour=tissue, shape=gene),size=3,alpha=.7) + 
		geom_abline(slope=1) +
		scale_colour_manual(values=cbPalette, name="Tissues") +
		theme_bw() +
		ggtitle("XIE") +
		geom_text(aes(label=dt.exp.xci[type=="XIE" & (gene=="XIST"|gene=="MAOA"|gene=="AVPR2"|gene=="ARSD"),gene]), hjust =-.5, nudge_x=0.05, size=3)
	print(p)

	# male-female expression dot plot for selected tissues only with gene names
	p<-ggplot(dt.exp.xci[type=="XIE" & (tissue=="PT"|tissue=="GA"|tissue=="SX")], aes(x=log2(female.exp),y=log2(male.exp))) + 
		geom_point(aes(col=tissue),size=3,alpha=.7) + 
		geom_abline(slope=1) + 
		scale_colour_manual(values=cbPalette, name="Tissues") +
		theme_bw() +
		ggtitle("XIE") +
		geom_text(aes(label=dt.exp.xci[type=="XIE" & (tissue=="PT"|tissue=="GA"|tissue=="SX"),gene]), hjust =-.5, nudge_x=0.05, size=3)
	print(p)

	#PT only
	p<-ggplot(dt.exp.xci[type=="XIE" & tissue=="PT"], aes(x=log2(female.exp),y=log2(male.exp))) + 
		geom_point(aes(col=tissue),size=3,alpha=.7) + 
		geom_abline(slope=1) + 
		scale_colour_manual(values=cbPalette, name="Tissues") +
		theme_bw() +
		ggtitle("XIE PT") +
		geom_text(aes(label=dt.exp.xci[type=="XIE" & tissue=="PT",gene]), hjust =-.5, nudge_x=-0.05, size=3)
	print(p)

	#########################################
	## Heatmap of XIE Expression log2ratio ##
	#########################################
	#0. female.fpkm-male.fpkm
	dt.xci.fpkm<-dt.exp.xci[type=="XIE",list(gene,tissue,exp.diff=female.exp-male.exp,female.exp,male.exp)]
	# replace 'NaN' and 'Inf' to 'NA'
	# http://stackoverflow.com/questions/12188509/cleaning-inf-values-from-an-r-dataframe
	set(dt.xci.fpkm, which(is.infinite(dt.xci.fpkm[["exp.diff"]]) | is.nan(dt.xci.fpkm[["exp.diff"]])), j = "exp.diff",value =NA)
	df.xci.fpkm=reshape2::dcast(dt.xci.fpkm, gene ~ tissue, value.var="exp.diff") # dcast: long(row-based) to wide(column-based)
	rownames(df.xci.fpkm)<-df.xci.fpkm$gene
	df.xci.fpkm$gene<-NULL

	heatmap.2(as.matrix(df.xci.fpkm), col=hmcol, 
		cellnote=round(as.matrix(df.xci.fpkm),2), notecol="black", 
		trace="none", dendrogram="both", Rowv=T, Colv=T,
		denscol="black", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, 
		key.xlab="fpkm(female-male)", main="Expression Difference of X-inactivation escapee (XIE)")


	#1. log2(female.fpkm)-log2(male.fpkm)
	#https://cran.r-project.org/web/packages/data.table/vignettes/datatable-reshape.html
	#reshape2::dcast(dt.exp.xci[type=="XIE",list(ensembl_gene_id,gene,tissue,log2ratio=log2(female.exp/male.exp))], ensembl_gene_id + gene ~ tissue, value.var="log2ratio")
	dt.xci.fpkm<-dt.exp.xci[type=="XIE",list(gene,tissue,log2ratio=log2((female.exp+0.01)/(male.exp+0.01)),female.exp,male.exp)]
	# replace 'NaN' and 'Inf' to 'NA'
	# http://stackoverflow.com/questions/12188509/cleaning-inf-values-from-an-r-dataframe
	set(dt.xci.fpkm, which(is.infinite(dt.xci.fpkm[["log2ratio"]]) | is.nan(dt.xci.fpkm[["log2ratio"]])), j = "log2ratio",value =NA)
	df.xci.fpkm=reshape2::dcast(dt.xci.fpkm, gene ~ tissue, value.var="log2ratio") # dcast: long(row-based) to wide(column-based)
	rownames(df.xci.fpkm)<-df.xci.fpkm$gene
	df.xci.fpkm$gene<-NULL

	heatmap.2(as.matrix(df.xci.fpkm), col=hmcol, 
		cellnote=round(as.matrix(df.xci.fpkm),2), notecol="black", 
		trace="none", dendrogram="both", Rowv=T, Colv=T,
		denscol="black", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, 
		key.xlab="log2(female/male)", main="Expression of X-inactivation escapee (XIE)")

	#2.female log2(fpkm) heatmap
	df.female.fpkm=reshape2::dcast(dt.xci.fpkm, gene ~ tissue, value.var="female.exp") # dcast: long(row-based) to wide(column-based)
	rownames(df.female.fpkm)<-df.female.fpkm$gene
	df.female.fpkm$gene<-NULL
	df.female.fpkm=log2(df.female.fpkm)
	df.female.fpkm[is.inf.data.frame(df.female.fpkm)]=NA

	heatmap.2(as.matrix(df.female.fpkm), col=hmcol, 
	cellnote=round(as.matrix(df.female.fpkm),2), notecol="black", 
	trace="none", dendrogram="none", Rowv=F, Colv=F,
	denscol="black", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, 
	key.xlab="log2(fpkm)", main="Female expression of X-inactivation escapee (XIE)")

	#3.male log2(fpkm) heatmap
	df.male.fpkm=reshape2::dcast(dt.xci.fpkm, gene ~ tissue, value.var="male.exp") # dcast: long(row-based) to wide(column-based)
	rownames(df.male.fpkm)<-df.male.fpkm$gene
	df.male.fpkm$gene<-NULL
	df.male.fpkm=log2(df.male.fpkm)
	df.male.fpkm[is.inf.data.frame(df.male.fpkm)]=NA

	heatmap.2(as.matrix(df.male.fpkm), col=hmcol, 
	cellnote=round(as.matrix(df.male.fpkm),2), notecol="black", 
	trace="none", dendrogram="none", Rowv=F, Colv=F,
	denscol="black", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, 
	key.xlab="log2(fpkm)", main="Male expression of X-inactivation escapee (XIE)")

	dev.off()
}

#################################
## 2. Heatmap of XIE Meth Diff ##
#################################
if(FALSE){
	for(my.xci in c("XIE","XIC","XI")){
		file.name<-file.path("~/results/RoadMap/X-inactivation/Heatmap",paste(my.xci,"meth.diff.heatmap",time.stamp,"pdf",sep="."))
		pdf(file=file.name, width=8.3, height=11.7, title=paste("Methylation difference by gender: ",my.xci))

		#cpg.contexts=c("Genes","Promo_15.05", "Promo_10.10","CPGi", "CPGi-shores", "Enhancer") # 1/Aug/2016
		for(my.region in dt.meth.xci[,unique(region)]){
			cat(paste("running ",my.xci,"at ",my.region,"...\n"))
			##
			## female-male methylation dot plot
			##
			if(FALSE){
				ggplot(dt.meth.xci[type==my.xci & region==my.region], aes(x=female.c/female.cov,y=male.c/male.cov, col=tissue)) + 
				geom_point(size=3,alpha=.6) + 
				ggtitle(paste("methylation level of ", my.region, "from ", my.xci, "genes")) + 
				geom_abline(slope=1) +
				theme_bw() +
				scale_colour_manual(values=cbPalette, name="Tissues") 
			}
			##################################
			## Heatmap of Methylation level ##
			##################################
			if(TRUE){
				dt.meth.diff<-dt.meth.xci[type==my.xci & region==my.region, list(meth.diff=sum(female.c)/sum(female.cov)-sum(male.c)/sum(male.cov)),by="gene,tissue"]
				# replace 'NaN' and 'Inf' to 'NA'
				# http://stackoverflow.com/questions/12188509/cleaning-inf-values-from-an-r-dataframe
				set(dt.meth.diff, which(is.infinite(dt.meth.diff[["meth.diff"]]) | is.nan(dt.meth.diff[["meth.diff"]])), j = "meth.diff",value =NA)
				# in case of promoter or cpgi, there could be multiple entries for the same gene
				df.meth.diff=reshape2::dcast(dt.meth.diff, gene ~ tissue, fun.aggregate=mean, value.var="meth.diff")
				rownames(df.meth.diff)<-df.meth.diff$gene
				df.meth.diff$gene<-NULL

				df.meth.diff[is.nan.data.frame(df.meth.diff)]=NA

				if(my.xci=="XIE"){
					heatmap.2(as.matrix(df.meth.diff), col=hmcol,
					cellnote=round(as.matrix(df.meth.diff),2), notecol="black", 
					trace="none", dendrogram="both", Rowv=T, Colv=T, 
					denscol="black", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, 
					key.xla="meth(female-male)", main=paste0("Methylation difference (by gender) of  ",my.xci," at ",my.region))
				}else{
					#it failed with hclust as 'NA' generaged by dist
					if((my.xci=="XI" | my.xci=="XIC") & !my.region=="Genes"){
						df.meth.diff<-df.meth.diff[order(rowSums(df.meth.diff, na.rm=T),decreasing=T),]
						heatmap.2(as.matrix(df.meth.diff),col=hmcol, 
						trace="none", dendrogram="none", Rowv=F, Colv=F, 
						denscol="black", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, 
						key.xlab="meth(female-male)", main=paste0("Methylation difference (by gender) of  ",my.xci," at ",my.region))
					}else{
						heatmap.2(as.matrix(df.meth.diff),col=hmcol, 
						trace="none", dendrogram="col", Rowv=T, Colv=T, 
						denscol="black", offsetCol=0.1, offsetRow=-0.2, srtCol=45, keysize=0.7, 
						key.xlab="meth(female-male)", main=paste0("Methylation difference (by gender) of  ",my.xci," at ",my.region))
					}
				}
			}
		}#end of my.region
		dev.off()
	}
}

#########################
## Met Diff Clustering ##
## for all chrX        ##
#########################
if(FALSE){
	file.name<-file.path("~/results/RoadMap/X-inactivation/Heatmap",paste("chrX.meth.diff.heatmap",time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=8.3, height=11.7, title="Clustering based on Methylation Difference @ chrX")

	## dt.meth.all.chrX: Note that ensembl_gene_id can be null <NA> because of left-join (all.x=TRUE)
	for(my.region in dt.meth.all.chrX[,unique(region)]){
		dt.meth.diff<-dt.meth.all.chrX[!is.na(ensembl_gene_id) & region==my.region, list(cpg.num=sum(num.sites), meth.diff=sum(female.c)/sum(female.cov)-sum(male.c)/sum(male.cov)),by="ensembl_gene_id,tissue"]
		df.meth.diff=reshape2::dcast(dt.meth.diff[!is.na(ensembl_gene_id)], ensembl_gene_id  ~ tissue, fun.aggregate=mean, value.var="meth.diff")
		rownames(df.meth.diff)<-df.meth.diff$ensembl_gene_id
		df.meth.diff$ensembl_gene_id<-NULL
		#df.meth.diff$SB<-NULL # remove 'SB' as it has lots of 'NA'
		df.meth.diff[is.nan.data.frame(df.meth.diff)]=NA # replace 'NaN' to 'NA'
		select<-rowSums(!is.na(df.meth.diff)) >= ncol(df.meth.diff)/2 # 'NA' less than half of the sample

		#make a heatmap
		heatmap.2(as.matrix(df.meth.diff[select,]),col=hmcol, 
		trace="none", dendrogram="both", Rowv=T, Colv=T, labRow="",
		denscol="black", offsetCol=0.1, srtCol=45, keysize=0.7, 
		key.xlab="meth(female-male)", main=paste0("Methylation difference (by gender) at ",my.region))
	}
	dev.off()
}

############################## 
## EXP and MET of Placenta  ##
## EXP/MET mearsued by ENSG ##
############################## 
# To match meth.diff and fpkm.diff in chrX
# have a look at very top of this script as well
if(TRUE){
	my.cpg.type="CG" # CG, CHH, CHG
	my.tissue="PT"

	## Load this RoadMap tissue (only for 'get.meth.for.xci')
	#dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)
	#xci.list<-split(dt.xci.anno, dt.xci.anno[,type]) # isa 'list'
	#cpg.contexts=c("Genes","Genes_10.10", "Promo_10.10", "CPGi_Genes_5K")

	##############
	## get Meth ##
	##############
	cpg.contexts=c("Genes", "Promo_10.10", "CPGi")
	## dt.meth.all.chrX: Note that ensembl_gene_id can be null <NA> because of left-join (all.x=TRUE)
	## ensembl_gene_id from RnBeads (Ensembl 73 or gencode 18)
	dt.meth.all.chrX.ensg.PT=dt.meth.all.chrX[tissue==my.tissue & !is.na(ensembl_gene_id), .(f.met=sum(female.c)/sum(female.cov), m.met=sum(male.c)/sum(male.cov), diff.met=sum(female.c)/sum(female.cov)-sum(male.c)/sum(male.cov), num.sites=sum(num.sites)),by="ensembl_gene_id,region"][region %in% cpg.contexts]

	############################
	## get PT expression      ##  
	## see dt.tissue.exp.list ##
	## based on ddsFpkm       ##
	## which is defined above ##
	############################
	# 1. FG.JD.AGA (GRCh37)
	myCaller='DESeq2'
	source("~/Pipelines/config/DEG.R") # load config
	deseq.RData <- paste0(deseq.dir, "/deseq.",sampleType,'.RData')
	load(deseq.RData) 
	cat("dds, res, rld, and (ddsExonicGene?) loaded\n")

	my.res<-res; my.res$ensembl_gene_id=rownames(my.res); dt.res=as.data.table(as.data.frame(my.res))
	#dt.exp.pops=merge(dt.tissue.exp.list[["PT"]], dt.res, by="ensembl_gene_id") # dt.tissue.exp.list based on ddsFpkm
	dt.exp.pops=dt.res # regardless of FPKM values

	# 2. JD DEG (GRCh38)
	#dt.jd=fread("~/Pipelines/data/JD/gender_DEG_JD.csv", header=T)

	#####################
	## Merge Met & Exp ##
	#####################
	## ensembl_gene_id from RnBeads (Ensembl 73 or gencode 18)
	dt.meth.exp=merge(dt.meth.all.chrX.ensg.PT, dt.exp.pops, by="ensembl_gene_id") # based on FG

	#dt.meth.exp.jd=merge(dt.meth.all.chrX.ensg.PT, dt.jd, by="ensembl_gene_id", all.x=TRUE) # based on JD
	#dt.meth.exp.jd[, `Within_JD_DEG?` := ifelse(is.na(pval),"No","Yes")]

	# No meth info
	#dt.jd[chr=="X" & ! ensembl_gene_id %in% dt.meth.all.chrX.ensg.PT[,unique(ensembl_gene_id)]]
	#   ensembl_gene_id    pval    padj  hgnc chr                                      desc
	#1: ENSG00000002586 1.5e-08 3.7e-06  CD99   X                            CD99 molecule 
	#2: ENSG00000124333 9.9e-06 2.0e-03 VAMP7   X    vesicle associated membrane protein 7 
	#3: ENSG00000169093 0.0e+00 7.7e-02 ASMTL   X acetylserotonin O-methyltransferase-like 

	#dt.meth.exp[,.N,.(region)]
	#dt.meth.exp[,.N,.(region,`is.na(padj)`=is.na(padj))][order(region)]
	#dt.meth.exp[padj<0.01,.N,.(region,`log2FoldChange>0`=log2FoldChange>0)][order(region)]

	##################################################
	## Plot Meth-Diff at CPGi (or promoters) & Gene ##
	##################################################
	if(TRUE){
		file.name<-file.path("~/results/RoadMap/X-inactivation/Heatmap",paste("PT.meth.exp.FG.JD.chrX",TR_PREFIX,time.stamp,"pdf",sep="."))
		pdf(file=file.name, width=12, height=8, title="% Diff Meth and Expression")

		min.num.sites=4 # at least this number of CG at a specificed region
		dt.foo<-as.data.table(reshape2::dcast(dt.meth.all.chrX.ensg.PT[num.sites>=min.num.sites], ensembl_gene_id ~ region, value.var="diff.met")) # dcast: long(row-based) to wide(column-based)
		dt.foo[,CPGi:=CPGi*100]
		dt.foo[,Genes:=Genes*100]
		dt.foo[,Promo_10.10:=Promo_10.10*100]
		## dt.foo[,ensembl_gene_id] from RnBeads (Ensembl 73 or gencode 18)
		## dt.exp.pops[,ensembl_gene_id] from DESeq2 (Ensembl 75 or gencode 19?)
		dt.meth.exp.fg<-merge(dt.foo, dt.exp.pops, by="ensembl_gene_id")
		dt.meth.exp.fg[, `Expression?` := ifelse(log2FoldChange<0, 'female>male', 'female<male')]
		dt.meth.exp.fg[,`new.padj` := p.adjust(pvalue, method="BH")]

		####################
		## CPGi vs. Genes ##
		####################
		# call the function 'PlotMethExp'
		PlotMethExp(my.padj=1,my.target=c("CPGi","Genes")) # chrX
		PlotMethExp(my.padj=0.5,my.target=c("CPGi","Genes")) # chrX
		PlotMethExp(my.padj=0.1,my.target=c("CPGi","Genes")) # chrX
		PlotMethExp(my.padj=0.01,my.target=c("CPGi","Genes")) # chrX

		PlotMethExp(my.padj=1,my.target=c("Promo_10.10","Genes")) # chrX
		PlotMethExp(my.padj=0.5,my.target=c("Promo_10.10","Genes")) # chrX
		PlotMethExp(my.padj=0.1,my.target=c("Promo_10.10","Genes")) # chrX
		PlotMethExp(my.padj=0.01,my.target=c("Promo_10.10","Genes")) # chrX
		dev.off()

		# or plot directly via ggplot (see below)
		if(FALSE){
			p<-ggplot(dt.meth.exp.fg, aes(CPGi, Genes)) + 
				geom_point(aes(col=padj),size=3,alpha=.9) + 
				ggtitle("% Diff Meth (chrX)") +
				labs(x="CPGi", y="Genes") +
				geom_hline(yintercept=0) + geom_vline(xintercept=0) +
				theme_Publication()
			print(p)

			# one-to-many: (1-ensg-to-2-hgnc)
			fields <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "gene_biotype", "description")
			# padj<0.5
			dt.meth.exp.fg.5<-merge(dt.meth.exp.fg[padj<0.5], getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.meth.exp.fg[padj<0.5,unique(ensembl_gene_id)], mart = myMart), by="ensembl_gene_id")
			p<-ggplot(dt.meth.exp.fg.5[!is.na(CPGi) & !is.na(Genes)], aes(CPGi, Genes)) + 
				geom_point(aes(col=padj,shape=`Expression?`),size=3,alpha=.9) + 
				ggtitle("% Diff Meth (chrX & padj<0.5)") +
				labs(x="CPGi", y="Genes") +
				geom_hline(yintercept=0) + geom_vline(xintercept=0) +
				geom_text(aes(label=dt.meth.exp.fg.5[!is.na(CPGi) & !is.na(Genes),hgnc_symbol]), hjust =1.2, nudge_x=0.05, size=4) +
				theme_Publication()
			print(p)

			# padj<0.1
			dt.meth.exp.fg.1<-merge(dt.meth.exp.fg[padj<0.1], getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.meth.exp.fg[padj<0.1,unique(ensembl_gene_id)], mart = myMart), by="ensembl_gene_id")
			p<-ggplot(dt.meth.exp.fg.1[!is.na(CPGi) & !is.na(Genes)], aes(CPGi, Genes)) + 
				geom_point(aes(col=padj,shape=`Expression?`),size=3,alpha=.9) + 
				ggtitle("% Diff Meth (chrX & padj<0.1)") +
				labs(x="CPGi", y="Genes") +
				geom_hline(yintercept=0) + geom_vline(xintercept=0) +
				geom_text(aes(label=dt.meth.exp.fg.1[!is.na(CPGi) & !is.na(Genes),hgnc_symbol]), hjust =1.2, nudge_x=0.05, size=4) +
				theme_Publication()
			print(p)

			# padj<0.01
			dt.meth.exp.fg.01<-merge(dt.meth.exp.fg[padj<0.01], getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.meth.exp.fg[padj<0.01,unique(ensembl_gene_id)], mart = myMart), by="ensembl_gene_id")
			p<-ggplot(dt.meth.exp.fg.01[!is.na(CPGi) & !is.na(Genes)], aes(CPGi, Genes)) + 
				geom_point(aes(col=padj,shape=`Expression?`),size=3,alpha=.9) + 
				ggtitle("% Diff Meth (chrX & padj<0.01)") +
				labs(x="CPGi", y="Genes") +
				geom_hline(yintercept=0) + geom_vline(xintercept=0) +
				geom_text(aes(label=dt.meth.exp.fg.01[!is.na(CPGi) & !is.na(Genes),hgnc_symbol]), hjust =1.2, nudge_x=0.05, size=4) +
				theme_Publication()
			print(p)

			#########################
			## Promoters vs. Genes ##
			#########################
			# chrX
			p<-ggplot(dt.meth.exp.fg, aes(Promo_10.10, Genes)) + 
				geom_point(aes(col=padj),size=3,alpha=.9) + 
				ggtitle("% Diff Meth (chrX)") +
				labs(x="Promo_10.10", y="Genes") +
				geom_hline(yintercept=0) + geom_vline(xintercept=0) +
				theme_Publication()
			print(p)

			# padj<0.5
			dt.meth.exp.fg.5<-merge(dt.meth.exp.fg[padj<0.5], getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.meth.exp.fg[padj<0.5,unique(ensembl_gene_id)], mart = myMart), by="ensembl_gene_id")
			p<-ggplot(dt.meth.exp.fg.5[!is.na(Promo_10.10) & !is.na(Genes)], aes(Promo_10.10, Genes)) + 
				geom_point(aes(col=padj,shape=`Expression?`,size=2^abs(log2FoldChange)),alpha=.9) + 
				ggtitle("% Diff Meth (chrX & padj<0.5)") +
				labs(x="Promo_10.10", y="Genes") +
				geom_hline(yintercept=0) + geom_vline(xintercept=0) +
				geom_text(aes(label=dt.meth.exp.fg.5[!is.na(Promo_10.10) & !is.na(Genes),hgnc_symbol]), hjust =1.2, nudge_x=0.05, size=4) +
				theme_Publication()
			print(p)

			# padj<0.1
			dt.meth.exp.fg.1<-merge(dt.meth.exp.fg[padj<0.1], getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.meth.exp.fg[padj<0.1,unique(ensembl_gene_id)], mart = myMart), by="ensembl_gene_id")
			p<-ggplot(dt.meth.exp.fg.1[!is.na(Promo_10.10) & !is.na(Genes)], aes(Promo_10.10, Genes)) + 
				#geom_point(aes(col=padj,shape=`Expression?`),size=3,alpha=.9) + 
				geom_point(aes(col=padj,shape=`Expression?`,size=2^abs(log2FoldChange)),alpha=.9) + 
				ggtitle("% Diff Meth (chrX & padj<0.1)") +
				labs(x="Promo_10.10", y="Genes") +
				geom_hline(yintercept=0) + geom_vline(xintercept=0) +
				geom_text(aes(label=dt.meth.exp.fg.1[!is.na(Promo_10.10) & !is.na(Genes),hgnc_symbol]), hjust =1.2, nudge_x=0.05, size=4) +
				theme_Publication()
			print(p)

			# padj<0.01
			dt.meth.exp.fg.01<-merge(dt.meth.exp.fg[padj<0.01], getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.meth.exp.fg[padj<0.01,unique(ensembl_gene_id)], mart = myMart), by="ensembl_gene_id")
			p<-ggplot(dt.meth.exp.fg.01[!is.na(Promo_10.10) & !is.na(Genes)], aes(Promo_10.10, Genes)) + 
				geom_point(aes(col=padj,shape=`Expression?`,size=2^abs(log2FoldChange)),alpha=.9) + 
				ggtitle("% Diff Meth (chrX & padj<0.01)") +
				labs(x="Promo_10.10", y="Genes") +
				geom_hline(yintercept=0) + geom_vline(xintercept=0) +
				geom_text(aes(label=dt.meth.exp.fg.01[!is.na(Promo_10.10) & !is.na(Genes),hgnc_symbol]), hjust =1.2, nudge_x=0.05, size=4) +
				theme_Publication()
			print(p)
		}
	}

	if(FALSE){
		###############################
		## Merge with XCI Annotation ##
		###############################
		dt.meth.exp.carrel<-merge(dt.meth.exp, dt.xci.anno)

		# XIE exp(female) < exp (male)
		dt.meth.exp.carrel[type=="XIE" & region=="Genes" & log2FoldChange>0] # 6 DEGs
		# XIE of rarely expressed genes
		dt.meth.exp.carrel[type=="XIE" & region=="Genes" & (log2FoldChange>0 | is.na(padj))] # not expressed enough or female less expressed than male

		##############################
		## Plot Meth-Exp of PT chrX ##
		##############################
		file.name<-file.path("~/results/RoadMap/X-inactivation/Heatmap",paste("PT.meth.exp",time.stamp,"pdf",sep="."))
		pdf(file=file.name, width=11.7, height=8.3, title="XCI Expression by Tissue")

		min.num.sites=3
		# all chrX genes
		p1<-ggplot(dt.meth.exp[num.sites>min.num.sites], aes(diff.met*100, -log2FoldChange)) + 
			geom_point(aes(col=padj),size=3) + 
			geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
			ggtitle("Methylation and Expression Difference of ChrX Genes by Gender") +
			labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
			theme_bw() +
			#facet_wrap(~region, nrow=2)
			facet_wrap(~region)
		#print(p1)

		# all chrX expressed genes (padj not null)
		p2<-ggplot(dt.meth.exp[num.sites>min.num.sites & !is.na(padj)], aes(diff.met*100, -log2FoldChange)) + 
			geom_point(aes(col=padj),size=3) + 
			geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
			ggtitle("Methylation and Expression Difference of ChrX (expressed) Genes by Gender") +
			labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
			theme_bw() +
			#facet_wrap(~region, nrow=2)
			facet_wrap(~region)
		#print(p2)

		multiplot(p1,p2)

		# chrX genes of padj less than 0.5
		p1<-ggplot(dt.meth.exp[num.sites>min.num.sites & padj<.5], aes(diff.met*100, -log2FoldChange)) + 
			geom_point(aes(col=padj),size=4) + 
			geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
			ggtitle("Methylation and Expression Difference of ChrX Genes by Gender (p-adj<0.5)") +
			labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
			theme_bw() +
			#facet_wrap(~region, nrow=2)
			facet_wrap(~region)
		#print(p1)

		# chrX genes of padj less than 0.1
		p2<-ggplot(dt.meth.exp[num.sites>min.num.sites & padj<.1], aes(diff.met*100, -log2FoldChange)) + 
			geom_point(aes(col=padj),size=4) + 
			geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
			ggtitle("Methylation and Expression Difference of ChrX Genes by Gender (p-adj<0.1)") +
			labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
			theme_bw() +
			#facet_wrap(~region, nrow=2)
			facet_wrap(~region)
		#print(p2)
		multiplot(p1,p2)

		# chrX genes of padj less than 0.05
		p1<-ggplot(dt.meth.exp[num.sites>min.num.sites & padj<.05], aes(diff.met*100, -log2FoldChange)) + 
			geom_point(aes(col=padj),size=4) + 
			geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
			ggtitle("Methylation and Expression Difference of ChrX Genes by Gender (p-adj<0.05)") +
			labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
			theme_bw() +
			#facet_wrap(~region, nrow=2)
			facet_wrap(~region)
		#print(p1)
		# chrX genes of padj less than 0.01
		p2<-ggplot(dt.meth.exp[num.sites>min.num.sites & padj<.01], aes(-log2FoldChange, diff.met*100)) + 
			geom_point(aes(col=padj),size=3) + 
			#geom_point(size=4,col="navy") + 
			geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
			#ggtitle("Methylation and Expression Difference of ChrX Genes by Gender (p-adj<0.01)") +
			labs(y="methylation difference (%)", x="log2(Fold-change)") +
			#theme_bw() +
			theme_Publication() +
			#facet_wrap(~region, nrow=2)
			facet_wrap(~region)
		#print(p2)
		multiplot(p1,p2)

		# chrX genes reported in Carrel et al. Nature  2005
		p<-ggplot(dt.meth.exp.carrel[num.sites>min.num.sites & type %in% c("XIC","XIE")], aes(diff.met*100, -log2FoldChange)) + 
			geom_point(aes(col=padj),size=3) + 
			geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
			ggtitle("Methylation and Expression Difference of Genes reported by Carrel et al.") +
			labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
			theme_bw() +
			facet_grid(type ~region)
		print(p)

		#dev.off()
		############################
		## Meth Boxplot by Gender ##
		############################
		# 1. FG
		#file.name<-file.path("~/results/RoadMap/X-inactivation/Heatmap",paste("meth.level.FG.DEG",time.stamp,"tiff",sep="."))
		#tiff(filename=file.name, width=11.7, height=8.27, units="in", res=300, compression='lzw')
		dt.meth.fg.gender=rbind(
							dt.meth.exp[log2FoldChange<0,.(`Gender`="Female",`met`=f.met, diff.met, num.sites, baseMean, log2FoldChange, pvalue, padj),by="ensembl_gene_id,region"],
							dt.meth.exp[log2FoldChange<0,.(`Gender`="Male",`met`=m.met, diff.met, num.sites, baseMean, log2FoldChange, pvalue, padj),by="ensembl_gene_id,region"]
						)
		dt.meth.fg.gender[, `IS_DEG?` := ifelse(baseMean<=1, 'NO_EXP', ifelse(padj<0.01 ,"DEG","NO_DEG"))]

		p1<-ggplot(dt.meth.fg.gender[num.sites>min.num.sites], aes(`IS_DEG?`,met*100)) + 
			geom_boxplot(aes(fill=Gender),show.legend = NA) + 
			scale_fill_manual(values=my.col[["Gender"]]) +	
			#ggtitle("Methylation Level of Differentially Expressed ChrX Genes by Gender (FG, p-adj<0.01 & f>m)") +
			labs(y="% methylation",x="") +
			#theme_bw() +
			theme_Publication() +
			facet_wrap(~region)
		print(p1)
		#dev.off()

		file.name<-file.path("~/results/RoadMap/X-inactivation/FG.DEG.gender.meth.exp.tiff")
		tiff(filename=file.name, width=11.7, height=8.27, units="in", res=300, compression='lzw')

		multiplot(p2,p1)

		dev.off()

		# 2. JD
		if(FALSE){
			#file.name<-file.path("~/results/RoadMap/X-inactivation/Heatmap",paste("meth.level.XIE",time.stamp,"tiff",sep="."))
			#tiff(filename=file.name, width=11.7, height=8.27, units="in", res=300, compression='lzw')
			dt.meth.jd.gender=rbind(
									dt.meth.exp.jd[,.(`Gender`="Female",`met`=f.met, diff.met, num.sites),by="ensembl_gene_id,region,Within_JD_DEG?"],
									dt.meth.exp.jd[,.(`Gender`="Male",`met`=m.met, diff.met, num.sites),by="ensembl_gene_id,region,Within_JD_DEG?"])

			p2<-ggplot(dt.meth.jd.gender[num.sites>min.num.sites], aes(`Within_JD_DEG?`,met*100)) + 
				geom_boxplot(aes(fill=Gender)) + 
				scale_fill_manual(values=my.col[["Gender"]]) +	
				ggtitle("Methylation Level of Differentially Expressed ChrX Genes by Gender (JD, f>m)") +
				theme_bw() +
				facet_wrap(~region)
			print(p2)

			multiplot(p1,p2)
			dev.off()
		}

		#> dt.meth.jd.gender[`Within_JD_DEG?`=="Yes",length(unique(ensembl_gene_id))]
		#[1] 32
		#> dt.meth.fg.gender[`IS_DEG?`=="Yes",length(unique(ensembl_gene_id))]
		#[1] 23
		#> length(intersect(dt.meth.jd.gender[`Within_JD_DEG?`=="Yes",unique(ensembl_gene_id)], dt.meth.fg.gender[`IS_DEG?`=="Yes",unique(ensembl_gene_id)]))
		#[1] 20

		###################
		## get gene name ##
		###################
		# one-to-many: (1-ensg-to-2-hgnc)
		fields <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "gene_biotype", "description")
		#fields <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "gene_biotype")
		df.gene<-getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.res[padj<0.01,unique(ensembl_gene_id)], mart = myMart) # is a data.frame
		#df.gene<-getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.meth.exp[,unique(ensembl_gene_id)], mart = myMart) # is a data.frame
		dt.deg.01<-merge(as.data.table(df.gene), dt.res, by="ensembl_gene_id")
		
		df.gene<-getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.res[padj<0.05,unique(ensembl_gene_id)], mart = myMart) # is a data.frame
		dt.deg.05<-merge(as.data.table(df.gene), dt.res, by="ensembl_gene_id")

		df.gene<-getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.res[padj<0.1,unique(ensembl_gene_id)], mart = myMart) # is a data.frame
		dt.deg.1<-merge(as.data.table(df.gene), dt.res, by="ensembl_gene_id")

		dt.deg.01[,.N,.(chromosome_name=="X")]
		#   chromosome_name  N
		#1:            TRUE 29
		#2:           FALSE 11
		dt.res[ensembl_gene_id %in% dt.deg.01[chromosome_name=="X",ensembl_gene_id],.N,.(`log2FoldChange>0`=log2FoldChange>0)]
		#   log2FoldChange>0  N
		#1:             TRUE  6
		#2:            FALSE 23
		dt.meth.all.chrX.ensg.PT[ensembl_gene_id %in% dt.deg.01[chromosome_name=="X",ensembl_gene_id],.N,region]
		#        region  N
		#1:       Genes 26
		#2: Promo_10.10 24
		#3:        CPGi 16
		dt.meth.exp[padj<0.01,.N,.(region,`log2FoldChange>0`=log2FoldChange>0)][order(region)]
		#        region log2FoldChange>0  N
		#1:        CPGi            FALSE 13
		#2:        CPGi             TRUE  3
		#3:       Genes            FALSE 22
		#4:       Genes             TRUE  4
		#5: Promo_10.10            FALSE 21
		#6: Promo_10.10             TRUE  3

		## missing 3 genes from 29 dt.deg.01[chromosome_name=="X"]
		## they do not have methyhlation info
		merge(dt.deg.01[chromosome_name=="X"], dt.meth.all.chrX.ensg.PT[region=="Genes"], by="ensembl_gene_id", all.x=TRUE)

		###########################
		## Save Meth & Exp Data  ##
		###########################
		# all chrX genes
		dt.foo<-as.data.table(reshape2::dcast(dt.meth.exp[, .(ensembl_gene_id, region, diff.met.pct=diff.met*100)], ensembl_gene_id ~ region, value.var="diff.met.pct"))
		setnames(dt.foo,c("ensembl_gene_id","wgbs.cg.diff.pct.cpgi","wgbs.cg.diff.pct.gene","wgbs.cg.diff.pct.promoter"))

		dt.fg.met<-merge(
			dt.foo,
			dt.meth.exp[, .(FG.AGA.log2fc=mean(-log2FoldChange), FG.AGA.pval=mean(pvalue), FG.AGA.padj=mean(padj)), ensembl_gene_id],
			by="ensembl_gene_id"
			)
		write.csv(dt.fg.met, file=file.path("~/results/RoadMap/X-inactivation/Heatmap",paste("PT.meth.exp",time.stamp,"csv",sep=".")))

		# DEG (padj<0.01)
		dt.exp.pops[padj<0.01]
		fields <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "gene_biotype", "description")
		df.gene<-getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.exp.pops[padj<0.01,ensembl_gene_id], mart = myMart) # is a data.frame
		dt.deg.fg<-merge(as.data.table(df.gene), dt.exp.pops[padj<0.01], by="ensembl_gene_id")

		# Meth female
		#dt.meth.female<-reshape2::dcast(dt.meth.exp[min.num.sites & padj<0.01, .(ensembl_gene_id, region, f.met.pct=round(f.met*100,2))], ensembl_gene_id ~ region, value.var="f.met.pct")
		dt.meth.female<-reshape2::dcast(dt.meth.exp[padj<0.01, .(ensembl_gene_id, region, f.met.pct=round(f.met*100,2))], ensembl_gene_id ~ region, value.var="f.met.pct")
		setnames(dt.meth.female,c("ensembl_gene_id","female.CPGi.meth","female.Genes.meth","female.Promo_10.10.meth"))

		# Meth male
		#dt.meth.male<-reshape2::dcast(dt.meth.exp[num.sites>min.num.sites & padj<0.01, .(ensembl_gene_id, region, m.met.pct=round(m.met*100,2))], ensembl_gene_id ~ region, value.var="m.met.pct")
		dt.meth.male<-reshape2::dcast(dt.meth.exp[padj<0.01, .(ensembl_gene_id, region, m.met.pct=round(m.met*100,2))], ensembl_gene_id ~ region, value.var="m.met.pct")
		setnames(dt.meth.male,c("ensembl_gene_id","male.CPGi.meth","male.Genes.meth","male.Promo_10.10.meth"))

		# num.sites (No. of CG) 
		dt.meth.num<-as.data.table(reshape2::dcast(dt.meth.exp[padj<0.01, .(ensembl_gene_id, region, num.sites)], ensembl_gene_id ~ region, value.var="num.sites"))
		setnames(dt.meth.num,c("ensembl_gene_id","CPGi.num.cg","Genes.num.cg","Promo_10.10.num.cg"))

		# p-value & adjustp
		dt.pval<-dt.exp.pops[padj<0.01 & ensembl_gene_id %in% dt.meth.num[,ensembl_gene_id],.(ensembl_gene_id,pvalue,padj)]

		#merge(dt.deg.fg[chromosome_name=="X",.(ensembl_gene_id,hgnc_symbol,gene_biotype,description,female.exp,male.exp,baseMean,log2FoldChange)], merge(dt.meth.female, dt.meth.male))
		dt.deg.meth<-merge(dt.deg.fg[chromosome_name=="X",.(ensembl_gene_id,hgnc_symbol,gene_biotype,baseMean,log2fc=-log2FoldChange, female.fpkm=female.exp,male.fpkm=male.exp)], merge(merge(merge(dt.meth.female, dt.meth.male), dt.meth.num), dt.pval), all.x=TRUE)
		write.csv(dt.deg.meth, file=file.path("~/results/RoadMap/X-inactivation/Heatmap",paste("PT.FG.DEG.meth.exp",time.stamp,"csv",sep=".")))


		####################################
		## Venn Diagram between DEG & XIE ##
		####################################
		dt.deg.list=list()
		#dt.deg.list=lapply(split(dt.xci.anno, dt.xci.anno[,type]), function(i) i[,ensembl_gene_id])
		dt.deg.list[["XIE"]]=dt.xci.anno[type=="XIE",ensembl_gene_id]
		dt.deg.list[["FG.DEG.01"]]=dt.deg.01[chromosome_name=="X",ensembl_gene_id]
		my.venn.file="~/results/RoadMap/X-inactivation/Heatmap/venn.fg.deg.xci.tiff"
		venn.diagram(
			x=dt.deg.list,
			filename = my.venn.file,
			col = "black",
			fill = c("dodgerblue", "goldenrod1"),
			alpha = 0.50,
			cat.col = c("dodgerblue", "goldenrod1"),
			cat.cex = 1.1,
			cat.fontface = "bold",
			margin = 0.05
		)

		#merge(dt.meth.exp, df.gene, by="ensembl_gene_id")

		#write.csv(merge(dt.fg.met, dt.gene, by="ensembl_gene_id")[order(FG.AGA.pval)], file=file.path("~/results/RoadMap/X-inactivation/Heatmap",paste("PT.meth.exp",time.stamp,"csv",sep=".")))
	}
}

cat("All done\n")
