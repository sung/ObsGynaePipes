#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# based on RoadMap/plot.met.exp.by.gender.R
# First created 30/Mar/2016
# Last modified 13/Jan/2017

# use bin/R/chrX/plot.meth.exp.diff.R

TR_PREFIX="GRCh37"
RNASEQ_TYPE="SMALL" # SMALL (miRNA)
METH_TYPE="SS.Bio"
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/chrX/local.R")# defines time.stamp
										# dt.xci.anno, dt.meth.chrX, dt.exp.pops, dt.exp.pops.chrX
file.name<-file.path("~/results/chrX",paste("PT.meth.exp.FG.JD.chrX",RNASEQ_TYPE,METH_TYPE,time.stamp,"pdf",sep="."))
pdf(file=file.name, width=12, height=12, title=paste("% Diff Meth and Expression",RNASEQ_TYPE))

##
## At least 10% difference and at least 10x coverage
#dt.balaton.deg[abs(log2FoldChange)>log2(1.1) & baseMean>10] 
#t(table(dt.balaton.deg[abs(log2FoldChange)>log2(1.1) & baseMean>10,list(XCI,`Expression?`)]))

##
## There is no meth info for the PAR genes
##dt.meth.chrX[ensembl_gene_id %in% dt.balaton.deg[grepl("PAR",Y_homology),ensembl_gene_id]]

if(RNASEQ_TYPE=="TOTAL"){
	## Methylation level of placenta chrX DEG (XCI annotation from Balaton)
	dt.balaton.deg.meth<-merge(dt.balaton.deg, dt.meth.chrX, by="ensembl_gene_id")
	#t(table(dt.balaton.deg.meth[tissue=="PT" & region=="Promo_10.10" & abs(log2FoldChange)>log2(1.1) & baseMean>10 & num.sites>=4,list(XCI,`Expression?`)]))

	##############################
	## X-axis: % Meth of female ##
	## Y-axis: % Meth of male   ##
	## DEG (FDR<0.01) only      ##
	## dt.balaton.deg.meth      ##
	## needs to be pre-defined  ##
	##############################
	#plotMethXY(METH_TYPE, "Promoter", my.padj=1, min.depth=0, min.num.sites)
	#plotMethXY(METH_TYPE, "Promoter", my.padj=0.1, min.depth=0, min.num.sites)
	plotMethXY(METH_TYPE, "Promoter", my.padj=0.01, min.depth=10, min.num.sites) # DEG only (with Balaton XCI)

	#plotMethXY(METH_TYPE, "CPGi", my.padj=1, min.depth=0, min.num.sites)
	#plotMethXY(METH_TYPE, "CPGi", my.padj=0.1, min.depth=10, min.num.sites)
	plotMethXY(METH_TYPE, "CPGi", my.padj=0.01, min.depth=10, min.num.sites) # DEG only

	#plotMethXY(METH_TYPE, "Gene-body", my.padj=1, min.depth=0, min.num.sites)
	#plotMethXY(METH_TYPE, "Gene-body", my.padj=0.01, min.depth=10, min.num.sites)
	plotMethXY(METH_TYPE, "Gene-body", my.padj=0.01, min.depth=10, min.num.sites) # DEG only

	##################################
	## X-axis: XCI                  ##
	## Y-axis: Boxplot of meth.diff ##
	## DEG (FDR<0.01) only          ##
	## dt.balaton.deg.meth          ##
	## needs to be pre-defined      ##
	##################################
	plotBoxBalatonMethDiff(METH_TYPE, "Promoter", min.depth=10, min.num.sites)
	plotBoxBalatonMethDiff(METH_TYPE, "CPGi", min.depth=10, min.num.sites)
	plotBoxBalatonMethDiff(METH_TYPE, "Gene-body", min.depth=10, min.num.sites)

	##################################
	## X-axis: XCI                  ##
	## Y-axis: % methylation by sex ##
	## DEG (FDR<0.01) only          ##
	## dt.balaton.deg.meth         ## 
	## needs to be pre-defined      ##
	##################################
	#plotBoxBalatonMethbySex("WGoxBS", "Promoter", min.depth=10, min.num.sites)
	#plotBoxBalatonMethbySex("WGoxBS", "CPGi", min.depth=10, min.num.sites)
	#plotBoxBalatonMethbySex("WGoxBS", "Gene-body", min.depth=10, min.num.sites)

	#####################################
	## X-axis: Diff-Meth (female-male) ##
	## Y-axis: Diff-Meth (female-male) ##
	## dt.meth.diff.exp                ## 
	## must be pre-defined             ##
	## dt.balaton used if padj<0.01    ##
	## min.num.site already applied >4 ##
	#####################################
	# 1. CPGi vs. Genes
	#plotMethExp(METH_TYPE, min.depth=0, my.padj=1, my.target=c("CPGi","Gene-body")) # chrX
	#plotMethExp(METH_TYPE, min.depth=0, my.padj=0.1, my.target=c("CPGi","Gene-body")) # chrX
	plotMethExp(METH_TYPE, min.depth=10, my.padj=0.01, my.target=c("CPGi","Gene-body")) # DEG only 
	# 2. Promo vs. Genes
	#plotMethExp(METH_TYPE, min.depth=0, my.padj=1, my.target=c("Promoter","Gene-body")) # chrX
	#plotMethExp(METH_TYPE,, min.depth=0, my.padj=0.1, my.target=c("Promoter","Gene-body")) # chrX
	plotMethExp(METH_TYPE, my.padj=0.01, my.target=c("Promoter","Gene-body")) # DEG only 
}else if(RNASEQ_TYPE=="SMALL"){
	##############################
	## X-axis: % Meth of female ##
	## Y-axis: % Meth of male   ##
	## DEG (FDR<0.01) only      ##
	## dt.balaton.deg.meth      ##
	## needs to be pre-defined  ##
	##############################
	plotMethXY(METH_TYPE, "miR Promoter", my.padj=1, min.depth=0, min.num.sites)
	plotMethXY(METH_TYPE, "miR Promoter", my.padj=0.6, min.depth=0, min.num.sites)
	plotMethXY(METH_TYPE, "miR Promoter", my.padj=0.4, min.depth=0, min.num.sites)
	plotMethXY(METH_TYPE, "miR Promoter", my.padj=0.1, min.depth=0, min.num.sites)

	plotMethXY(METH_TYPE, "miR Gene-body", my.padj=1, min.depth=0, min.num.sites)
	plotMethXY(METH_TYPE, "miR Gene-body", my.padj=0.6, min.depth=0, min.num.sites)
	plotMethXY(METH_TYPE, "miR Gene-body", my.padj=0.4, min.depth=0, min.num.sites)
	plotMethXY(METH_TYPE, "miR Gene-body", my.padj=0.1, min.depth=0, min.num.sites)

	plotMethXY(METH_TYPE, "miR CPGi", my.padj=1, min.depth=0, min.num.sites)
	plotMethXY(METH_TYPE, "miR CPGi", my.padj=0.6, min.depth=0, min.num.sites)
	plotMethXY(METH_TYPE, "miR CPGi", my.padj=0.4, min.depth=0, min.num.sites)
	plotMethXY(METH_TYPE, "miR CPGi", my.padj=0.1, min.depth=0, min.num.sites)

	#####################################
	## X-axis: Diff-Meth (female-male) ##
	## Y-axis: Diff-Meth (female-male) ##
	## dt.meth.diff.exp                ## 
	## must be pre-defined             ##
	## dt.balaton used if padj<0.01    ##
	## min.num.site already applied >4 ##
	#####################################
	# 1. CPGi vs. Genes
	plotMethExp(METH_TYPE, min.depth=0, my.padj=1, my.target=c("miR CPGi","miR Gene-body")) # chrX
	# 2. Promo vs. Genes
	plotMethExp(METH_TYPE, min.depth=0, my.padj=1, my.target=c("miR Promoter","miR Gene-body")) # chrX
}else{
	stop(paste(RNASEQ_TYPE, "Not supported\nEither TOTAL or SMALL"))
}


######################################################################################
## ChrX Methylation Difference by Tissue                                            ##
## see also "bin/R/RoadMap/plot.roadmap.gender.meth.diff.boxplot.tissue.cpgi.etc.R" ##
######################################################################################
if(FALSE){
	load("~/results/RoadMap/X-inactivation/RData/dt.meth.chrX.per.tissue.CG.Genes_20.20.RData")
	dt.meth.chrX.all<-rbindlist(dt.meth.chrX.per.tissue)[!is.na(ensembl_gene_id), .(f.met=sum(female.c)/sum(female.cov)*100, m.met=sum(male.c)/sum(male.cov)*100, diff.met=sum(female.c)/sum(female.cov)*100-sum(male.c)/sum(male.cov)*100, num.sites=sum(num.sites)), by="tissue,ensembl_gene_id,region"]
	p<-ggplot(dt.meth.chrX.all[region!="Promo_15.05"], aes(tissue, diff.met)) + 
		geom_boxplot(aes(fill=region),outlier.shape=NA,alpha=0.7) +
		coord_cartesian(ylim=c(-50,80)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
		scale_y_continuous(breaks=seq(from=-50,to=80,by=20)) + # Ticks every .2
		geom_hline(yintercept=0) +
		#scale_fill_manual(name="CpG\nContext", labels=c("CPGi","CPGi-shore","Enhancer","Gene-body","Promoter"),values=cbPalette) +	
		scale_fill_Publication(name="CpG\nContext", labels=c("CPGi","CPGi-shore","Enhancer","Gene-body","Promoter")) +	
		ggtitle("Methylation Difference of Chromosome X") +
		labs(y="% Methylation Difference\n(Female - Male)", x="Tissues") +
		theme_Publication()

	tiff(filename="~/results/chrX/meth.diff.by.tissue.tiff", width=11.7, height=8.27, units="in", res=300, compression='lzw')
	print(p)
	dev.off()
}

##############################
## Plot Meth-Exp of PT chrX ##
##############################
if(FALSE){
	# all chrX genes
	p1<-ggplot(dt.meth.exp[region %in% cpg.contexts & num.sites>=min.num.sites], aes(diff.met, log2FoldChange)) + 
		geom_point(aes(col=padj),size=3) + 
		geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
		ggtitle("Methylation and Expression Difference of ChrX Genes by Gender") +
		labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
		theme_Publication() +
		facet_wrap(~region)
	# all chrX expressed genes (padj not null)
	p2<-ggplot(dt.meth.exp[region %in% cpg.contexts & num.sites>=min.num.sites & !is.na(padj)], aes(diff.met, log2FoldChange)) + 
		geom_point(aes(col=padj),size=3) + 
		geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
		ggtitle("Methylation and Expression Difference of ChrX (expressed) Genes by Gender") +
		labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
		theme_Publication() +
		facet_wrap(~region)
	multiplot(p1,p2)

	##################
	# chrX genes of new.padj less than 0.05
	p1<-ggplot(dt.meth.exp[region %in% cpg.contexts & num.sites>=min.num.sites & new.padj<.1], aes(diff.met, log2FoldChange)) + 
		geom_point(aes(col=new.padj),size=4) + 
		geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
		ggtitle("Methylation and Expression Difference of ChrX Genes by Gender (new.p-adj<0.05)") +
		labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
		theme_Publication() +
		facet_wrap(~region)
	# chrX genes of new.padj less than 0.01
	p2<-ggplot(dt.meth.exp[region %in% cpg.contexts & num.sites>=min.num.sites & new.padj<.01], aes(diff.met, log2FoldChange)) + 
		geom_point(aes(col=new.padj),size=4) + 
		geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
		ggtitle("Methylation and Expression Difference of ChrX Genes by Gender (new.p-adj<0.01)") +
		labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
		theme_Publication() +
		#facet_wrap(~region, nrow=2)
		facet_wrap(~region)
	multiplot(p1,p2)

	##################
	# chrX genes of new.padj less than 0.5
	p1<-ggplot(dt.meth.exp[region %in% cpg.contexts & num.sites>=min.num.sites & new.padj<.5], aes(diff.met, log2FoldChange)) + 
		geom_point(aes(col=new.padj),size=4) + 
		geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
		ggtitle("Methylation and Expression Difference of ChrX Genes by Gender (new.p-adj<0.5)") +
		labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
		theme_Publication() +
		facet_wrap(~region)
	# chrX genes of new.padj less than 0.1
	p2<-ggplot(dt.meth.exp[region %in% cpg.contexts & num.sites>=min.num.sites & new.padj<.1], aes(diff.met, log2FoldChange)) + 
		geom_point(aes(col=new.padj),size=4) + 
		geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
		ggtitle("Methylation and Expression Difference of ChrX Genes by Gender (new.p-adj<0.1)") +
		labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
		theme_Publication() +
		facet_wrap(~region)
	multiplot(p1,p2)

	dev.off()
}

############################
## Meth Boxplot by Gender ##
############################
if(FALSE){
	dt.meth.exp=merge(dt.meth.chrX, dt.exp.pops.chrX, by="ensembl_gene_id")
	# log2FoldChange>0: F-M>0: F>M
	dt.meth.gender=rbind(
						dt.meth.exp[log2FoldChange>0,.(`Gender`="Female",`met`=f.met, num.sites, baseMean, log2FoldChange, pvalue, padj, new.padj),by="ensembl_gene_id,region"],
						dt.meth.exp[log2FoldChange>0,.(`Gender`="Male",`met`=m.met, num.sites, baseMean, log2FoldChange, pvalue, padj, new.padj),by="ensembl_gene_id,region"]
					)
	dt.meth.gender[, `IS_DEG?` := ifelse(is.na(padj), 'NO_EXP', ifelse(new.padj<0.01 ,"DEG","NO_DEG"))]

	p1<-ggplot(dt.meth.gender[region %in% cpg.contexts & num.sites>=min.num.sites], aes(`IS_DEG?`,met)) + 
		geom_boxplot(aes(fill=Gender),show.legend = NA) + 
		scale_fill_manual(values=my.col[["Gender"]]) +	
		#ggtitle("Methylation Level of Differentially Expressed ChrX Genes by Gender (FG, p-adj<0.01 & f>m)") +
		labs(y="% methylation",x="") +
		#theme_bw() +
		theme_Publication() +
		facet_wrap(~region)

	multiplot(p2,p1)
}

	#> dt.meth.gender[`IS_DEG?`=="DEG",length(unique(ensembl_gene_id))]
	#[1] 7
	#> dt.meth.gender[`IS_DEG?`=="Yes",length(unique(ensembl_gene_id))]
	#[1] 23

if(FALSE){
	####################################
	## Venn Diagram between DEG & XIE ##
	####################################
	dt.xci=rbindlist(dt.xci.ensg)
	library(VennDiagram)
	dt.deg.list=list()
	dt.deg.list[["XIE"]]=dt.xci[score==9,ensembl_gene_id]
	dt.deg.list[["chrX.DEG.01"]]=dt.deg.chrX[new.padj<0.01,ensembl_gene_id]
	my.venn.file="~/results/chrX/venn.pops.deg.xci.new.padj.tiff"
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
}

if(FALSE){
	###################
	## get gene name ##
	###################
	# one-to-many: (1-ensg-to-2-hgnc)
	fields <- c("ensembl_gene_id", "hgnc_symbol","gene_biotype", "description")
	dt.gene<-as.data.table(getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.exp.pops[padj<0.1,unique(ensembl_gene_id)], mart = myMart)) # is a data.table
	dt.deg.1<-merge(dt.gene, dt.exp.pops, by="ensembl_gene_id")

	dt.deg.1[padj<0.01 & chromosome_name=="X",.N,(`log2FoldChange>0`=log2FoldChange>0)]
	dt.deg.1[padj<0.05 & chromosome_name=="X",.N,(`log2FoldChange>0`=log2FoldChange>0)]
	dt.deg.1[chromosome_name=="X",.N,(`log2FoldChange>0`=log2FoldChange>0)]
	dt.deg.1[padj<0.01 & chromosome_name=="X",.(ensembl_gene_id,hgnc_symbol,gene_biotype,description,baseMean,log2FoldChange,padj)][order(log2FoldChange)]

	write.csv(dt.deg.1[padj<0.01 & chromosome_name=="X"], file="~/results/chrX/dt.deg.chrX.01.csv")

	##############
	# chisq.test #
	##############
	#sum.ens=dt.ens.chr[chromosome_name %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,'X'),.N,chromosome_name][order(chromosome_name)][,sum(N)]
	sum.ens=dt.ens.chr[chromosome_name %in% dt.deg.1[padj<0.01,unique(chromosome_name)],.N,chromosome_name][order(chromosome_name)][,sum(N)]
	chi1<-chisq.test(
		# observed 
		dt.deg.1[padj<0.01,.N/nrow(dt.exp.pops[padj<0.01]),chromosome_name][order(chromosome_name)][,V1],
		# expected
		p=dt.ens.chr[chromosome_name %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,'X'),.N/sum.ens,chromosome_name][chromosome_name %in% dt.deg.1[padj<0.01,.N,chromosome_name][,chromosome_name]][order(chromosome_name)][,V1]
	)

	chi2<-chisq.test(
		# observed 
		dt.deg.1[padj<0.01,.N,chromosome_name][order(chromosome_name)][,N],
		# expected
		p=dt.ens.chr[chromosome_name %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,'X'),.N,chromosome_name][chromosome_name %in% dt.deg.1[padj<0.01,.N,chromosome_name][,chromosome_name]][order(chromosome_name)][,N],
		rescale.p = TRUE
	)

	chi3<-prop.test(
		# observed 
		dt.deg.1[padj<0.01,.N,chromosome_name][order(chromosome_name)][,N],
		# expected
		dt.ens.chr[chromosome_name %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,'X'),.N,chromosome_name][chromosome_name %in% dt.deg.1[padj<0.01,.N,chromosome_name][,chromosome_name]][order(chromosome_name)][,N]
	)
}# if FALSE

#SAT1 (spermine acetyltransferase 1): ENSG00000130066
#SMS (spermine synthase): ENSG00000102172
#XIST:	ENSG00000229807

