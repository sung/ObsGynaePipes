#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# based on RoadMap/plot.met.exp.by.gender.R
# First created 30/Mar/2017
# Last modified 30/Mar/2017

TR_PREFIX="GRCh37"
RNASEQ_TYPE="TOTAL" # TOTAL|SMALL (miRNA)
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/chrX/local.R")# defines time.stamp dt.xci.anno, dt.meth.chrX, dt.exp.pops
source("~/Pipelines/bin/R/GTEx/local.R")# defines dt.chrX.gene 

my.tissue="PT.SS.Tech" # PT, PT.SS.Bio, PT.SS.Tech
METH_TYPE=meth.type.rev[[my.tissue]] # WGoxBS, SS.Bio, SS.Tech

if(FALSE){
	##############################
	## Prepare a data structure ##
	##############################
	dt.foo<-dt.meth.exp[region %in% cpg.contexts & num.sites>=min.num.sites] # at least 4 CG within the region 
	dt.bar=rbind(
			#dt.foo[baseMean<10][,`type`:=factor("Rarely expressed")], # any genes from chrX
			#dt.foo[TRUE][,`type`:=factor("All chrX")], # any genes from chrX
			dt.foo[baseMean>=10 & new.padj<0.01 & log2FoldChange > log2(1.1)][,`type`:=factor("Escaped")], # Escaped from XCI
			dt.foo[baseMean>=10 & new.padj>0.3  & log2FoldChange <= log2(1.1) & log2FoldChange>0][,`type`:=factor("Inactivated")], # Subject to XCI
			dt.foo[baseMean>=10 & new.padj<0.01 & log2FoldChange < -log2(1.1)][,`type`:=factor("Male-biased")] # Genes significantly expressed more in male
			)
	dt.bar$region<-factor(dt.bar$region, levels=c('CPGi','Promo_10.10','Genes'))

	#####################################
	## Boxplot of Meth Diff for by XCI ##
	#####################################
	# with outlier: All chrX, Escaped, Inactivated, Male-biased
	p1<-ggplot(dt.bar[tissue==my.tissue], aes(type, diff.met)) + 
		geom_boxplot(aes(fill=region),alpha=0.8, outlier.shape=1, outlier.size=5, size=1, width=.5) +
		geom_hline(yintercept=0) +
		scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		scale_y_continuous(breaks=seq(from=-50,to=50,by=10)) + 
		ggtitle(paste0("Human Placenta (",METH_TYPE,")")) +
		labs(x="Chromosome X Genes", y="% Methylation Difference (Female - Male)") +
		theme_Publication() 

	file.name<-file.path("~/results/chrX/Figures",paste("meth.diff.chrX",RNASEQ_TYPE,METH_TYPE,"tiff",sep="."))
	tiff(filename=file.name,width=12, height=8.5,units="in",res=300, compression = 'lzw') #A4 size
	print(p1)
	dev.off()

	# without outlier: All chrX, Escaped, Inactivated, Male-biased
	p2<-ggplot(dt.bar[tissue==my.tissue], aes(type, diff.met)) + 
		geom_boxplot(aes(fill=region),alpha=0.8, outlier.shape=NA, size=1, width=.5) +
		geom_hline(yintercept=0) +
		scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		scale_y_continuous(breaks=seq(from=-50,to=50,by=10)) + 
		ggtitle(paste0("Human Placenta (",METH_TYPE,")")) +
		labs(x="Chromosome X Genes", y="% Methylation Difference (Female - Male)") +
		theme_Publication() 

	file.name<-file.path("~/results/chrX/Figures",paste("meth.diff.chrX",RNASEQ_TYPE,METH_TYPE,"no.outlier","tiff",sep="."))
	tiff(filename=file.name,width=12, height=8.5,units="in",res=300, compression = 'lzw') #A4 size
	print(p2)
	dev.off()

	###############################################
	## Scatterplot of meth-diff & log2FoldChange ##
	###############################################
	p.metH.fc.grid<-ggplot(merge(dt.bar[tissue==my.tissue], dt.chrX.gene[,.(ensembl_gene_id,hgnc_symbol)]), aes(diff.met, log2FoldChange, label=hgnc_symbol)) +
		geom_point(aes(col=region),size=6, alpha=.5) + 
		geom_text(check_overlap = TRUE, hjust =1.2, nudge_x=0.05, size=4.5) +
		coord_cartesian(ylim=c(-1, 1)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
		geom_hline(yintercept=0) +
		geom_vline(xintercept=0) +
		scale_color_manual(values=cbPalette2, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		labs(x="% Methylation Difference (Female - Male)", y="Log2 Fold-Change") +
		facet_grid(region~type, labeller=labeller(region=unlist(my.region.rev))) + 
		theme_Publication() +
		theme(legend.position="none")
	print(p.meth.fc.grid)

	file.name<-file.path("~/results/chrX/Figures/MethFC",paste("log2FC.diffMeth.chrX.grid",time.stamp,"tiff",sep="."))
	tiff(filename=file.name,width=12, height=11,units="in",res=300, compression = 'lzw') #A4 size
	print(p.meth.fc.grid)
	dev.off()

	####################
	## % Meth XY plot ##
	####################
	p.xy<-ggplot(dt.bar[tissue==my.tissue], aes(f.met, m.met)) + 
		geom_point(aes(col=region),size=5,alpha=.6) +
		scale_color_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		labs(x="% Methylation (Female)", y="% Methylation (Male)") +
		geom_abline(slope=1, lty=3) +
		facet_grid(region~type, labeller=labeller(region=unlist(my.region.rev))) + 
		theme_Publication() +
		theme(legend.position="none")
	print(p.xy)

	file.name<-file.path("~/results/chrX/Figures",paste("meth.female.male.scatter.chrX",RNASEQ_TYPE,METH_TYPE,"tiff",sep="."))
	tiff(filename=file.name,width=12, height=11,units="in",res=300, compression = 'lzw') #A4 size
	print(p.xy)
	dev.off()

	###################################
	## % Meth XY plot with gene name ##
	###################################
	p.xy.gene<-ggplot(merge(dt.bar[tissue==my.tissue], dt.chrX.gene[,.(ensembl_gene_id,hgnc_symbol)]), aes(f.met, m.met, label=hgnc_symbol)) +
				geom_point(aes(col=region),size=5,alpha=.6) +
				geom_text(check_overlap = TRUE, hjust =1.2, nudge_x=0.05, size=4.5) +
				scale_color_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
				labs(x="% Methylation (Female)", y="% Methylation (Male)") +
				geom_abline(slope=1, lty=3) +
				facet_grid(region~type, labeller=labeller(region=unlist(my.region.rev))) + 
				theme_Publication() +
				theme(legend.position="none")
	file.name<-file.path("~/results/chrX/Figures",paste("meth.female.male.scatter.chrX.gene.name",RNASEQ_TYPE,METH_TYPE,"tiff",sep="."))
	tiff(filename=file.name,width=12, height=11,units="in",res=300, compression = 'lzw') #A4 size
	print(p.xy.gene)
	dev.off()

	my.type="Escaped"
	my.regions=c('CPGi','Promo_10.10','Genes')
	for(i in seq_along(my.regions)){
		p.xy.gene1<-ggplot(merge(dt.bar[tissue==my.tissue & region==my.regions[i] & type==my.type], dt.chrX.gene[,.(ensembl_gene_id,hgnc_symbol)]), aes(f.met, m.met, label=hgnc_symbol)) +
					geom_point(aes(col=region),size=5, alpha=.7,col=cbPalette[i]) + 
					geom_text(check_overlap = TRUE, hjust =1.2, nudge_x=0.05, size=4.5) +
					labs(x="% Methylation (Female)", y="% Methylation (Male)") +
					geom_abline(slope=1, lty=3) +
					ggtitle(paste(my.region.rev[[my.regions[i]]],"region of chromosome X")) +
					theme_Publication() +
					theme(legend.position="none")
		file.name<-file.path("~/results/chrX/Figures",paste("meth.female.male.scatter.chrX.gene.name",my.regions[i],my.type,RNASEQ_TYPE,METH_TYPE,"tiff",sep="."))
		tiff(filename=file.name,width=12, height=11,units="in",res=300, compression = 'lzw') #A4 size
		print(p.xy.gene1)
		dev.off()
	}

	#################################
	## Boxplot of Abs Meth by grid ##
	#################################
	dt.boxplot.meth.abs=rbind(
		dt.bar[,.(tissue,ensembl_gene_id,region,type,`Gender`="Female",`met`=f.met, num.sites, baseMean, log2FoldChange, pvalue, padj, new.padj)],
		dt.bar[,.(tissue,ensembl_gene_id,region,type,`Gender`="Male",`met`=m.met, num.sites, baseMean, log2FoldChange, pvalue, padj, new.padj)]
	)

	p.abs.meth2<-ggplot(dt.boxplot.meth.abs[tissue==my.tissue], aes(type, met)) + 
		geom_boxplot(aes(fill=region,col=Gender),alpha=0.7,size=.85) + 
		scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		scale_colour_manual(values=my.col[["Gender"]]) +	
		labs(x="", y="% Methylation") +
		facet_grid(type~Gender, labeller=labeller(region=unlist(my.region.rev))) + 
		theme_Publication() +
		theme(axis.text.x = element_blank())
	print(p.abs.meth2)

	p.abs.meth3<-ggplot(dt.boxplot.meth.abs[tissue==my.tissue], aes(region,met)) + 
		geom_boxplot(aes(fill=region,col=Gender),alpha=0.7,size=.85, outlier.shape=1, outlier.size=5, width=.5) +
		scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		scale_x_discrete(labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		scale_colour_manual(values=my.col[["Gender"]]) +	
		labs(x="", y="% Methylation") +
		facet_grid(type~Gender, labeller=labeller(region=unlist(my.region.rev))) + 
		theme_Publication() +
		theme(legend.position="none")
	print(p.abs.meth3)

	file.name<-file.path("~/results/chrX/Figures",paste("boxplot.abs.meth.by.sex.chrX.",RNASEQ_TYPE,METH_TYPE,"tiff",sep="."))
	tiff(filename=file.name,width=8.5, height=11,units="in",res=300, compression = 'lzw') #A4 size
	print(p.abs.meth3)
	dev.off()

	###################
	## All chrX only ##
	###################
	file.name<-file.path("~/results/chrX/Figures",paste("meth.exp.diff.chrX",RNASEQ_TYPE,METH_TYPE,time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=12, height=8, title=paste("% Meth Diff: ", METH_TYPE, RNASEQ_TYPE, "RNA-Seq"))

	my.regions=c('CPGi','Promo_10.10','Genes')
	for(i in seq_along(my.regions)){
		p4<-ggplot(dt.bar[tissue==my.tissue & type=="All chrX" & region==my.regions[i]], aes(diff.met, log2FoldChange)) + 
			geom_point(aes(col=region),size=5, alpha=.7,col=cbPalette[i]) + 
			geom_hline(yintercept=0) +
			geom_vline(xintercept=0) +
			labs(x="% Methylation Difference (Female - Male)", y="Log2 Fold-Change") +
			ggtitle(paste(my.region.rev[[my.regions[i]]],"region of chromosome X")) +
			theme_Publication()
		print(p4)
	}

	p5<-ggplot(dt.bar[tissue==my.tissue & type=="All chrX"], aes(diff.met, log2FoldChange, group=type)) + 
		geom_point(aes(col=region),size=5, alpha=.7) + 
		geom_hline(yintercept=0) +
		geom_vline(xintercept=0) +
		scale_color_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		labs(x="% Methylation Difference (Female - Male)", y="Log2 Fold-Change") +
		facet_wrap(~region, labeller=labeller(region=unlist(my.region.rev))) + 
		theme_Publication()
	print(p5)

	dev.off()


	##################
	## Fisher Test ##
	##################
	# 1. total-RNA
	fisher.test(matrix(c(2392, 51, 54849, 50),nrow=2))$p.value # based on total number of genes
	fisher.test(matrix(c(639, 51, 18801, 50),nrow=2))$p.value # based on total number of expressed (depth>10) genes

	# 2. small RNA
	fisher.test(matrix(c(133, 4, 1943, 3),nrow=2))$p.value # based on total number of miRNA
	fisher.test(matrix(c(109, 4, 1592, 3),nrow=2))$p.value # based on total number of expressed (depth>10) miRNA

	############################################################
	## % Methylation Difference vs. Boxplot of Log2FoldChange ##
	############################################################
	dt.foo<-dt.meth.exp[tissue==my.tissue & region %in% cpg.contexts & num.sites>=min.num.sites] # at least 4 CG within the region 
	#table(dt.foo[,cut(diff.met,c(-100,-20,-10,-5,0,5,10,20,100),include.lowest =T, right=F)])
	dt.foo[,meth.bins:=cut(diff.met,c(-100,-20,-10,-5,0,5,10,20,100),include.lowest =T, right=F)]
	dt.foo$region<-factor(dt.foo$region, levels=c('CPGi','Promo_10.10','Genes'))

	#http://stackoverflow.com/questions/21310609/ggplot2-box-whisker-plot-show-95-confidence-intervals-remove-outliers
	quantiles_95 <- function(x) {
		r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
		names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
		r
	}
	#https://rpubs.com/dgolicher/median_boot
	median_cl_boot <- function(x, conf = 0.95) {
		lconf <- (1 - conf)/2
		uconf <- 1 - lconf
		require(boot)
		bmedian <- function(x, ind) median(x[ind])
		bt <- boot(x, bmedian, 1000)
		bb <- boot.ci(bt, type = "perc")
		data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t, uconf))
	}

	plotBoxBin<-function(my.padj=1){
		p.bin<-ggplot(dt.foo[baseMean>=10 & new.padj <= my.padj], aes(meth.bins,log2FoldChange)) + 
			geom_boxplot(aes(fill=region), alpha=0.8, outlier.shape=NA, size=.7, width=.8) +
			geom_hline(yintercept=0) +
			scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
			coord_cartesian(ylim=c(-.6, .6)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
			scale_y_continuous(breaks=seq(from=-1,to=1,by=.2)) + 
			ggtitle(paste("Chromosome X (DEG P<",my.padj,")")) +
			labs(x="Bins of % Methylation Difference (Female-Male)", y="Log2FoldChange (Female - Male)") +
			theme_Publication() 
		return(p.bin)
	}
	plotBoxBin1<-function(my.padj=1){
		p.bin<-ggplot(dt.foo[baseMean>=10 & new.padj <= my.padj], aes(meth.bins,log2FoldChange)) + 
			geom_boxplot(aes(col=region), alpha=0.8, outlier.shape=NA, size=.9, width=.8) +
			geom_hline(yintercept=0) +
			scale_colour_manual(values=cbPalette2, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
			coord_cartesian(ylim=c(-.6, .6)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
			scale_y_continuous(breaks=seq(from=-1,to=1,by=.2)) + 
			ggtitle(paste("Chromosome X (DEG P<",my.padj,")")) +
			labs(x="Bins of % Methylation Difference (Female-Male)", y="Log2FoldChange (Female - Male)") +
			theme_Publication() 
		return(p.bin)
	}
	plotBoxBin2<-function(my.padj=1){
		p.bin<-ggplot(dt.foo[baseMean>=10 & new.padj <= my.padj], aes(meth.bins,log2FoldChange)) + 
			stat_summary(aes(fill=region), alpha=.8, position="dodge", fun.data = quantiles_95, geom="boxplot", size=.7) +
			geom_hline(yintercept=0) +
			scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
			coord_cartesian(ylim=c(-.6, .6)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
			scale_y_continuous(breaks=seq(from=-1,to=1,by=.2)) + 
			ggtitle(paste("Chromosome X (DEG P<",my.padj,")")) +
			labs(x="Bins of % Methylation Difference (Female-Male)", y="Log2FoldChange (Female - Male)") +
			theme_Publication() 
		print(p.bin)
	}
	plotBoxBin3<-function(my.padj=1){
		p.bin<-ggplot(dt.foo[baseMean>=10 & new.padj <= my.padj], aes(meth.bins,log2FoldChange,group=region)) + 
			stat_summary(aes(col=region), position="dodge", fun.data = median_cl_boot, geom = "errorbar") + 
			stat_summary(aes(col=region), position="dodge", fun.y = median, geom = "point") +
			#position_dodge(width = .7) +
			geom_hline(yintercept=0) +
			scale_colour_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
			coord_cartesian(ylim=c(-.6, .6)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
			scale_y_continuous(breaks=seq(from=-1,to=1,by=.2)) + 
			ggtitle(paste("Chromosome X (DEG P<",my.padj,")")) +
			labs(x="Bins of % Methylation Difference (Female-Male)", y="Log2FoldChange (Female - Male)") +
			theme_Publication() 
		print(p.bin)
	}

	p1<-plotBoxBin(1)
	p2<-plotBoxBin(.5)
	p3<-plotBoxBin(.05)
	p4<-plotBoxBin(.01)

	p1.1<-plotBoxBin1(1)
	p1.2<-plotBoxBin1(.5)
	p1.3<-plotBoxBin1(.05)
	p1.4<-plotBoxBin1(.01)

	file.name<-file.path("~/results/chrX/Figures/MethFC",paste("log2FC.diffmethBin.chrX",time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=12, height=10, title=paste("% Diff Meth and Expression"))
	#multiplot(p1,p2)
	#multiplot(p3,p4)
	multiplot(p1.1,p1.2)
	multiplot(p1.3,p1.4)
	dev.off()

	###############################################
	## Plot Meth Diff and Log2FoldChange of chrX ##
	###############################################
	dt.foo<-dt.meth.exp[tissue==my.tissue & region %in% cpg.contexts & num.sites>=min.num.sites] # at least 4 CG within the region 
	dt.foo[,meth.bins:=cut(diff.met,c(-100,-20,-10,-5,0,5,10,20,100),include.lowest =T, right=F)]
	dt.foo[,padj.bins:=cut(new.padj,c(0,1e-2,1e-2*5,1e-1*5,1),include.lowest =T, right=F)]
	dt.foo$region<-factor(dt.foo$region, levels=c('CPGi','Promo_10.10','Genes'))

	plotDiff<-function(my.padj=1){
		p.diff<-ggplot(dt.foo[baseMean>=10 & new.padj<my.padj], aes(diff.met, log2FoldChange)) + 
			geom_point(aes(col=new.padj),size=4,alpha=.7) + 
			geom_smooth(se=FALSE,col="red") +
			coord_cartesian(ylim=c(-1, 1)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
			geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
			ggtitle(paste("Methylation and Expression Difference of ChrX Genes by Fetal Sex (new.p-adj<",my.padj,")")) +
			labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
			theme_Publication() +
			facet_grid(~region, labeller=labeller(region=unlist(my.region.rev)))
		return(p.diff)
	}

	p.diff1<-plotDiff(1)
	p.diff2<-plotDiff(.5)
	p.diff3<-plotDiff(.05)
	p.diff4<-plotDiff(.01)

	time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p') #2015-01-07_01_30PM
	file.name<-file.path("~/results/chrX/Figures/MethFC",paste("log2FC.diffmeth.chrX",time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=12, height=10, title=paste("% Diff Meth and Log2FoldChange"))
	multiplot(p.diff1,p.diff2)
	multiplot(p.diff3,p.diff4)
	dev.off()
}

##################
## All together ##
##################
if(FALSE){
	## Prepare a data structure ##
	## abs(log2FoldChange)<1 to remove XIST from plotting 
	dt.foo<-dt.meth.exp[tissue==my.tissue & region %in% cpg.contexts & num.sites>=min.num.sites & baseMean>=10 & abs(log2FoldChange)<1] # at least 4 CG within the region 
	dt.foo[,meth.bins:=cut(diff.met,c(-100,-20,-10,-5,0,5,10,20,100),include.lowest =T, right=F)]
	dt.foo[,padj.bins:=cut(new.padj,c(0,1e-2,1e-2*5,1e-1*5,1),include.lowest =T, right=F)]
	dt.foo$region<-factor(dt.foo$region, levels=c('CPGi','Promo_10.10','Genes'))

	# for meth.fc.grid
	dt.bar=rbind(
			dt.foo[new.padj<0.01 & log2FoldChange > log2(1.1)][,`type`:=factor("Escaped")], # Escaped from XCI
			dt.foo[new.padj>0.3  & log2FoldChange <= log2(1.1) & log2FoldChange>0][,`type`:=factor("Inactivated")], # Subject to XCI
			dt.foo[new.padj<0.01 & log2FoldChange < -log2(1.1)][,`type`:=factor("Male-biased")] # Genes significantly expressed more in male
			)
	dt.bar$region<-factor(dt.bar$region, levels=c('CPGi','Promo_10.10','Genes'))

	## Correlation Coefficient (diff.meth vs. log2FC)
	dt.cf<-dt.foo[baseMean>=10,.("pearson.cf"=cor(diff.met,log2FoldChange),"spearman.cf"=cor(diff.met,log2FoldChange,method="spearman")),.(region,padj.bins)][order(region,padj.bins)]
	#dcast( dt.cf, region~padj.bins, value.var="pearson.cf")
	#dcast( dt.cf, region~padj.bins, value.var="spearman.cf")

	p.diff<-ggplot(dt.foo, aes(diff.met, log2FoldChange)) + 
		geom_point(aes(col=new.padj),size=4,alpha=.7) + 
		geom_smooth(se=FALSE,col="red") +
		coord_cartesian(ylim=c(-1, 1)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
		geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
		ggtitle("Methylation and Expression Difference of ChrX Genes by P-value Ranges & Genomic Regions") +
		labs(x="% Methylation Difference (Female - Male)", y="log2FC") +
		theme_Publication() +
		facet_grid(region~padj.bins, labeller=labeller(region=unlist(my.region.rev)))

	p.bin<-ggplot(dt.foo, aes(meth.bins,log2FoldChange)) + 
		geom_boxplot(aes(col=region), alpha=0.8, outlier.shape=NA, size=.9, width=.8) +
		geom_hline(yintercept=0) +
		scale_colour_manual(values=cbPalette2, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		coord_cartesian(ylim=c(-.6, .6)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
		scale_y_continuous(breaks=seq(from=-1,to=1,by=.2)) + 
		ggtitle("Methylation and Expression Difference of ChrX Genes by P-value Ranges & Genomic Regions") +
		labs(x="Bins of % Methylation Difference (Female-Male)", y="Log2FoldChange (Female - Male)") +
		theme_Publication() +
		facet_wrap(~padj.bins,nrow=2)

	p.bin2<-ggplot(dt.foo, aes(meth.bins,log2FoldChange)) + 
		geom_boxplot(aes(col=region), alpha=0.8, outlier.shape=NA, size=.7, width=.8) +
		geom_hline(yintercept=0) +
		scale_colour_manual(values=cbPalette2, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		coord_cartesian(ylim=c(-.6, .6)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
		scale_y_continuous(breaks=seq(from=-1,to=1,by=.2)) + 
		ggtitle("Methylation and Expression Difference of ChrX Genes by P-value Ranges & Genomic Regions") +
		labs(x="Bins of % Methylation Difference (Female-Male)", y="Log2FoldChange (Female - Male)") +
		theme_Publication() +
		theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
		facet_grid(region~padj.bins, labeller=labeller(region=unlist(my.region.rev)))

	p.cf1<-ggplot(dt.cf, aes(padj.bins, pearson.cf)) + 
		geom_bar(stat="identity",aes(fill=region),position="dodge") +
		scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		labs(x="Bins of P-value range", y="Correlation Coefficient (Pearson)") +
		theme_Publication() 
	#print(p.cf1)

	p.cf2<-ggplot(dt.cf, aes(padj.bins, spearman.cf)) + 
		geom_bar(stat="identity",aes(fill=region),position="dodge") +
		scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		labs(x="Bins of P-value range", y="Correlation Coefficient (Spearman)") +
		theme_Publication() 
	#print(p.cf2)

	p.meth.fc.grid<-ggplot(merge(dt.bar, dt.chrX.gene[,.(ensembl_gene_id,hgnc_symbol)]), aes(diff.met, log2FoldChange, label=hgnc_symbol)) +
		geom_point(aes(col=region),size=6, alpha=.5) + 
		geom_text(check_overlap = TRUE, hjust =1.2, nudge_x=0.05, size=4.5) +
		coord_cartesian(ylim=c(-1, 1)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
		geom_hline(yintercept=0) +
		GEOM_vline(xintercept=0) +
		scale_color_manual(values=cbPalette2, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		labs(x="% Methylation Difference (Female - Male)", y="Log2 Fold-Change") +
		facet_grid(region~type, labeller=labeller(region=unlist(my.region.rev))) + 
		theme_Publication() +
		theme(legend.position="none")

	p.box<-ggplot(dt.bar[tissue==my.tissue], aes(type, diff.met)) + 
		geom_boxplot(aes(fill=region),alpha=0.8, outlier.shape=NA, size=1, width=.5) +
		geom_hline(yintercept=0) +
		scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Promo_10.10"]],my.region.rev[["Genes"]])) +	
		scale_y_continuous(breaks=seq(from=-50,to=50,by=10)) + 
		coord_cartesian(ylim=c(-30, 30)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
		ggtitle(paste0("Human Placenta (",METH_TYPE,")")) +
		labs(x="Chromosome X Genes", y="% Methylation Difference (Female - Male)") +
		theme_Publication() 

	time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p') #2015-01-07_01_30PM
	file.name<-file.path("~/results/chrX/Figures/MethFC",paste(METH_TYPE,"log2FC.diffmeth.chrX.all",time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=12, height=10, title=paste("% Diff Meth and Log2FoldChange"))

	print(p.diff)
	print(p.bin2)
	#multiplot(p.cf1, p.cf2)
	#print(p.xy.gene)
	print(p.meth.fc.grid)
	print(p.box)

	dev.off()
}

##########################################################
## Fisher-test between %Meth-Differecne and Fold-Change ## 
##########################################################
if(FALSE){
	# 5% diff
	dt.foo<-dt.meth.exp[tissue==my.tissue & region %in% cpg.contexts & num.sites>=min.num.sites & baseMean>=10] # at least 4 CG within the region 
	dt.foo[,meth.bins:=cut(diff.met,c(-100,-5,5,100),include.lowest =T, right=F)]
	dt.foo$region<-factor(dt.foo$region, levels=c('CPGi','Promo_10.10','Genes'))

	dt.bar=rbind(
			dt.foo[new.padj<0.01 & log2FoldChange > log2(1.1)][,`type`:=factor("Female-biased")], # Escaped from XCI
			dt.foo[new.padj>0.3  & log2FoldChange <= log2(1.1) & log2FoldChange>0][,`type`:=factor("Inactivated")], # Subject to XCI
			dt.foo[new.padj<0.01 & log2FoldChange < -log2(1.1)][,`type`:=factor("Male-biased")] # Genes significantly expressed more in male
			)
	dt.bar$region<-factor(dt.bar$region, levels=c('CPGi','Promo_10.10','Genes'))

	dcast(dt.bar[,.N,.(type,region,meth.bins)], region+meth.bins~type,value.var="N")
	write.csv(dcast(dt.bar[,.N,.(type,region,meth.bins)], region+meth.bins~type,value.var="N"), file=file.path("~/results/chrX/Figures/MethFC",paste(METH_TYPE,"meth.diff.count.chrX.5.percent.csv",sep=".")),row.names=F)

	# 10% diff
	dt.foo<-dt.meth.exp[tissue==my.tissue & region %in% cpg.contexts & num.sites>=min.num.sites & baseMean>=10] # at least 4 CG within the region 
	dt.foo[,meth.bins:=cut(diff.met,c(-100,-10,10,100),include.lowest =T, right=F)]
	dt.foo$region<-factor(dt.foo$region, levels=c('CPGi','Promo_10.10','Genes'))

	dt.bar=rbind(
			dt.foo[new.padj<0.01 & log2FoldChange > log2(1.1)][,`type`:=factor("Female-biased")], # Escaped from XCI
			dt.foo[new.padj>0.3  & log2FoldChange <= log2(1.1) & log2FoldChange>0][,`type`:=factor("Inactivated")], # Subject to XCI
			dt.foo[new.padj<0.01 & log2FoldChange < -log2(1.1)][,`type`:=factor("Male-biased")] # Genes significantly expressed more in male
			)
	dt.bar$region<-factor(dt.bar$region, levels=c('CPGi','Promo_10.10','Genes'))

	dcast(dt.bar[,.N,.(type,region,meth.bins)], region+meth.bins~type,value.var="N")
	write.csv(dcast(dt.bar[,.N,.(type,region,meth.bins)], region+meth.bins~type,value.var="N"), file=file.path("~/results/chrX/Figures/MethFC",paste(METH_TYPE,"meth.diff.count.chrX.10.percent.csv",sep=".")),row.names=F)

	#fisher.test(t(matrix(c(2,11,0,22,53,0,5,108,5), nrow=3)))
}
