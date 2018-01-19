#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# based on RoadMap/plot.met.exp.by.gender.R
# First created 30/Mar/2016
# Last modified 13/Jan/2017

TR_PREFIX="GRCh37"
RNASEQ_TYPE="TOTAL" # TOTAL|SMALL (miRNA)
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/chrX/local.R")# defines time.stamp
										# dt.xci.anno, dt.meth.chrX, dt.exp.pops
my.tissue="PT" # PT, PT.SS.Bio, PT.SS.Tech
file.name<-file.path("~/results/chrX",paste("boxplot.meth.diff.by.exp.chrX",RNASEQ_TYPE,my.tissue,time.stamp,"pdf",sep="."))
pdf(file=file.name, width=12, height=8, title=paste("% Meth Diff: ", my.tissue, RNASEQ_TYPE, "RNA-Seq"))

#####################################################
## Boxplot Meth-Diff at CPGi (or promoters) & Gene ##
## for WGoxBS, SS.Bio, SS.Tech                     ##
#####################################################
p1<-ggplot(dt.meth.exp[region %in% cpg.contexts & num.sites>=min.num.sites], aes(tissue, diff.met)) + 
	geom_boxplot(aes(fill=region),outlier.shape=NA,alpha=0.9) +
	geom_hline(yintercept=0) +
	scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Genes"]],my.region.rev[["Promo_10.10"]])) +	
	scale_x_discrete(labels=c("WGoxBS","Bio Rep","Tech Rep")) +
	ggtitle("Chromosome X of human placenta") +
	labs(x="Sequencing Type", y="% Methylation Difference\n(Female - Male)") +
	coord_cartesian(ylim=c(-50,50)) +
	theme_Publication() 
print(p1)

#############################################################
## Boxplot of Meth Diff for WGoxBS, SS.Bio, SS.Tech by XCI ##
#############################################################
dt.foo<-dt.meth.exp[region %in% cpg.contexts & num.sites>=min.num.sites] # at least 4 CG within the region 
dt.boxplot.meth.diff=rbind(
		dt.foo[][,`type`:=factor("All")], # any genes from chrX
		dt.foo[baseMean<10][,`type`:=factor("Rarely expressed")], # any genes from chrX
		#dt.foo[baseMean>=10][,`type`:="Coverage>10"], # any genes from chrX
		dt.foo[baseMean>=10 & new.padj>0.3  & log2FoldChange <= log2(1.1) & log2FoldChange>0][,`type`:=factor("Inactivated")], # Subject to XCI
		dt.foo[baseMean>=10 & new.padj<0.01 & log2FoldChange > log2(1.1)][,`type`:=factor("Escaped")], # Escaped from XCI
		dt.foo[baseMean>=10 & new.padj<0.01 & log2FoldChange < -log2(1.1)][,`type`:=factor("F<<M")] # Genes significantly expressed more in male
		#dt.foo[baseMean>=10 & new.padj>=0.01 & log2FoldChange >= -log2(1.1) & log2FoldChange < 0][,`type`:="F<M"] # Genes expressed more in male
		)
p3<-ggplot(dt.boxplot.meth.diff, aes(type, diff.met)) + 
	geom_boxplot(aes(fill=region),outlier.shape=NA,alpha=0.9) +
	geom_hline(yintercept=0) +
	scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Genes"]],my.region.rev[["Promo_10.10"]])) +	
	scale_x_discrete(limit=c("All","Rarely expressed","Inactivated","Escaped","F<<M")) +
	coord_cartesian(ylim=c(-50,50)) +
	ggtitle("Chromosome X of human placenta") +
	labs(x="Chromosome X Genes", y="% Methylation Difference\n(Female - Male)") +
	facet_wrap(~tissue,nrow=2) + 
	theme_Publication() 
print(p3)

############################################
## Boxplot of Meth Diff for WGoxBS by XCI ##
############################################
p2<-ggplot(dt.boxplot.meth.diff[tissue==my.tissue], aes(type, diff.met)) + 
	geom_boxplot(aes(fill=region),outlier.shape=NA,alpha=0.9) +
	geom_hline(yintercept=0) +
	scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Genes"]],my.region.rev[["Promo_10.10"]])) +	
	scale_x_discrete(limit=c("All","Rarely expressed","Inactivated","Escaped","F<<M")) +
	scale_y_continuous(breaks=seq(from=-50,to=50,by=10)) + 
	ggtitle("Chromosome X of human placenta (WGoxBS)") +
	labs(x="Chromosome X Genes", y="% Methylation Difference\n(Female - Male)") +
	theme_Publication() 
print(p2)

###########################################
## Boxplot of Abs Meth for WGoxBS by XCI ##
###########################################
dt.boxplot.meth.abs=rbind(
	dt.boxplot.meth.diff[,.(tissue,ensembl_gene_id,region,type,`Gender`="Female",`met`=f.met, num.sites, baseMean, log2FoldChange, pvalue, padj, new.padj)],
	dt.boxplot.meth.diff[,.(tissue,ensembl_gene_id,region,type,`Gender`="Male",`met`=m.met, num.sites, baseMean, log2FoldChange, pvalue, padj, new.padj)]
)

# % Methylation
p.abs.meth<-ggplot(dt.boxplot.meth.abs[tissue==my.tissue], aes(type, met)) + 
	geom_boxplot(aes(fill=region,col=Gender,shape=Gender),outlier.shape=NA,alpha=0.7,size=.85) + 
	scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Genes"]],my.region.rev[["Promo_10.10"]])) +	
	scale_colour_manual(values=my.col[["Gender"]]) +	
	scale_x_discrete(limit=c("All","Rarely expressed","Inactivated","Escaped","F<<M")) +
	labs(x="Chromosome X Genes", y="% Methylation") +
	theme_Publication()
print(p.abs.meth)

# % Methylation by grid 1
p.abs.meth2<-ggplot(dt.boxplot.meth.abs[tissue==my.tissue], aes(type, met)) + 
	geom_boxplot(aes(fill=region,col=Gender),alpha=0.7,size=.85) + 
	scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Genes"]],my.region.rev[["Promo_10.10"]])) +	
	scale_colour_manual(values=my.col[["Gender"]]) +	
	scale_x_discrete(limit=c("All","Rarely expressed","Inactivated","Escaped","F<<M")) +
	labs(x="", y="% Methylation") +
	facet_grid(type~Gender) +
	theme_Publication() +
	theme(axis.text.x = element_blank())
print(p.abs.meth2)

# % Methylation by grid 2
if(FALSE){
	file.name<-file.path("~/results/chrX",paste("boxplot.abs.meth.diff.by.XCI.Sex",RNASEQ_TYPE,my.tissue,time.stamp,"tiff",sep="."))
	tiff(filename=file.name,width=8.5, height=12,units="in",res=300, compression = 'lzw') #A4 size
	p.abs.meth3<-ggplot(dt.boxplot.meth.abs[tissue==my.tissue], aes(region,met)) + 
		geom_boxplot(aes(fill=region,col=Gender),alpha=0.7,size=.85) + 
		scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Genes"]],my.region.rev[["Promo_10.10"]])) +	
		scale_colour_manual(values=my.col[["Gender"]]) +	
		labs(x="", y="% Methylation") +
		facet_grid(type~Gender) +
		theme_Publication() +
		theme(axis.text.x = element_blank())
	print(p.abs.meth3)
	dev.off()
}

####################
## % Meth XY plot ##
####################
dt.meth.xy<-dt.meth.exp[region %in% cpg.contexts & num.sites>=min.num.sites] # at least 4 CG within the region 
dt.meth.xy[baseMean<10,`type`:=factor("Rarely expressed")] # any genes from chrX
dt.meth.xy[baseMean>=10 & new.padj>0.3  & log2FoldChange <= log2(1.1) & log2FoldChange>0,`type`:=factor("Inactivated")] # Subject to XCI
dt.meth.xy[baseMean>=10 & new.padj<0.01 & log2FoldChange > log2(1.1),`type`:=factor("Escaped")] # Escaped from XCI
dt.meth.xy[baseMean>=10 & new.padj<0.01 & log2FoldChange < -log2(1.1),`type`:=factor("F<<M")] # Genes significantly expressed more in male
dt.meth.xy[,`region`:=ifelse(region=="Genes","Gene-body",ifelse(region=="Promo_10.10","Promoter","CPGi"))] # update values

p.xy<-ggplot(dt.meth.xy[!is.na(type) & tissue==my.tissue], aes(f.met, m.met)) + 
	geom_point(aes(col=region),size=5,alpha=.6) +
	scale_colour_manual(values=cbPalette) +	
	labs(x="% Methylation (Female)", y="% Methylation (Male)") +
	geom_abline(slope=1, lty=3) +
	facet_grid(region~type) + 
	theme_Publication() +
	theme(legend.position="none")
print(p.xy)

dev.off()


if(FALSE){
	#########################
	# Female Expressed More #
	#########################
	dt.boxplot.meth.diff<-list()
	dt.meth.gender<-list()
	# 1. all chrX genes
	dt.more.female<-dt.meth.exp[region %in% cpg.contexts & num.sites>=min.num.sites & log2FoldChange>0]
	my.ensg.num.more.female<-c(dt.more.female[,length(unique(ensembl_gene_id))],
					dt.more.female[is.na(padj),length(unique(ensembl_gene_id))],
					dt.more.female[new.padj<my.pval[1],length(unique(ensembl_gene_id))],
					dt.more.female[new.padj<my.pval[2],length(unique(ensembl_gene_id))],
					dt.more.female[new.padj<my.pval[3],length(unique(ensembl_gene_id))]
					)

	dt.boxplot.meth.diff[["expressed.more.in.female"]]=rbind(
			dt.more.female[][,`type`:=my.exp.sig[1]],
			# 2. non or rarely expressed chrX gene 
			dt.more.female[is.na(padj)][,`type`:=my.exp.sig[2]],
			dt.more.female[new.padj<my.pval[1]][,`type`:=my.exp.sig[3]],
			dt.more.female[new.padj<my.pval[2]][,`type`:=my.exp.sig[4]],
			dt.more.female[new.padj<my.pval[3]][,`type`:=my.exp.sig[5]]
			)

	dt.meth.gender[["expressed.more.in.female"]]=rbind(
						dt.boxplot.meth.diff[["expressed.more.in.female"]][,.(`Gender`="Female",`met`=f.met, num.sites, baseMean, log2FoldChange, pvalue, padj, new.padj),by="tissue,ensembl_gene_id,region,type"],
						dt.boxplot.meth.diff[["expressed.more.in.female"]][,.(`Gender`="Male",`met`=m.met, num.sites, baseMean, log2FoldChange, pvalue, padj, new.padj),by="tissue,ensembl_gene_id,region,type"]
					)
	# Meth Diff (Female - Male)
	p.meth.diff.f<-ggplot(dt.boxplot.meth.diff[["expressed.more.in.female"]][tissue==my.tissue], aes(type, diff.met)) + 
		geom_boxplot(aes(col=region),alpha=0.7) +
		geom_hline(yintercept=0) +
		ggtitle("Methylation Difference of ChrX Genes (females are more expressed)") +
		labs(x="", y="% Methylation Difference\n(Female - Male)") +
		scale_x_discrete(limits=my.exp.sig, labels=paste0(my.exp.sig, "(n=", my.ensg.num.more.female, ")")) +
		theme_Publication() 
	# % Methylation
	p.abs.meth.f<-ggplot(dt.meth.gender[["expressed.more.in.female"]][tissue==my.tissue], aes(type, met)) + 
		geom_boxplot(aes(col=region,fill=Gender),alpha=0.7) + 
		scale_fill_manual(values=my.col[["Gender"]]) +	
		scale_x_discrete(limits=my.exp.sig, labels=paste0(my.exp.sig, "(n=", my.ensg.num.more.female, ")")) +
		labs(x="Significance Level of Expression Difference (Female - Male)", y="% Methylation") +
		theme_Publication()

	#multiplot(p.meth.diff.f,p.abs.meth.f)

	#########################
	# Male Expressed More #
	#########################
	dt.more.male<-dt.meth.exp[region %in% cpg.contexts & num.sites>=min.num.sites & log2FoldChange<0]
	my.ensg.num.more.male<-c(dt.more.male[,length(unique(ensembl_gene_id))],
					dt.more.male[is.na(padj),length(unique(ensembl_gene_id))],
					dt.more.male[new.padj<my.pval[1],length(unique(ensembl_gene_id))],
					dt.more.male[new.padj<my.pval[2],length(unique(ensembl_gene_id))],
					dt.more.male[new.padj<my.pval[3],length(unique(ensembl_gene_id))]
					)

	dt.boxplot.meth.diff[["expressed.more.in.male"]]=rbind(
		dt.more.male[][,`type`:=my.exp.sig[1]],
		# 2. non or rarely expressed chrX gene 
		dt.more.male[is.na(padj)][,`type`:=my.exp.sig[2]],
		dt.more.male[new.padj<my.pval[1]][,`type`:=my.exp.sig[3]],
		dt.more.male[new.padj<my.pval[2]][,`type`:=my.exp.sig[4]],
		dt.more.male[new.padj<my.pval[3]][,`type`:=my.exp.sig[5]]
		)

	dt.meth.gender[["expressed.more.in.male"]]=rbind(
						dt.boxplot.meth.diff[["expressed.more.in.male"]][,.(`Gender`="Female",`met`=f.met, num.sites, baseMean, log2FoldChange, pvalue, padj, new.padj),by="tissue,ensembl_gene_id,region,type"],
						dt.boxplot.meth.diff[["expressed.more.in.male"]][,.(`Gender`="Male",`met`=m.met, num.sites, baseMean, log2FoldChange, pvalue, padj, new.padj),by="tissue,ensembl_gene_id,region,type"]
					)

	# Meth Diff (Female - Male)
	p.meth.diff.m<-ggplot(dt.boxplot.meth.diff[["expressed.more.in.male"]][tissue==my.tissue], aes(type, diff.met)) + 
		geom_boxplot(aes(col=region),alpha=0.7) +
		geom_hline(yintercept=0) +
		ggtitle("Methylation Difference of ChrX Genes (males are more expressed)") +
		labs(x="", y="% Methylation Difference\n(Female - Male)") +
		scale_x_discrete(limits=my.exp.sig, labels=paste0(my.exp.sig, "(n=", my.ensg.num.more.male, ")")) +
		theme_Publication() 
	# % Methylation
	p.abs.meth.m<-ggplot(dt.meth.gender[["expressed.more.in.male"]][tissue==my.tissue], aes(type, met)) + 
		#geom_boxplot(aes(col=Gender,fill=region),alpha=0.7) + 
		geom_boxplot(aes(col=region,fill=Gender),alpha=0.7) + 
		scale_fill_manual(values=my.col[["Gender"]]) +	
		scale_x_discrete(limits=my.exp.sig, labels=paste0(my.exp.sig, "(n=", my.ensg.num.more.male, ")")) +
		labs(x="Significance Level of Expression Difference (Female - Male)", y="% Methylation") +
		theme_Publication()

	multiplot(p.meth.diff.f,p.meth.diff.m)
	multiplot(p.abs.meth.f,p.abs.meth.m)
	#multiplot(p.meth.diff.m,p.abs.meth.m)
	dev.off()
}


cat("All done\n")
