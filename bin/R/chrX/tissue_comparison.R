#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 20/Mar/2017
# Last modified 20/Mar/2017

library(pheatmap)
library(data.table)

TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/GTEx/local.R")# defines my.tissue dt.tissue
										# dl.exp.gtex, dt.chrX.gene, dt.exp.all.chrX, dt.deg.all.chrX, dt.balaton.all.deg
########################
## Genes of Interests ##
########################
target.ensg<-list(
	# 1. Escapee defined by Balaton (n=26+29) + PAR (n=22) + XIST
	`Balaton`=dt.chrX.gene[hgnc_symbol %in% c("XIST",dt.balaton[Balaton %in% c("Mostly E","E","PAR"),hgnc_symbol])]$ensembl_gene_id,
	# 2. DEG from GTEx tissues + Placenta (n=157) at least 30% fold change log2(1.3)
	`DEG10`=dt.deg.all.chrX[abs(log2FoldChange)>log2(1.1),unique(ensembl_gene_id)],
	`DEG20`=dt.deg.all.chrX[abs(log2FoldChange)>log2(1.2),unique(ensembl_gene_id)],
	`DEG30`=dt.deg.all.chrX[abs(log2FoldChange)>log2(1.3),unique(ensembl_gene_id)],
	# 3. DEG from Placenta
	# remove ENSG00000252113 (RNU6-523P) where GTEx data not available
	`PT.DEG`=dt.deg.all.chrX[Tissue=="Placenta", ensembl_gene_id],
	# 4. Female-biased DEG from Placenta (F>M)
	# remove ENSG00000252113 (RNU6-523P) where GTEx data not available
	`PT.DEG.Female`=dt.deg.all.chrX[Tissue=="Placenta" & log2FoldChange>0, ensembl_gene_id],
	# 5. Male-biased DEG from Placenta (F<M)
	`PT.DEG.Male`=dt.deg.all.chrX[Tissue=="Placenta" & log2FoldChange<0, ensembl_gene_id],
	# 6. Placenta-specific DEG (DEG only found in Placenta)
	# remove ENSG00000252113 (RNU6-523P) where GTEx data not available
	`PT.DEG.UNIQ`=dt.deg.all.chrX[Tissue=="Placenta" & !ensembl_gene_id %in% dt.deg.all.chrX[Tissue!="Placenta",unique(ensembl_gene_id)], ensembl_gene_id],
	# 7. Placenta DEG also found in other tissues
	`PT.DEG.COMMON`=dt.deg.all.chrX[Tissue=="Placenta" & ensembl_gene_id %in% dt.deg.all.chrX[Tissue!="Placenta",unique(ensembl_gene_id)], ensembl_gene_id],
	# 8. Placenta-specific escapee (DEG, F>M, not reported as escapee)
	# remove ENSG00000252113 (RNU6-523P) where GTEx data not available
	`PT.Novel.Escaped`=dt.balaton.all.deg[Tissue=="Placenta" & `Expression?`=="female>male" & XCI!="Escaped", ensembl_gene_id],
	# Tissue-wide escapee
	`TWE`=dt.chrX.gene[ensembl_gene_id %in% dt.deg.all.chrX[,.N,ensembl_gene_id][N>=12,ensembl_gene_id],ensembl_gene_id],
	# Tissue-specific escapee
	`TSE`=dt.chrX.gene[ensembl_gene_id %in% dt.deg.all.chrX[,.N,ensembl_gene_id][N<=4,ensembl_gene_id],ensembl_gene_id]
	)
#dcast(dt.exp.all.chrX[ensembl_gene_id %in% target.ensg[["PT.Novel.Escaped"]]], ensembl_gene_id~Tissue, value.var="log2FoldChange")

#############
## Heatmap ##
#############
plotHeatMaplog2FC(my.target="Balaton", is.save=T, my.fontsize_row=7)
plotHeatMaplog2FC(my.target="DEG10", is.save=T, my.fontsize_row=7)
plotHeatMaplog2FC(my.target="PT.DEG", is.order.position=F, is.sort.row=T,is.save=T)
plotHeatMaplog2FC(my.target="PT.DEG.UNIQ", is.order.position=F, is.sort.row=T, is.save=T)
plotHeatMaplog2FC(my.target="PT.DEG.COMMON", is.save=T)

plotHeatMaplog2FC(my.target="PT.DEG.Male", is.save=T)
plotHeatMaplog2FC(my.target="PT.DEG.Female", is.save=T)
plotHeatMaplog2FC(my.target="TWE", is.save=T)
plotHeatMaplog2FC(my.target="TSE", is.save=T, my.fontsize_row=7)
#plotHeatMaplog2FC(my.target="PT.Novel.Escaped", is.order.position=F)

#############
## Boxplot ##
#############
#dt.target<-dt.exp.all.chrX[ensembl_gene_id %in% target.ensg[["PT.DEG.Female"]]]
dt.target<-dt.exp.all.chrX[ensembl_gene_id %in% dt.balaton.all.deg[Tissue=="Placenta" & ensembl_gene_id %in% target.ensg[["PT.DEG.UNIQ"]] & `Expression?`=="female>male",ensembl_gene_id]]
dt.target[,`Source`:=ifelse(Tissue=="Placenta","This Study","GTEx")]
p1<-ggplot(dt.target, aes(Tissue, log2FoldChange)) + 
	geom_boxplot(aes(fill=Source),alpha=.7) +
	scale_x_discrete(limits=c(dt.tissue$SMTS,"Placenta")) +
	scale_fill_manual(values=c(`GTEx`=cbPalette[1],`This Study`=cbPalette[2])) +
	#scale_fill_grey() +
	labs(x="") +
	theme_Publication() +
	theme(axis.text.x=element_text(angle=45, hjust=1),plot.title = element_text(hjust = 0), legend.position="none")

#dt.target<-dt.exp.all.chrX[ensembl_gene_id %in% target.ensg[["PT.DEG.Male"]]]
dt.target<-dt.exp.all.chrX[ensembl_gene_id %in% dt.balaton.all.deg[Tissue=="Placenta" & ensembl_gene_id %in% target.ensg[["PT.DEG.UNIQ"]] & `Expression?`=="female<male",ensembl_gene_id]]
dt.target[,`Source`:=ifelse(Tissue=="Placenta","This Study","GTEx")]
p2<-ggplot(dt.target, aes(Tissue, log2FoldChange)) + 
	geom_boxplot(aes(fill=Source),alpha=.7) +
	scale_x_discrete(limits=c(dt.tissue$SMTS,"Placenta")) +
	scale_fill_manual(values=c(`GTEx`=cbPalette[1],`This Study`=cbPalette[2])) +
	theme_Publication() +
	theme(axis.text.x=element_text(angle=45, hjust=1),plot.title = element_text(hjust = 0), legend.position=c(.92,1))

#file.name<-file.path("~/results/chrX/Figures/Boxplot/Boxplot.PT.chrX.DEG.GTEx.tiff")
file.name<-file.path("~/results/chrX/Figures/Boxplot/Boxplot.PT.chrX.DEG.UNIQ.GTEx.tiff")
tiff(filename=file.name,width=8, height=10,units="in",res=300, compression = 'lzw') #A4 size (2400px x 1500px)
multiplot(p1, p2)
dev.off()

