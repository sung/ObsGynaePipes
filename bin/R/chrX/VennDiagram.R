#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 12/Apr/2017
# Last modified 12/Apr/2017

library(data.table)
library(VennDiagram)
TR_PREFIX="GRCh37"

load(paste0("~/data/Annotation/Ensembl/",TR_PREFIX,"/dt.ensg.RData")) # load 'gr.ensg, dt.ensg'
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/GTEx/local.R")# defines my.tissue dt.tissue
										# dl.exp.gtex, dt.chrX.gene, dt.exp.all.chrX, dt.deg.all.chrX, dt.balaton.all.deg
##################
## Balaton 2015 ##
##################
# dt.balaton from GTEx/local.R

###################
## Venn Diagram  ##
###################
venn.list=list()
venn.list[["Tukiainen.Escaped.Female-biased.tuki"]][["Placenta"]]=dt.deg.all.chrX[Tissue=="Placenta" & log2FoldChange>0, ensembl_gene_id] # 47 ENSG
venn.list[["Tukiainen.Escaped.Female-biased.tuki"]][["GTEx"]]=dt.deg.all.chrX[Tissue!="Placenta" & log2FoldChange>0, unique(ensembl_gene_id)] # 94 ENSG (female-biased)
venn.list[["Tukiainen.Escaped.Female-biased.tuki"]][["Known Escaped"]]=dt.tuki[reported_XCI=="Escape",ensembl_gene_id]
venn.list[["Tukiainen.Escaped.Female-biased.tuki"]][["Known Inact"]]=dt.tuki[reported_XCI=="Inactive",ensembl_gene_id]

venn.list[["Escaped.Female-biased"]][["Placenta"]]=dt.deg.all.chrX[Tissue=="Placenta" & log2FoldChange>0, unique(as.character(ensembl_gene_id))] # 47 ENSG
venn.list[["Escaped.Female-biased"]][["GTEx"]]=dt.deg.all.chrX[Tissue!="Placenta" & log2FoldChange>0, unique(as.character(ensembl_gene_id))] # 94 ENSG (female-biased)

venn.list[["PAR"]][["Placenta"]]=dt.deg.all.chrX[Tissue=="Placenta", ensembl_gene_id] # 60 ENSG
venn.list[["PAR"]][["GTEx"]]=dt.deg.all.chrX[Tissue!="Placenta",unique(ensembl_gene_id)] # 176 ENSG
venn.list[["PAR"]][["PAR"]]=dt.ensg[chromosome_name=="X" & (end_position<2699520 | start_position>154931044), as.character(ensembl_gene_id) ] # 47 ENSG

my.type="Tukiainen.Escaped.Female-biased"
file.name<-file.path(paste0("~/results/chrX/Figures/Venn/Venn.PT.GTEx.",my.type,".tiff"))
venn.diagram(
	x=venn.list[[my.type]],
	filename = file.name,
	height=1500,
	width=2400,
	#main="C",
	#main.cex = 1.5,
	main.fontfamily = "sans",
	main.fontface = "bold",
	col = "black",
	fill = c(cbPalette[2],cbPalette[1],"white","white"),
	alpha = 0.50,
	#cat.col = c("dodgerblue", "goldenrod1", "darkorange1"),
	cex = 1.4,
	fontfamily = "sans",
	#cat.cex = 1.3,
	cat.fontfamily = "sans",
	cat.fontface = "bold",
	margin = 0.05
)

my.type="Escaped.Female-biased"
file.name<-file.path(paste0("~/results/chrX/Figures/Venn/Venn.PT.GTEx.",my.type,".tiff"))
venn.diagram(
	x=venn.list[[my.type]],
	filename = file.name,
	height=1500,
	width=2400,
	main.fontfamily = "sans",
	main.fontface = "bold",
	col = "black",
	fill = c(cbPalette[2],cbPalette[1]),
	alpha = 0.50,
	cex = 1.4,
	fontfamily = "sans",
	cat.fontfamily = "sans",
	cat.fontface = "bold",
	margin = 0.05
)

