library(data.table)
library(VennDiagram)
TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R")
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/GTEx/local.R")

#mergeDEG(flag="withChrY")
load("~/results/RNA-Seq/GTEx/DESeq2/res.all.tissues.RData")
dt.exp.all<-rbindlist(dl.exp.gtex)
dt.exp.all[,`baseFemale`:=2*baseMean/(1+1/2^log2FoldChange)]
dt.exp.all[,`baseMale`:=2*baseMean/(1+2^log2FoldChange)]

# ChrY only
dt.exp.all.chrY<-dt.exp.all[chromosome_name=="Y"]
# adjust p-value based on chr only
dt.exp.all.chrY[,`new.padj` := p.adjust(pvalue, method="BH"),by="Tissue"]

dt.exp.all.chrX<-dt.exp.all[chromosome_name=="X"]
dt.exp.all.chrX[,`new.padj` := p.adjust(pvalue, method="BH"),by="Tissue"]

################## 
## Venn diagram ##
################## 
my.baseMean=1
dt.deg.all<-dt.exp.all.chrY[new.padj<0.01 & abs(log2FoldChange)>log2(1.1) & baseMean>=my.baseMean]
venn.list=list()
venn.list[["Placenta"]]<-dt.deg.all[chromosome_name=="Y" & Tissue=="Placenta",as.character(ensembl_gene_id)] # chrY in Placenta 
venn.list[["GTEx"]]<-dt.deg.all[chromosome_name=="Y" & Tissue!="Placenta",as.character(unique(ensembl_gene_id))] # chrY in GTEx

if(FALSE){
file.name<-file.path("~/results/RNA-Seq/GTEx/Figures/Venn",paste0("Venn.GTEx.chrY.",my.baseMean,".tiff"))
VennDiagram::venn.diagram(
	x=venn.list,
	filename = file.name,
	height=900,
	width=1400,
	main.fontfamily = "sans",
	main.fontface = "bold",
	col = "black",
	fill = c(cbPalette[2],cbPalette[1]),
	alpha = 0.50,
	cex = 1,
	cat.cex = 0.7,
	fontfamily = "sans",
	cat.fontfamily = "sans",
	cat.fontface = "bold",
	cat.dist=c(0.05,0.05),
	margin = 0.05
)
}

#############
## Heatmap ##
#############
# placenta-specific DE chrY genes
#dt.target<-dt.exp.all[Tissue!="Breast" & ensembl_gene_id %in% unique(unlist(venn.list))]
dt.target<-dt.exp.all.chrY[Tissue!="Breast" & ensembl_gene_id %in% unique(unlist(venn.list))]
dt.target[,`Tissue`:=gsub("_"," ",Tissue)] # replace "_" with " "
dt.target[baseMean<my.baseMean,log2FoldChange:=NA]

dt.bar<-data.table::dcast(merge(dt.target,dt.ensg[,.(ensembl_gene_id,hgnc_symbol)]), ensembl_gene_id+hgnc_symbol~Tissue, value.var="log2FoldChange")
dt.bar$rowMean=rowMeans(dt.bar[,3:22],na.rm=T); dt.bar<-dt.bar[order(rowMean)]
dt.bar[,`Type`:=ifelse(ensembl_gene_id %in% venn.list$Placenta[!venn.list$Placenta %in% venn.list$GTEx], "Placenta-specific",
					   ifelse(ensembl_gene_id %in% venn.list$GTEx[!venn.list$GTEx%in% venn.list$Placenta], "GTEx-specific", "Common")
					   )]
merged.log2fc<-as.matrix(dt.bar[,3:22])
rownames(merged.log2fc)<-dt.bar$hgnc_symbol

my.min.scale=floor(min(as.vector(merged.log2fc),na.rm=T))
my.max.scale=ceiling(max(as.vector(merged.log2fc),na.rm=T))

# breaks: a sequence of numbers that covers the range of values in mat and is one element longer than color vector
my.break=seq(my.min.scale,my.max.scale,by=0.1)
# colors: one less than breaks
my.col.scale=colorRampPalette(rev(brewer.pal(9, "Blues")))(length(my.break)-1)

ann_genes<-data.frame(`Type`=factor(dt.bar$Type),row.names=dt.bar$hgnc_symbol)
ann_colors<-list(`Type`=c(`Common`=cbPalette[4],`GTEx-specific`='black', `Placenta-specific`=cbPalette[2]))

#pheatmap::pheatmap(merged.log2fc, cluster_rows=F, cluster_cols=T, breaks=my.break, annotation_row=ann_genes, annotation_colors=ann_colors, color =my.col.scale, width=10, height=15, fontsize_col=17, fontsize_row=17)
file.name<-file.path("~/results/RNA-Seq/GTEx/Figures/Heatmap",paste0("Heatmap.GTEx.chrY.",my.baseMean,".tiff"))
pheatmap::pheatmap(merged.log2fc, cluster_rows=F, cluster_cols=T, breaks=my.break, annotation_row=ann_genes, annotation_colors=ann_colors, color =my.col.scale, width=10, height=15, fontsize_col=17, fontsize_row=14, filename=file.name)
