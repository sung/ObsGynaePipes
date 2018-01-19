#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 2/Nov/2017
# Last modified 2/Nov/2017

library(data.table)
TR_PREFIX="GRCh37"
load(paste0("~/data/Annotation/Ensembl/",TR_PREFIX,"/dt.ensg.RData")) # load 'gr.ensg, dt.ensg'
source("~/Pipelines/bin/R/GTEx/local.R")# defines my.tissue dt.tissue
source("~/Pipelines/config/graphic.R")

# Female-biased genes
dt.gtex=dt.deg.all.chrX[Tissue!="Placenta" & log2FoldChange>0,.(`source`="GTEx",`No_gtex_tissue`=.N),ensembl_gene_id]
dt.pt=dt.deg.all.chrX[Tissue=="Placenta" & log2FoldChange>0,.(`source`="Placenta"),ensembl_gene_id]
dt.merged<-merge(dt.gtex, dt.pt, all.x=T, all.y=T)

dt.target.ensg<-rbind(
					data.table(`Type`="Escaped, Common",`ensembl_gene_id`=dt.merged[!is.na(source.x) & !is.na(source.y),ensembl_gene_id]),
					data.table(`Type`="Escaped, GTEx",`ensembl_gene_id`=dt.merged[!is.na(source.x) & is.na(source.y),ensembl_gene_id]),
					data.table(`Type`="Escaped, Placenta",`ensembl_gene_id`=dt.merged[is.na(source.x) & !is.na(source.y),ensembl_gene_id])
					)

##############################################
## Boxplot of log2FoldChange across tissues ##
##############################################
dt.target=merge(dt.exp.all.chrX, dt.target.ensg)
#dt.target[,`Tissue`:=gsub("_"," ",Tissue)] # replace "_" with " "
dt.target[,`Source`:=ifelse(Tissue=="Placenta","Placenta","GTEx")]
p1<-ggplot(dt.target, aes(Tissue, log2FoldChange)) + 
	geom_boxplot(aes(fill=Source),outlier.shape=NA,alpha=.7) +
	scale_x_discrete(limits=c(dt.tissue$SMTS,"Placenta")) +
	scale_fill_manual(values=c(`GTEx`=cbPalette[1],`This Study`=cbPalette[2])) +
	ylim(-.5,1) +
	facet_wrap(~Type,nrow=3) +
	theme_Publication() +
	theme(axis.text.x=element_text(angle=45, hjust=1),plot.title = element_text(hjust = 0), legend.position="none")
##
dt.fc<-rbind(
		dt.target[,.(`Avg`="Mean",`log2FC`=mean(log2FoldChange),.N),.(Type,Source,ensembl_gene_id)],
		dt.target[,.(`Avg`="Median",`log2FC`=median(log2FoldChange),.N),.(Type,Source,ensembl_gene_id)]
		)
p2<-ggplot(dt.fc, aes(Type, log2FC)) + 
	geom_boxplot(aes(fill=Source),outlier.shape=NA,alpha=.7) +
	scale_fill_manual(values=c(`GTEx`=cbPalette[1],`Placenta`=cbPalette[2])) +
	ylim(-.5,1) +
	labs(x="") +
	facet_wrap(~Avg) +
	theme_Publication() +
	theme(axis.text.x=element_text(angle=45, hjust=1),plot.title = element_text(hjust = 0))

file.name<-file.path("~/results/chrX/Figures/Boxplot",paste("chrX.escapee.log2FC",my.type,time.stamp,"pdf",sep="."))
pdf(file=file.name, width=12, height=8, title="Escapee: Log2FC")
print(p1)
print(p2)
dev.off()

write.csv(dt.target[order(Type,Source,Tissue,ensembl_gene_id)], file=gzfile(file.path("~/results/chrX/CSV/chrX.escapee.GTEx.Placenta.log2FoldChange.csv.gz")),row.names=F, quote=F)

