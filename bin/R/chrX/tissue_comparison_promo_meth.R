#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 3/Jul/2017
# Last modified 3/Jul/2017

library(data.table)
TR_PREFIX="GRCh37"
RNASEQ_TYPE="TOTAL" # TOTAL|SMALL (miRNA)
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/chrX/local.R")# defines time.stamp dt.xci.anno, dt.meth.chrX, dt.exp.pops
source("~/Pipelines/bin/R/GTEx/local.R")# defines my.tissue dt.tissue
										# dl.exp.gtex, dt.chrX.gene, dt.exp.all.chrX, dt.deg.all.chrX, dt.balaton.all.deg
########################################
## 3/Jul/2017 as per Gordon's request ##
########################################

########################
## Genes of Interests ##
########################
target.ensg<-list(
	`Inactivated`=dt.exp.all.chrX[Tissue=="Placenta" & baseMean>10 & new.padj>0.3  & abs(log2FoldChange) <= log2(1.1) & 
								  ensembl_gene_id %in% dt.exp.all.chrX[Tissue!="Placenta" & baseMean>10 & new.padj>0.3  & abs(log2FoldChange) <= log2(1.1), .N, ensembl_gene_id][N>=14, ensembl_gene_id],as.character(ensembl_gene_id)], 
	`Known.Escaped`=dt.deg.all.chrX[Tissue=="Placenta" & log2FoldChange>0 & ensembl_gene_id %in% dt.deg.all.chrX[Tissue!="Placenta",unique(ensembl_gene_id)], ensembl_gene_id],
	`PT.Novel.Escaped`=dt.deg.all.chrX[Tissue=="Placenta" & log2FoldChange>0 & !ensembl_gene_id %in% dt.deg.all.chrX[Tissue!="Placenta",unique(ensembl_gene_id)], ensembl_gene_id],
	`PT.Male.biased`=dt.deg.all.chrX[Tissue=="Placenta" & log2FoldChange<0,ensembl_gene_id],
	`PT.Novel.Male.biased`= dt.deg.all.chrX[Tissue=="Placenta" & log2FoldChange<0 & !ensembl_gene_id %in% dt.deg.all.chrX[Tissue!="Placenta",unique(ensembl_gene_id)], ensembl_gene_id]
)
# remove 5 genes below from class3 (low-expressed genes (baeMean<10) in other tissues, but not in placenta)
#target.ensg[["PT.Novel.Escaped"]]=target.ensg[["PT.Novel.Escaped"]][!target.ensg[["PT.Novel.Escaped"]] %in% dt.exp.all.chrX[ensembl_gene_id %in% target.ensg[["PT.Novel.Escaped"]] & Tissue!="Placenta",median(baseMean),ensembl_gene_id][V1<10,ensembl_gene_id]]
#target.ensg[["PT.Novel.Male.biased"]]=target.ensg[["PT.Novel.Male.biased"]][!target.ensg[["PT.Novel.Male.biased"]] %in% dt.exp.all.chrX[ensembl_gene_id %in% target.ensg[["PT.Novel.Male.biased"]] & Tissue!="Placenta",median(baseMean),ensembl_gene_id][V1<10,ensembl_gene_id]]
target.ensg[["PT.DEG.UNIQ"]]=c(target.ensg[["PT.Novel.Escaped"]], target.ensg[["PT.Novel.Male.biased"]])

########################## 
## Prepare a data table ##
##########################
dt.chrX.gene.class<-dt.ensg[chromosome_name=="X",.(ensembl_gene_id,hgnc_symbol)]
dt.chrX.gene.class[ensembl_gene_id %in% target.ensg[["Inactivated"]],`type`:=1] # N=282
dt.chrX.gene.class[ensembl_gene_id %in% target.ensg[["Known.Escaped"]],`type`:=2] # N=25
dt.chrX.gene.class[ensembl_gene_id %in% target.ensg[["PT.Novel.Escaped"]],`type`:=3] # N=22 (or 22-5)
dt.chrX.gene.class[ensembl_gene_id %in% target.ensg[["PT.Male.biased"]],`type`:=4] # N=13

#################
## methlatyion ##
#################
# this R object from the following script: 'bin/R/PT.WGoxBS/local.R'
min.doc=10 # 1, 5, 10
load(file.path("~/results/RoadMap/X-inactivation/RData",paste0("dt.meth.chrX.per.tissue.CG.",min.doc,".RData"))) 
# this overwrites 'dt.meth.chrX' from chrX/local.R
dt.meth.chrX=dt.meth.chrX.per.tissue[,.(tissue,region,ensembl_gene_id,hgnc_symbol,female.met,male.met,diff.met=female.met-male.met,num.sites)]

###################################################
## % Methylation Difference at the ChrX Promoter ##
###################################################
foo<-dcast(dt.meth.chrX[region=="Promo_10.10",.(ensembl_gene_id, diff.met, num.sites, tissue)], ensembl_gene_id~tissue, value.var="diff.met")
colnames(foo)<-c("ensembl_gene_id","WGoxBS.diff.met","SS.Bio.diff.met","SS.Tech.diff.met")
bar<-dcast(dt.meth.chrX[region=="Promo_10.10",.(ensembl_gene_id, diff.met, num.sites, tissue)], ensembl_gene_id~tissue, value.var="num.sites")
colnames(bar)<-c("ensembl_gene_id","WGoxBS.num.CpG","SS.Bio.num.CpG","SS.Tech.num.CpG")
dt.chrX.meth.promo<-merge(dt.chrX.gene.class, merge(foo,bar),by="ensembl_gene_id", all.x=TRUE)
write.csv(dt.chrX.meth.promo, file=file.path("~/results/chrX/CSV",paste0("chrX.promoter.meth.diff.WGoxBS.SureSelect.DoC.",min.doc,".csv")),row.names=F)

dt.chrX.meth.promo[,`Category`:=ifelse(type==1,"Inactivated",ifelse(type==2,"Escaped common",ifelse(type==3,"Escaped placenta","Male-biased")))]
dt.chrX.meth.promo[WGoxBS.num.CpG>=10 & SS.Bio.num.CpG>=10,.(ensembl_gene_id,hgnc_symbol,Category,WGoxBS.diff.met,SS.Bio.diff.met,WGoxBS.num.CpG,SS.Bio.num.CpG)]
write.csv(dt.chrX.meth.promo[(WGoxBS.num.CpG>=10 | SS.Bio.num.CpG>=10) & !is.na(Category),.(ensembl_gene_id,hgnc_symbol,Category,WGoxBS.diff.met,SS.Bio.diff.met,WGoxBS.num.CpG,SS.Bio.num.CpG)][order(Category)], file=file.path("~/results/chrX/CSV",paste0("chrX.promoter.meth.diff.WGoxBS.SureSelect.DoC.",10,".clean.csv")),row.names=F, quote=F)

####################################################
## % Methylation Difference at the ChrX gene-body ##
####################################################
foo<-dcast(dt.meth.chrX[region=="Genes",.(ensembl_gene_id, diff.met, num.sites, tissue)], ensembl_gene_id~tissue, value.var="diff.met")
colnames(foo)<-c("ensembl_gene_id","WGoxBS.diff.met","SS.Bio.diff.met","SS.Tech.diff.met")
bar<-dcast(dt.meth.chrX[region=="Genes",.(ensembl_gene_id, diff.met, num.sites, tissue)], ensembl_gene_id~tissue, value.var="num.sites")
colnames(bar)<-c("ensembl_gene_id","WGoxBS.num.CpG","SS.Bio.num.CpG","SS.Tech.num.CpG")
dt.chrX.meth.gene.body<-merge(dt.chrX.gene.class, merge(foo,bar), by="ensembl_gene_id",all.x=TRUE)
write.csv(dt.chrX.meth.gene.body, file=file.path("~/results/chrX/CSV",paste0("chrX.gene.body.meth.diff.WGoxBS.SureSelect.DoC.",min.doc,".csv")),row.names=F)

##############################################
## Boxplot of log2FoldChange across tissues ##
##############################################
dt.target<-dt.exp.all.chrX[ensembl_gene_id %in% target.ensg[["PT.Novel.Escaped"]]]
#dt.target[,`Tissue`:=gsub("_"," ",Tissue)] # replace "_" with " "
dt.target[,`Source`:=ifelse(Tissue=="Placenta","This Study","GTEx")]
p1<-ggplot(dt.target, aes(Tissue, log2FoldChange)) + 
	geom_boxplot(aes(fill=Source),alpha=.7) +
	scale_x_discrete(limits=c(dt.tissue$SMTS,"Placenta")) +
	scale_fill_manual(values=c(`GTEx`=cbPalette[1],`This Study`=cbPalette[2])) +
	#labs(x="") +
	theme_Publication() +
	theme(axis.text.x=element_text(angle=45, hjust=1),plot.title = element_text(hjust = 0), legend.position="none")

file.name<-file.path("~/results/chrX/Figures/Boxplot/Boxplot.PT.Novel.Escaped.GTEx.tiff")
tiff(filename=file.name,width=8, height=5,units="in",res=300, compression = 'lzw')
print(p1)
dev.off()

write.csv(merge(dt.target[,.(ensembl_gene_id, Tissue, baseMean, log2FoldChange, new.padj)], dt.ensg[,.(ensembl_gene_id,hgnc_symbol)]), file=gzfile(file.path("~/results/chrX/CSV/Placenta.specific.escaped.genes.log2FoldChange.csv.gz")),row.names=F, quote=F)

#############
## Heatmap ##
#############
# defined from GTEx/local.R
#plotHeatMaplog2FC("PT.DEG.UNIQ", is.order.position=F, is.sort.row=T, is.save=T)

dt.target<-dt.exp.all.chrX[ensembl_gene_id %in% target.ensg[["PT.DEG.UNIQ"]]]
dt.target[,`Tissue`:=gsub("_"," ",Tissue)] # replace "_" with " "
dt.target[baseMean<10,log2FoldChange:=NA]

dt.foo<-dcast(dt.target, ensembl_gene_id~Tissue, value.var="log2FoldChange")
dt.bar<-data.table(merge(dt.foo, dt.ensg[,.(ensembl_gene_id, hgnc_symbol)]))
#dt.bar[,gene_name:=ifelse(hgnc_symbol=="",as.character(ensembl_gene_id),as.character(hgnc_symbol))] # replace empy hgnc_symbol with ensembl_gene_id
#save(dt.bar, file="~/results/chrX/RData/log2FoldChange.PT.DEG.UNIQ.R")

	my.end<-ncol(dt.bar)-1 # remove last one or two column ('hgnc_symbol', 'gene_name')
	merged.log2fc<-as.matrix(dt.bar[,2:my.end])
	rownames(merged.log2fc)<-dt.bar$hgnc_symbol

	my.min.scale=floor(min(as.vector(merged.log2fc),na.rm=T))
	my.max.scale=ceiling(max(as.vector(merged.log2fc),na.rm=T))
	# breaks: a sequence of numbers that covers the range of values in mat and is one element longer than color vector
	if(my.min.scale < -1){
		if(my.max.scale>1){
			my.break=c(my.min.scale,seq(-1,1,by=0.1),my.max.scale)
		}else{
			my.break=c(my.min.scale,seq(-1,1,by=0.1))
		}
	}else{
		# assume min.scale < 1
		if(my.max.scale>1){
			my.break=c(seq(my.min.scale, 1, by=0.1), my.max.scale)
		}else{
			my.break=c(seq(my.min.scale, 1, by=0.1))
		}
	}

	my.num.less.than.zero=sum(as.numeric(my.break<0))
	# colors: one less than breaks
	my.col.blue=colorRampPalette(colors = c("blue","#f2f2ff"))(my.num.less.than.zero)
	my.col.red=colorRampPalette(colors = c("#f7eaea","red"))(length(my.break)-my.num.less.than.zero-1)
	my.col.scale=c(my.col.blue, my.col.red)

	my.filename="~/results/chrX/Figures/Heatmap/pheatmap.PT.DEG.UNIQ.gene.clustered.new.tiff"
	pheatmap::pheatmap(merged.log2fc, cluster_rows=T, cluster_cols=T, breaks=my.break, color =my.col.scale, width=10, height=15, fontsize_col=17, fontsize_row=17, filename=my.filename)

	my.filename="~/results/chrX/Figures/Heatmap/pheatmap.PT.DEG.UNIQ.gene.clustered.new.pdf"
	pheatmap::pheatmap(merged.log2fc, cluster_rows=T, cluster_cols=T, breaks=my.break, color =my.col.scale, width=10, height=15, fontsize_col=17, fontsize_row=17, filename=my.filename)
	#pheatmap::pheatmap(merged.log2fc, cluster_rows=T, cluster_cols=T, breaks=my.break, color =my.col.scale, na_col="black", width=10, height=15, fontsize_col=17, fontsize_row=17, filename=my.filename)


#######################
## Log2FC - MethDiff ##
#######################
my.type="PT"
file.name<-file.path("~/results/chrX/Figures/Boxplot",paste("chrX.log2FC.MethDiff",my.type,time.stamp,"pdf",sep="."))
pdf(file=file.name, width=12, height=8, title="ChrX: Log2FC and % Meth Diff-Seq")
dt.foo<-merge(merge(dt.exp.all.chrX[Tissue=="Placenta"], dt.chrX.gene.class), dt.meth.chrX[tissue==my.type & region=="Promo_10.10"], by=c("ensembl_gene_id","hgnc_symbol"))
p0<-ggplot(dt.foo, aes(log2FoldChange, diff.met,label=hgnc_symbol)) + 
	geom_point(color='grey10',alpha=0.5,size=6) + 
	geom_text(check_overlap = TRUE, hjust =1.2, nudge_x=0.05, size=8) +
	labs(y="% Methylation Difference at Promoter\n(Female - Male)") +
	theme_Publication()

plot.meth.diff<-function(min.num=10){
	dt.foo<-merge(merge(dt.exp.all.chrX[Tissue=="Placenta"], dt.chrX.gene.class[!is.na(type)]), dt.meth.chrX[tissue==my.type & region=="Promo_10.10" & num.sites>=min.num], by=c("ensembl_gene_id","hgnc_symbol"))
	p<-ggplot(dt.foo, aes(log2FoldChange, diff.met,label=hgnc_symbol)) + 
		coord_cartesian(xlim=c(-1,1)) +
		geom_point(aes(col=as.factor(type)),alpha=0.5,size=6) + 
		geom_text(check_overlap = TRUE, hjust =1.2, nudge_x=0.05, size=5) +
		scale_y_continuous(breaks=seq(-40,40,5)) +
		scale_colour_discrete(name="Types",labels=c(`1`="Inactivated",`2`="Escaped, common",`3`="Escaped, placenta",`4`="Male-biased")) +
		geom_hline(yintercept=0) + 
		labs(y="% Methylation Difference at Promoter\n(Female - Male)") +
		ggtitle(paste("No. CpG>=",min.num)) +
		theme_Publication() +
		theme(legend.position=c(0.85,0.9))

	p.box<-ggplot(dt.foo, aes(as.factor(type), diff.met)) + 
		geom_boxplot(aes(fill=as.factor(type)),alpha=0.5, size=1.2, outlier.size=6) + 
		scale_x_discrete(labels=c(`1`="Inactivated",`2`="Escaped, common",`3`="Escaped, placenta",`4`="Male-biased")) +
		scale_y_continuous(breaks=seq(-40,40,5)) +
		geom_hline(yintercept=0) + 
		ggtitle(paste("No. CpG>=",min.num)) + 
		labs(x='', y="") +
		theme_Publication() +
		theme(legend.position="none", axis.text.x=element_blank())
	multiplot(p,p.box,cols=2)
}

print(p0)
plot.meth.diff(5)
plot.meth.diff(10)
plot.meth.diff(50)
dev.off()
