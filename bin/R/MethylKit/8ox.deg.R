#!/usr/bin/Rscript --vanilla

pdf ( file=paste0(methylkit_dir,'/',my_prefix,'.deg.',format(Sys.time(), "%Y-%m-%d_%I%p"), ".pdf") ) # /scratch/ssg29/results/SLX-8074.SLX-8080.v1/MethylKit/8oxBS.DATE.pdf

#############
## Set DEG ##
#############
deg=read.csv(file=toptags) # DEG from the edgeR result
top.deg=deg[(deg$logFC > 1 | deg$logFC < -1) & deg$PValue < deg_PValue , ] # show genes of fold change > 2 & P<0.001

# annotation
fields <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "strand", "start_position", "end_position", "description", "gene_biotype")
top.deg.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = top.deg$ensg, mart = ensembl)
top.deg.anno=merge(top.deg, top.deg.anno, by.x="ensg", by.y="ensembl_gene_id")
top.deg.anno$chromosome_name[top.deg.anno$chromosome_name=='MT']<-'M' # replace 'MT' with 'M'
top.deg.anno$chromosome_name<-paste("chr",top.deg.anno$chromosome_name,sep='') # add 'chr' prefix
write.csv(top.deg.anno, file='/home/ssg29/scratch/results/RNA-Seq/SGA.AGA/edgeR/8oxBS.toptags.pair.edgeR.csv')

## get the transcripts of DEG 
fields <- c("hgnc_symbol","ensembl_gene_id","ensembl_transcript_id")
top.deg.enst=getBM(attributes = fields, filters = "ensembl_gene_id", values = top.deg$ensg, mart = ensembl) # is a data.frame

top.deg.range=IRanges( start=top.deg.anno$start_position,end=top.deg.anno$end_position,names=top.deg.anno$ensembl_gene_id)
top.deg.gr <- GRanges(seqnames<-Rle(top.deg.anno$chromosome_name) , ranges<-top.deg.range, strand<-Rle(top.deg.anno$strand), gene<-top.deg.anno$ensembl_gene_id)
length(top.deg.gr) # should be the same with dim(top.deg)

# gene-wise
#g=getGene(id="ENSG00000037280", type="ensembl_gene_id", mart=ensembl)
#show(g)

####################################
# Make a bed file of DEG           #
# to feed read.transcript.features #
####################################
gene.bed=read.delim(gene.file, header=F) # is a data-frame
top.deg.bed = gene.bed[gene.bed$V4 %in% top.deg.enst$ensembl_transcript_id,]
write.table(top.deg.bed, file=paste0("~/data/Annotation/",my_prefix,".top.deg.bed"), row.names=F, col.names=F, quote = F, sep = "\t")
top.deg.feature.gr=read.transcript.features(location=paste0("~/data/Annotation/",my_prefix,".top.deg.bed"), up.flank=promoter.up, down.flank=promoter.down) # a ‘GRangesList’ containing locations of exons/introns/promoters/TSSes

## intersect
#top.deg.feature.gr.promoter=intersect(gene.obj$promoters, top.deg.gr)

################################################
# Annotate differentially methylated CpG bases #
# over exons/introns/promoters/TSSes of DEG    #
# Resolution: base                             #
################################################
# annotate differentially methylated Cs with promoter/exon/intron using annotation data
deg.met.anno.b=annotate.WithGenicParts(top.dmc, top.deg.feature.gr) # isa 'annotationByGenicParts'
deg.met.anno.b
deg.met.enst=getAssociationWithTSS(deg.met.anno.b) # isa data.frame (todo: print out gene names and compared with DEG)
# get percentage/number of differently methylated regions that overlap with intron/exon/promoters
getTargetAnnotationStats(deg.met.anno.b, percentage=TRUE, precedence=TRUE)
getFeatsWithTargetsStats(deg.met.anno.b, percentage=TRUE)
# plot the percentage of diffentially methylated bases overlapping with exon/intron/promoters
plotTargetAnnotation(deg.met.anno.b, precedence=TRUE, main="differential methylation annotation over DEG")

#######################################
# Regional analysis                   #
# DEG as a unit of analysis           # 
# rather than individual methylated-C #
# Resolution: region (top.deg.enst.gr)#
#######################################
top.deg.met=regionCounts(filtered.myobj, top.deg.gr) # 'top.deg.met' is a methylRawList object with nrow(top.deg.anno) methylRaw objects
top.deg.met.unite=unite(top.deg.met,destrand=TRUE) # 'top.deg.met.unite' is a 'methylBase'
getCorrelation(top.deg.met.unite,method="pearson",plot=T)
clusterSamples(top.deg.met.unite, dist="correlation", method="ward", plot=TRUE)
PCASamples(top.deg.met.unite, screeplot=TRUE)
PCASamples(top.deg.met.unite)

deg.met.diff=calculateDiffMeth(top.deg.met.unite,num.cores=np) # is a 'methylDiff'
getData(deg.met.diff)
top.deg.dmr.anno=merge(getData(deg.met.diff), top.deg.anno, by.x=c("chr", "start", "end"), by.y=c("chromosome_name", "start_position", "end_position") )
write.csv(top.deg.dmr.anno, file='/home/ssg29/scratch/results/RNA-Seq/SGA.AGA/edgeR/8oxBS.toptags.dmr.csv')
# plot 1) methylation difference and 2) the expression fold Change per each DEG
plot(top.deg.dmr.anno$meth.diff, top.deg.dmr.anno$logFC, main="Fold Change vs. % methylation difference", ylab='logFC', xlab='% met diff')
cor(top.deg.dmr.anno$meth.diff, top.deg.dmr.anno$logFC) # correlation between 1) and 2) above

top.deg.dmr=get.methylDiff(deg.met.diff, difference=myDegDiffvalue,qvalue=myDegQvalue) # isa methylDiff (resulution: base)
getData(top.deg.dmr)

deg.met.anno=annotate.WithGenicParts(top.deg.dmr, gene.obj) # isa 'annotationByGenicParts'
deg.met.anno
getTargetAnnotationStats(deg.met.anno, percentage=TRUE, precedence=TRUE)
getFeatsWithTargetsStats(deg.met.anno, percentage=TRUE)
plotTargetAnnotation(deg.met.anno, precedence=TRUE, main="differential methylation annotation")

########################################################
## methylation profiles at the promoter regions of DEG
########################################################
top.deg.met.pro=regionCounts(filtered.myobj, top.deg.feature.gr$promoters) # 'top.deg.met' is a methylRawList object with nrow(top.deg.anno) methylRaw objects
top.deg.met.pro.unite=unite(top.deg.met.pro,destrand=TRUE) # 'top.deg.met.unite' is a 'methylBase'
getCorrelation(top.deg.met.pro.unite,method="pearson",plot=T)
clusterSamples(top.deg.met.pro.unite, dist="correlation", method="ward", plot=TRUE)
PCASamples(top.deg.met.pro.unite, screeplot=TRUE)
PCASamples(top.deg.met.pro.unite)

deg.met.diff.pro=calculateDiffMeth(top.deg.met.pro.unite,num.cores=np) # is a 'methylDiff'
getData(deg.met.diff.pro)
#top.deg.dmr.pro.anno=merge(getData(deg.met.diff.pro), top.deg.anno, by.x=c("chr", "start", "end"), by.y=c("chromosome_name", "start_position", "end_position") )
#write.csv(top.deg.dmr.pro.anno, file='/home/ssg29/scratch/results/RNA-Seq/SGA.AGA/edgeR/8oxBS.toptags.dmr.promoter.csv')
#plot(top.deg.dmr.pro.anno$meth.diff, top.deg.dmr.pro.anno$logFC) # plot 1) methylation difference and 2) the expression fold Change per each DEG
#cor(top.deg.dmr.pro.anno$meth.diff, top.deg.dmr.pro.anno$logFC) # correlation between 1) and 2) above

top.deg.dmr.pro=get.methylDiff(deg.met.diff.pro, difference=myDegDiffvalue,qvalue=myDegQvalue) # isa methylDiff (resulution: base)
getData(top.deg.dmr.pro)

#####################################################
## methylation profiles at the exon regions of DEG  #
#####################################################
top.deg.met.exon=regionCounts(filtered.myobj, unique(top.deg.feature.gr$exons)) # 'top.deg.met' is a methylRawList object with nrow(top.deg.anno) methylRaw objects
top.deg.met.exon.unite=unite(top.deg.met.exon,destrand=TRUE) # 'top.deg.met.unite' is a 'methylBase'
getCorrelation(top.deg.met.exon.unite,method="pearson",plot=T)
clusterSamples(top.deg.met.exon.unite, dist="correlation", method="ward", plot=TRUE)
PCASamples(top.deg.met.exon.unite, screeplot=TRUE)
PCASamples(top.deg.met.exon.unite)

top.deg.met.diff.exon=calculateDiffMeth(top.deg.met.exon.unite,num.cores=np) # is a 'methylDiff'
getData(top.deg.met.diff.exon)
top.deg.dmr.exon=get.methylDiff(top.deg.met.diff.exon, difference=myDegDiffvalue,qvalue=myDegQvalue) # isa methylDiff (resulution: base)
getData(top.deg.dmr.exon)

## methylation profiles at the intron regions of DEG
top.deg.met.intron=regionCounts(filtered.myobj, unique(top.deg.feature.gr$introns)) # 'top.deg.met' is a methylRawList object with nrow(top.deg.anno) methylRaw objects
top.deg.met.intron.unite=unite(top.deg.met.intron,destrand=TRUE) # 'top.deg.met.unite' is a 'methylBase'
getCorrelation(top.deg.met.intron.unite,method="pearson",plot=T)
clusterSamples(top.deg.met.intron.unite, dist="correlation", method="ward", plot=TRUE)
PCASamples(top.deg.met.intron.unite, screeplot=TRUE)
PCASamples(top.deg.met.intron.unite)

top.deg.met.diff.intron=calculateDiffMeth(top.deg.met.intron.unite,num.cores=np) # is a 'methylDiff'
getData(top.deg.met.diff.intron)
top.deg.dmr.intron=get.methylDiff(top.deg.met.diff.intron, difference=myDegDiffvalue,qvalue=myDegQvalue) # isa methylDiff (resulution: base)
getData(top.deg.dmr.intron)

## methylation profiles at the TSS of DEG
top.deg.met.tss=regionCounts(filtered.myobj, unique(top.deg.feature.gr$TSSes)) # 'top.deg.met' is a methylRawList object with nrow(top.deg.anno) methylRaw objects
top.deg.met.tss.unite=unite(top.deg.met.tss,destrand=TRUE) # 'top.deg.met.unite' is a 'methylBase'
getCorrelation(top.deg.met.tss.unite,method="pearson",plot=T)
clusterSamples(top.deg.met.tss.unite, dist="correlation", method="ward", plot=TRUE)
PCASamples(top.deg.met.tss.unite, screeplot=TRUE)
PCASamples(top.deg.met.tss.unite)

top.deg.met.diff.tss=calculateDiffMeth(top.deg.met.tss.unite,num.cores=np) # is a 'methylDiff'
getData(top.deg.met.diff.tss)
top.deg.dmr.tss=get.methylDiff(top.deg.met.diff.tss, difference=myDegDiffvalue,qvalue=myDegQvalue) # isa methylDiff (resulution: base)
getData(top.deg.dmr.tss)

			
#########################
## DMR analysis for DEG #
#########################
# get DMC within DEG
diffMetDEG=annotate.WithFeature(top.dmc, top.deg.gr, feature.name="DEG") # isa 'annotationByFeature'
top.dmc[getMembers(diffMetDEG)==1,] # getMembers returns 'matrix'

dev.off()
