#!/usr/bin/Rscript --vanilla

pdf ( file=paste0(methylkit_dir,'/',my_prefix,'.genome.wide.',format(Sys.time(), "%Y-%m-%d_%I%p"), ".pdf") ) # /scratch/ssg29/results/SLX-8075.SLX-8077.SLX-8081.v1/MethylKit/8oxBS.DATE.pdf

########################################
# Get descriptive stats on methylation #
########################################
print ("counting records within filtered.myobj")
for (i in 1:length(filtered.myobj)) {
	print(nrow(filtered.myobj[[i]]))

	getMethylationStats(filtered.myobj[[i]],plot=F,both.strands=F)
	getCoverageStats(filtered.myobj[[i]],plot=F,both.strands=F)

	# plot histogram
	getMethylationStats(filtered.myobj[[i]],plot=T,both.strands=F)
	getCoverageStats(filtered.myobj[[i]],plot=T,both.strands=F)
}


print ("counting records of meth")
nrow(meth)

getCorrelation(meth, method="pearson",plot=T) #default:"pearson", other options are "kendall" and "spearman")
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
PCASamples(meth, screeplot=TRUE)
PCASamples(meth)

##
## Differentially Methylated Cytosine
# visualize the distribution of hypo/hyper-methylated bases/regions per chromosome
diffMethPerChr(top.dmc,plot=F, qvalue.cutoff=myQvalue, meth.cutoff=myDiffvalue)
diffMethPerChr(top.dmc,plot=T, qvalue.cutoff=myQvalue, meth.cutoff=myDiffvalue)
write.csv(diffMethPerChr(top.dmc,plot=FALSE, qvalue.cutoff=myQvalue, meth.cutoff=myDiffvalue), file=paste0(methylkit_dir,'/',my_prefix,".dmc.per.chr.txt"))

################################################
# Annotate differentially methylated CpG bases #
# over exons/introns/promoters/TSSes
# Resolution: base
################################################
# gene.obj is defined from 'bin/R/run_methylkit.8ox.R'
# promoter/exon/intron/TSS
met.anno=annotate.WithGenicParts(top.dmc, gene.obj) # isa 'annotationByGenicParts'
met.anno

getTargetAnnotationStats(met.anno, percentage=TRUE, precedence=TRUE)
getFeatsWithTargetsStats(met.anno, percentage=TRUE)
plotTargetAnnotation(met.anno, precedence=TRUE, main="differential methylation annotation")

## get the transcripts 
## of DMC
dmc<-getData(top.dmc) # isa data.frame
dmc.enst=getAssociationWithTSS(met.anno) # isa data.frame (todo: print out gene names and compared with DEG)
dmc$enst<-dmc.enst$feature.name # add a new column (assuming the order or record is the same between dmc & dmc.enst)
dmc$dist.to.feature<-dmc.enst$dist.to.feature

# get the gene names
# of DMC
fields <- c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol","description", "gene_biotype", "transcript_biotype")
dummy.enst=getBM(attributes = fields, filters = "ensembl_transcript_id", values = dmc$enst, mart = ensembl) # is a data.frame
dmc=merge(dmc, dummy.enst, by.x="enst", by.y="ensembl_transcript_id")
# order by meth.diff
write.csv(head(dmc[order(dmc$ensembl_gene_id, dmc$meth.diff),]), file=paste0(methylkit_dir,'/',my_prefix,".dmc",myDiffvalue,"p.txt"))

################################################
# Annotate differentially methylated CpG bases #
# over CPGi
# Resolution: base
################################################
cpg.obj=read.feature.flank(location=cgi.file, feature.flank.name=c("CpGi","shores")) # a ‘GenomicRangesList’ contatining one GRanges object for flanks and one for GRanges object for the main feature.
met.cgi.anno=annotate.WithFeature.Flank(top.dmc, cpg.obj$CpGi, cpg.obj$shores, feature.name="CpGi",flank.name="shores")
met.cgi.anno

getTargetAnnotationStats(met.cgi.anno, percentage=TRUE, precedence=TRUE)
getFeatsWithTargetsStats(met.cgi.anno, percentage=TRUE)
plotTargetAnnotation(met.cgi.anno,col=c("green","gray","white"), main="differential methylation annotation")

# user-defined annotations here
#diffmyAnn=annotate.WithFeature(top.dmc, my.anno.obj, feature.name="my_feature")

#######################################
# Regional Analysis                   #
# promoter as a unit of analysis      # 
# rather than individual methylated-C #
# Resolution: region (promoter)       #
#######################################
promoters=regionCounts(filtered.myobj, gene.obj$promoters) # is a methylRawList object with length(filtered.myobj) methylRaw objects
promoters.unite=unite(promoters,destrand=TRUE) # is a 'methylBase'
getCorrelation(promoters.unite,method="pearson",plot=T)
clusterSamples(promoters.unite, dist="correlation", method="ward", plot=TRUE)
PCASamples(promoters.unite, screeplot=TRUE)
PCASamples(promoters.unite)

promoters.diff=calculateDiffMeth(promoters.unite, num.cores=np) # is a 'methylDiff'
promoters.top.dmr=get.methylDiff(promoters.diff,difference=myDiffvalue,qvalue=myQvalue) # isa methylDiff
promoters.dmr=getData(promoters.top.dmr) 

# visualize the distribution of hypo/hyper-methylated bases/regions per chromosome
diffMethPerChr(promoters.top.dmr,plot=F, qvalue.cutoff=myQvalue, meth.cutoff=myDiffvalue)
diffMethPerChr(promoters.top.dmr,plot=T, qvalue.cutoff=myQvalue, meth.cutoff=myDiffvalue)

# promoter/exon/intron/TSS
promoters.anno=annotate.WithGenicParts(promoters.top.dmr, gene.obj) # isa 'annotationByGenicParts'
promoters.anno
promoters.enst=getAssociationWithTSS(promoters.anno)   #todo: compared with DEG 
														# where does the row name (index) come from?
promoters.dmr$enst<-promoters.enst$feature.name # add a new column (assuming the order or record is the same between promoters.dmr & promoters.enst)

# get gene name
fields <- c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol","description", "gene_biotype", "transcript_biotype")
dummy.enst=getBM(attributes = fields, filters = "ensembl_transcript_id", values = promoters.dmr$enst, mart = ensembl) # is a data.frame
promoters.dmr=merge(promoters.dmr, dummy.enst, by.x="enst", by.y="ensembl_transcript_id")
# order by meth.diff
write.csv(promoters.dmr[order(promoters.dmr$meth.diff, decreasing=T),], file=paste0(methylkit_dir,'/',my_prefix,".promoter.dmr",myDiffvalue,"p.txt"))

# CPGi
promoters.cgi.anno=annotate.WithFeature.Flank(promoters.top.dmr, cpg.obj$CpGi, cpg.obj$shores, feature.name="CpGi",flank.name="shores")

dev.off()

if(FALSE){
##########################
# Tiling Window Analysis
##########################
tiles=tileMethylCounts(filtered.myobj,win.size=1000,step.size=1000)
head(tiles[[1]],3)

##
## Coerce
##
class(meth)
as(meth,"GRanges")

class(myDiff)
as(myDiff,"GRanges")

##
## Select
##
select(meth,1:5) # get first 10 rows of a methylBase object
myDiff[21:25,] # get 5 rows of a methylDiff object

# Percent methylation values can be extracted from methylBase object by using percMethylation function.
# creates a matrix containing percent methylation values
#perc.meth=percMethylation(meth)

}
