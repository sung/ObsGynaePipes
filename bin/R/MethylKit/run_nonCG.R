#!/usr/bin/Rscript --vanilla
#http://code.google.com/p/methylkit

#source("http://methylkit.googlecode.com/files/install.methylKit.R")
#install.methylKit(ver="0.9.2",dependencies=TRUE)

#
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("GenomicFeatures")
#biocLite("GenomicRanges")
#biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(methylKit)
#library("GenomicFeatures")
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)

top_result_dir<-"/scratch/ssg29/results"

#result
#result_dir=paste0(methylkit_dir, "/Replicates") # /scratch/ssg29/results/SLX-8074.SLX-8080.v1/MethylKit/Replicates 
pdf ( file=paste0(result_dir,'.',format(Sys.time(), "%Y-%m-%d_%I%p"), ".pdf") ) # /scratch/ssg29/results/SLX-8074.SLX-8080.v1/MethylKit/Replicates/A010.A001.DATE.pdf

#############################
## 1. Read methylation files
#############################
file.list=list( 
"/home/ssg29/scratch/results/SLX-8074.SLX-8080.v1/BisSNP/A001/SLX-8074.SLX-8080.A001.cpg.filtered.non.CG.methylkit.txt",
"/home/ssg29/scratch/results/SLX-8074.SLX-8080.v1/BisSNP/A003/SLX-8074.SLX-8080.A003.cpg.filtered.non.CG.methylkit.txt",
"/home/ssg29/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/BisSNP/A005/SLX-8075.SLX-8077.SLX-8081.A005.cpg.filtered.non.CG.methylkit.txt",
"/home/ssg29/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/BisSNP/A007/SLX-8075.SLX-8077.SLX-8081.A007.cpg.filtered.non.CG.methylkit.txt",
"/home/ssg29/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/BisSNP/A002/SLX-8075.SLX-8077.SLX-8081.A002.cpg.filtered.non.CG.methylkit.txt",
"/home/ssg29/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/BisSNP/A004/SLX-8075.SLX-8077.SLX-8081.A004.cpg.filtered.non.CG.methylkit.txt",
"/home/ssg29/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/BisSNP/A008/SLX-8075.SLX-8077.SLX-8081.A008.cpg.filtered.non.CG.methylkit.txt",
"/home/ssg29/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/BisSNP/A006/SLX-8075.SLX-8077.SLX-8081.A006.cpg.filtered.non.CG.methylkit.txt"
)
sample.id=list("AGA2F","SGA2F","AGA3F","SGA3F","AGA5M","SGA5M","AGA8M","SGA8M")
#myobj is a class of 'methylRaw' ('methylRawList')
#myobj=read( file.list, sample.id=list(paste0(project1,".A001"),paste0(project2,".A001"),paste0(project1,".A010"),paste(project1,".A003")),
myobj=read( file.list, sample.id=sample.id,
			assembly="hg19",treatment=c(0,1,0,1,0,1,0,1),context="CpH") # is a 'methylRawList'

# Filtering samples based on read coverage
# filter if read count<5
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=NULL) # is a 'methylRawList' 

##########################################
## 2. Get descriptive stats on methylation
##########################################
# plot methylation statistics 

#> getMethylationStats(myobj[[1]],plot=F,both.strands=F)
#methylation statistics per base
#summary:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.00    0.00   66.67   56.74  100.00  100.00 
#percentiles:
#       0%       10%       20%       30%       40%       50%       60%       70% 
#  0.00000   0.00000   0.00000  25.00000  50.00000  66.66667  80.00000 100.00000 
#      80%       90%       95%       99%     99.5%     99.9%      100% 
#100.00000 100.00000 100.00000 100.00000 100.00000 100.00000 100.00000 
par(mfrow=c(4,2))
for(i in 1:length(filtered.myobj)){
	#getMethylationStats(filtered.myobj[[1]],plot=T,both.strands=F)
	getMethylationStats(filtered.myobj[[i]],plot=T,both.strands=F)
}
par(mfrow=c(1,1))

# see what the data looks like for sample 2 in myobj methylRawList
head(filtered.myobj[[1]])
#methylRaw object with 6 rows
#--------------
#    chr start   end strand coverage numCs numTs
#11 chr9 10648 10648      -        7     0     7
#12 chr9 10689 10689      +        5     0     5
#13 chr9 10690 10690      -        8     0     8
#14 chr9 10708 10708      +        5     0     5
#15 chr9 10709 10709      -        8     0     8
#16 chr9 10716 10716      +        5     0     5
#--------------
#sample.id: A002
#assembly: hg19
#context: CpG
#resolution: base


###########################################################
## 3. Get bases covered by all samples and cluster samples
###########################################################
# merge all samples to one table by using base-pair locations that are covered in all samples
# setting destrand=TRUE, will merge reads on both strans of a CpG dinucleotide. This provides better 
# coverage, but only advised when looking at CpG methylation
# 
meth.destr=unite(filtered.myobj,destrand=TRUE) # is a 'methylBase'
meth=unite(filtered.myobj,destrand=FALSE) # is a 'methylBase'
head(meth)

getCorrelation(meth,method="pearson",plot=T) #default:"pearson", other options are "kendall" and "spearman")

# Percent methylation values can be extracted from methylBase object by using percMethylation function.
# creates a matrix containing percent methylation values
perc.meth=percMethylation(meth)

# cluster all samples using correlation distance and return a tree object for plclust
#hc = clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

# cluster all samples using correlation distance and plot hiarachical clustering
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

# screeplot of principal component analysis.
PCASamples(meth, screeplot=TRUE)

# principal component anlaysis of all samples.
PCASamples(meth)

if(FALSE){
	##########################################
	## 4. Calculate differential methylation
	##########################################
	# calculate differential methylation p-values and q-values
	myDiff=calculateDiffMeth(meth,num.cores=4) # is a 'methylDiff'

	# check how data part of methylDiff object looks like
	head( myDiff )

	# get differentially methylated regions with 25% difference and qvalue<0.01
	myDiffvalue=25
	myQvalue=0.01
	myDiff25p=get.methylDiff(myDiff,difference=myDiffvalue,qvalue=myQvalue)
		# difference: cutoff for absolute value of methylation change between test and control (default:25)
		# qvalue: cutoff for qvalue of differential methylation statistic (default:0.01)
		# type: one of the "hyper","hypo" or "all" strings.  Specifies what
		#	type of differentially menthylated bases/regions should be
		#	returned.  For retrieving Hyper-methylated regions/bases
		#	type="hyper", for hypo-methylated type="hypo" (default:"all")
	write.csv(getData(myDiff25p), file=paste0(result_dir,".diff.25p.txt"))

	# get differentially hypo methylated regions with 25% difference and qvalue<0.01
	myDiff25pHypo =get.methylDiff(myDiff,difference=myDiffvalue,qvalue=myQvalue,type="hypo") 
	write.csv(getData(myDiff25pHypo), file=paste0(result_dir,".diff.25p.Hypo.txt"))

	# get differentially hyper methylated regions with 25% difference and qvalue<0.01
	myDiff25pHyper=get.methylDiff(myDiff,difference=myDiffvalue,qvalue=myQvalue,type="hyper")
	write.csv(getData(myDiff25pHyper), file=paste0(result_dir,".diff.25p.Hyper.txt"))

	# visualize the distribution of hypo/hyper-methylated bases/regions per chromosome
	diffMethPerChr(myDiff25p,plot=TRUE, qvalue.cutoff=myQvalue, meth.cutoff=myDiffvalue)
	write.csv(diffMethPerChr(myDiff25p,plot=FALSE, qvalue.cutoff=myQvalue, meth.cutoff=25), file=paste0(result_dir,".diff.per.chr.txt"))

	#######################################################
	## 5. Annotate differentially methylated bases/regions
	#######################################################
	# read-in transcript locations to be used in annotation
	# IMPORTANT: annotation files that come with the package (version >=0.5) are a subset of full annotation
	# files. Download appropriate annotation files from UCSC (or other sources) in BED format
	#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
	gene.obj=read.transcript.features(system.file("extdata", "refseq.hg19.bed.txt", package = "methylKit")) # see ~/R/x86_64-pc-linux-gnu-library/3.1/methylKit/extdata/
	cpg.obj=read.feature.flank(system.file("extdata", "cpgi.hg19.bed.txt",package = "methylKit"),feature.flank.name=c("CpGi","shores"))
	# read my own bed file
	my.anno.obj=read.bed("my_annotation.bed")

	# annotate differentially methylated Cs with promoter/exon/intron using annotation data
	diffAnn=annotate.WithGenicParts(myDiff25p, gene.obj)
	diffCpGann=annotate.WithFeature.Flank(myDiff25p,cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")
	diffmyAnn=annotate.WithFeature(myDiff25p, my.anno.obj, feature.name="my_feature")
}

#######################
## Regional analysis ##
#######################
promoters=regionCounts(filtered.myobj, gene.obj$promoters)
head(promoters[[1]])

################################################
# Convenience functions for annotation objects
################################################
# After getting the annotation of diㄦentially methylated regions, we can get the distance to TSS and nearest gene name using the getAssociationWithTSS function target.row is the row number in myDiff25p
head(getAssociationWithTSS(diffAnn))
# get percentage/number of differently methylated regions that overlap with intron/exon/promoters
getTargetAnnotationStats(diffAnn, percentage=TRUE, precedence=TRUE)
# plot the percentage of diffentially methylated bases overlapping with exon/intron/promoters
plotTargetAnnotation(diffAnn, precedence=TRUE, main="differential methylation annotation")

head(getAssociationWithTSS(diffCpGann))
# also plot the CpG island annotation the same way
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"), main="differential methylation annotation")
# get percentage of intron/exon/promoters that overlap with diㄦentially methylated bases.
getFeatsWithTargetsStats(diffCpGann,percentage=TRUE)


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

dev.off()
