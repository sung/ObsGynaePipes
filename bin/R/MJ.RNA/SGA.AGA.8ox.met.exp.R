#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

myCaller="8ox"
source ("~/Pipelines/config/DEG.R") # load config

# pdf output filename 
pdf_file_name<-"~/scratch/results/RNA-Seq/SGA.AGA/MetExp/plot"
pdf(file=paste(pdf_file_name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

#load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/top.diff.meth.Tss2.RData") # top <1% Tss (pair2,pair3,pair5,pair8)
#load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/mRNA/All.Tss.regions.common.sig.meth.diff.expressed.RData") # meanCount>=20 (comon meth.diff Tss2000, Tss1000, Tss500)
load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/mRNA/All.Tss.regions.common.sig.meth.diff.expressed2.RData") # meanCount>=1 (common meth.diff Tss2000, Tss1000, Tss500)
# /scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/mRNA/pair2.tss.regions.methdiff.expressed.RData (pair-wise, not common across pairs)

#######################################
# DEG (unified from edgeR and DESeq2) #
#######################################
#union.top.deg<-read.csv(file="~/scratch/results/RNA-Seq/SGA.AGA/DEG/OXBS/unified.top.deg.top.100.OXBS.csv", stringsAsFactors=FALSE)
#union.top.deg<-read.csv(file="~/scratch/results/RNA-Seq/SGA.AGA/DEG/OXBS/unified.top.deg.top.200.OXBS.csv", stringsAsFactors=FALSE)
#union.top.deg<-read.csv(file="~/scratch/results/RNA-Seq/SGA.AGA/DEG/LOW.PAPPA/unified.top.deg.top.100.LOW.PAPPA.csv", stringsAsFactors=FALSE)
union.top.deg<-read.csv(file="~/scratch/results/RNA-Seq/SGA.AGA/DEG/LOW.PAPPA/unified.top.deg.top.200.LOW.PAPPA.csv", stringsAsFactors=FALSE)

############################
## Overlap between DM & DE #
############################
#Tss2000
dmr.tss.deg<-merge(Tss2000, union.top.deg, by.x="Ensembl.Gene.ID", by.y="ensg")
for(i in c("Ensembl.Gene.ID","feature.name","chr","strand","Associated.Gene.Name"))(dmr.tss.deg[i]<-droplevels(dmr.tss.deg[i])) #re-level
lapply(split(dmr.tss.deg, dmr.tss.deg$Associated.Gene.Name), function(i) round(mean(c(i$meth.diff.x.x, i$meth.diff.y.x, i$meth.diff.x.y, i$meth.diff.y.y)),1) )

#Tss1000
dmr.tss.deg<-merge(Tss1000, union.top.deg, by.x="Ensembl.Gene.ID", by.y="ensg")
for(i in c("Ensembl.Gene.ID","feature.name","chr","strand","Associated.Gene.Name"))(dmr.tss.deg[i]<-droplevels(dmr.tss.deg[i])) #re-level
lapply(split(dmr.tss.deg, dmr.tss.deg$Associated.Gene.Name), function(i) round(mean(c(i$meth.diff.x.x, i$meth.diff.y.x, i$meth.diff.x.y, i$meth.diff.y.y)),1) )

#Tss500
dmr.tss.deg<-merge(Tss500, union.top.deg, by.x="Ensembl.Gene.ID", by.y="ensg")
for(i in c("Ensembl.Gene.ID","feature.name","chr","strand","Associated.Gene.Name"))(dmr.tss.deg[i]<-droplevels(dmr.tss.deg[i])) #re-level
lapply(split(dmr.tss.deg, dmr.tss.deg$Associated.Gene.Name), function(i) round(mean(c(i$meth.diff.x.x, i$meth.diff.y.x, i$meth.diff.x.y, i$meth.diff.y.y)),1) )

##########################
## Read MJ RNA-Seq data ##
##########################
#DESeq2 RData for this sampleType
deseq.RData<-"/home/ssg29/scratch/results/RNA-Seq/SGA.AGA/DESeq2/OXBS/deseq.OXBS.RData" # dds, ddsFlt, res, resFlt, rld, rld.nb, vsd 
load(deseq.RData) #dds, ddsFlt, res, resFlt, rld, rld.nb, vsd

readCount<-as.data.frame(counts(dds))
# raw count
Tss2000.readCount<-readCount[rownames(readCount) %in% Tss2000$Ensembl.Gene.ID,] # a list of data.frame
Tss2000.readCount$Ensembl.Gene.ID<-rownames(Tss2000.readCount)
Tss1000.readCount<-readCount[rownames(readCount) %in% Tss1000$Ensembl.Gene.ID,] # a list of data.frame
Tss1000.readCount$Ensembl.Gene.ID<-rownames(Tss1000.readCount)
Tss500.readCount<-readCount[rownames(readCount) %in% Tss500$Ensembl.Gene.ID,] # a list of data.frame
Tss500.readCount$Ensembl.Gene.ID<-rownames(Tss500.readCount)

# Merge methylation & expression 
Tss2000<-merge(Tss2000, Tss2000.readCount, by="Ensembl.Gene.ID")
Tss1000<-merge(Tss1000, Tss1000.readCount, by="Ensembl.Gene.ID")
Tss500<-merge(Tss500, Tss500.readCount, by="Ensembl.Gene.ID")

# group by gene (aggregate transcripts)
Tss500.agg<-aggregate(cbind(AGA2F, SGA2F, AGA3F, SGA3F, AGA5M, SGA5M, AGA8M, SGA8M, `66`,`68`,`72`,`74`,`75`,`77`,`73`,`65`) ~ Associated.Gene.Name, Tss500, mean)
Tss1000.agg<-aggregate(cbind(AGA2F, SGA2F, AGA3F, SGA3F, AGA5M, SGA5M, AGA8M, SGA8M, `66`,`68`,`72`,`74`,`75`,`77`,`73`,`65`) ~ Associated.Gene.Name, Tss1000, mean)
Tss2000.agg<-aggregate(cbind(AGA2F, SGA2F, AGA3F, SGA3F, AGA5M, SGA5M, AGA8M, SGA8M, `66`,`68`,`72`,`74`,`75`,`77`,`73`,`65`) ~ Associated.Gene.Name, Tss2000, mean)

#Tss500 (read-count vs. % met)
if(FALSE){
	par(mfrow=c(2,2))
	plot(Tss500.agg$AGA2F, log2(Tss500.agg$`66`+1), xlab='%met', ylab='log2(read count+1)', main="AGA2F", cex=0.3, col="#00000020", pch=20)
	text(Tss500.agg$AGA2F, log2(Tss500.agg$`66`+1), labels=Tss500.agg$Associated.Gene.Name)

	plot(Tss500.agg$SGA2F, log2(Tss500.agg$`68`+1), xlab='%met', ylab='log2(read count+1)', main="SGA2F", cex=0.3, col="#00000020", pch=20)
	text(Tss500.agg$SGA2F, log2(Tss500.agg$`68`+1), labels=Tss500.agg$Associated.Gene.Name)

	plot(Tss500.agg$AGA3F, log2(Tss500.agg$`72`+1), xlab='%met', ylab='log2(read count+1)', main="AGA3F", cex=0.3, col="#00000020", pch=20)
	text(Tss500.agg$AGA3F, log2(Tss500.agg$`72`+1), labels=Tss500.agg$Associated.Gene.Name)

	plot(Tss500.agg$SGA3F, log2(Tss500.agg$`74`+1), xlab='%met', ylab='log2(read count+1)', main="SGA3F", cex=0.3, col="#00000020", pch=20)
	text(Tss500.agg$SGA3F, log2(Tss500.agg$`74`+1), labels=Tss500.agg$Associated.Gene.Name)

	par(mfrow=c(2,2))
	plot(Tss500.agg$AGA5M, log2(Tss500.agg$`75`+1), xlab='%met', ylab='log2(read count+1)', main="AGA5M")
	plot(Tss500.agg$SGA5M, log2(Tss500.agg$`77`+1), xlab='%met', ylab='log2(read count+1)', main="SGA5M")
	plot(Tss500.agg$AGA8M, log2(Tss500.agg$`73`+1), xlab='%met', ylab='log2(read count+1)', main="AGA8M")
	plot(Tss500.agg$SGA8M, log2(Tss500.agg$`65`+1), xlab='%met', ylab='log2(read count+1)', main="SGA8M")

	#Tss1000 (read-count vs. % met)
	par(mfrow=c(4,2))
	plot(Tss1000.agg[,c("AGA2F", "66")], xlab='%met', ylab='read count', main="AGA2F")
	plot(Tss1000.agg[,c("SGA2F", "68")], xlab='%met', ylab='read count', main="SGA2F")
	plot(Tss1000.agg[,c("AGA3F", "72")], xlab='%met', ylab='read count', main="AGA3F")
	plot(Tss1000.agg[,c("SGA3F", "74")], xlab='%met', ylab='read count', main="SGA3F")
	plot(Tss1000.agg[,c("AGA5M", "75")], xlab='%met', ylab='read count', main="AGA5M")
	plot(Tss1000.agg[,c("SGA5M", "77")], xlab='%met', ylab='read count', main="SGA5M")
	plot(Tss1000.agg[,c("AGA8M", "73")], xlab='%met', ylab='read count', main="AGA8M")
	plot(Tss1000.agg[,c("SGA8M", "65")], xlab='%met', ylab='read count', main="SGA8M")
	par(mfrow=c(1,1))

	#Tss2000 (read-count vs. % met)
	par(mfrow=c(2,4))
	plot(Tss2000.agg[,c("AGA2F", "66")], xlab='%met', ylab='read count', main="AGA2F")
	plot(Tss2000.agg[,c("SGA2F", "68")], xlab='%met', ylab='read count', main="SGA2F")
	plot(Tss2000.agg[,c("AGA3F", "72")], xlab='%met', ylab='read count', main="AGA3F")
	plot(Tss2000.agg[,c("SGA3F", "74")], xlab='%met', ylab='read count', main="SGA3F")
	plot(Tss2000.agg[,c("AGA5M", "75")], xlab='%met', ylab='read count', main="AGA5M")
	plot(Tss2000.agg[,c("SGA5M", "77")], xlab='%met', ylab='read count', main="SGA5M")
	plot(Tss2000.agg[,c("AGA8M", "73")], xlab='%met', ylab='read count', main="AGA8M")
	plot(Tss2000.agg[,c("SGA8M", "65")], xlab='%met', ylab='read count', main="SGA8M")
	par(mfrow=c(1,1))
}

#####################
# %met & expression #
#####################
par(mfrow=c(2,2))
# Tss500
#Pair2
dummy<-rbind(
	as.matrix(with(Tss500.agg, data.frame(met=AGA2F, cnt=`66`, logCnt=log2(`66`+1), gene=Associated.Gene.Name, condition="AGA", col="blue"))),
	as.matrix(with(Tss500.agg, data.frame(met=SGA2F, cnt=`68`, logCnt=log2(`68`+1), gene=Associated.Gene.Name, condition="SGA", col="red")))
)
plot(dummy[,c("met","logCnt")], main="Tss500 Pair2 (female)", cex=0.3, col="#00000020", pch=20)
text(dummy[,c("met","logCnt")], labels=dummy[,c("gene")], col=dummy[,c('col')])
#Pair3
dummy<-rbind(
	as.matrix(with(Tss500.agg, data.frame(met=AGA3F, cnt=`72`, logCnt=log2(`72`+1), gene=Associated.Gene.Name, condition="AGA", col="blue"))),
	as.matrix(with(Tss500.agg, data.frame(met=SGA3F, cnt=`74`, logCnt=log2(`74`+1), gene=Associated.Gene.Name, condition="SGA", col="red")))
)
plot(dummy[,c("met","logCnt")], main="Tss500 Pair3 (female)", cex=0.3, col="#00000020", pch=20)
text(dummy[,c("met","logCnt")], labels=dummy[,c("gene")], col=dummy[,c('col')])
#Pair5
dummy<-rbind(
	as.matrix(with(Tss500.agg, data.frame(met=AGA5M, cnt=`75`, logCnt=log2(`75`+1), gene=Associated.Gene.Name, condition="AGA", col="blue"))),
	as.matrix(with(Tss500.agg, data.frame(met=SGA5M, cnt=`77`, logCnt=log2(`77`+1), gene=Associated.Gene.Name, condition="SGA", col="red")))
)
plot(dummy[,c("met","logCnt")], main="Tss500 Pair5 (male)", cex=0.3, col="#00000020", pch=20)
text(dummy[,c("met","logCnt")], labels=dummy[,c("gene")], col=dummy[,c('col')])
#Pair8
dummy<-rbind(
	as.matrix(with(Tss500.agg, data.frame(met=AGA8M, cnt=`73`, logCnt=log2(`73`+1), gene=Associated.Gene.Name, condition="AGA", col="blue"))),
	as.matrix(with(Tss500.agg, data.frame(met=SGA8M, cnt=`65`, logCnt=log2(`65`+1), gene=Associated.Gene.Name, condition="SGA", col="red")))
)
plot(dummy[,c("met","logCnt")], main="Tss500 Pair8 (male)", cex=0.3, col="#00000020", pch=20)
text(dummy[,c("met","logCnt")], labels=dummy[,c("gene")], col=dummy[,c('col')])

# Tss1000
#Pair2
dummy<-rbind(
	as.matrix(with(Tss1000.agg, data.frame(met=AGA2F, cnt=`66`, logCnt=log2(`66`+1), gene=Associated.Gene.Name, condition="AGA", col="blue"))),
	as.matrix(with(Tss1000.agg, data.frame(met=SGA2F, cnt=`68`, logCnt=log2(`68`+1), gene=Associated.Gene.Name, condition="SGA", col="red")))
)
plot(dummy[,c("met","logCnt")], main="Tss1000 Pair2 (female)", cex=0.3, col="#00000020", pch=20)
text(dummy[,c("met","logCnt")], labels=dummy[,c("gene")], col=dummy[,c('col')])
#Pair3
dummy<-rbind(
	as.matrix(with(Tss1000.agg, data.frame(met=AGA3F, cnt=`72`, logCnt=log2(`72`+1), gene=Associated.Gene.Name, condition="AGA", col="blue"))),
	as.matrix(with(Tss1000.agg, data.frame(met=SGA3F, cnt=`74`, logCnt=log2(`74`+1), gene=Associated.Gene.Name, condition="SGA", col="red")))
)
plot(dummy[,c("met","logCnt")], main="Tss1000 Pair3 (female)", cex=0.3, col="#00000020", pch=20)
text(dummy[,c("met","logCnt")], labels=dummy[,c("gene")], col=dummy[,c('col')])
#Pair5
dummy<-rbind(
	as.matrix(with(Tss1000.agg, data.frame(met=AGA5M, cnt=`75`, logCnt=log2(`75`+1), gene=Associated.Gene.Name, condition="AGA", col="blue"))),
	as.matrix(with(Tss1000.agg, data.frame(met=SGA5M, cnt=`77`, logCnt=log2(`77`+1), gene=Associated.Gene.Name, condition="SGA", col="red")))
)
plot(dummy[,c("met","logCnt")], main="Tss1000 Pair5 (male)", cex=0.3, col="#00000020", pch=20)
text(dummy[,c("met","logCnt")], labels=dummy[,c("gene")], col=dummy[,c('col')])
#Pair8
dummy<-rbind(
	as.matrix(with(Tss1000.agg, data.frame(met=AGA8M, cnt=`73`, logCnt=log2(`73`+1), gene=Associated.Gene.Name, condition="AGA", col="blue"))),
	as.matrix(with(Tss1000.agg, data.frame(met=SGA8M, cnt=`65`, logCnt=log2(`65`+1), gene=Associated.Gene.Name, condition="SGA", col="red")))
)
plot(dummy[,c("met","logCnt")], main="Tss1000 Pair8 (male)", cex=0.3, col="#00000020", pch=20)
text(dummy[,c("met","logCnt")], labels=dummy[,c("gene")], col=dummy[,c('col')])

# Tss2000
#Pair2
dummy<-rbind(
	as.matrix(with(Tss2000.agg, data.frame(met=AGA2F, cnt=`66`, logCnt=log2(`66`+1), gene=Associated.Gene.Name, condition="AGA", col="blue"))),
	as.matrix(with(Tss2000.agg, data.frame(met=SGA2F, cnt=`68`, logCnt=log2(`68`+1), gene=Associated.Gene.Name, condition="SGA", col="red")))
)
plot(dummy[,c("met","logCnt")], main="Tss2000 Pair2 (female)", cex=0.3, col="#00000020", pch=20)
text(dummy[,c("met","logCnt")], labels=dummy[,c("gene")], col=dummy[,c('col')])
#Pair3
dummy<-rbind(
	as.matrix(with(Tss2000.agg, data.frame(met=AGA3F, cnt=`72`, logCnt=log2(`72`+1), gene=Associated.Gene.Name, condition="AGA", col="blue"))),
	as.matrix(with(Tss2000.agg, data.frame(met=SGA3F, cnt=`74`, logCnt=log2(`74`+1), gene=Associated.Gene.Name, condition="SGA", col="red")))
)
plot(dummy[,c("met","logCnt")], main="Tss2000 Pair3 (female)", cex=0.3, col="#00000020", pch=20)
text(dummy[,c("met","logCnt")], labels=dummy[,c("gene")], col=dummy[,c('col')])
#Pair5
dummy<-rbind(
	as.matrix(with(Tss2000.agg, data.frame(met=AGA5M, cnt=`75`, logCnt=log2(`75`+1), gene=Associated.Gene.Name, condition="AGA", col="blue"))),
	as.matrix(with(Tss2000.agg, data.frame(met=SGA5M, cnt=`77`, logCnt=log2(`77`+1), gene=Associated.Gene.Name, condition="SGA", col="red")))
)
plot(dummy[,c("met","logCnt")], main="Tss2000 Pair5 (male)", cex=0.3, col="#00000020", pch=20)
text(dummy[,c("met","logCnt")], labels=dummy[,c("gene")], col=dummy[,c('col')])
#Pair8
dummy<-rbind(
	as.matrix(with(Tss2000.agg, data.frame(met=AGA8M, cnt=`73`, logCnt=log2(`73`+1), gene=Associated.Gene.Name, condition="AGA", col="blue"))),
	as.matrix(with(Tss2000.agg, data.frame(met=SGA8M, cnt=`65`, logCnt=log2(`65`+1), gene=Associated.Gene.Name, condition="SGA", col="red")))
)
plot(dummy[,c("met","logCnt")], main="Tss2000 Pair8 (male)", cex=0.3, col="#00000020", pch=20)
text(dummy[,c("met","logCnt")], labels=dummy[,c("gene")], col=dummy[,c('col')])

###############################
## diff.Met & diff.expression #
###############################
par(mfrow=c(2,2))
#Tss500
dummy<-rbind(
	as.matrix(with(Tss500.agg, data.frame(diffMet=SGA2F-AGA2F, diffExp=log2(`68`+1)-log2(`66`+1), gene=Associated.Gene.Name, pair="Pair2", col='orange'))),
	as.matrix(with(Tss500.agg, data.frame(diffMet=SGA3F-AGA3F, diffExp=log2(`74`+1)-log2(`72`+1), gene=Associated.Gene.Name, pair="Pair3", col='grey'))),
	as.matrix(with(Tss500.agg, data.frame(diffMet=SGA5M-AGA5M, diffExp=log2(`77`+1)-log2(`75`+1), gene=Associated.Gene.Name, pair="Pair5", col='maroon'))),
	as.matrix(with(Tss500.agg, data.frame(diffMet=SGA8M-AGA8M, diffExp=log2(`65`+1)-log2(`73`+1), gene=Associated.Gene.Name, pair="Pair8", col='black')))
)
plot(dummy[,c("diffMet","diffExp")], main="Tss500 diff(methylation) & diff(expression)", col=dummy[,c('col')])
abline(h=c(0),v=c(0), col="blue")

#Tss1000
dummy<-rbind(
	as.matrix(with(Tss1000.agg, data.frame(diffMet=SGA2F-AGA2F, diffExp=log2(`68`+1)-log2(`66`+1), gene=Associated.Gene.Name, pair="Pair2", col='orange'))),
	as.matrix(with(Tss1000.agg, data.frame(diffMet=SGA3F-AGA3F, diffExp=log2(`74`+1)-log2(`72`+1), gene=Associated.Gene.Name, pair="Pair3", col='grey'))),
	as.matrix(with(Tss1000.agg, data.frame(diffMet=SGA5M-AGA5M, diffExp=log2(`77`+1)-log2(`75`+1), gene=Associated.Gene.Name, pair="Pair5", col='maroon'))),
	as.matrix(with(Tss1000.agg, data.frame(diffMet=SGA8M-AGA8M, diffExp=log2(`65`+1)-log2(`73`+1), gene=Associated.Gene.Name, pair="Pair8", col='black')))
)
plot(dummy[,c("diffMet","diffExp")], main="Tss1000 diff(methylation) & diff(expression)", col=dummy[,c('col')])
abline(h=c(0),v=c(0), col="blue")

#Tss2000
dummy<-rbind(
	as.matrix(with(Tss2000.agg, data.frame(diffMet=SGA2F-AGA2F, diffExp=log2(`68`+1)-log2(`66`+1), gene=Associated.Gene.Name, pair="Pair2", col='orange'))),
	as.matrix(with(Tss2000.agg, data.frame(diffMet=SGA3F-AGA3F, diffExp=log2(`74`+1)-log2(`72`+1), gene=Associated.Gene.Name, pair="Pair3", col='grey'))),
	as.matrix(with(Tss2000.agg, data.frame(diffMet=SGA5M-AGA5M, diffExp=log2(`77`+1)-log2(`75`+1), gene=Associated.Gene.Name, pair="Pair5", col='maroon'))),
	as.matrix(with(Tss2000.agg, data.frame(diffMet=SGA8M-AGA8M, diffExp=log2(`65`+1)-log2(`73`+1), gene=Associated.Gene.Name, pair="Pair8", col='black')))
)
plot(dummy[,c("diffMet","diffExp")], main="Tss2000 diff(methylation) & diff(expression)", col=dummy[,c('col')])
abline(h=c(0),v=c(0), col="blue")


###################################
# Sample Clustering Based on rlog #
###################################
rlogMat <- assay(rld) #
sampleDists<- dist(t(rlogMat))
mat <- as.matrix(sampleDists) # row: samples, col: samples
colnames(mat) <- rownames(mat)<- paste(rownames(mat), paste(samples$Code,samples$Pair,samples$Sex,sep="."), sep=":")

# All Genes
my.filename=paste0(pdf_file_name,".rld.heatmap")
tiff(filename=paste0(my.filename,".tiff"),width=12,height=9,units="in",res=300, compression = 'lzw')
heatmap.2(mat, trace="none", col = hmcol, margin=c(13, 13))
dev.off()

# Only from Selected Genes
select<-cell.types$ensg

if(FALSE){
	############################################
	## oxBS Significantly DM promoter regions ##
	############################################
	tss.RData<-"~/scratch/results/Methyl-Seq/SGA.AGA/methylKit/dmr.tss.RData"
	if(file.exists(tss.RData)){
		cat("loading TSS RData...\n")
		load (tss.RData) # 
	}else{
		cat("creating TSS RData...\n")
		#load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/top.diff.meth.Tss2.RData") # top <1% Tss (pair2,pair3,pair5,pair8)
		#load("~/scratch/results/Methyl-Seq/SGA.AGA/methylKit/top.diff.meth.Tss.RData") # old version (duplicated ENST)

		tss.5k<-list(); tss1k<-list(); tss2k<-list()
		# read from the source
		# DMR Tss2000 from Pair2
		cat("Reading pair2...\n")
		load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair2.annotate.Tss.regions3.RData") # d1:Tss500, d2:Tss1000, d3:Tss2000, e1:Tss500(q<0.05), e1:Tss1000(q<0.05), e3:Tss2000(q<0.05) ...
		rm(d1,d2,d3)
		tss.5k$pair2<-e1  # tss500
		tss1k$pair2<-e2  # tss1000
		tss2k$pair2<-e3  # tss2000
		rm(e1,e2,e3)

		# DMR Tss2000 from Pair3
		cat("Reading pair3...\n")
		load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair3.annotate.Tss.regions3.RData")
		rm(d1,d2,d3)
		tss.5k$pair3<-e1 
		tss1k$pair3<-e2 
		tss2k$pair3<-e3 
		rm(e1,e2,e3)

		# DMR Tss2000 from Pair5
		cat("Reading pair5...\n")
		load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair5.annotate.Tss.regions3.RData")
		rm(d1,d2,d3)
		tss.5k$pair5<-e1 
		tss1k$pair5<-e2 
		tss2k$pair5<-e3 
		rm(e1,e2,e3)

		# DMR Tss2000 from Pair8
		cat("Reading pair8...\n")
		load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair8.annotate.Tss.regions3.RData")
		rm(d1,d2,d3)
		tss.5k$pair8<-e1 
		tss1k$pair8<-e2 
		tss2k$pair8<-e3 
		rm(e1,e2,e3)

		##########################
		# common directional DMR #
		##########################
		cat("creating dmr.tss...\n")
		dmr.tss<-list() # common DMR across all pairs
		# Tss2000
		cat("creating dmr.tss$tss2000...\n")
		dummy.pair23<-merge(tss2k$pair2, tss2k$pair3, by=c("feature.name", "chr","start","end","strand", "Ensembl.Gene.ID"))
		dummy.pair23<-dummy.pair23[dummy.pair23$meth.diff.x * dummy.pair23$meth.diff.y > 0,]
		dummy.pair58<-merge(tss2k$pair5, tss2k$pair8, by=c("feature.name", "chr","start","end","strand", "Ensembl.Gene.ID"))
		dummy.pair58<-dummy.pair58[dummy.pair58$meth.diff.x * dummy.pair58$meth.diff.y > 0,]
		dmr.tss$tss2000<-merge(dummy.pair23, dummy.pair58, by=c("feature.name", "chr","start","end","strand", "Ensembl.Gene.ID"))
		dmr.tss$tss2000<-dmr.tss$tss2000[dmr.tss$tss2000$meth.diff.x.x * dmr.tss$tss2000$meth.diff.x.y > 0,]
		# Tss1000
		cat("creating dmr.tss$tss1000...\n")
		dummy.pair23<-merge(tss1k$pair2, tss1k$pair3, by=c("feature.name", "chr","start","end","strand", "Ensembl.Gene.ID"))
		dummy.pair23<-dummy.pair23[dummy.pair23$meth.diff.x * dummy.pair23$meth.diff.y > 0,]
		dummy.pair58<-merge(tss1k$pair5, tss1k$pair8, by=c("feature.name", "chr","start","end","strand", "Ensembl.Gene.ID"))
		dummy.pair58<-dummy.pair58[dummy.pair58$meth.diff.x * dummy.pair58$meth.diff.y > 0,]
		dmr.tss$tss1000<-merge(dummy.pair23, dummy.pair58, by=c("feature.name", "chr","start","end","strand", "Ensembl.Gene.ID"))
		dmr.tss$tss1000<-dmr.tss$tss1000[dmr.tss$tss1000$meth.diff.x.x * dmr.tss$tss1000$meth.diff.x.y > 0,]
		# Tss5000
		cat("creating dmr.tss$tss500...\n")
		dummy.pair23<-merge(tss.5k$pair2, tss.5k$pair3, by=c("feature.name", "chr","start","end","strand", "Ensembl.Gene.ID"))
		dummy.pair23<-dummy.pair23[dummy.pair23$meth.diff.x * dummy.pair23$meth.diff.y > 0,]
		dummy.pair58<-merge(tss.5k$pair5, tss.5k$pair8, by=c("feature.name", "chr","start","end","strand", "Ensembl.Gene.ID"))
		dummy.pair58<-dummy.pair58[dummy.pair58$meth.diff.x * dummy.pair58$meth.diff.y > 0,]
		dmr.tss$tss500<-merge(dummy.pair23, dummy.pair58, by=c("feature.name", "chr","start","end","strand", "Ensembl.Gene.ID"))
		dmr.tss$tss500<-dmr.tss$tss500[dmr.tss$tss500$meth.diff.x.x * dmr.tss$tss500$meth.diff.x.y > 0,]

		rm(dummy.pair23,dummy.pair58)
		save(tss.5k, tss1k, tss2k, dmr.tss, file=tss.RData) 

		cat("dmr.tss RData crated...\n")
	}

	########################
	# Top 10/20% meth.diff #
	########################
	tss2k.top<-lapply(tss2k, function(i) i[abs(i$meth.diff) >= quantile( abs(i$meth.diff), 0.90 ),] ) 
	tss2k.top<-lapply(tss2k.top, function(i) unique(merge(union.top.deg, i, by.x="ensg", by.y="Ensembl.Gene.ID"))) # a list of data.frame
	tss2k.top<-lapply(tss2k.top, 
		function(i) with(i,
			data.frame(
				ensg=ensg, 
				hgnc=ifelse(is.na(hgnc_symbol.y),hgnc_symbol.x,hgnc_symbol.y),
				description=ifelse(is.na(description.y),description.x,description.y),
				gene_biotype=ifelse(is.na(gene_biotype.y),gene_biotype.x,gene_biotype.y),
				logFC=ifelse(is.na(log2FoldChange),logFC,log2FoldChange),
				logCNT=ifelse(is.na(baseMean),log(2^logCPM*20),log(baseMean)),
				FDR=ifelse(is.na(padj),FDR,padj),
				caller=ifelse(is.na(X.y),"edgeR","DESeq2"),
				meth.diff=meth.diff,
				feature.name=feature.name
			)
		)
	)

	# avg met.diff per hgnc
	sapply(split(tss2k.top$pair2, tss2k.top$pair2$hgnc), function(i) round(mean(i$meth.diff),1))
	sapply(split(tss2k.top$pair3, tss2k.top$pair3$hgnc), function(i) round(mean(i$meth.diff),1))
	sapply(split(tss2k.top$pair5, tss2k.top$pair5$hgnc), function(i) round(mean(i$meth.diff),1))
	sapply(split(tss2k.top$pair8, tss2k.top$pair8$hgnc), function(i) round(mean(i$meth.diff),1))

	sapply(split(tss2k.top$pair2, tss2k.top$pair2$hgnc), function(i) round(mean(i$logFC),2))

	# show me
	tss2k.top$pair5[order(tss2k.top$pair5$hgnc),]

	# Pair2: ENST00000504003 ENSG00000214944 ARHGEF28
	# Pair8: ENST00000510131 ENSG00000214944 ARHGEF28
	tss2k$pair3[tss2k$pair3$Ensembl.Gene.ID=="ENSG00000214944",]
	tss2k$pair8[tss2k$pair8$Ensembl.Gene.ID=="ENSG00000214944",]

	######################
	## Common meth.diff ##
	######################
	# same direction in all pairs
	dmr.tss.deg<-lapply(dmr.tss, function(i) merge(union.top.deg, i, by.x="ensg", by.y="Ensembl.Gene.ID")) # a list of data.frame
	dmr.tss.deg<-lapply(dmr.tss.deg, 
		function(i) with(i,
			data.frame(
				ensg=ensg, 
				hgnc=ifelse(is.na(hgnc_symbol.y),hgnc_symbol.x,hgnc_symbol.y),
				description=ifelse(is.na(description.y),description.x,description.y),
				gene_biotype=ifelse(is.na(gene_biotype.y),gene_biotype.x,gene_biotype.y),
				logFC=ifelse(is.na(log2FoldChange),logFC,log2FoldChange),
				logCNT=ifelse(is.na(baseMean),log(2^logCPM*20),log(baseMean)),
				FDR=ifelse(is.na(padj),FDR,padj),
				caller=ifelse(is.na(X.y),"edgeR","DESeq2"),
				#caller=ifelse(is.na(log2FoldChange),"edgeR","DESeq2"),
				col=ifelse(is.na(padj),col.x,col.y),
				enst=feature.name,
				chr=chr,
				start=start,
				end=end,
				strand=strand,
				AGA2F=AGA2F,
				SGA2F=SGA2F,
				AGA3F=AGA3F,
				SGA3F=SGA3F,
				AGA5M=AGA5M,
				SGA5M=SGA5M,
				AGA8M=AGA8M,
				SGA8M=SGA8M
			)
		)
	)

	lapply(split(dmr.tss.deg$tss2000, dmr.tss.deg$tss2000$hgnc), function(i) mean(c(i$SGA2F-i$AGA2F,i$SGA3F-i$AGA3F,i$SGA5M-i$AGA5M,i$SGA8M-i$AGA8M)))
	lapply(split(dmr.tss.deg$tss1000, dmr.tss.deg$tss1000$hgnc), function(i) mean(c(i$SGA2F-i$AGA2F,i$SGA3F-i$AGA3F,i$SGA5M-i$AGA5M,i$SGA8M-i$AGA8M)))
	lapply(split(dmr.tss.deg$tss500, dmr.tss.deg$tss500$hgnc), function(i) mean(c(i$SGA2F-i$AGA2F,i$SGA3F-i$AGA3F,i$SGA5M-i$AGA5M,i$SGA8M-i$AGA8M)))

	######################################
	# Average Read Counts of common DMRs #
	######################################
	library(DESeq2)
	#DESeq2 RData for this sampleType
	deseq.RData<-"/home/ssg29/scratch/results/RNA-Seq/SGA.AGA/DESeq2/OXBS/deseq.OXBS.RData" # dds, ddsFlt, res, resFlt, rld, rld.nb, vsd 
	load(deseq.RData)
	rlogMat <- as.data.frame(assay(rld)) #
	readCount<-as.data.frame(counts(dds))
	gene.cnt<-lapply(dmr.tss, function(i) readCount[rownames(readCount) %in% i$Ensembl.Gene.ID,]) # a list of data.frame

	#bins<-cut( rowMeans(gene.cnt$tss2000), c(0,5,10,20,50,100,max(rowMeans(gene.cnt$tss2000))), include.lowest =T, right=F)
	#bins<-cut( rowMeans(gene.cnt$tss2000), c(0,20,max(rowMeans(gene.cnt$tss2000))), include.lowest =T, right=F)
	#pie( round(table(bins)/length(bins)*100,2), labels=paste(levels(bins), ":", round(table(bins)/length(bins)*100,2), "%" ), clockwise=T)
	#pie( table(bins), labels=paste(levels(bins), ":", table(bins), " genes"), clockwise=T, main="Tss2000")
	#barplot(table(bins))

	makePie<-function(i){
		#bins<-cut( rowMeans(i), c(0,5,10,20,50,100,max(rowMeans(i))), include.lowest =T, right=F)
		bins<-cut( rowMeans(i), c(0,20,max(rowMeans(i))), include.lowest =T, right=F)
		pie( table(bins), labels=paste(levels(bins), ":", table(bins), " genes"), clockwise=T, cex=1.5)
	}
	par(mfrow=c(1,3))
	lapply(gene.cnt, makePie)
	par(mfrow=c(1,1))

	gene.cnt$tss2000[,"ensg"]<-rownames(gene.cnt$tss2000)
	gene.cnt$tss1000[,"ensg"]<-rownames(gene.cnt$tss1000)
	gene.cnt$tss500[,"ensg"]<-rownames(gene.cnt$tss500)

	readCount$ensg<-rownames(readCount)
	rlogMat$ensg<-rownames(rlogMat)
	#dmr.tss.deseq<-lapply(dmr.tss, function(i) merge(readCount, i, by.x="ensg", by.y="Ensembl.Gene.ID")) # a list of data.frame
	dmr.tss.deseq<-lapply(dmr.tss, function(i) merge(rlogMat, i, by.x="ensg", by.y="Ensembl.Gene.ID")) # a list of data.frame

	plot(x=dmr.tss.deseq$tss2000[,"66"],y=dmr.tss.deseq$tss2000[,"AGA2F"], xlab="log(cnt)", ylab="% met")

	par(mfrow=c(2,2))
	#Pair2
	plot(y=dmr.tss.deseq$tss2000[,"68"]-dmr.tss.deseq$tss2000[,"66"],x=dmr.tss.deseq$tss2000[,"meth.diff.x.x"], ylab="cnt(SGA-AGA)", xlab="%met(SGA-AGA)", main="Pair2")
	abline(h=c(0),v=c(0), col="blue")
	#Pair3
	plot(y=dmr.tss.deseq$tss2000[,"74"]-dmr.tss.deseq$tss2000[,"72"],x=dmr.tss.deseq$tss2000[,"meth.diff.y.x"], ylab="cnt(SGA-AGA)", xlab="%met(SGA-AGA)", main="Pair3")
	abline(h=c(0),v=c(0), col="blue")
	#Pair5
	plot(y=dmr.tss.deseq$tss2000[,"77"]-dmr.tss.deseq$tss2000[,"75"],x=dmr.tss.deseq$tss2000[,"meth.diff.x.y"], ylab="cnt(SGA-AGA)", xlab="%met(SGA-AGA)", main="Pair5")
	abline(h=c(0),v=c(0), col="blue")
	#Pair8
	plot(y=dmr.tss.deseq$tss2000[,"65"]-dmr.tss.deseq$tss2000[,"73"],x=dmr.tss.deseq$tss2000[,"meth.diff.y.y"], ylab="cnt(SGA-AGA)", xlab="%met(SGA-AGA)", main="Pair8")
	abline(h=c(0),v=c(0), col="blue")
	par(mfrow=c(1,1))
}



if(FALSE){
	fields <- c("ensembl_gene_id", "hgnc_symbol","description", "gene_biotype")
	p2.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = pair2$Ensembl.Gene.ID, mart = ensembl)
	p3.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = pair3$Ensembl.Gene.ID, mart = ensembl)
	p5.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = pair5$Ensembl.Gene.ID, mart = ensembl)
	p8.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = pair8$Ensembl.Gene.ID, mart = ensembl)
	##############################################################
	## correlation between diff(expression) & diff(methylation) ##
	##############################################################
	#ddsCpms from DESeq2
	#load("~/scratch/results/RNA-Seq/SGA.AGA/DESeq2/8oxBS/8oxBS.paired.dds.RData") # dds/ddsFlt
	mj.cpm.diff=data.frame(p2.cpm.diff=log2(ddsCpms[,"68"]/ddsCpms[,"66"]), p3.cpm.diff=log2(ddsCpms[,"74"]/ddsCpms[,"72"]), p5.cpm.diff=log2(ddsCpms[,"77"]/ddsCpms[,"75"]), p8.cpm.diff=log2(ddsCpms[,"65"]/ddsCpms[,"73"]))
	mj.cpm.diff$ensg<-rownames(mj.cpm.diff)
	ddsKeep = rowSums(ddsCpms >1) >= nrow(samples)/2
	mj.cpm.diff.flt=mj.cpm.diff[ddsKeep,]

	fields <- c("ensembl_gene_id", "hgnc_symbol","description", "gene_biotype", "chromosome_name","start_position","end_position","strand")
	dummy.ensg=getBM(attributes = fields, filters = "ensembl_gene_id", values = rownames(mj.cpm.diff.flt) , mart = ensembl)
	dummy.ensg$chromosome_name[dummy.ensg$chromosome_name=='MT']<-'M' # replace 'MT' with 'M'
	dummy.ensg$chromosome_name<-paste("chr",dummy.ensg$chromosome_name,sep='') # add 'chr' prefix
	dummy.range=IRanges( start=dummy.ensg$start_position,end=dummy.ensg$end_position,names=dummy.ensg$ensembl_gene_id)
	dummy.gr=GRanges(seqnames<-Rle(dummy.ensg$chromosome_name) , ranges<-dummy.range, strand<-Rle(dummy.ensg$strand), mcols<-DataFrame(data.frame(ensg=dummy.ensg$ensembl_gene_id, hgnc=dummy.ensg$hgnc_symbol, biotype=dummy.ensg$gene_biotype)))

	# DMR Tss2000 from Pair2
	load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair2.diffMeth.Tss.RData") # old version
	rm(myDiff1,myDiff2,myDiffp1,myDiffp2,myDiffp3)
	#head(myDiff3) # Tss2000
	p2.tss2k.gr=as(myDiff3,"GRanges")

	# find overlap between tss2000 and deg
	#findOverlaps(p2.tss2k.gr, dummy.gr, minoverlap=99L)
	#Hits of length 113181
	#queryLength: 159715
	#subjectLength: 13977

	dummy.ovl=as.matrix(findOverlaps(p2.tss2k.gr, dummy.gr, minoverlap=99L))
	p2.tss2k.gr.flt<-p2.tss2k.gr[dummy.ovl[,"queryHits"]]
	mcols(p2.tss2k.gr.flt)<-DataFrame(mcols(p2.tss2k.gr.flt), mcols(dummy.gr[dummy.ovl[,"subjectHits"]])) # attach gene info (ensg, desc, biotype) to p2.tss2k.gr.flt
	mcols(p2.tss2k.gr.flt)<-DataFrame(mcols(p2.tss2k.gr.flt), p2.cpm.diff=mj.cpm.diff.flt[mcols(p2.tss2k.gr.flt)$ensg,]$p2.cpm.diff) # attach log2(CPM(SGA)/CPM(AGA)) 
	# plot
	par(mfrow=c(2,2))
	plot(as.matrix(mcols(p2.tss2k.gr.flt)[,c("meth.diff","p2.cpm.diff")]), col="#00000020", pch=20, cex=0.5, main="Pair2", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(as.matrix(mcols(p3.tss2k.gr.flt)[,c("meth.diff","p3.cpm.diff")]), col="#00000020", pch=20, cex=0.5, main="Pair3", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(as.matrix(mcols(p5.tss2k.gr.flt)[,c("meth.diff","p5.cpm.diff")]), col="#00000020", pch=20, cex=0.5, main="Pair5", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(as.matrix(mcols(p8.tss2k.gr.flt)[,c("meth.diff","p8.cpm.diff")]), col="#00000020", pch=20, cex=0.5, main="Pair8", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	# correlation
	cor(as.matrix(mcols(p2.tss2k.gr.flt)[,c("meth.diff","p2.cpm.diff")]), method='spearman')
	cor(as.matrix(mcols(p3.tss2k.gr.flt)[,c("meth.diff","p3.cpm.diff")]), method='spearman')
	cor(as.matrix(mcols(p5.tss2k.gr.flt)[,c("meth.diff","p5.cpm.diff")]), method='spearman')
	cor(as.matrix(mcols(p8.tss2k.gr.flt)[,c("meth.diff","p8.cpm.diff")]), method='spearman')

	# DMR Tss2000 from Pair2
	load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair2.annotate.Tss.regions3.RData") # d1:Tss500, d2:Tss1000, d3:Tss2000, e1:Tss500(q<0.05) ...
	rm(d1,d2,e1,e2,e3)
	p2.tss2k=merge(d3, with(mj.cpm.diff, data.frame(ensg, p2.cpm.diff)), by.x="Ensembl.Gene.ID", by.y="ensg") # attach log2(CPM(SGA)/CPM(AGA))

	# DMR Tss2000 from Pair3
	load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair3.annotate.Tss.regions3.RData")
	rm(d1,d2,e1,e2,e3)
	p3.tss2k=merge(d3, with(mj.cpm.diff, data.frame(ensg, p3.cpm.diff)), by.x="Ensembl.Gene.ID", by.y="ensg") # attach log2(CPM(SGA)/CPM(AGA))

	# DMR Tss2000 from Pair5
	load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair5.annotate.Tss.regions3.RData")
	rm(d1,d2,e1,e2,e3)
	p5.tss2k=merge(d3, with(mj.cpm.diff, data.frame(ensg, p5.cpm.diff)), by.x="Ensembl.Gene.ID", by.y="ensg") # attach log2(CPM(SGA)/CPM(AGA))

	# DMR Tss2000 from Pair8
	load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair8.annotate.Tss.regions3.RData")
	rm(d1,d2,e1,e2,e3)
	p8.tss2k=merge(d3, with(mj.cpm.diff, data.frame(ensg, p8.cpm.diff)), by.x="Ensembl.Gene.ID", by.y="ensg") # attach log2(CPM(SGA)/CPM(AGA))

	#################################
	# Plot diff(cpm) vs. diff(%met) #
	#################################
	par(mfrow=c(2,2))
	plot(p2.tss2k[,c("meth.diff", "p2.cpm.diff")], col="#00000020", pch=20, cex=0.5, main="Pair2", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(p3.tss2k[,c("meth.diff", "p3.cpm.diff")], col="#00000020", pch=20, cex=0.5, main="Pair3", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(p5.tss2k[,c("meth.diff", "p5.cpm.diff")], col="#00000020", pch=20, cex=0.5, main="Pair3", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(p8.tss2k[,c("meth.diff", "p8.cpm.diff")], col="#00000020", pch=20, cex=0.5, main="Pair3", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")

	#################################################
	# Plot diff(cpm) vs. diff(%met) of p-value<-.05 #
	#################################################
	par(mfrow=c(2,2))
	plot(p2.tss2k[p2.tss2k$qvalue<=0.05,c("meth.diff", "p2.cpm.diff")], col="#00000020", pch=20, cex=0.5, main="Pair2 (%met q-value<=0.05)", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(p3.tss2k[p3.tss2k$qvalue<=0.05,c("meth.diff", "p3.cpm.diff")], col="#00000020", pch=20, cex=0.5, main="Pair3 (%met q-value<=0.05)", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(p5.tss2k[p5.tss2k$qvalue<=0.05,c("meth.diff", "p5.cpm.diff")], col="#00000020", pch=20, cex=0.5, main="Pair5 (%met q-value<=0.05)", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(p8.tss2k[p8.tss2k$qvalue<=0.05,c("meth.diff", "p8.cpm.diff")], col="#00000020", pch=20, cex=0.5, main="Pair8 (%met q-value<=0.05)", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")

	###############################################
	# Plot diff(cpm) vs. diff(%met) of top 1% DMR #
	###############################################
	top.p2=merge(pair2[pair2$percent==1,], p2.anno, by.x="Ensembl.Gene.ID", by.y="ensembl_gene_id")
	top.p3=merge(pair3[pair3$percent==1,], p3.anno, by.x="Ensembl.Gene.ID", by.y="ensembl_gene_id")
	top.p5=merge(pair5[pair5$percent==1,], p5.anno, by.x="Ensembl.Gene.ID", by.y="ensembl_gene_id")
	top.p8=merge(pair8[pair8$percent==1,], p8.anno, by.x="Ensembl.Gene.ID", by.y="ensembl_gene_id")

	top.p2=merge(top.p2, with(mj.cpm.diff, data.frame(ensg, p2.cpm.diff)), by.x="Ensembl.Gene.ID", by.y="ensg") # attach log2(CPM(SGA)/CPM(AGA))
	top.p3=merge(top.p3, with(mj.cpm.diff, data.frame(ensg, p3.cpm.diff)), by.x="Ensembl.Gene.ID", by.y="ensg")
	top.p5=merge(top.p5, with(mj.cpm.diff, data.frame(ensg, p5.cpm.diff)), by.x="Ensembl.Gene.ID", by.y="ensg")
	top.p8=merge(top.p8, with(mj.cpm.diff, data.frame(ensg, p8.cpm.diff)), by.x="Ensembl.Gene.ID", by.y="ensg")

	par(mfrow=c(2,2))
	plot(top.p2[top.p2$region=="Tss2000",c("meth.diff", "p2.cpm.diff")], col="#00000020", pch=19,  main="Pair2 (top 1% Tss2000 promoters)", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(top.p3[top.p3$region=="Tss2000",c("meth.diff", "p3.cpm.diff")], col="#00000020", pch=19,  main="Pair3 (top 1% Tss2000 promoters)", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(top.p5[top.p5$region=="Tss2000",c("meth.diff", "p5.cpm.diff")], col="#00000020", pch=19,  main="Pair5 (top 1% Tss2000 promoters)", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(top.p8[top.p8$region=="Tss2000",c("meth.diff", "p8.cpm.diff")], col="#00000020", pch=19,  main="Pair8 (top 1% Tss2000 promoters)", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")

	top.p2 <- with(top.p2, GRanges(chr, IRanges(start, end), strand, qvalue, meth.diff, met.AGA=AGA2F, met.SGA=SGA2F, region, enst=feature.name, ensg=Ensembl.Gene.ID, hgnc_symbol, description, gene_biotype))
	top.p3 <- with(top.p3, GRanges(chr, IRanges(start, end), strand, qvalue, meth.diff, met.AGA=AGA3F, met.SGA=SGA3F, region, enst=feature.name, ensg=Ensembl.Gene.ID, hgnc_symbol, description, gene_biotype))
	top.p5 <- with(top.p5, GRanges(chr, IRanges(start, end), strand, qvalue, meth.diff, met.AGA=AGA5M, met.SGA=SGA5M, region, enst=feature.name, ensg=Ensembl.Gene.ID, hgnc_symbol, description, gene_biotype))
	top.p8 <- with(top.p8, GRanges(chr, IRanges(start, end), strand, qvalue, meth.diff, met.AGA=AGA8M, met.SGA=SGA8M, region, enst=feature.name, ensg=Ensembl.Gene.ID, hgnc_symbol, description, gene_biotype))

	plot(as.matrix(mcols(top.p2)[,c("meth.diff", "p2.cpm.diff")]), col="#00000020", pch=19,  main="Pair2 (top 1% promoters)", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(as.matrix(mcols(top.p3)[,c("meth.diff", "p3.cpm.diff")]), col="#00000020", pch=19,  main="Pair3 (top 1% promoters)", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(as.matrix(mcols(top.p5)[,c("meth.diff", "p5.cpm.diff")]), col="#00000020", pch=19,  main="Pair2 (top 1% promoters)", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	plot(as.matrix(mcols(top.p8)[,c("meth.diff", "p8.cpm.diff")]), col="#00000020", pch=19,  main="Pair2 (top 1% promoters)", xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")

############################################
# Plot diff(cpm) vs. diff(%met) of DEG(MJ) # 
############################################
	# 1. DEG & DMRs
	par(mfrow=c(2,2))
	num=length(unique(p2.tss2k[p2.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,]$Ensembl.Gene.ID)) # NO. of DEGs having diff.met
	plot(p2.tss2k[p2.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,c("meth.diff", "p2.cpm.diff")], col="#00000020", pch=19,  main=paste0("Pair2 (MJ ", num, " DEG)"), xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	num=length(unique(p3.tss2k[p3.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,]$Ensembl.Gene.ID)) # NO. of DEGs having diff.met
	plot(p3.tss2k[p3.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,c("meth.diff", "p3.cpm.diff")], col="#00000020", pch=19,  main=paste0("Pair3 (MJ ", num, " DEG)"), xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	num=length(unique(p5.tss2k[p5.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,]$Ensembl.Gene.ID)) # NO. of DEGs having diff.met
	plot(p5.tss2k[p5.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,c("meth.diff", "p5.cpm.diff")], col="#00000020", pch=19,  main=paste0("Pair5 (MJ ", num, " DEG)"), xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	num=length(unique(p8.tss2k[p8.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,]$Ensembl.Gene.ID)) # NO. of DEGs having diff.met
	plot(p8.tss2k[p8.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,c("meth.diff", "p8.cpm.diff")], col="#00000020", pch=19,  main=paste0("Pair8 (MJ ", num, " DEG)"), xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")

	# 2. DEG & sigDMRs
	par(mfrow=c(2,2))
	num=length(unique(p2.tss2k[p2.tss2k$qvalue<=0.05 & p2.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,]$Ensembl.Gene.ID)) # NO. of DEGs having diff.met
	plot(p2.tss2k[p2.tss2k$qvalue<=0.05 & p2.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,c("meth.diff", "p2.cpm.diff")], col="#00000020", pch=19,  main=paste0("Pair2 (%met q-value<=0.05 & MJ ", num, " DEG)"), xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	num=length(unique(p3.tss2k[p3.tss2k$qvalue<=0.05 & p3.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,]$Ensembl.Gene.ID)) # NO. of DEGs having diff.met
	plot(p3.tss2k[p3.tss2k$qvalue<=0.05 & p3.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,c("meth.diff", "p3.cpm.diff")], col="#00000020", pch=19,  main=paste0("Pair3 (%met q-value<=0.05 & MJ ", num, " DEG)"), xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	num=length(unique(p5.tss2k[p5.tss2k$qvalue<=0.05 & p5.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,]$Ensembl.Gene.ID)) # NO. of DEGs having diff.met
	plot(p5.tss2k[p5.tss2k$qvalue<=0.05 & p5.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,c("meth.diff", "p5.cpm.diff")], col="#00000020", pch=19,  main=paste0("Pair5 (%met q-value<=0.05 & MJ ", num, " DEG)"), xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")
	num=length(unique(p8.tss2k[p8.tss2k$qvalue<=0.05 & p8.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,]$Ensembl.Gene.ID)) # NO. of DEGs having diff.met
	plot(p8.tss2k[p8.tss2k$qvalue<=0.05 & p8.tss2k$Ensembl.Gene.ID %in% union.top.deg$ensg,c("meth.diff", "p8.cpm.diff")], col="#00000020", pch=19,  main=paste0("Pair8 (%met q-value<=0.05 & MJ ", num, " DEG)"), xlab="% methylation(SGA-AGA)", ylab="log2(CPM(SGA)/CPM(AGA))")
	abline(h=c(0),v=c(0), col="blue")

	##########################
	# check bio-type of DMRs #
	##########################
	# 1. bio-type of top DMRs (top.p2, top.p3, top.p5 and top.p8)
	dummy.p2=with(top.p2[top.p2$region=="Tss2000",], data.frame(ensg=Ensembl.Gene.ID, biotype=gene_biotype))
	dummy.p3=with(top.p3[top.p3$region=="Tss2000",], data.frame(ensg=Ensembl.Gene.ID, biotype=gene_biotype))
	dummy.p5=with(top.p5[top.p5$region=="Tss2000",], data.frame(ensg=Ensembl.Gene.ID, biotype=gene_biotype))
	dummy.p8=with(top.p8[top.p8$region=="Tss2000",], data.frame(ensg=Ensembl.Gene.ID, biotype=gene_biotype))

	get_biotype<-function() {
		dummy.p2=dummy.p2[!duplicated(dummy.p2),] # remove duplicated entries
		dummy.p3=dummy.p3[!duplicated(dummy.p3),]
		dummy.p5=dummy.p5[!duplicated(dummy.p5),]
		dummy.p8=dummy.p8[!duplicated(dummy.p8),]

		with(dummy.p2,table(biotype))
		with(dummy.p3,table(biotype))
		with(dummy.p5,table(biotype))
		with(dummy.p8,table(biotype))

		#mymerge = function(x,y) merge.data.frame(x,y,all=TRUE)
		#Reduce(mymerge,list(with(dummy.p2,table(biotype)),with(dummy.p3,table(biotype)),with(dummy.p5,table(biotype)),with(dummy.p8,table(biotype))))

		dummy.p=rbind(dummy.p2,dummy.p3,dummy.p5,dummy.p8)
		dummy.p=dummy.p[!duplicated(dummy.p[,c('ensg','biotype')]),] # remove duplicated entries
		length(unique(dummy.p$ensg)) # NO. of gene
		with(dummy.p,table(biotype)) # Count number of genes by gene_biotype
		dummy.p.freq=as.data.frame(with(dummy.p,table(biotype))) # Count number of genes by gene_biotype as a matrix
		dummy.p.freq[order(dummy.p.freq$Freq, decreasing=T),]
	}

	get_biotype()

	# bio-type of significant DMRs (q-value<=0.05)
	fields <- c("ensembl_gene_id", "hgnc_symbol","gene_biotype")
	p2.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = p2.tss2k[p2.tss2k$qvalue<=0.05,]$Ensembl.Gene.ID, mart = ensembl)
	p3.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = p3.tss2k[p2.tss2k$qvalue<=0.05,]$Ensembl.Gene.ID, mart = ensembl)
	p5.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = p5.tss2k[p2.tss2k$qvalue<=0.05,]$Ensembl.Gene.ID, mart = ensembl)
	p8.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = p8.tss2k[p2.tss2k$qvalue<=0.05,]$Ensembl.Gene.ID, mart = ensembl)
	dummy.p2=with(merge(p2.tss2k[p2.tss2k$qvalue<=0.05,], p2.anno, by.x="Ensembl.Gene.ID", by.y="ensembl_gene_id"), data.frame(ensg=Ensembl.Gene.ID, biotype=gene_biotype))
	dummy.p3=with(merge(p3.tss2k[p3.tss2k$qvalue<=0.05,], p3.anno, by.x="Ensembl.Gene.ID", by.y="ensembl_gene_id"), data.frame(ensg=Ensembl.Gene.ID, biotype=gene_biotype))
	dummy.p5=with(merge(p5.tss2k[p5.tss2k$qvalue<=0.05,], p5.anno, by.x="Ensembl.Gene.ID", by.y="ensembl_gene_id"), data.frame(ensg=Ensembl.Gene.ID, biotype=gene_biotype))
	dummy.p8=with(merge(p8.tss2k[p8.tss2k$qvalue<=0.05,], p8.anno, by.x="Ensembl.Gene.ID", by.y="ensembl_gene_id"), data.frame(ensg=Ensembl.Gene.ID, biotype=gene_biotype))

	# bio-type of any DMRs
	fields <- c("ensembl_gene_id", "hgnc_symbol","gene_biotype")
	p2.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = p2.tss2k$Ensembl.Gene.ID, mart = ensembl)
	p3.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = p3.tss2k$Ensembl.Gene.ID, mart = ensembl)
	p5.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = p5.tss2k$Ensembl.Gene.ID, mart = ensembl)
	p8.anno=getBM(attributes = fields, filters = "ensembl_gene_id", values = p8.tss2k$Ensembl.Gene.ID, mart = ensembl)
	dummy.p2=with(merge(p2.tss2k, p2.anno, by.x="Ensembl.Gene.ID", by.y="ensembl_gene_id"), data.frame(ensg=Ensembl.Gene.ID, biotype=gene_biotype))
	dummy.p3=with(merge(p3.tss2k, p3.anno, by.x="Ensembl.Gene.ID", by.y="ensembl_gene_id"), data.frame(ensg=Ensembl.Gene.ID, biotype=gene_biotype))
	dummy.p5=with(merge(p5.tss2k, p5.anno, by.x="Ensembl.Gene.ID", by.y="ensembl_gene_id"), data.frame(ensg=Ensembl.Gene.ID, biotype=gene_biotype))
	dummy.p8=with(merge(p8.tss2k, p8.anno, by.x="Ensembl.Gene.ID", by.y="ensembl_gene_id"), data.frame(ensg=Ensembl.Gene.ID, biotype=gene_biotype))
}


############
## DAVID  ##
############
if(FALSE){
	# DAVID result (default option) using union.gene
	david.fa.chart.file=paste0(DavidHome,"/deg.unified.david.fa.chart.txt")
	david.fa.cluster.file=paste0(DavidHome,"/deg.unified.david.fa.cluster.txt")
	david.fa.table.file=paste0(DavidHome,"/deg.unified.david.fa.table.txt")

	david.fa.chart=read.table(file=david.fa.chart.file, header=TRUE, sep="\t")

	#GOTERM_BP_FAT, GOTERM_CC_FAT, GOTERM_MF_FAT
	dummy.david=david.fa.chart[david.fa.chart$Category=="GOTERM_BP_FAT",]
	write.csv(dummy.david[order(dummy.david$Count, decreasing=T), c("Category", "Term", "Count", "Genes", "PValue", "Bonferroni", "Benjamini", "FDR")], file=paste0(DavidHome,"/deg.unified.david.go.bp.fat.csv"))

	dummy.david=david.fa.chart[david.fa.chart$Category=="GOTERM_MF_FAT",]
	write.csv(dummy.david[order(dummy.david$Count, decreasing=T), c("Category", "Term", "Count", "Genes", "PValue", "Bonferroni", "Benjamini", "FDR")], file=paste0(DavidHome,"/deg.unified.david.go.mf.fat.csv"))

	dummy.david=david.fa.chart[david.fa.chart$Category=="GOTERM_CC_FAT",]
	write.csv(dummy.david[order(dummy.david$Count, decreasing=T), c("Category", "Term", "Count", "Genes", "PValue", "Bonferroni", "Benjamini", "FDR")], file=paste0(DavidHome,"/deg.unified.david.go.cc.fat.csv"))

	# LEVEL5
	dummy.david=david.fa.chart[david.fa.chart$Category=="GOTERM_BP_5",]
	write.csv(dummy.david[order(dummy.david$Count, decreasing=T), c("Category", "Term", "Count", "Genes", "PValue", "Bonferroni", "Benjamini", "FDR")], file=paste0(DavidHome,"/deg.unified.david.go.bp.level5.csv"))

	dummy.david=david.fa.chart[david.fa.chart$Category=="GOTERM_MF_5",]
	write.csv(dummy.david[order(dummy.david$Count, decreasing=T), c("Category", "Term", "Count", "Genes", "PValue", "Bonferroni", "Benjamini", "FDR")], file=paste0(DavidHome,"/deg.unified.david.go.mf.level5.csv"))

	dummy.david=david.fa.chart[david.fa.chart$Category=="GOTERM_CC_5",]
	write.csv(dummy.david[order(dummy.david$Count, decreasing=T), c("Category", "Term", "Count", "Genes", "PValue", "Bonferroni", "Benjamini", "FDR")], file=paste0(DavidHome,"/deg.unified.david.go.cc.level5.csv"))

	## OMIM_DISEASE
	dummy.david=david.fa.chart[david.fa.chart$Category=="OMIM_DISEASE",]
	write.csv(dummy.david[order(dummy.david$Count, decreasing=T), c("Category", "Term", "Count", "Genes", "PValue", "Bonferroni", "Benjamini", "FDR")], file=paste0(DavidHome,"/deg.unified.david.omim.disease.csv"))

	## GENETIC_ASSOCIATION_DB_DISEASE
	dummy.david=david.fa.chart[david.fa.chart$Category=="GENETIC_ASSOCIATION_DB_DISEASE",]
	write.csv(dummy.david[order(dummy.david$Count, decreasing=T), c("Category", "Term", "Count", "Genes", "PValue", "Bonferroni", "Benjamini", "FDR")], file=paste0(DavidHome,"/deg.unified.david.gene.asso.db.disease.csv"))

	## KEGG_PATHWAY
	dummy.david=david.fa.chart[david.fa.chart$Category=="KEGG_PATHWAY",]
	write.csv(dummy.david[order(dummy.david$Count, decreasing=T), c("Category", "Term", "Count", "Genes", "PValue", "Bonferroni", "Benjamini", "FDR")], file=paste0(DavidHome,"/deg.unified.kegg.pathway.csv"))

	## SP_PIR_KEYWORDS
	dummy.david=david.fa.chart[david.fa.chart$Category=="SP_PIR_KEYWORDS",]
	write.csv(dummy.david[order(dummy.david$Count, decreasing=T), c("Category", "Term", "Count", "Genes", "PValue", "Bonferroni", "Benjamini", "FDR")], file=paste0(DavidHome,"/deg.unified.pir.keywords.csv"))
}

dev.off()
cat("All is done", "\n")
