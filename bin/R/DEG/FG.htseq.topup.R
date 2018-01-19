#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

myCaller="FG.topup"
source ("~/Pipelines/config/DEG.R") # load config

cat("Final top-up for FG HTSEQ...\n")
targets<-c("78","97","142","155")

# SLX-8551.SLX-8552
samples<-samples[samples$SampleName %in% targets,]
counts = readDGE(samples$HtseqFile, header=FALSE, sep="")$counts

# SLX-9170.SLX-9171 (final mRNA top-up with total RNA assay)
samples.topup= read.csv(file="~/scratch/results/RNA-Seq/SGA.AGA/Meta/meta.topup.csv", stringsAsFactors=FALSE)
counts.topup = readDGE(samples.topup$HtseqFile, header=FALSE, sep="")$counts

counts.new<- counts + counts.topup
colnames(counts.new)<-targets

#ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = samples, design= formula(~ Condition)) # isa 'DESeqDataSet'

for(i in colnames(counts.new)){ 
	dummy<-as.data.frame(counts.new[,i])
	colnames(dummy)<-i
	write.table(dummy, file=paste0("~/scratch/results/RNA-Seq/Topup/",i,".igenome.HTSeq.count.txt"), quote=FALSE, col.names = FALSE)
}

dev.off()
cat("All is done\n")
