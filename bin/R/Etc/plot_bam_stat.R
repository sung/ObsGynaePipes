#!/usr/bin/Rscript --vanilla

# output filename 

result_top_dir="/whale-data/ssg29/Methyl-Seq/results/SGA_v3/Bismark/Alignment"

if(FALSE){
	plot_stat<-function(x) { #x should be a bardcode
		barcode<-x
		cat(barcode,"\n")

		result_dir=paste0(result_top_dir, "/", barcode, "/", barcode) #/whale-data/ssg29/Methyl-Seq/results/SGA_v2/Bismark/Alignment/A001.old/A001

		pdf ( file=paste0(result_dir,'.',format(Sys.time(), "%Y-%m-%d_%I%p"), ".pdf") )

		mapq <- read.table(paste0(result_dir,'.r_1_val_1.fq.gz_bismark_bt2_pe.bam.mapq.txt'))

		quantile(mapq$V1, seq(0, 1, by=.05))

		hist(mapq$V1, main=paste("MQ distribution for ",barcode), probability=FALSE, xlab="MQ", ylab="NO. of reads")

		boxplot(mapq$V1, names=c("MQ"), horizontal=T, main='Boxplot of MQ', xlab='MQ', las=1)

		write(quantile(mapq$V1, seq(0, 1, by=.05)), file=paste0(result_dir,".quantile.txt"))

		# taking too long within the PDF
		#n <- length(mapq$V1)
		#plot((1:n - 1)/(n - 1), sort(mapq$V1), type="l", main = "MQ Quantiles for A003", xlab = "Sample Fraction", ylab = "Sample Quantile")
	}

	barcodes<-c('A003')
	for(i in barcodes)(plot_stat(i))
}


barcode='A003'
result_dir=paste0(result_top_dir, "/", barcode, "/", barcode) #/whale-data/ssg29/Methyl-Seq/results/SGA_v2/Bismark/Alignment/A003/A003
pdf ( file=paste0(result_dir,'.mapq_mismatch',format(Sys.time(), "%Y-%m-%d_%I%p"), ".pdf") )

## 1. MQ distribution Before / After dedup 
if(FALSE){
# before dedup
input1=paste0(result_dir,'.r_1_val_1.fq.gz_bismark_bt2_pe.bam.mapq.txt')


	mapq <- read.table(input1)
	quantile(mapq$V1, seq(0, 1, by=.05))
	hist(mapq$V1, main="MQ distribution before de-dup", probability=FALSE, xlab="MQ", ylab="NO. of reads")
	boxplot(mapq$V1, names=c("MQ"), horizontal=T, main='Boxplot of MQ', xlab='MQ', las=1)
	write(quantile(mapq$V1, seq(0, 1, by=.05)), file=paste0(result_dir,".quantile.before.dedup.txt"))

# after dedup
input2=paste0(result_dir,'.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.bam.mapq.txt')

	mapq <- read.table(input2)
	quantile(mapq$V1, seq(0, 1, by=.05))
	hist(mapq$V1, main="MQ distribution after de-dup", probability=FALSE, xlab="MQ", ylab="NO. of reads")
	boxplot(mapq$V1, names=c("MQ"), horizontal=T, main='Boxplot of MQ', xlab='MQ', las=1)
	write(quantile(mapq$V1, seq(0, 1, by=.05)), file=paste0(result_dir,".quantile.after.dedup.txt"))
}

## 2. MQ and NO of [I/D/X]

# before dedup
input1=paste0(result_dir,'.r_1_val_1.fq.gz_bismark_bt2_pe.bam.mapq.mismatch.txt') 
mapq <- read.table(input1)
cor(mapq)
#plot(mapq$V1, mapq$V2, main="Scatterplot before de-dup", xlab="MAPQ", ylab="sum(I,D,X)")

zero.mapq<-mapq[mapq$V1<1,]
hist(zero.mapq$V2, main="Sum of [I|D|X]", probability=FALSE, xlab="sum(IDX)", ylab="NO. of reads MAPQ<1")

three.mapq<-mapq[mapq$V1<3,]
hist(three.mapq$V2, main="Sum of [I|D|X]", probability=FALSE, xlab="sum(IDX)", ylab="NO. of reads MAPQ<3")

five.mapq<-mapq[mapq$V1<5,]
hist(five.mapq$V2, main="Sum of [I|D|X]", probability=FALSE, xlab="sum(IDX)", ylab="NO. of reads MAPQ<5")

ten.mapq<-mapq[mapq$V1<10,]
hist(ten.mapq$V2, main="Sum of [I|D|X]", probability=FALSE, xlab="sum(IDX)", ylab="NO. of reads MAPQ<10")

dev.off()
