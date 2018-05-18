#
# TO DO:
# make this as a standardised reporting template 
# using R-makrdown
#
library(data.table)
library(ggExtra) # to use ggMarginal
source("~/Pipelines/config/graphic.R")
time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p')

SLX="SLX-SML1"
SPECIES='Homo_sapiens' # 'Homo_sapiens' or 'Mus_musculus' or 'PhiX'
VERSION=paste(SPECIES,"v4",sep=".")
PROJECT=paste(SLX,VERSION,sep=".")
resultDir="/home/ssg29/results"
BARCODES=dir(file.path(resultDir,PROJECT,"TopHat"))

cntFile=file.path(resultDir,PROJECT,paste0(SLX,".fq.read.base.cnt.txt"))
if(!file.exists(cntFile)){stop(paste(cntFile,"not found"))}

foo=list()
for(i in BARCODES){
    cat(i,"\n")
    # No. of raw reads
    readIn1=sum(as.integer(simplify2array(strsplit(system(paste(paste("grep ",i, cntFile,sep=" "), "| grep r_1"),intern=T), " "))[2,])) # a list of char v

    myFile=file.path(resultDir,PROJECT,"TopHat",i,"align_summary.txt")
    if(!file.exists(myFile)){stop(paste(myFile,"not found"))}
    # No. of reads (r1+r2) for TopHat
    readIn2=as.numeric(system(paste0("grep Input ",myFile, "| cut -d: -f2 | awk '{sum+=$1}END{print sum}'"),intern=T))
    # No. of 1st mapped reads (r1+r2) by TopHat
    readOut1=as.numeric(system(paste0("grep Mapped ",myFile, "| cut -d: -f2 |  cut -d'(' -f1 | awk '{sum+=$1}END{print sum}'"),intern=T))
    # No. of 2nd mapped reads (r1+r2) by TopHat
    myFile=file.path(resultDir,PROJECT,"TopHat",i,"unmapped","align_summary.txt")
    if(!file.exists(myFile)){stop(paste(myFile,"not found"))}
    readOut2=as.numeric(system(paste0("grep Mapped ",myFile, "| cut -d: -f2 |  cut -d'(' -f1 | awk '{sum+=$1}END{print sum}'"),intern=T))

    readOut=readOut1+readOut2
    foo[[i]]=c(readIn1,readIn2,readOut1,readOut2,readOut) # a vector
}
bar=data.frame(t(data.frame(foo))); colnames(bar) <- c("readIn1","readIn2","readOut1","readOut2","readOut")
p <- ggplot(bar, aes(readIn1/10^6, readOut/readIn2*100)) + 
    geom_point(size=6,alpha=.5) + 
    labs(x="No. of Reads (millions)", y="Mapping Efficiency (%)") + 
    ggtitle(paste0(SLX, "(n=",nrow(bar),")")) + 
    theme_Publication()

my.filename <- file.path(resultDir,PROJECT,paste0(SLX,".read.count.mapping.efficiency"))
my.filename <- file.path(resultDir,paste0("RNA-Seq/Cholestasis.2017/DESeq2/ALL/",SLX,".read.count.mapping.efficiency"))
if(capabilities()[["tiff"]]){
    tiff(filename=paste0(my.filename,".tiff"),width=9,height=9,units="in",res=300, compression = 'lzw')
}else{
    jpeg(filename=paste0(my.filename,".jpeg"),width=9,height=9,units="in",res=300)
}
ggMarginal(p, type = "boxplot", size=8)
dev.off()
