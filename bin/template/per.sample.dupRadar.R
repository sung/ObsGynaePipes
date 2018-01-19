#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

library(data.table)
library(dupRadar) # R >= 3.2

time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p') #2015-01-07_01_30PM
cols <- colorRampPalette(c("black","blue","green","yellow","red"))
##########
# Config #
##########
args <- commandArgs(TRUE)
SLX <- args[1] #SLX="SLX-9342" # e.g. SLX-8080 
Barcode <- args[2] #Barcode="D701_D502" # e.g. A001
bamDuprm <- args[3]
#bamDuprm <- "/home/ssg29/results/SLX-9342.Homo_sapiens.PE75.v1/TopHat/D701_D502/merged_accepted_hits.uniq.dedup.bam"
threads <- as.numeric(args[4]) #threads  <- 4 

gtf <- "/home/ssg29/data/genome/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.82.gtf"  # the gene model
stranded <- 1       # '0' (unstranded), '1' (stranded) and '2' (reversely stranded)
paired   <- TRUE # is the library paired end?

##########
## Main ##
##########
ProjectName=strsplit(bamDuprm,"/")[[1]][5] # SLX-9342.Homo_sapiens.PE75.v1 
BamfileName=strsplit(bamDuprm,"/")[[1]][8] # merged_accepted_hits.uniq.dedup.bam
file.name<-file.path("/home/ssg29/results",paste0(ProjectName,"/dupRadar/",Barcode,"/",BamfileName,".RData"))
dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads,verbose=T)
save(dm, file=file.name)

file.name<-file.path("/home/ssg29/results",paste0(ProjectName,"/dupRadar/",Barcode,"/",BamfileName,".pdf"))
pdf(file=file.name, width=12, height=8, title=paste("dupRadar:",SLX, Barcode, BamfileName))

cat("duprateExpDensPlot\n")
duprateExpDensPlot(DupMat=dm)
cat("duprateExpBoxplot\n")
duprateExpBoxplot(DupMat=dm)    
cat("duprateExpPlot\n")
duprateExpPlot(DupMat=dm)
cat("readcountExpBoxplot\n")
readcountExpBoxplot(DupMat=dm)
cat("expressionHist\n")
expressionHist(DupMat=dm)

dev.off()
