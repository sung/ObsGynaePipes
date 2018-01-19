#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

library(data.table)
library(dupRadar) # R >= 3.2

##########
# Config #
##########
time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p') #2015-01-07_01_30PM
cols <- colorRampPalette(c("black","blue","green","yellow","red"))

SLX="SLX-9342" # e.g. SLX-8080 
#PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
#PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="D701_D502" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

#bamDuprm <- "/home/ssg29/results/SLX-9342.Homo_sapiens.PE75.v1/TopHat/D701_D502/merged_accepted_hits.dedup.bam"
bamDuprm <- "/home/ssg29/results/SLX-9342.Homo_sapiens.PE75.v1/TopHat/D701_D502/merged_accepted_hits.uniq.dedup.bam"
gtf <- "/home/ssg29/data/genome/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.82.gtf"  # the gene model
stranded <- 1       # '0' (unstranded), '1' (stranded) and '2' (reversely stranded)
paired   <- TRUE # is the library paired end?
threads  <- 4 

dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads,verbose=T)

#file.name<-file.path("/home/ssg29/results/SLX-9342.Homo_sapiens.PE75.v1/dupRadar/D701_D502","merged_accepted_hits.dedup.bam.pdf")
file.name<-file.path("/home/ssg29/results/SLX-9342.Homo_sapiens.PE75.v1/dupRadar/D701_D502","merged_accepted_hits.uniq.dedup.bam.pdf")
pdf(file=file.name, width=12, height=8, title=paste0("dupRadar:",SLX, Barcode))

duprateExpDensPlot(DupMat=dm)
duprateExpBoxplot(DupMat=dm)    
duprateExpPlot(DupMat=dm)
readcountExpBoxplot(DupMat=dm)
expressionHist(DupMat=dm)

dev.off()
