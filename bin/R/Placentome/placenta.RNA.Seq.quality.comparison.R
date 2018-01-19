#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

library(data.table) 
library(scales) # for 'lables=comma'
source("~/Pipelines/config/graphic.R")

mySource="Placentome" # "FG"
myCohort="POPS" # "POPS", "CLT", "PET", "SGA"
myProject=paste(mySource,myCohort,sep=".") # e.g. Placentome.CTL
top.dir=file.path("~/results/RNA-Seq",mySource) #top.dir="~/results/RNA-Seq/Placentome"

########################
# pdf output filename ##
########################
my.file.name<- file.path(top.dir, paste0("ReadCounts/placenta.rna.seq.read.count.comparision.ggplot"))
pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title="Placenta RNA-Seq Comparision") # A4 size

#############################
## Lib Size vs. NO. Sample ##
#############################
dummy<-fread(file.path(top.dir, "Meta/placenta.rna.seq.comparision.txt"), header=T)
p<-ggplot2::ggplot(dummy, aes(x=NO.Sample,y=Lib.Size,colour=paste0(Lib.Type,Read.Len),label=paste(Author, Year, sep=", "))) + 
	geom_point(size=8) + 
	geom_text(colour="black") + 
	xlab("Number of Samples") + ylab("Library Size (per million per sample)") + 
	ggtitle("Placenta Transcriptome Studies")
print(p)

p<-ggplot2::ggplot(dummy, aes(x=NO.Sample,y=ifelse(Lib.Type=="PE", Lib.Size*2*Read.Len, Lib.Size*Read.Len),colour=Lib.Type,label=paste(Author, Year, sep=", "))) +
	geom_point(size=8) + 
	geom_text(colour="black") +
	xlab("Number of Samples") + ylab("Avg. Library Size (per million bp)") + 
	ggtitle("Placenta Transcriptome Studies")
print(p)

dev.off()
