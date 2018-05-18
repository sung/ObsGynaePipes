#!/usr/bin/RsCRIPT --vanilla
# Sung Gong <sung@bio.cc>

TR_PREFIX='GRCh38' # GRCh37|GRCh38
source("~/Pipelines/bin/R/Placentome/local.R") # load 'dt.pops 'dt.read.cnt'
source("~/Pipelines/config/graphic.R")
mySource="Placentome" # "FG"
top.dir=file.path("~/results/RNA-Seq",mySource) #top.dir="~/results/RNA-Seq/Placentome"

############################
## 1. No Sample by Cohort ##
############################
library(scales)
p.samples<-ggplot(dt.pops, aes(x = factor(1), fill = Cohort)) + 
	geom_bar(width = 1) + 
	coord_polar(theta = "y") +
	scale_fill_manual(values=cbPalette) +
	#geom_text(data=dt.pops[,.N,Cohort], aes(y = N/3 + c(0, cumsum(N)[-length(N)]), label = percent(N/100)), size=5) +
	geom_text(data=dt.pops[,.N,Cohort], aes(y = N/3 + c(0, cumsum(N)[-length(N)]), label = N), size=5) +
	ggtitle(paste("Samples (n=", nrow(dt.pops), ")")) +
	theme_Publication() +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks = element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank()
		)

tiff(filename=file.path(top.dir,"Figures/samples.pie.tiff"), width=10, height=8, units="in", res=300, compression='lzw')
print(p.samples)
dev.off()

#########################
## 2. Read Count Data  ##
#########################
dt.pops.read.cnt<-merge(dt.pops, dt.read.cnt, by=c("Library","BarCode"))
dummy<-dt.pops.read.cnt[,list(Read=round(sum(as.numeric(Read))/10^9,2)),"Category"]
p.read<-ggplot(dummy, aes(x="",y=Read, fill=Category)) + 
	geom_bar(stat = "identity", width=1) + 
	coord_polar(theta = "y") + 
	scale_fill_grey() + 
	theme(axis.text.x=element_blank()) +
	geom_text(aes(y = Read/3 + c(0, cumsum(Read)[-length(Read)]), label = Read), size=5) +
	ggtitle(paste("Reads (", dummy[,sum(Read)], " Billion (10^9))")) +
	theme_Publication() +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks = element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank()
		)

tiff(filename=file.path(top.dir,"Figures/reads.pie.tiff"), width=10, height=8, units="in", res=300, compression='lzw')
print(p.read)
dev.off()

############################
## 3. No Read per cohort  ##
############################
dt.dummy=rbind(dt.pops.read.cnt[Category!="Unmapped",list(SampleName,Cohort,Category,Read)],

dt.foo<-rbind(dt.pops.read.cnt[,list(Read=sum(Read)),"SampleName,Cohort"][,`Category`:="Sequenced Read"],
			  dt.pops.read.cnt[Category!="Unmapped",list(Read=sum(Read)),"SampleName,Cohort"][,`Category`:="Total Mapped Read"]
			  )

p.read2<-ggplot(dt.foo, aes(Cohort,Read/10^6)) + 
	geom_boxplot(aes(colour=Category),width=.5,size=1) +
	scale_colour_manual(values=cbPalette) +
	theme_Publication()

tiff(filename=file.path(top.dir,"Figures/reads.per.sample.box.tiff"), width=12, height=8, units="in", res=300, compression='lzw')
print(p.read2)
dev.off()

