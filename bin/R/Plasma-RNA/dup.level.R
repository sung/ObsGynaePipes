library(data.table)
source("~/Pipelines/config/graphic.R")

foo<-fread("~/Pipelines/bin/R/Plasma-RNA/CSV/SLX-11368.Homo_sapiens.PE150.v1.flagstat.txt") # Ultra-deep PE150 (344M reads)
tar<-fread("~/Pipelines/bin/R/Plasma-RNA/CSV/SLX-9342.Homo_sapiens.PE75.v2.flagstat.txt") # Moderate PE75 (32M reads)
yar<-fread("~/Pipelines/bin/R/Plasma-RNA/CSV/SLX-10285.Homo_sapiens.SE125.D708_D503.50pct.sampled.flagstat.txt") # Placenta sample SE125 (106M reads)
#bar<-fread("~/Pipelines/bin/R/Plasma-RNA/CSV/SLX-9342.Homo_sapiens.PE150.v1.flagstat.txt") # Shallow PE150 (4M reads)
#var<-fread("~/Pipelines/bin/R/Plasma-RNA/CSV/SLX-9342.Homo_sapiens.SE50.v1.flagstat.txt") # Shallow SE50 (3M reads)

# Shallow PE150 (4M reads)
#echo -e "Type Sampled Total Secondary Duplicated Paired Properly_paired Singleton" > SLX-9342.Homo_sapiens.PE150.v1.flagstat.txt
#foo="01 1 2 3 4 5 6 7 8 9 10"
#for i in $foo;do FILE=/home/ssg29/results/SLX-9342.Homo_sapiens.PE150.v1/TopHat/D701_D508/accepted_hits.$i.dedup.flagstat.txt ; printf "Plasma PE150_4M $i "; awk '/total/{t=$1}/secondary/{s=$1}/duplicates/{d=$1}/paired in sequencing/{p=$1}/properly paired/{pp=$1}/singletons/{st=$1}END{print t,s,d,p,pp,st}' $FILE; done >> SLX-9342.Homo_sapiens.PE150.v1.flagstat.txt

# Moderate PE75 (32M reads)
#echo -e "Type Sampled Total Secondary Duplicated Paired Properly_paired Singleton" > SLX-9342.Homo_sapiens.PE75.v2.flagstat.txt
#foo="01 1 2 3 4 5 6 7 8 9 10"
#for i in $foo;do FILE=/home/ssg29/results/SLX-9342.Homo_sapiens.PE75.v2/TopHat/D701_D508/accepted_hits.$i.dedup.flagstat.txt ; printf "Plasma PE75_32M $i "; awk '/total/{t=$1}/secondary/{s=$1}/duplicates/{d=$1}/paired in sequencing/{p=$1}/properly paired/{pp=$1}/singletons/{st=$1}END{print t,s,d,p,pp,st}' $FILE; done >> SLX-9342.Homo_sapiens.PE75.v2.flagstat.txt

# Shallow SE50 (3M reads)
#echo -e "Type Sampled Total Secondary Duplicated Paired Properly_paired Singleton" > SLX-9342.Homo_sapiens.SE50.v1.flagstat.txt
#foo="01 1 2 3 4 5 6 7 8 9 10"
#for i in $foo;do FILE=/scratch/obsgynae/POPS/results/SLX-9342.Homo_sapiens.SE50.v1/TopHat/D701_D508/accepted_hits.$i.dedup.flagstat.txt ; printf "Plasma SE50_3M $i "; awk '/total/{t=$1}/secondary/{s=$1}/duplicates/{d=$1}/paired in sequencing/{p=$1}/properly paired/{pp=$1}/singletons/{st=$1}END{print t,s,d,p,pp,st}' $FILE; done >> SLX-9342.Homo_sapiens.SE50.v1.flagstat.txt

#########
## ALL ##
#########
#dt.dup<-rbind(foo,bar,tar,var,yar)
dt.dup<-rbind(foo,tar,yar)
p150.h<-ggplot(dt.dup[Type=="PE150_344M"], aes(Sampled, Duplicated/(Total-Secondary)*100)) + 
	geom_bar(stat="identity",fill="darkgrey") + 
	geom_point(size=5) + 
	geom_line(size=1.5) + 
	ggtitle("SLX-11368: PE150 (344M reads)") + 
	labs(x="Sampled Reads (%)", y="Duplicated Reads (%)") + 
	scale_x_continuous(breaks=bar$Sampled, limit=c(-5,105)) +
	scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) +
	theme_Publication()

p150.l<-ggplot(dt.dup[Type=="PE150_4M"], aes(Sampled, Duplicated/(Total-Secondary)*100)) + 
	geom_bar(stat="identity",fill="darkgrey") + 
	geom_point(size=5) + 
	geom_line(size=1.5) + 
	ggtitle("SLX-9342: PE150 (4M reads)") + 
	labs(x="Sampled Reads (%)", y="Duplicated Reads (%)") + 
	scale_x_continuous(breaks=bar$Sampled, limit=c(-5,105)) +
	scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) +
	theme_Publication()

p75<-ggplot(dt.dup[Type=="PE75_32M"], aes(Sampled, Duplicated/(Total-Secondary)*100)) + 
	geom_bar(stat="identity",fill="darkgrey") + 
	geom_point(size=5) + 
	geom_line(size=1.5) + 
	ggtitle("SLX-9342: PE75 (32M reads)") + 
	labs(x="Sampled Reads (%)", y="Duplicated Reads (%)") + 
	scale_x_continuous(breaks=bar$Sampled, limit=c(-5,105)) +
	scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) +
	theme_Publication()

s50<-ggplot(dt.dup[Type=="SE50_3M"], aes(Sampled, Duplicated/(Total-Secondary)*100)) + 
	geom_bar(stat="identity",fill="darkgrey") + 
	geom_point(size=5) + 
	geom_line(size=1.5) + 
	ggtitle("SLX-9342: SE50 (3M reads)") + 
	labs(x="Sampled Reads (%)", y="Duplicated Reads (%)") + 
	scale_x_continuous(breaks=bar$Sampled, limit=c(-5,105)) +
	scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) +
	theme_Publication()

p150<-ggplot(dt.dup[grepl("PE150",Type)], aes(Sampled, Duplicated/(Total-Secondary)*100)) + 
	geom_bar(stat="identity",fill="darkgrey") + 
	geom_point(size=5) + 
	geom_line(size=1.5) + 
	ggtitle("Duplication Ratio") + 
	labs(x="Sampled Reads (%)", y="Duplicated Reads (%)") + 
	scale_x_continuous(breaks=bar$Sampled, limit=c(-5,105)) +
	scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) +
	facet_wrap(~Type) +
	theme_Publication()

p75.s50<-ggplot(dt.dup[grepl("PE75",Type) | grepl("SE50",Type)], aes(Sampled, Duplicated/(Total-Secondary)*100)) + 
	geom_bar(stat="identity",fill="darkgrey") + 
	geom_point(size=5) + 
	geom_line(size=1.5) + 
	ggtitle("Duplication Ratio") + 
	labs(x="Sampled Reads (%)", y="Duplicated Reads (%)") + 
	scale_x_continuous(breaks=bar$Sampled, limit=c(-5,105)) +
	scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) +
	facet_wrap(~Type) +
	theme_Publication()
# x-axis: % sampled, y-axis: %
p.all.dup1<-ggplot(dt.dup, aes(Sampled, Duplicated/(Total-Secondary)*100)) + 
	geom_bar(stat="identity",fill="darkgrey") + 
	geom_point(size=5) + 
	geom_line(size=1.5) + 
	ggtitle("Duplication Ratio") + 
	labs(x="Sampled Reads (%)", y="Duplicated Reads (%)") + 
	scale_x_continuous(breaks=bar$Sampled, limit=c(-5,105)) +
	scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) +
	facet_wrap(~Type, nrow=2) +
	theme_Publication()
# x-axis: % sampled, y-axis: %
p.all.dup2<-ggplot(dt.dup, aes(Sampled, Duplicated/(Total-Secondary)*100,group=Type)) + 
	geom_point(aes(color=Type),size=5,alpha=.7) + 
	#geom_line(linetype=2,size=1.5,color='darkgrey') + 
	geom_line(aes(color=Type,linetype=Tissue),size=1.5,alpha=.6) + 
	ggtitle("Duplication Ratio") + 
	labs(x="Sampled Reads (%)", y="Duplicated Reads (%)") + 
	scale_x_continuous(breaks=bar$Sampled) +
	scale_y_continuous(breaks=seq(0,100,10),limits=c(20,100)) +
	theme_Publication()
# x-axis: % sampled, y-axis: no.
p.all.dup3<-ggplot(dt.dup, aes(Sampled, Duplicated, group=Type)) + 
	geom_point(aes(color=Type),size=5,alpha=.7) + 
	#geom_line(linetype=2,size=1.5,color='darkgrey') + 
	geom_line(aes(color=Type,linetype=Tissue),size=1.5,alpha=.6) + 
	ggtitle("Number of Duplicated Reads") + 
	labs(x="Sampled Reads (%)", y="Number") + 
	scale_x_continuous(breaks=bar$Sampled) +
	theme_Publication()
# x-axis: no read, y-axis: no.
p.all.dup4<-ggplot(dt.dup, aes(No_read*(Sampled/100)/10^6, Duplicated/10^6, group=Type)) + 
	geom_point(aes(color=Type),size=5,alpha=.7) + 
	#geom_line(linetype=2,size=1.5,color='darkgrey') + 
	geom_line(aes(color=Type,linetype=Tissue),size=1.5,alpha=.6) + 
	ggtitle("Number of Duplicated Reads") + 
	labs(x="No of Sampled Reads (million)", y="Number of Duplicated Reads (million)") + 
	#scale_x_continuous(breaks=bar$Sampled) +
	theme_Publication()

p.all.uniq1<-ggplot(dt.dup, aes(Sampled, (Total-Secondary-Duplicated)/(Total-Secondary)*100, group=Type)) + 
	geom_point(aes(color=Type),size=7,alpha=.7) + 
	geom_line(aes(color=Type,linetype=Tissue),size=1.5,alpha=.6) + 
	ggtitle("Mapped Unique Reads") + 
	labs(x="Sampled Reads (%)", y="Percentage of Uniq Reads (%)") + 
	scale_x_continuous(breaks=bar$Sampled) +
	scale_y_continuous(breaks=seq(0,100,10),limits=c(0,80)) +
	theme_Publication()

p.all.uniq2<-ggplot(dt.dup, aes(Sampled, Total-Secondary-Duplicated, group=Type)) + 
	geom_point(aes(color=Type),size=7,alpha=.7) + 
	#geom_line(aes(color=Type),linetype=2,size=1.5,alpha=.5) + 
	geom_line(aes(color=Type,linetype=Tissue),size=1.5,alpha=.6) + 
	ggtitle("Mapped Unique Reads") + 
	labs(x="Sampled Reads (%)", y="Number of Uniq Reads") + 
	scale_x_continuous(breaks=bar$Sampled) +
	theme_Publication()
# x-axis: No. of sampled read; y-axis: No. of uniq read
p.all.uniq3<-ggplot(dt.dup, aes(No_read*(Sampled/100)/10^6, (Total-Secondary-Duplicated)/10^6, group=Type)) + 
	geom_point(aes(color=Type),size=7,alpha=.7) + 
	#geom_line(aes(color=Type),linetype=2,size=1.5,alpha=.5) + 
	geom_line(aes(color=Type,linetype=Tissue),size=1.5,alpha=.6) + 
	ggtitle("Mapped Unique Reads") + 
	labs(x="No. of Sampled Reads (million)", y="Number of Uniq Reads (milliion") + 
	#scale_x_continuous(breaks=bar$Sampled) +
	theme_Publication()

##
## FeatureCount
##
if(FALSE){
# 1. Plasma PE150_344M: Ultra-deep PE150 (344M reads)
sampling.levels=list("01"=1,"1"=10,"3"=30,"5"=50,"7"=70,"10"=100)
#sampling.levels=list("01"=1,"1"=10,"2"=20,"3"=30,"4"=40,"5"=50,"7"=70,"9"=90,"10"=100)
dl.count.dup<-list()
dl.count.dedup<-list()
for(i in names(sampling.levels)){
	# run at the OBPC66?
	myRData=file.path("~/results/SLX-11368.Homo_sapiens.v1/dupRadar/NoIndex",paste("accepted_hits",i,"dedup.bam.RData",sep="."))
	if(!file.exists(myRData)){
		stop(paste(myRData, "not found\n"))
	}
	cat(myRData,"\n")
	load(myRData) # 
	dt.dm<-data.table(dm)
	dl.count.dup[["PE150_344M"]][[as.character(sampling.levels[[i]])]]=dt.dm[allCounts>=10,.N] # NO. of genes having at least 10 reads (including duplicates)
	dl.count.dedup[["PE150_344M"]][[as.character(sampling.levels[[i]])]]=dt.dm[filteredCounts>=10,.N] # NO. of genes having at least 10 reads (including duplicates)
}

# 2. Plasma PE75_32M: Moderate PE75 (32M reads)
#sampling.levels=list("01"=1,"1"=10,"2"=20,"3"=30,"4"=40,"5"=50,"6"=60,"7"=70,"8"=80,"9"=90,"10"=100)
sampling.levels=list("01"=1,"1"=10,"3"=30,"5"=50,"7"=70,"10"=100)
for(i in names(sampling.levels)){
	# run at the HPC
	myRData=file.path("~/results/SLX-9342.Homo_sapiens.PE75.v2/dupRadar/D701_D508",paste("accepted_hits",i,"dedup.bam.RData",sep="."))
	if(!file.exists(myRData)){
		stop(paste(myRData, "not found\n"))
	}
	cat(myRData,"\n")
	load(myRData) # 
	dt.dm<-data.table(dm)
	dl.count.dup[["PE75_32M"]][[as.character(sampling.levels[[i]])]]=dt.dm[allCounts>=10,.N] # NO. of genes having at least 10 reads (including duplicates)
	dl.count.dedup[["PE75_32M"]][[as.character(sampling.levels[[i]])]]=dt.dm[filteredCounts>=10,.N] # NO. of genes having at least 10 reads (including duplicates)
}

# 3. Placenta_SE125_106M: Placenta sample SE125 (106M reads)
sampling.levels=list("01"=1,"1"=10,"3"=30,"5"=50,"7"=70,"10"=100)
for(i in names(sampling.levels)){
	# run at the OBPC66
	myRData=file.path("~/results/SLX-10285.Homo_sapiens.v1/dupRadar/D708_D503",paste("accepted_hits",i,"dedup.bam.RData",sep="."))
	if(!file.exists(myRData)){
		stop(paste(myRData, "not found\n"))
	}
	cat(myRData,"\n")
	load(myRData) # 
	dt.dm<-data.table(dm)
	dl.count.dup[["SE125_106M"]][[as.character(sampling.levels[[i]])]]=dt.dm[allCounts>=10,.N] # NO. of genes having at least 10 reads (including duplicates)
	dl.count.dedup[["SE125_106M"]][[as.character(sampling.levels[[i]])]]=dt.dm[filteredCounts>=10,.N] # NO. of genes having at least 10 reads (including duplicates)
}

dt.count.dup<-rbind(
	data.table(`Tissue`="Plasma",`Class`="Duplicate-included",`Type`=melt(dl.count.dup[1])$L1, `Sampled`=as.numeric(rownames(melt(dl.count.dup[1]))), `Count`=melt(dl.count.dup[1])$value),
	data.table(`Tissue`="Plasma",`Class`="Duplicate-included",`Type`=melt(dl.count.dup[2])$L1, `Sampled`=as.numeric(rownames(melt(dl.count.dup[2]))), `Count`=melt(dl.count.dup[2])$value),
	data.table(`Tissue`="Placenta",`Class`="Duplicate-included",`Type`=melt(dl.count.dup[3])$L1, `Sampled`=as.numeric(rownames(melt(dl.count.dup[3]))), `Count`=melt(dl.count.dup[3])$value)
	)
dt.count.dedup<-rbind(
	data.table(`Tissue`="Plasma",`Class`="De-duplicated",`Type`=melt(dl.count.dedup[1])$L1, `Sampled`=as.numeric(rownames(melt(dl.count.dedup[1]))), `Count`=melt(dl.count.dedup[1])$value),
	data.table(`Tissue`="Plasma",`Class`="De-duplicated",`Type`=melt(dl.count.dedup[2])$L1, `Sampled`=as.numeric(rownames(melt(dl.count.dedup[2]))), `Count`=melt(dl.count.dedup[2])$value),
	data.table(`Tissue`="Placenta",`Class`="De-duplicated",`Type`=melt(dl.count.dedup[3])$L1, `Sampled`=as.numeric(rownames(melt(dl.count.dedup[3]))), `Count`=melt(dl.count.dedup[3])$value)
	)

p.cnt.gene<-ggplot(rbind(dt.count.dup,dt.count.dedup), aes(Sampled, Count, group=Type)) + 
	geom_point(aes(color=Type),size=7,alpha=.7) + 
	#geom_line(aes(color=Type),linetype=2,size=1.5,alpha=.5) + 
	geom_line(aes(color=Type,linetype=Tissue),size=1.5,alpha=.6) + 
	ggtitle("NO. of Genes") + 
	labs(x="Sampled Reads (%)", y="Number of Genes (at least 10 reads)") + 
	scale_x_continuous(breaks=bar$Sampled) +
	facet_wrap(~Class) +
	theme_Publication()
} # end of if FALSE

##
## Save to FILE
##
my.file.name<- "~/results/RNA-Seq/Plasma.2017/Cluster/low.PSG7.sample.duplication.ratio"
pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11, height=8, title="Low PSG7 Sample: Plasma and Tissue") # A4 size
#print(p150.h)
#print(p150.l)
#print(p150)
#print(p75)
#print(s50)
#print(p75.s50)
#multiplot(p1,p2,cols=2)
#multiplot(p3,p4,cols=2)
#print(p.all.dup1)
#print(p.all.dup2)
#print(p.all.dup3)

#print(p.all.uniq1)
#print(p.all.uniq2)
print(p.all.uniq3)
print(p.all.dup4)
#multiplot(p.all.uniq1, p.all.uniq2, cols=2)
#print(p.cnt.gene)
dev.off()


############################################################
# 10K Down-sampling from the ultra-deep (300M) sequencing ##
############################################################
#dt.downsample<-fread("~/results/RNA-Seq/Plasma.2017/Meta/PE150.10K.sampling.flagstat.txt")

