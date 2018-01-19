#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

library(data.table) 
library(scales) # for 'lables=comma'
source("~/Pipelines/bin/R/Placentome/local.R")
source("~/Pipelines/config/graphic.R")

mySource="Placentome" # "FG"
myCohort="POPS" # "POPS", "CLT", "PET", "SGA"
myProject=paste(mySource,myCohort,sep=".") # e.g. Placentome.CTL
top.dir=file.path("~/results/RNA-Seq",mySource) #top.dir="~/results/RNA-Seq/Placentome"

########################
# pdf output filename ##
########################
cat("making PDF file...\n")
my.file.name<- file.path(top.dir, paste0("ReadCounts/",myProject,".read.count.ggplot")) # per transcript (TCON_)
pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title=paste0("Read Count of ",myProject)) # A4 size

###################################################
## Read Count Data                               ##
## see ~/results/RNA-Seq/Placentome/Meta/README  ##
## dt.cnt.cohort defined within local.R          ##
###################################################

cat("Making some graph...\n")
#http://www.sthda.com/english/wiki/ggplot2-pie-chart-quick-start-guide-r-software-and-data-visualization
## 2. No of Read by Mapping Category 
dummy<-dt.cnt.cohort[,list(Read=round(sum(as.numeric(Read))/10^9,2)),"Category"]
p<-ggplot(dummy, aes(x="",y=Read, fill=Category)) + 
	geom_bar(stat = "identity", width=1) + 
	coord_polar(theta = "y") + 
	scale_fill_grey() + 
	theme(axis.text.x=element_blank()) +
	geom_text(aes(y = Read/3 + c(0, cumsum(Read)[-length(Read)]), label = Read), size=5) +
	ggtitle(paste("Reads (", dummy[,sum(Read)], " Gb (10^9))"))
print(p)

## 3. No Read per cohort 
dummy<-dt.cnt.cohort[,list(Read=sum(Read)),"Cohort,Library,Barcode"]
p<-ggplot(dummy, aes(Cohort,Read/10^6,col=Cohort)) + 
	geom_violin(width=.4, trim=FALSE) + 
	geom_boxplot(width=.1) + 
	#geom_jitter(position=position_jitter(0.2)) +
	geom_dotplot(binaxis='y', stackdir='center', dotsize=.2, position=position_dodge(1) # same as above
	scale_colour_manual(values=cbPalette) 
print(p)

dev.off()


#SLX-9169.SE50
if(FALSE){

#SLX=SLX-9169
#BARCODE=$(for i in `ls /home/ssg29/scratch/data/fastq/$SLX/$SLX*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f8 | cut -d'.' -f2 ; done | uniq)
#FILE="~/results/$SLX.Homo_sapiens.SE125.v2/$SLX.SE125.mapping.barchart.grch38.txt"
#echo "Barcode Category Read" > $FILE 
#for j in $BARCODE; do printf "$j Mapped "; awk '/Mapped/{print $3}' ~/results/$SLX.Homo_sapiens.v4/TopHat/$j/align_summary.txt; done >> $FILE 
#for j in $BARCODE; do printf "$j Rescued "; awk '/Mapped/{print $3}' ~/results/$SLX.Homo_sapiens.v4/TopHat/$j/unmapped/align_summary.txt; done >> $FILE 
#for j in $BARCODE; do printf "$j Unmapped "; awk '/Input/{input=$3};/Mapped/{mapped=$3}END{print input-mapped}' ~/results/$SLX.Homo_sapiens.v4/TopHat/$j/unmapped/align_summary.txt; done >> $FILE 

	cat("SE50\n")
	my.input<-read.delim("~/results/SLX-9169.Homo_sapiens.v4/SLX-9169.mapping.barchart.txt", sep=" ")

	out.file.name<-"~/results/RNA-Seq/SGA.AGA.total.SE50/SLX-9169.SE50.mapping.barchart.grch37"
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw')

	ggplot(my.input, aes(x = Barcode, y = Read, fill=Category)) + 
	geom_bar(stat = "identity") + 
	ggtitle("Read Counts of SLX-9169 (SE50)") + 
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	scale_y_continuous(labels = comma) + 
	scale_fill_grey()

	dev.off()

#SLX-9169.SE125 (GRCh37)
	cat("SE125\n")
	my.input<-read.delim("~/results/SLX-9169.Homo_sapiens.SE125.v1/SLX-9169.SE125.mapping.barchart.txt", sep=" ")

	out.file.name<-"~/results/RNA-Seq/SGA.AGA.total.SE125/SLX-9169.SE125.mapping.barchart.grch37"
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw')

	ggplot(my.input, aes(x = Barcode, y = Read, fill=Category)) + 
	geom_bar(stat = "identity") + 
	ggtitle("Read Counts of SLX-9169 (SE125)") + 
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	scale_y_continuous(labels = comma) + 
	scale_fill_grey()
	dev.off()

	#SLX-9169.SE125 (GRCh38)
	cat("SE125\n")
	my.input<-read.delim("~/results/SLX-9169.Homo_sapiens.SE125.v2/SLX-9169.SE125.mapping.barchart.txt", sep=" ")

	out.file.name<-"~/results/RNA-Seq/SGA.AGA.total.SE125/SLX-9169.SE125.mapping.barchart.grch38"
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw')

	ggplot(my.input, aes(x = Barcode, y = Read, fill=Category)) + 
	geom_bar(stat = "identity") + 
	ggtitle("Read Counts of SLX-9169 (SE125, GRCH38)") + 
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	scale_y_continuous(labels = comma) + 
	scale_fill_grey()

	dev.off()

	#SLX-9169.SE125 unmapped reads (GRCh37 vs. GRCh38)
	cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	my.input<-read.delim("~/results/SLX-9169.Homo_sapiens.SE125.v2/SLX-9169.SE125.mapping.barchart.grch37.grch38.txt", sep=" ")

	out.file.name<-"~/results/RNA-Seq/SGA.AGA.total.SE125/SLX-9169.SE125.unmapped.barchart.grch37.grch38"
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw')

	ggplot(my.input, aes(x = Barcode, y = Read, fill=Category)) + 
	geom_bar(stat="identity", position="dodge") +
	ggtitle("Final Unmapped Read Counts of SLX-9169") + 
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	scale_y_continuous(labels = comma) +
	scale_fill_manual(values=cbPalette)
	dev.off()
	# histogram
	ggplot(my.input, aes(Read, fill=Category)) + 
	geom_density(alpha=0.2) + 
	ggtitle("Final Unmapped Read Counts of SLX-9169") +
	scale_x_continuous(labels = comma)
	# as a box-plot
	out.file.name<-"~/results/RNA-Seq/SGA.AGA.total.SE125/SLX-9169.SE125.unmapped.boxplot.grch37.grch38"
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw')
	ggplot(my.input, aes(x = Category, y = Read, fill=Category)) + 
	geom_boxplot() +
	ggtitle("Final Unmapped Read Counts of SLX-9169") + 
	scale_y_continuous(labels = comma) +
	scale_fill_manual(values=cbPalette)
	dev.off()

	#small RNA-Seq
	#SLX-9176.v2
	cat("small RNA\n")
	my.input<-read.delim("~/results/SLX-9176.Homo_sapiens.v2/SLX-9176.v2.mapping.barchart.txt", sep=" ")
	out.file.name<-"~/results/RNA-Seq/SGA.AGA.small.v2/SLX-9176.v2.mapping.barchart"
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw')

	ggplot(my.input, aes(x = Barcode, y = Read, fill=Category)) + 
	geom_bar(stat = "identity") + 
	ggtitle("Read Counts of SLX-9176 (SE50)") + 
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	scale_y_continuous(labels = comma) + 
	scale_fill_grey(breaks=c("Trimmed","Filtered","Mapped","Unmapped"))

	dev.off()

}
