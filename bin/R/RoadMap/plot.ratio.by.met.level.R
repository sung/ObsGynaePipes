#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# Last modified 18/Aug/2016

TR_PREFIX="GRCh37"
library(data.table)
library(reshape2)
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/RoadMap/local.R")
source("~/lib/theme_publish.R")

################################
## RoadMap                    ##
## Schultz et al. Nature 2015 ##
## Tissue comparison          ##
## Female-Male %Met Difference##
################################
my.RData="~/results/RoadMap/CH/RData/CX.freq.RData"
if(file.exists(my.RData)){
	load(my.RData)
	cat("CX.freq loaded\n")
}else{
	my.bin=.05
	all.types=c("CG","CHH","CHG")
	CX.freq<-list()
	for(my.cpg.type in all.types){ 
		dummy1<-list()
		for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){
			## Load this RoadMap tissue
			dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)

			# averaged by gender
			dummy1[[my.tissue]]=data.table(
					tissue=my.tissue,
					context=my.cpg.type,
					min.met=seq(0,1,by=my.bin)*100,
					ratio=sapply(seq(0,1,by=my.bin), function(i){nrow(dt.query[(V5.x+V5.y)/(V6.x+V6.y)>=i])})/nrow(dt.query),
					ratio.female=sapply(seq(0,1,by=my.bin), function(i){nrow(dt.query[V5.x/V6.x>=i])})/nrow(dt.query),
					ratio.male=sapply(seq(0,1,by=my.bin), function(i){nrow(dt.query[V5.y/V6.y>=i])})/nrow(dt.query)
					)
		}# end of my.tissue 
		CX.freq[[my.cpg.type]]<-data.table::rbindlist(dummy1) # aggregate tissue
	}
	save(CX.freq, file=my.RData)
}

file.name<-file.path("~/results/RoadMap/CH/Figures",paste("percentage.by.meth.level",time.stamp,"pdf",sep="."))
pdf(file=file.name, width=11.7, height=8.3, title="Percentage of Occurrence by Methylation Level")

if(TRUE){
	##############################################
	# 1. Fraction of CX by min-methylation level #
	##############################################
	# 1.1. CH   
	p<-ggplot(rbindlist(CX.freq)[context!="CG" & min.met>0 & min.met<40], aes(min.met, ratio*100, col=tissue)) + 
		geom_point() + 
		geom_line(aes(linetype=context)) +  
		scale_colour_manual(values=cbPalette, name="Tissues") + 
		labs(y="Fraction of CH (%)", x="Mimimum methylation level") + 
		ggtitle("Percentage of CH by methylation level") 
		#theme_bw() +
		#theme(axis.text=element_text(size=14), legend.text=element_text(size=14), legend.title=element_text(size=14), axis.title=element_text(size=16), plot.title=element_text(size=20), title.margin=unit(.5, "cm")) +
	print(p+theme_Publication())
	# 1.2. CHG
	p<-ggplot(rbindlist(CX.freq)[context=="CHG" & min.met>0 & min.met<40], aes(min.met, ratio*100, col=tissue)) + 
		geom_point() + 
		geom_line(aes(linetype=context)) +  
		scale_colour_manual(values=cbPalette, name="Tissues") + 
		labs(y="Fraction of CHG (%)", x="Mimimum methylation level") + 
		ggtitle("Percentage of CHG by methylation level")
	print(p+theme_Publication())
	# 1.3. CHH
	p<-ggplot(rbindlist(CX.freq)[context=="CHH" & min.met>0 & min.met<40], aes(min.met, ratio*100, col=tissue)) + 
		geom_point() + 
		geom_line(aes(linetype=context)) +  
		scale_colour_manual(values=cbPalette, name="Tissues") + 
		labs(y="Fraction of CHH (%)", x="Mimimum methylation level") + 
		ggtitle("Percentage of CHH by methylation level") + 
		scale_linetype_manual(values=2)
	print(p+theme_Publication())
	# 1.4. CG
	p<-ggplot(rbindlist(CX.freq)[context=="CG" & min.met>0 & min.met<40], aes(min.met, ratio*100, col=tissue)) + 
		geom_point() + 
		geom_line(aes(linetype=context)) +  
		scale_colour_manual(values=cbPalette, name="Tissues") + 
		labs(y="Fraction of CG (%)", x="Mimimum methylation level") + 
		ggtitle("Percentage of CG by methylation level")
	print(p+theme_Publication())
}

	###############################
	# 2. Fraction of CX by Gender #
	###############################
	# 3-a. All tissues - CHG by gender
if(TRUE){
	p<-ggplot(melt(rbindlist(CX.freq)[context=="CHG"], id=1:3,measure=5:6,variable.name="gender",value.name="ratio")[min.met>0 & min.met<30], aes(min.met, ratio*100, col=tissue)) + 
		geom_point() + 
		geom_line(aes(linetype=gender)) + 
		scale_colour_manual(values=cbPalette, name="Tissues") + 
		labs(y="Fraction of CHG (%)", x="Mimimum methylation level") + 
		ggtitle("CHG")
	print(p+theme_Publication())
	# 3-a. All tissues - CHH by gender
	p<-ggplot(melt(rbindlist(CX.freq)[context=="CHH"], id=1:3,measure=5:6,variable.name="gender",value.name="ratio")[min.met>0 & min.met<30], aes(min.met, ratio*100, col=tissue)) + 
		geom_point() + 
		geom_line(aes(linetype=gender)) + 
		scale_colour_manual(values=cbPalette, name="Tissues") + 
		labs(y="Fraction of CHH (%)", x="Mimimum methylation level") + 
		ggtitle("CHH")
	print(p+theme_Publication())
}


if(FALSE){
	############################
	# 2. PT only - CHH and CHG #
	############################
	p<-ggplot(rbindlist(CX.freq)[tissue=="PT" & context!="CG" & min.met>0 & min.met<50], aes(min.met, ratio*100, col=tissue)) + 
		geom_point() + 
		geom_line(aes(linetype=context)) + 
		scale_colour_manual(values=c(`PT`=cbPalette[7]), name="Tissues") + 
		labs(y="Fraction of CH", x="Mimimum methylation level") + 
		ggtitle("Percentage of CH by methylation level")
	print(p+theme_Publication())

	# 3-a. PT only - CH by gender
	p<-ggplot(melt(rbindlist(CX.freq)[tissue=="PT" & context!="CG"], id=1:3,measure=5:6,variable.name="gender",value.name="ratio")[min.met>0 & min.met<40], aes(min.met, ratio*100, col=context)) + 
		geom_point() + 
		geom_line(aes(linetype=gender)) + 
		ggtitle("Placenta CH by Gender")
	print(p+theme_Publication())
	# 3-a. PT only - CHG by gender
	p<-ggplot(melt(rbindlist(CX.freq)[tissue=="PT" & context=="CHG"], id=1:3,measure=5:6,variable.name="gender",value.name="ratio")[min.met>0 & min.met<40], aes(min.met, ratio*100, col=tissue)) + 
		geom_point() + 
		geom_line(aes(linetype=gender)) + 
		scale_colour_manual(values=c(`PT`=cbPalette[7]), name="Tissues") + 
		ggtitle("Placenta CHG")
	print(p+theme_Publication())
	# 3-a. PT only - CHG by gender
	p<-ggplot(melt(rbindlist(CX.freq)[tissue=="PT" & context=="CHG"], id=1:3,measure=5:6,variable.name="gender",value.name="ratio")[min.met>5 & min.met<40], aes(min.met, ratio*100, col=tissue)) + 
		geom_point() + 
		geom_line(aes(linetype=gender)) + 
		scale_colour_manual(values=c(`PT`=cbPalette[7]), name="Tissues") + 
		ggtitle("Placenta CHG")
	print(p+theme_Publication())
	# 3-c. PT only - CHH by gender
	p<-ggplot(melt(rbindlist(CX.freq)[tissue=="PT" & context=="CHH"], id=1:3,measure=5:6,variable.name="gender",value.name="ratio")[min.met>0 & min.met<40], aes(min.met, ratio*100, col=tissue)) + 
		geom_point() + 
		geom_line(aes(linetype=gender)) + 
		scale_colour_manual(values=c(`PT`=cbPalette[7]), name="Tissues") + 
		ggtitle("Placenta CHH")
	print(p+theme_Publication())
	# 3-c. PT only - CHH by gender
	p<-ggplot(melt(rbindlist(CX.freq)[tissue=="PT" & context=="CHH"], id=1:3,measure=5:6,variable.name="gender",value.name="ratio")[min.met>5 & min.met<40], aes(min.met, ratio*100, col=tissue)) + 
		geom_point() + 
		geom_line(aes(linetype=gender)) + 
		scale_colour_manual(values=c(`PT`=cbPalette[7]), name="Tissues") + 
		ggtitle("Placenta CHH")
	print(p+theme_Publication())
}

dev.off()

cat("All done\n")
