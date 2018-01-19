#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

source("~/Pipelines/config/graphic.R")
library(data.table)
library(DESeq2)

dt.target=data.table(
	Gene=c("SLC38A1","SLC38A2","SLC38A4","IGF2"),
	ENSG=c("ENSG00000111371","ENSG00000134294","ENSG00000139209","ENSG00000167244"))

dt.bw<-fread("~/Pipelines/data/RNA-Seq/rna.seq.sample.bw.pw.csv")

#myProject="Boy.Girl.FG.JD.GRCh38" # Health Placenta Samples
#sampleType="AGA"

myProject="SGA.AGA.SE125.GRCh38" # SGA & AGA Samples
sampleType="ALL"

pdf(file=file.path(paste0("~/results/RNA-Seq/",myProject),paste0(myProject,".Slc.bw.fpkm.pdf")), width=11.7, height=11, title="FPKM vs BW (PW) for Slc")

myRData=paste("~/results/RNA-Seq",myProject,paste0("DESeq2/",sampleType,"/deseq.",sampleType,".RData"),sep="/")
load(myRData)
cat("dds and res loaded\n")

ddsFpkm <-fpkm(dds) # isa 'matrix'

dt.fpkm<-data.table(melt(t(ddsFpkm[dt.target$ENSG,]), varnames=c("SampleName","ENSG"),value.name="FPKM"))
dt.fpkm<-merge(dt.fpkm, data.table(SampleName=rownames(colData(dds)), data.table(as.data.frame(colData(dds)))), by="SampleName")
dt.fpkm<-merge(dt.fpkm, dt.target, by="ENSG")

dt.fpkm.bw<-merge(dt.fpkm, dt.bw, by="SampleName")

# PW
#ggplot(dt.fpkm.bw, aes(FPKM, PW_Zscore, col=Gene)) + geom_point(alpha=.7,size=2) + geom_smooth(method="lm",se=FALSE) + scale_colour_manual(values=cbPalette2) + theme_Publication()
#ggplot(dt.fpkm.bw, aes(FPKM, PW_Zscore, col=Gene)) + geom_point(alpha=.7,size=2) + geom_smooth() + scale_colour_manual(values=cbPalette2) + facet_wrap(~Gene,nrow=2) + theme_Publication()

# Healthy
if(grepl("Boy.Girl",myProject)){
	for(i in c("BW","BW_Zscore", "PW","PW_Zscore")){
		r_squared=lapply(split(dt.fpkm.bw, dt.fpkm.bw[,unique(Gene)]), function(j) round(summary(lm(j[,.(FPKM,get(i))]))$r.squared,4))

		my_labeller<- function(variable,value){
			  return(paste(value,"(R^2=",r_squared[value],")"))
		}

		p<-ggplot(dt.fpkm.bw, aes_string("FPKM", i, col="Gene")) + 
			geom_point(alpha=.7,size=2) + 
			geom_smooth(method="lm",se=FALSE) + 
			scale_colour_manual(values=cbPalette2) + 
			theme_Publication()
		print(p)

		p<-ggplot(dt.fpkm.bw, aes_string("FPKM", i, col="Gene")) + 
			geom_point(alpha=.7,size=2) + 
			geom_smooth(method="lm",se=FALSE) + 
			scale_colour_manual(values=cbPalette2) + 
			facet_wrap(~Gene,nrow=2, labeller=my_labeller) + 
			theme_Publication()
		print(p)
	}
}else if(grepl("SGA.AGA",myProject)){
	# Boxplot of FPKM by SGA/AGA
	p<-ggplot(dt.fpkm.bw, aes(Gene, FPKM, fill=Condition)) + 
		geom_boxplot(alpha=.8) + 
		annotate("text",x=dt.target[order(Gene),Gene],y=210,label=paste0("P=",round(res[dt.target[order(Gene)]$ENSG,"padj"],4)),size=6) +
		scale_fill_manual(values=cbPalette, labels=c("AGA","SGA")) + 
		theme_Publication()
	print(p)

	# FPKM vs. BW etc. (xy-plot) by SGA/AGA
	for(i in c("BW","BW_Zscore", "PW","PW_Zscore")){
		r_squared_aga=lapply(split(dt.fpkm.bw[Condition==0], dt.fpkm.bw[Condition==0,unique(Gene)]), function(j) round(summary(lm(j[,.(FPKM,get(i))]))$r.squared,4))
		r_squared_sga=lapply(split(dt.fpkm.bw[Condition==1], dt.fpkm.bw[Condition==1,unique(Gene)]), function(j) round(summary(lm(j[,.(FPKM,get(i))]))$r.squared,4))

		my_labeller<- function(variable,value){
			  return(paste(value,"(AGA=",r_squared_aga[value],"SGA=",r_squared_sga[value],")"))
		}

		p<-ggplot(dt.fpkm.bw, aes_string("FPKM", i, col="Condition")) + 
			geom_point(alpha=.7,size=2) + 
			geom_smooth(method="lm",se=FALSE) + 
			scale_colour_manual(values=cbPalette2, labels=c("AGA","SGA")) + 
			facet_wrap(~Gene,nrow=2, labeller=my_labeller) + 
			theme_Publication()
		print(p)
	}

}else{
	stop("Not supported")
}

dev.off()

cat("All done\n")
