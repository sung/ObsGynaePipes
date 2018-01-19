library(data.table)
source("~/Pipelines/config/graphic.R")

foo<-fread("~/Pipelines/bin/R/Plasma-RNA/CSV/SLX-11369-PE150.PhiX.csv")
foo[,`Unmapped`:=2*Sampled_reads-(PhiX+Human)]

bar<-melt(foo, id.vars=c("DownSampled","Sampled_reads"), measure.vars=c("PhiX","Human","Unmapped"), variable.name="Reference", value.name="Number") # column-wise to row-wise
bar$DownSampled<-factor(bar$DownSampled, levels=c('10K','60K','100K'))

p.bar<-ggplot(bar, aes(DownSampled, Number)) + 
	geom_bar(aes(fill=`Reference`),color='black',stat="identity") + 
	ggtitle("No. of mapped reads for SLX-11369 (PE150: 365M reads)") +
	scale_fill_manual(values=cbPalette) +
	theme_Publication()

var<-bar[,.(DownSampled,Reference,`Ratio`=Number/Sampled_reads*.5*100)]
p.pie1<-ggplot(var[DownSampled=="10K"], aes(x = "", y=Ratio, fill=Reference)) + 
	geom_bar(stat="identity",color='black',width = 1) + 
	coord_polar(theta = "y") +
	scale_fill_manual(values=cbPalette) +
	theme(axis.text.x=element_blank()) +
	geom_text(aes(y=cumsum(Ratio) - Ratio/ 2,label = paste0(round(Ratio,1),"%")), size=5) +
	theme_Publication() +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks = element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		legend.position="none"
		) +
	facet_wrap(~DownSampled) 
p.pie2<-ggplot(var[DownSampled=="60K"], aes(x = "", y=Ratio, fill=Reference)) + 
	geom_bar(stat="identity",color='black',width = 1) + 
	coord_polar(theta = "y") +
	scale_fill_manual(values=cbPalette) +
	theme(axis.text.x=element_blank()) +
	geom_text(aes(y=cumsum(Ratio) - Ratio/ 2,label = paste0(round(Ratio,1),"%")), size=5) +
	theme_Publication() +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks = element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		legend.position="none"
		) +
	facet_wrap(~DownSampled) 
p.pie3<-ggplot(var[DownSampled=="100K"], aes(x = "", y=Ratio, fill=Reference)) + 
	geom_bar(stat="identity",color='black',width = 1) + 
	coord_polar(theta = "y") +
	scale_fill_manual(values=cbPalette) +
	theme(axis.text.x=element_blank()) +
	geom_text(aes(y=cumsum(Ratio) - Ratio/ 2,label = paste0(round(Ratio,1),"%")), size=5) +
	theme_Publication() +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks = element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		legend.position="none"
		) +
	facet_wrap(~DownSampled) 

my.file.name<- "~/results/RNA-Seq/Plasma.2017/Cluster/SLX-11369.phiX.mapping.ratio"
pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11, height=8, title="Low PSG7 Sample: Plasma and Tissue") # A4 size
print(p.bar)
multiplot(p.pie1,p.pie2,p.pie3,cols=3)
dev.off()

mapply()
