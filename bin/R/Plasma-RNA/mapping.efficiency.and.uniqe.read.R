library(data.table)
library(ggplot2)
source("~/Pipelines/config/graphic.R")
time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p')
		########################
		## Mapping efficiency ##
        ## Low-PSG7 sample    ##
		########################
		dt.map<-fread("~/results/RNA-Seq/Plasma.2017/Meta/PE75.PE150.read.cnt.txt",stringsAsFactors=TRUE); dt.map$Type<-factor(dt.map$Type,levels=rev(levels(dt.map$Type)))
		p0<-ggplot(dt.map[Library %in% c("SLX-9342","SLX-9345","SLX-11368","SLX-11369") & Type %in% c("PE150","PE75") & Yield>10^7], aes(Yield/10^6,Mapping_Efficiency)) + 
				geom_point(aes(col=Type),size=6,alpha=.7) + 
				labs(x="No. of Reads (Millions)", y="Mapping Efficiency(%)") + 
				ggtitle("Plasma RNA-Seq") + 
				scale_colour_manual(values=rev(cbPalette)) + 
				scale_x_continuous(breaks=seq(0,400,50)) +
				theme_Publication() + theme(legend.position=c(0.8,.8))

		p0.1<-ggplot(dt.map[Library %in% c("SLX-9342","SLX-9345") & Type=="PE75"], aes(Yield/10^6,Mapping_Efficiency)) + 
				geom_point(aes(col=Library),size=6,alpha=.5) + 
				labs(x="No. of Reads (millions)", y="Mapping Efficiency (%)") + 
				#ggtitle("Plasma RNA-Seq") + 
				scale_colour_manual(values=cbPalette, label=c(paste0(c("SLX-9342 (N=", eval(dt.map[Library=="SLX-9342" & Type=="PE75",.N]),")"), collapse=""), paste0(c("SLX-9345 (N=", eval(dt.map[Library=="SLX-9345" & Type=="PE75",.N]),")"), collapse=""))) + 
				theme_Publication()

		p0.2<-ggplot(dt.map[Library %in% c("SLX-9342","SLX-9345") & Type=="PE75"], aes(Library,Yield/10^6)) + 
				geom_boxplot(aes(fill=Library),alpha=0.5, outlier.shape=NA, width=.3) + 
				#coord_cartesian(ylim=c(20,50)) +
				#labs(x="Library", y="No. of reads (millions)") + 
				scale_fill_manual(values=cbPalette,label=c(paste0(c("SLX-9342 (N=", eval(dt.map[Library=="SLX-9342" & Type=="PE75",.N]),")"), collapse=""), paste0(c("SLX-9345 (N=", eval(dt.map[Library=="SLX-9345" & Type=="PE75",.N]),")"), collapse=""))) + 
				theme_Publication() 

		p0.3<-ggplot(dt.map[Library %in% c("SLX-9342","SLX-9345") & Type=="PE75"], aes(Library,Mapping_Efficiency)) + 
				geom_boxplot(aes(fill=Library),alpha=0.5, outlier.shape=NA, width=.3) + 
				#coord_cartesian(ylim=c(40,90)) +
				#ggtitle("Plasma RNA-Seq (PE75)") + 
				labs(x="Library", y="Mapping Efficiency(%)") + 
				scale_fill_manual(values=cbPalette,label=c(paste0(c("SLX-9342 (N=", eval(dt.map[Library=="SLX-9342" & Type=="PE75",.N]),")"), collapse=""), paste0(c("SLX-9345 (N=", eval(dt.map[Library=="SLX-9345" & Type=="PE75",.N]),")"), collapse=""))) + 
				theme_Publication()

		p1<-ggplot(dt.map, aes(Yield/10^6,Mapping_Efficiency)) + 
				geom_point(aes(col=Type),size=5,alpha=.5) + 
				#ylim(c(0,90)) + 
				#xlim(c(0,4e+8)) + 
				labs(x="No. of Reads (Millions)", y="Mapping Efficiency(%)") + 
				ggtitle("Plasma RNA-Seq") + 
				scale_colour_manual(values=cbPalette) + 
				#scale_x_log10() +
				theme_Publication() + theme(legend.position=c(0.8,.8))

		# XY-plot
		p2<-ggplot(dcast(dt.map[Library=="SLX-9342" & Type %in% c("PE75","PE150")], Barcode~Type, value.var="Mapping_Efficiency"), aes(PE150, PE75)) + 
				geom_point(size=6,alpha=.5) + 
				ylim(c(0,100)) + 
				xlim(c(0,100)) + 
				ggtitle("Plasma RNA-Seq (SLX-9342)") + 
				theme_Publication()

		# XY-plot
		p3<-ggplot(dcast(dt.map[Library=="SLX-9342" & Type %in% c("PE75","SE50")], Barcode~Type, value.var="Mapping_Efficiency"), aes(SE50, PE75)) + 
				geom_point(size=6,alpha=.5) + 
				ylim(c(0,100)) + 
				xlim(c(0,100)) + 
				ggtitle("Plasma RNA-Seq (SLX-9342)") + 
				theme_Publication()
		#multiplot(p2,p3)

		dt.dummy<-dt.map[Library %in% c("SLX-11368","SLX-11369") | Barcode %in% c("D701_D508","D711_D508")]; dt.dummy$Individual<-c("A","B","A","B","A","A","B")
		p4<-ggplot(dt.dummy[Yield>10^7], aes(Yield/10^6,Mapping_Efficiency)) + 
				geom_point(aes(col=Type,shape=Individual),size=5) + 
				ylim(c(0,100)) + 
				labs(x="No. of Reads (millions)", y="Mapping Efficiency (%)") + 
				ggtitle("Two Individuals of Low (A) and High (B) PSG7 Level") + 
				scale_colour_manual(values=rev(cbPalette)) + 
				theme_Publication()

		###################################################################################
		# 10K Down-sampling from the ultra-deep (300M) sequencing (SLX-11369, High PSG7) ##
		# graph of 1) mapping ratio (%) and 2) uniqe/multi reads                         ##
		###################################################################################
		dt.downsample<-fread("~/results/RNA-Seq/Plasma.2017/Meta/PE150.10K.sampling.read.cnt.txt")
		dt.downsample$Length<-factor(dt.downsample$Length, levels=c('PE50','PE75','PE100','PE125','PE150'))

		p5<-ggplot(dt.downsample[Source!="SLX-11369"], aes(Length, `Mapping Efficiency`,group=BaseQ)) + 
				geom_point(aes(col=as.factor(BaseQ)),size=5) +
				geom_line(aes(col=as.factor(BaseQ)),linetype=2,size=1.2) +
				labs(y="Mapping Efficiency (%)") +
				scale_colour_brewer(name="Base-quality\nThreshold",palette="Set1") +
				ggtitle("Down-Sampled Reads (n=10K from 365M)") + 
				theme_Publication()

		dt.dummy<-melt(dt.downsample[Source!="SLX-11369"], id.vars=c("Length","BaseQ"), measure.vars=c("After Trimming","Mapped Read"), variable.name="Type", value.name="Number")
		dt.dummy[,Number:=ifelse(Type=="After Trimming",Number*2,Number)]
		p6<-ggplot(dt.dummy, aes(Length, Number)) +
				geom_bar(aes(fill=Type),stat="identity",position="dodge") +
				labs(y="No. of Read (R1+R2)") +
				scale_fill_grey() +
				ggtitle("Down-Sampled Reads (n=10K from 365M)\nTrimming at different base-quality threshold") + 
				facet_grid(~BaseQ) +
				theme_Publication() +
				theme(axis.text.x = element_text(angle = 45, hjust = 1))

		dt.dummy<-melt(dt.downsample[Source!="SLX-11369"], id.vars=c("Length","BaseQ"), measure.vars=c("Multi Mapping","Assigned Read","Ambiguity"), variable.name="Type", value.name="Number")
		#dt.dummy<-melt(dt.downsample[Source!="SLX-11369"], id.vars=c("Length","BaseQ"), measure.vars=c("Multi Mapping","Assigned Read"), variable.name="Type", value.name="Number")
		p7<-ggplot(dt.dummy, aes(Length, Number,group=Type)) +
				geom_point(aes(col=Type),size=5) +
				geom_line(aes(col=Type)) +
				labs(y="Number of Quantified Reads (by featureCount)") +
				scale_colour_manual(values=cbPalette2) +
				ggtitle("Down-Sampled Reads (n=10K from 365M)\nTrimming at different base-quality threshold") + 
				facet_grid(~BaseQ) + 
				theme_Publication() +
				theme(axis.text.x = element_text(angle = 45, hjust = 1))

		dt.dummy<-melt(dt.downsample[Source!="SLX-11369",.(Length,BaseQ,`Uniq Reads`=Total-Secondary-Duplicated)], id.vars=c("Length","BaseQ"), measure.vars=c("Uniq Reads"), variable.name="Type", value.name="Number")
		p8<-ggplot(dt.dummy, aes(Length, Number,group=Type)) +
				geom_point(aes(col=Type),size=5) +
				geom_line(aes(col=Type)) +
				scale_colour_manual(values=cbPalette2) +
				ggtitle("Down-Sampled Reads (n=10K from 365M)\nTrimming at different base-quality threshold") + 
				facet_grid(~BaseQ) + 
				theme_Publication() +
				theme(axis.text.x = element_text(angle = 45, hjust = 1))

		my.filename <- paste("~/results/RNA-Seq/Plasma.2017/Cluster/PE150.10K.downsample.mapping.efficiency",time.stamp,"pdf",sep=".")
		pdf(file=my.filename, width=12, height=10, title="Mapping Efficiency of PE150")
		#print(p0)
		#print(p1)
		#print(p2)
		#print(p3)
		#print(p4)
		#print(p5)
		#print(p6)
		#print(p7)

		print(p0)
		print(
			gridExtra::grid.arrange(p0.3+labs(x="", y="")+theme(axis.text.x = element_blank(),legend.position=""),  
									p0.1+theme(legend.position=c(.9,.9), legend.background = element_rect(colour = 'black', linetype='dashed')), 
									blankPlot, # defined within graphic.R
									p0.2+coord_flip()+labs(x="", y="")+theme(axis.text.y = element_blank(),legend.position=""), 
									ncol=2, nrow=2, widths=c(1, 4), heights=c(4,1))
			)
		print(p4)
		dev.off()
