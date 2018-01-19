
	# see bin/R/RoadMap/get.dt.meth.level.per.region.R for dt.meth.region.RData
	my.RData=file.path("~/results/RoadMap/BS-Seq/RData/10X",paste(my.tissue,my.cpg.type,"dt.meth.region.RData",sep="."))
	if(file.exists(my.RData)){
		cat("loading dt.meth.region.RData...\n")
		load(my.RData)
		cat("dt.meth.region loaded\n")
	}else{
		stop("RData not exist")
	}
	#> dcast(dt.meth.region[,.N,"Region,Gender"], Region ~ Gender, value.var="N")
	#        Region   Female     Male
	#1:        CPGi    26620    26620
	#2: CPGi-shores    48208    48208
	#3:    Enhancer    27653    27653
	#4:       Genes    44270    44270
	#5: Promo_15.05    50781    50781
	#6:       Sites 14960649 14960649
	#7:    Tiling5K   555787   555787

	#######################
	## Figure 2          ##
	## Meth Diff Boxplot ##
	#######################
	dt.meth<-dt.meth.region[Tissue==my.tissue & Region %in% c("Sites","Genes","Promo_15.05","CPGi")]
	dt.meth[,`Chromosome Type`:=ifelse(chr=="X",'X-chromosome','Autosomes')]
	dt.meth.gender=data.table(dcast(dt.meth[Tissue==my.tissue], `Chromosome Type`+Region+chr+strand+start+end~Gender, value.var="meth")) 
	p.box<-ggplot(dt.meth.gender, aes(Region, (Female-Male)*100)) + 
				geom_boxplot(aes(fill=`Chromosome Type`), outlier.shape=NA, alpha=.7, size=.7, width=.5) + 
				labs(x="Genomic Regions",y="% Methylation Difference\n   (more in male)                     (more in female)") + 
				geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
				#ylim(c(-50,50)) +
				coord_cartesian(ylim=c(-60, 40)) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
				scale_y_continuous(breaks=seq(from=-60,to=40,by=10)) +
				scale_x_discrete(limits=c("Sites","Genes","Promo_15.05","CPGi"),labels=c("CpG sites","Gene-bodies","Promoters","CPG islands")) +
				scale_fill_manual(values=my.col[["Chromosome Types"]]) +
				theme_Publication() +
				theme(legend.position=c(0.87,0.91), plot.title = element_text(hjust = 0, size = 15,face="bold"))

	file.name<-file.path("~/results/CSMD1.2016.paper/Figures/Meth.diff",paste("Fig2",my.tissue,my.cpg.type,time.stamp,"tiff",sep="."))
	tiff(filename=file.name,width=14, height=10,units="in",res=300, compression = 'lzw') #A4 size 
	print(p.box)
	dev.off()


