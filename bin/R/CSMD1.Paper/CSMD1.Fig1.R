
	# see bin/R/RoadMap/get.dt.meth.level.per.region.R for dt.meth.region.RData
	my.RData=file.path("~/results/RoadMap/BS-Seq/RData/10X",paste(my.tissue,my.cpg.type,"dt.meth.region.RData",sep="."))
	if(file.exists(my.RData)){
		cat("loading dt.meth.region.RData...\n")
		load(my.RData)
		cat("dt.meth.region loaded\n")
	}else{
		stop("RData not exist")
	}

	########################
	## Figure 1           ##
	## Beta Distribution  ##
	## via plot_grid      ##
	########################
	library("cowplot") # only for R >= 3.2.1
	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="Sites"]
	dt.meth[,chr.type:=ifelse(chr=="X","X-chromosome","Autosomes")]
	p.site.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_line(stat="density", size=1) +
			ggtitle("CpG sites") +
			labs(x="",y="Density") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			facet_grid(.~chr.type) +
			theme_Publication() +
			theme(legend.position="none")

	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="Genes"]
	dt.meth[,chr.type:=ifelse(chr=="X","X-chromosome","Autosomes")]
	min.cpg.gene<-quantile(dt.meth[,num.sites], .25) # No. of CpG at 25%-percentile # min.cpg: 6
	#p.gene.beta<-ggplot(dt.meth[num.sites>=min.cpg.gene], aes(meth*100, col=Gender)) + 
	p.gene.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			#geom_density(size=1) + 
			geom_line(stat="density", size=1) +
			ggtitle("Gene-bodies") +
			labs(x="",y="") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			#ylim(c(0,0.062)) +
			facet_grid(.~chr.type) +
			theme_Publication() +
			theme(legend.position="none")
	
	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="Promo_15.05"]
	dt.meth[,chr.type:=ifelse(chr=="X","X-chromosome","Autosomes")]
	p.promoter.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_line(stat="density", size=1) +
			ggtitle("Promoters") +
			labs(x="% Methylation Level", y="Density") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			facet_grid(.~chr.type) +
			theme_Publication() +
			theme(legend.position="none")

	dt.meth=dt.meth.region[Tissue==my.tissue & Region=="CPGi"]
	dt.meth[,chr.type:=ifelse(chr=="X","X-chromosome","Autosomes")]
	p.cpgi.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_line(stat="density", size=1) +
			ggtitle("CpG islands") +
			labs(x="% Methylation Level", y="") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			facet_grid(.~chr.type) +
			theme_Publication() + 
			theme(legend.position=c(0.91,0.87)) 

	#file.name<-file.path("~/results/CSMD1.2016.paper/Figures/Beta.dist",paste("Fig1",my.tissue,my.cpg.type,time.stamp,"line.ylim.tiff",sep="."))
	file.name<-file.path("~/results/CSMD1.2016.paper/Figures/Beta.dist",paste("Fig1",my.tissue,my.cpg.type,time.stamp,"line.chrX.autosomes.tiff",sep="."))
	tiff(filename=file.name, width=12, height=10, units="in",res=300, compression = 'lzw') #A4 size 
	cowplot::plot_grid(p.site.beta, p.gene.beta, p.promoter.beta, p.cpgi.beta, labels=c("A","B","C","D"), label_size=20, ncol = 2, nrow = 2)
	dev.off()

