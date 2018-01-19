
	# see bin/R/RoadMap/get.dt.meth.level.per.region.R for dt.meth.region.RData
	my.RData=file.path("~/results/RoadMap/BS-Seq/RData",paste(my.tissue,my.cpg.type,"dt.meth.region.RData",sep="."))
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
	dt.meth=dt.meth.region[chr=="X" & Tissue==my.tissue & Region=="Sites"]
	p.site.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			#geom_density(size=1) + 
			geom_line(stat="density", size=1) +
			labs(x="",y="Density") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			#ylim(c(0,0.062)) +
			theme_Publication() +
			theme(legend.position="none")

	dt.meth=dt.meth.region[chr=="X" & Tissue==my.tissue & Region=="Genes"]
	min.cpg.gene<-quantile(dt.meth[,num.sites], .25) # No. of CpG at 25%-percentile # min.cpg: 6
	#p.gene.beta<-ggplot(dt.meth[num.sites>=min.cpg.gene], aes(meth*100, col=Gender)) + 
	p.gene.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			#geom_density(size=1) + 
			geom_line(stat="density", size=1) +
			labs(x="",y="") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			#ylim(c(0,0.062)) +
			theme_Publication() +
			theme(legend.position=c(0.91,0.91)) 
	
	dt.meth=dt.meth.region[chr=="X" & Tissue==my.tissue & Region=="Promo_15.05"]
	min.cpg.promo<-quantile(dt.meth[,num.sites], .25) # No. of CpG at 25%-percentile # min.cpg: 4
	#p.promoter.beta<-ggplot(dt.meth[num.sites>=min.cpg.promo], aes(meth*100, col=Gender)) + 
	p.promoter.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			#geom_density(size=1) + 
			geom_line(stat="density", size=1) +
			labs(x="% Methylation Level", y="Density") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			#ylim(c(0,0.062)) +
			theme_Publication() +
			theme(legend.position="none")

	dt.meth=dt.meth.region[chr=="X" & Tissue==my.tissue & Region=="CPGi"]
	min.cpg.cpgi<-quantile(dt.meth[,num.sites], .25) # No. of CpG at 25%-percentile # min.cpg: 28
	#p.cpgi.beta<-ggplot(dt.meth[num.sites>=min.cpg.cpgi], aes(meth*100, col=Gender)) + 
	p.cpgi.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			#geom_density(size=1) + 
			geom_line(stat="density", size=1) +
			labs(x="% Methylation Level", y="") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			#ylim(c(0,0.062)) +
			theme_Publication() + 
			theme(legend.position="none")

	#file.name<-file.path("~/results/CSMD1.2016.paper/Figures/Beta.dist",paste("Fig1",my.tissue,my.cpg.type,time.stamp,"line.ylim.tiff",sep="."))
	file.name<-file.path("~/results/CSMD1.2016.paper/Figures/Beta.dist",paste("Fig1",my.tissue,my.cpg.type,time.stamp,"line.chrX.tiff",sep="."))
	#file.name<-file.path("~/results/CSMD1.2016.paper/Figures/Beta.dist",paste("Fig1",my.tissue,my.cpg.type,time.stamp,"line.autosome.tiff",sep="."))
	tiff(filename=file.name,width=3000, height=3000,units="px",res=300, compression = 'lzw') 
	cowplot::plot_grid(p.site.beta, p.gene.beta, p.promoter.beta, p.cpgi.beta, labels=c("B","C","D","E"), label_size=15, ncol = 2, nrow = 2)
	dev.off()

