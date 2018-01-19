# part of ~/Pipelines/bin/R/CSMD1.Paper/CSMD1.paper.figures.main.R

	library(ggbio)
	# R/3.1.1 (this ggbio (autoplot specifically) conflicts with ggplot which may be >2.0)
	# use R/3.2.1 (this version of ggbio works, but color/fill does does not work :()

	dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type) # ida 'data.table' for PT

	## CSMD1 Gene +- 1kb (Grch37) ##
	####################
	# transcript track #
	####################
	# hg19ensGene defined in 'Annotation.R'
	#p.tr.reduced<-autoplot(hg19ensGene, which=gr.csmd.ext, stat="reduce", color="brown", fill="brown") + theme_alignment() 
	p.tr.reduced<-ggbio::autoplot(hg19ensGene, which=gr.csmd.ext, stat="reduce") + ggbio::theme_alignment() 

	####################
	# WGBS methylation #
	####################
	my.point.size=1.5; my.alpha=0.5
	# gr.csmd, gr.csmd.ext defined in ~/Pipelines/bin/R/CSMD1.Paper/local.R 
	# Ensembl (or NCBI) style coordinate (e.g. 8, not chr8)
	dt.csmd.meth<-dt.query[V1==mapSeqlevels( seqlevels(gr.csmd.ext),style="NCBI" ) & V2>=start(gr.csmd.ext) & V2<=end(gr.csmd.ext)
						, .(`chr`=paste("chr",V1), `start`=V2,`strand`=V3,`Female`=round(V5.x/V6.x*100,2), `Male`=round(V5.y/V6.y*100,2))]

	p.pt.meth<-ggplot(melt(dt.csmd.meth, id=1:3,measure=4:5,variable.name="Gender",value.name="Meth"), aes(start, Meth, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) +
				ylim(0,100) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				labs(y="% Methylation") +
				theme_Publication() + 
				theme(legend.position="none")
	####################
	# RNA-Seq Coverage #
	####################
	# UCSC-style (e.g. chr7) 
	rna.cov.list<-list(
		`Female`=as.data.table(as.data.frame(import.bedGraph("~/results/RNA-Seq/Placentome/BedGraph/CTLF/CTLF.GRCh37.75.CSMD1.bedgraph.gz"))),
		`Male`=as.data.table(as.data.frame(import.bedGraph("~/results/RNA-Seq/Placentome/BedGraph/CTLM/CTLM.GRCh37.75.CSMD1.bedgraph.gz")))
	)
	dt.rna.cov<-rbindlist(lapply(names(rna.cov.list), function(i) rna.cov.list[[i]][,`Gender`:=i]))

	p.rna.cov<-ggplot(dt.rna.cov[score!=0], aes(xmin=start, xmax=end, ymin=0, ymax=score, col=Gender, fill=Gender),alpha=0.1) + 
				geom_rect(alpha=0.1) +
				scale_fill_manual(name="Sex", values=my.col[["Gender"]]) +
				scale_color_manual(name="Sex", values=my.col[["Gender"]]) +
				labs(y="Average Coverage") +
				theme_Publication() + 
				theme(legend.position="none") + 
				facet_grid(Gender~ .)

	p.rna.cov2<-ggplot(dt.rna.cov[score!=0]) + 
				geom_rect(aes(xmin=start,xmax=end,ymin=0,ymax=score, col=Gender, fill=Gender),alpha=0.1) +
				scale_fill_manual(name="Sex", values=my.col[["Gender"]]) +
				scale_color_manual(name="Sex", values=my.col[["Gender"]]) +
				labs(y="Average Coverage") +
				theme_Publication() + 
				theme(legend.position="none") + 
				facet_grid(Gender~ .)

	my.file.name<-"~/results/CSMD1.2016.paper/Figures/Fig.CSMD1.methy.exp.profile.tiff"
	tiff(filename=my.file.name,width=15, height=8.5 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(ggbio::tracks(`CSMD1`=p.tr.reduced, `WGoxBS (n=4)`=p.pt.meth, `RNA-Seq (n=131)`=p.rna.cov, heights = c(1.2,3.5,7), scale.height=2.5, label.text.cex = 1.3, xlim=gr.csmd.ext))
	dev.off()

	my.file.name<-"~/results/CSMD1.2016.paper/Figures/Fig.CSMD1.methy.exp.profile.4800.1500.tiff"
	tiff(filename=my.file.name,width=4800, height=1500 ,units="px",res=300, compression = 'lzw') #A4 size 
	print(ggbio::tracks(`Exon`=p.tr.reduced, `WGoxBS (n=4)`=p.pt.meth, `RNA-Seq (n=131)`=p.rna.cov, heights = c(1.2,3.5,4.6), scale.height=2.5, label.text.cex = 1.3, xlim=gr.csmd.ext))
	dev.off()

