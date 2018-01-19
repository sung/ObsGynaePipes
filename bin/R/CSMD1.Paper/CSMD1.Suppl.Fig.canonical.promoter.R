# part of ~/Pipelines/bin/R/CSMD1.Paper/CSMD1.paper.figures.main.R

	library(ggbio)
	source("~/Pipelines/bin/R/RnBeads/local.R") # defines 'target.region', 'prep,wgb' and 'prep.sureselect'
	################################
	## Figure 5.                  ##
	## CSMD1 Putative Promoter    ##
	## Methylation Profiles       ##
	################################
	# R/3.1.1 (this ggbio (autoplot specifically) conflicts with ggplot which may be >2.0)
	# use R/3.2.1 (this version of ggbio works, but color/fill does does not work :()
	## see bin/R/RnBeads/SureSelect.ggplot.R
	my.project="SureSelect"
	my.slx<-"SLX-10409.v3"

	my.upflank=50 # up from TSS
	my.downflank=50 # down form TES

	min.missing.sample.ratio <- 0   # similar with filtering.missing.value.quantile of RnBeads

	wgbs.locus<-prep.wgbs()
	sureselect.v3.set1<-prep.sureselect("SureSelect.v3","Set1") # technical replicates
	sureselect.v3.set2<-prep.sureselect("SureSelect.v3","Set2") # biological replicates

	my.target<-"CSMD1.legacy.promo" # CSMD1.DMR, CSMD1.legacy.promo, CSMD1.putative.promo
	my.contrast<-"Gender"
	mycol<-my.col[[my.contrast]]
	my.gr.ensg <- rtracklayer::import.bed(target.region[[my.target]]) # TSS (3,040,668) +2K and -2K
	my.gr.ensg <- reduce(my.gr.ensg)
	names(my.gr.ensg)<-my.target
	my.gr.ext.ensg<-promoters(my.gr.ensg, upstream=my.upflank, downstream=end(my.gr.ensg)-start(my.gr.ensg)+1+my.downflank) # includes 2000 +TSS and 100 +TES 

	####################
	# transcript track #
	####################
	#p.tr.reduced<-autoplot(hg19ensGene, which=my.gr.ext.ensg, stat="reduce", color="brown", fill="brown") + theme_alignment() # does not work for uknown reason
	p.tr.reduced<-autoplot(hg19ensGene, which=my.gr.ext.ensg, color="brown", fill="brown") + theme_alignment() # hg19ensGene defined in 'Annotation.R'

	##########################
	# WGBS methylation GGBIO #
	##########################
	select<-mcols(wgbs.locus[["gl"]][[my.target]])$Condition=="AGA"; my.gr.wgbs<-as.data.frame(wgbs.locus[["gl"]][[my.target]][select,])
	my.gr.wgbs$type='WGoxBS (N=4)' # isa 'data.frame'

	################ 
	## Validation ##
	################ 
	# PT technical validation
	filter.tech<-start(sureselect.v3.set1[["gl"]][[my.target]]) %in% sureselect.v3.set1[["gender"]][[my.target]][["AGA"]] & mcols(sureselect.v3.set1[["gl"]][[my.target]])$Condition=="AGA"
	my.gr.tech<-as.data.frame(sureselect.v3.set1[["gl"]][[my.target]][filter.tech,])
	my.gr.tech$type="Technical (n=2)"
	# PT biological validation
	filter.bio<-start(sureselect.v3.set2[["gl"]][[my.target]]) %in% sureselect.v3.set2[["gender"]][[my.target]][["AGA"]] & mcols(sureselect.v3.set2[["gl"]][[my.target]])$Condition=="AGA"
	my.gr.bio<-as.data.frame(sureselect.v3.set2[["gl"]][[my.target]][filter.bio,])
	my.gr.bio$type="Biological (n=6)"
	# merger of the two above
	df.target<-rbind(my.gr.tech, my.gr.bio); 
	df.target$type<-factor(df.target$type,levels = c("Technical (n=2)", "Biological (n=6)")) # re-order this facet
	dt.target<-as.data.table(df.target)

	##############################################
	# PLOT WGBS & SureSelect methylation profile #
	##############################################
	my.point.size=1.7; my.alpha=0.5
	p.wgbs<-ggplot(my.gr.wgbs, aes(start, score, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) +
				ylim(0,100) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="", values=my.col[["Gender"]]) + 
				labs(y="% Methylation") +
				theme_Publication() + 
				theme(legend.position="top")
	p.validation<-ggplot(dt.target[type=="Technical (n=2)"|type=="Biological (n=6)"], aes(start, score, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) +
				ylim(0,100) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="", values=my.col[["Gender"]]) + 
				labs(y="% Methylation") +
				theme_Publication() + 
				#theme(legend.position=c(0.85,0.6)) +
				theme(legend.position="none") +
				facet_grid(type~ .)
	#p.pt.meth<-tracks(`CSMD1`=p.tr.reduced, `WGoxBS (n=4)`=p.wgbs, `In-solution target capture`=p.validation, heights = c(1.2,4,7), scale.height=2.5, label.text.cex = 1.4, xlim=resize(my.gr.ensg, fix="center",1700))
	p.pt.meth<-tracks(`CSMD1`=p.tr.reduced, `WGoxBS (n=4)`=p.wgbs, `In-solution target capture`=p.validation, heights = c(3,4,7), scale.height=2.5, label.text.cex = 1.4, xlim=my.gr.ensg)
	my.file.name<-"~/results/CSMD1.2016.paper/Figures/Suppl/Suppl.Fig.CSMD1.canonical.promoter.POPS.tiff"
	tiff(filename=my.file.name,width=15, height=8.5 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(p.pt.meth)
	dev.off()
