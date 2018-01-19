	# part of ~/Pipelines/bin/R/CSMD1.Paper/CSMD1.paper.figures.R

	library(ggbio)
	source("~/Pipelines/bin/R/RnBeads/local.R") # defines 'target.region', 'prep,wgb' and 'prep.sureselect'
	################################
	## Figure 4.                  ##
	## CSMD1 DMR + 2 exons        ## 
	## Methylation Profiles       ##
	################################
	# R/3.1.1 (this ggbio (autoplot specifically) conflicts with ggplot which may be >2.0)
	# use R/3.2.1 (this version of ggbio works, but color/fill does does not work :()
	## see bin/R/RnBeads/SureSelect.ggplot.R
	my.project="SureSelect"
	my.slx<-"SLX-10409.v3"

	my.upflank=100 # up from TSS
	my.downflank=100 # down form TES

	min.missing.sample.ratio <- 0   # similar with filtering.missing.value.quantile of RnBeads

	wgbs.locus<-prep.wgbs()
	sureselect.v3.set1<-prep.sureselect("SureSelect.v3","Set1") # technical replicates
	sureselect.v3.set2<-prep.sureselect("SureSelect.v3","Set2") # biological replicates

	my.target<-"CSMD1.DMR2" # CSMD1.DMR, CSMD1.legacy.promo, CSMD1.putative.promo
	my.contrast<-"Gender"
	mycol<-my.col[[my.contrast]]
	my.gr.ensg <- rtracklayer::import.bed(target.region[[my.target]])
	my.gr.ensg <- reduce(my.gr.ensg)
	names(my.gr.ensg)<-my.target
	my.gr.ext.ensg<-promoters(my.gr.ensg, upstream=my.upflank, downstream=end(my.gr.ensg)-start(my.gr.ensg)+1+my.downflank) # includes 2000 +TSS and 100 +TES 

	# overlap between cpgi and my.target 
	#cpg.obj=methylKit::read.feature.flank(location=my.cgi.file, feature.flank.name=c("CpGi","shores")) # a ‘GenomicRangesList’ contatining one GRanges object for flanks and one for GRanges object for the main feature.
	#my.overlap=as.matrix( findOverlaps(cpg.obj[["CpGi"]], my.gr.ext.ensg) )
	#my.cpgi <-cpg.obj$CpGi[my.overlap[,"queryHits"]]

	##############
	# chromosome #
	##############
	#p.chr <- Ideogram(genome = "hg19", xlabel=TRUE) + xlim(ranges(my.gr.ext.ensg)) + theme_Publication() 

	############################
	# cpgi track for my.target #
	############################
	#p.cpgi <-autoplot(my.cpgi, fill = "darkgreen", color = "grey") + theme_Publication() 

	####################
	# transcript track #
	####################
	#p.tr.reduced<-autoplot(hg19ensGene, which=gr.csmd, stat="reduce", color="brown", fill="brown") + theme_alignment() # hg19ensGene defined in 'Annotation.R'
	p.tr.reduced<-autoplot(hg19ensGene, which=my.gr.ensg, stat="reduce", color="brown", fill="brown") + theme_alignment() # hg19ensGene defined in 'Annotation.R'

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
	# merger of two above
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
	p.pt.meth<-tracks(`CSMD1`=p.tr.reduced, `WGoxBS (n=4)`=p.wgbs, `In-solution target capture`=p.validation, heights = c(1.2,4,7), scale.height=2.5, label.text.cex = 1.4, xlim=my.gr.ensg)
	my.file.name<-"~/results/CSMD1.2016.paper/Figures/Fig4.CSMD1.DMR+2exons.methylation.profiles.POPS.ver1.tiff"
	tiff(filename=my.file.name,width=15, height=8.5 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(p.pt.meth)
	dev.off()

	#########################################################
	# Supplementary Figure
	# PLOT WGBS, SureSelect and RoadMap methylation profile #
	#########################################################
	# PT WGBS
	my.point.size=.8; my.alpha=0.4
	p.wgbs2<-ggplot(dt.target[type=="WGBS"], aes(start, score, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) +
				ylim(0,100) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				labs(y="") +
				theme_Publication() + 
				theme(legend.position="top")
	# RoadMap
	subject.keys=c("seqnames","start","end")
	query.key=c("V1","V2","End")
	subject<-as.data.frame(my.gr.ensg)
	cat("\tRemoving leading 'chr'...\n")
	#subject$seqnames<-simplify2array(strsplit(as.character(subject$seqnames), "chr"))[2,] # chrX = > X
	subject$seqnames<-substr(subject$seqnames, 4, 5) #chrX=>X
	cat("\tconvert to data.table...\n")
	dt.subject<-as.data.table(subject) # isa data.table
	setkeyv(dt.subject,subject.keys) # regions of interests 
	rm(subject)
	cat("Subject loaded\n")

	## Load RoadMap Tissue
	dt.tissue.locus<-list()
	for(my.tissue in c(avail.tissues[!avail.tissues %in% "PA"])){ # except 'PT': placenta
		## Load this RoadMap tissue
		dt.roadmap=load.my.tissue.dt.merged(my.tissue,my.cpg.type)

		cat("\tFinding CpGs overlap...\n")
		system.time(dt.overlap<-foverlaps(dt.roadmap, dt.subject, by.x=query.key, type="any", nomatch=0L))
		#> dt.overlap
		#   V1    Start      End Strand       V2 V3  V4 V7 V5.x V6.x V5.y V6.y    i.End
		#1:  X  2746554  2748235      *  2746590  - CGA  1    6   17    2   10  2746590
		#2:  X  2835950  2836236      *  2836235  + CGG  1   14   14    9   10  2836235
		if(nrow(dt.overlap)){
			dt.meth<-rbind(dt.overlap[,.(V1,V2,V5.x/V6.x*100,"Female",my.tissue)], dt.overlap[,.(V1,V2,V5.y/V6.y*100,"Male",my.tissue)]) # isa 'data.table'
			setnames(dt.meth,c("chr","start","methylation_level","Gender","Tissue"))
			dt.tissue.locus[[my.tissue]]<-dt.meth
		}
	}# end of my.tissue 
	# Gender-specific
	p.roadmap<-ggplot(rbindlist(dt.tissue.locus), aes(start, methylation_level, col=Gender)) + geom_point(size=my.point.size,alpha=my.alpha) + ylim(0, 100) + 
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				labs(y="% Methylation") +
				theme_Publication() + 
				theme(legend.position="none") + 
				facet_grid(Tissue ~ .)

	# wgbs, tech, bio, roadmap
	my.file.name<-"~/results/CSMD1.2016.paper/Figures/Suppl/SupplFig.CSMD1.DMR+2exons.methylation.profiles.RoadMap.tiff"
	tiff(filename=my.file.name,width=10, height=15 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(tracks(` `=p.tr.reduced, `Placenta`=p.wgbs2, `Schultz et al`=p.roadmap, heights = c(1,4,28), scale.height=2.5, label.text.cex = 1.5, xlim=my.gr.ensg))
	dev.off()

	################################
	## Figure 4-ver2.             ##
	## CSMD1 Methylation Profiles ##
	## 2 RoadMap tissues          ##
	################################
	my.point.size=1.2; my.alpha=0.4
	p.placenta3<-ggplot(dt.target, aes(start, score, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) +
				ylim(0,100) +
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				labs(y="% Methylation") +
				theme_Publication() + 
				theme(legend.position=c(0.85,0.12)) +
				facet_grid(type~ .)
	# RoadMap 2 tissues
	p.roadmap2<-ggplot(rbindlist(list(`FT`=dt.tissue.locus[["FT"]],`SX`=dt.tissue.locus[["SX"]])), aes(start, methylation_level, col=Gender)) + 
				geom_point(size=my.point.size,alpha=my.alpha) + 
				ylim(0, 100) + 
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				labs(y="") +
				theme_Publication() + 
				theme(legend.position="none") + 
				facet_grid(Tissue ~ .)
	

	my.file.name<-"~/results/CSMD1.2016.paper/Figures/Fig4.CSMD1.DMR+2exons.methylation.profiles.POPS.ver2.tiff"
	tiff(filename=my.file.name,width=15, height=8.5 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(tracks(` `=p.tr.reduced, `Placenta`=p.placenta3, `Schultz et al`=p.roadmap2, heights = c(1,10.5,7), scale.height=2.5, label.text.cex = 1.5, xlim=my.gr.ext.ensg))
	dev.off()
