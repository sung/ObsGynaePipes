	# part of ~/Pipelines/bin/R/CSMD1.Paper/CSMD1.paper.figures.R
	library(data.table)

	#my.exon.RData<-paste0("~/results/CSMD1.2016.paper/RData/CSMD1.",TR_PREFIX,".POPS.Schultz.RoadMap.exon.RData")
	my.exon.RData<-paste0("~/results/CSMD1.2016.paper/RData/CSMD1.",TR_PREFIX,".POPS.RoadMap.Fetal.exon.RData")
	if(file.exists(my.exon.RData)){
		load(my.exon.RData)
	}else{
		if(FALSE){
			## Schultz data
			load("~/results/RNA-Seq/RoadMap.exon.GRCh37/DESeq2/Schultz/deseq.Schultz.RData")
			Schultz.count=counts(dds)
			Schultz.ddsFpkm=fpkm(dds)
			Schultz.colData=as.data.frame(colData(dds)) # isa 'data.frame'

			## RoadMap Consortiumn Placenta data
			load("~/results/RNA-Seq/RoadMap.exon.GRCh37/DESeq2/Consortium/deseq.Consortium.RData")
			RoadMap.counts=counts(dds)
			RoadMap.ddsFpkm=fpkm(dds)

			## RoadMap Consortiumn Fetal data
			load("~/results/RNA-Seq/RoadMap.exon.GRCh37/DESeq2/FetalEmbryo/deseq.FetalEmbryo.RData")
			RoadMap.counts=counts(dds)
			RoadMap.ddsFpkm=fpkm(dds)
			RoadMap.colData=as.data.frame(colData(dds)) # isa 'data.frame'
		}

		## RoadMap Schultz + Consortiumn
		load("~/results/RNA-Seq/RoadMap.exon.GRCh37/DESeq2/ALL/deseq.ALL.RData")
		RoadMap.counts=counts(dds) # isa 'matrix'
		RoadMap.ddsFpkm=fpkm(dds) # isa 'matrix'
		RoadMap.colData=as.data.frame(colData(dds)) # isa 'data.frame'

		## POPS (Boy vs. Girl)
		my.exon.POPS.RData<-paste0("~/results/CSMD1.2016.paper/RData/CSMD1.",TR_PREFIX,".POPS.boy.girl.exon.RData")
		if(file.exists(my.exon.POPS.RData)){
			load(my.exon.POPS.RData) # dt.res, rn. rn.control, rn.case, ddsFpkm, ddsFpm
		}else{
			## at this point, R version should be the same when deseq.RData was compiled.
			## probably 3.1
			deseq.RData <- file.path("~/results/RNA-Seq",paste0("Boy.Girl.exon.",TR_PREFIX,"/DESeq2/AGA/deseq.AGA.RData"))
			load (deseq.RData) # 

			my.contrast="Sex"
			dt.res<-as.data.frame(res); dt.res$ensg<-rownames(dt.res); dt.res<-as.data.table(dt.res)
			rn=rownames(colData(dds)) #as.character(samples[["SampleName"]]) # samples$SampleName (or samples[,c("SampleName")])
			rn.control=rownames(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1==0,]) 
			rn.case=rownames(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1!=0,])
			ddsFpkm <-fpkm(dds) # isa 'matrix'
			ddsFpm <-fpm(dds) # isa 'matrix'
			save(dt.res, rn, rn.control, rn.case, ddsFpkm, ddsFpm, file=my.exon.POPS.RData) 
		}
		POPS.ddsFpkm=ddsFpkm

		}
		#save(POPS.ddsFpkm, Schultz.ddsFpkm, RoadMap.ddsFpkm, file=my.exon.RData) 
		save(Fetal.colData, Fetal.ddsFpkm, POPS.ddsFpkm, RoadMap.counts, RoadMap.ddsFpkm, RoadMap.colData, file=my.exon.RData) 
	}

	#########################
	# tissue-level pheatmap #
	#########################
	library(pheatmap)
	# Schultz.ddsFpkm is 'ddsFpkm' of RoadMap Schultz tissues (no placenta)
	# RoadMap.ddsFpkm is 'ddsFpkm' of RoadMap Consortium placenta data
	# POPS.ddsFpkm is 'ddsFpkm' of Boy.Girl.exon.GRCh37 AGA samples
	#merged.fpkm<-cbind(Schultz.ddsFpkm, RoadMap.ddsFpkm, `Placenta`=rowMeans(POPS.ddsFpkm)) # isa 'matrix'
	#merged.fpkm[rev(rownames(merged.fpkm)),] # 5' exon to 3' exon
	merged.fpkm<-cbind(RoadMap.ddsFpkm, `Placenta`=rowMeans(POPS.ddsFpkm)) # isa 'matrix'
	colnames(merged.fpkm)<-gsub("_"," ",colnames(merged.fpkm))

	ann_tissues<-data.frame(
					#`Tissue`=factor(c(as.character(RoadMap.colData$Origin), "Placenta-derived")),
					`Tissue`=factor(c(ifelse(grepl("Brain",rownames(RoadMap.colData)), "Brain", ifelse(RoadMap.colData$Origin=="Non-Placenta","Other-somatic",as.character(RoadMap.colData$Origin))),"Placenta-derived")),
					`Source`=factor(c(ifelse(RoadMap.colData$Source=="Schultz","Schultz et al.",as.character(RoadMap.colData$Source)), "This Study")),
					row.names=colnames(merged.fpkm)
				)
	ann_exons<-data.frame(
					`TSS`=ifelse(rownames(merged.fpkm)=="ENSE00002108991","Canonical TSS",ifelse(rownames(merged.fpkm)=="ENSE00002127078","Placental TSS",NA)),
					row.names=rownames(merged.fpkm)
						  )
	ann_colors<-list(
					 `Tissue`=c(`Other-somatic`=brewer.pal(3, "Dark2")[1], `Placenta-derived`=brewer.pal(3, "Dark2")[2], `Brain`=brewer.pal(3, "Dark2")[3]),
					 `Source`=c(`Schultz et al.`="grey0", `UCSF-UBC`=cbPalette2[2], `BI`=cbPalette2[5], `This Study`=cbPalette2[3]),
					 `TSS`=c(`Canonical TSS`="#377EB8",`Placental TSS`="#E41A1C")
				)


	# portrait (row; tissues, col: exons)
	my.filename<-"~/results/CSMD1.2016.paper/Figures/Pheatmap.exon/Fig.pheatmap.exon.pops.roadmap.GRCh37.portrait.tiff"
	pheatmap(log(merged.fpkm[rev(rownames(merged.fpkm)),]/10^3+1), cluster_rows=F, cluster_cols=T, annotation_col=ann_tissues, annotation_colors=ann_colors, width=10, height=15, fontsize_col=10, fontsize_row=7, cellheight=9, gaps_col=T, cutree_col=2, filename=my.filename)

	######################
	## Fetal and Embyro ##
	######################
	ann_colors<-list(`Tissue`=c(`Placenta-derived`='darkmagenta', `Embryo-Somatic`='yellow4', `Fetal-Somatic`='yellowgreen'),
					 `Source`=c(`RoadMap Consortium`=cbPalette2[2], `This Study`=cbPalette2[3], `Gerrad et al.`=cbPalette2[4]))
	my.filename<-"~/results/CSMD1.2016.paper/Figures/Pheatmap.exon/Fig.pheatmap.exon.pops.roadmap.fetal.GRCh37.portrait.tiff"
	pheatmap(log(merged.fpkm[rev(rownames(merged.fpkm)),]/10^3+1), cluster_rows=F, cluster_cols=T, annotation_col=ann_tissues, annotation_colors=ann_colors, width=10, height=15, fontsize_col=10, fontsize_row=7, cellheight=9, filename=my.filename)

	# horizontal (row: exons, col: tissues)
	my.filename<-"~/results/CSMD1.2016.paper/Figures/Pheatmap.exon/Fig.pheatmap.exon.pops.roadmap.GRCh37.portrait.new.TSS.tiff"
	#pheatmap(t(log(merged.fpkm/10^3+1)), cluster_rows=T, cluster_cols=F, annotation_row=ann_tissues, annotation_colors=ann_colors, width=16, height=5, fontsize_col=7, fontsize_row=10, cellheight=10, cutree_row=2, filename=my.filename)
	pheatmap(t(log(merged.fpkm/10^3+1)), cluster_rows=T, cluster_cols=F, annotation_row=ann_tissues, annotation_colors=ann_colors, width=16, height=5, fontsize_col=7, fontsize_row=10, cellheight=10, filename=my.filename)

