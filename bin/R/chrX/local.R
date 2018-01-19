library(data.table)
#options(scipen=999) # diable scientific notation
##########
# Config #
##########
time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p') #2015-01-07_01_30PM

#########
## PAR ##
#########
gr.par<-GenomicRanges::GRanges(Rle(c("chrX","chrX")),IRanges(start=c(60001,154931044),end=c(2699520,155260560),names=c("PAR1","PAR2")))
genome(gr.par)<-"hg19"

#############
## miR  #####
#############
# Promoter (PROmiRNA)
gr.mir.pro<-rtracklayer::import.bed("~/data/Annotation/PROmiRNA/PROmiRNA.grch37.bed")
# Gene-body (GRCh37)
gr.mir.gene<-rtracklayer::import.gff("~/data/Annotation/miRbase.20/hsa.gff3")


##########################################
## Genes mentioned by Barbara R. Migeon ##
##########################################
bm.gene=c(`G6PD`="ENSG00000160211",`HPRT1`="ENSG00000165704",`PGK1`="ENSG00000102144",`TIMP1`="ENSG00000102265")

if(RNASEQ_TYPE=="TOTAL"){
	cpg.contexts=c("Genes", "Promo_10.10", "CPGi")
	my.region=list(`Promoter`="Promo_10.10", `CPGi`="CPGi", `Gene-body`="Genes", `CPGi-shore`="CPGi-shores", `Enhancer`="Enhancer")
	my.region.rev=list(`Promo_10.10`="Promoter", `CPGi`="CPGi", `Genes`="Gene-body", `CPGi-shores`="CPGi-shore", `Enhancer`="Enhancer")
	min.num.sites<- 4
	my.id="ensembl_gene_id"
	my.formula<-formula(tissue + ensembl_gene_id ~ region)
	my.abline=c(log10(0.01),0,-log10(0.01))
	my.beta=1.e-100
	my.label="A"

	my.pval=c(0.5,0.1,0.01)
	my.exp.sig<-c(`all chrX genes`='ALL', `rarely expressed`='LOW_EXP', `adjusted p-val<0.5`='P<0.5', `adjusted p-val<0.1`='P<0.1', `adjusted p-val<0.01`='P<0.01')
}else if(RNASEQ_TYPE=="SMALL"){
	cpg.contexts=c("miR_Genes_01.01", "miR_Promo_01.01", "CPGi")
	my.region=list(`miR Promoter`="miR_Promo_01.01", `miR CPGi`="CPGi", `miR Gene-body`="miR_Genes_01.01")
	my.region.rev=list(`miR_Promo_10.10`="Promoter", `CPGi`="CPGi", `miR_Genes_01.01`="Gene-body")
	min.num.sites<- 1 
	my.id="mirbase_id"
	my.formula<-formula(tissue + mirbase_id ~ region)
	my.abline=c(log10(0.1),0,-log10(0.1))
	my.beta=0
	my.label="C"

	my.pval=c(0.6,0.4,0.1)
	my.exp.sig<-c(`all chrX genes`='ALL', `rarely expressed`='LOW_EXP', `adjusted p-val<0.6`='P<0.6', `adjusted p-val<0.4`='P<0.4', `adjusted p-val<0.1`='P<0.1')
}else{
	stop(paste(RNASEQ_TYPE, "Not supported\nEither TOTAL or SMALL"))
}

##################
## Balaton 2015 ##
##################
if(!exists("dt.balaton")){
	dt.balaton<-fread("~/data/X-inactivation/Balaton.2015/Balaton.2015.csv", na.strings="")
	dt.balaton[,XCI:=ifelse(grepl("PAR",Y_homology),"PAR", ifelse(is.na(Balaton)|Balaton=="No call",'Unknown',ifelse(grepl("VE$",Balaton),"Variable", ifelse(grepl("E$",Balaton),"Escaped",ifelse(grepl("S$",Balaton),"Inactivated",Balaton)))))]
	dt.balaton[hgnc_symbol=="XIST",XCI:="Escaped"]
	dt.balaton[XCI=="Discordant", XCI:=ifelse(is.na(Carrel), XCI , ifelse(Carrel==1,"Escaped",ifelse(Carrel==0,"Inactivated","Variable")))] # Classify Discondordant case
}

############################
## PT WGoxBS & SureSelect ##
############################
meth.type=list(`WGoxBS`="PT", `SS.Tech`="PT.SS.Tech", `SS.Bio`="PT.SS.Bio")
meth.type.rev=list(`PT`="WGoxBS", `PT.SS.Tech`="SS.Tech", `PT.SS.Bio`="SS.Bio")
#load("/home/ssg29/data/RoadMap/BS-Seq/Schultz2015/PT.CG.dt.merged.RData")
#cat("PT.CG.dt.merged loaded\n") # PT WGoxBS (merged by sex)

#load("/home/ssg29/data/RoadMap/BS-Seq/Schultz2015/PT.CG.dt.merged.sureselect.RData")
#cat("dl.sureselect loaded\n") # PT SureSelect (tech and bio)

# R object will be saved into this file
my.chrX.RData<-paste0("~/results/chrX/RData/pops.meth.exp.chrX.",RNASEQ_TYPE,".RData")
if(file.exists(my.chrX.RData)){
	load(my.chrX.RData)
	cat("dt.meth.chrX, dt.exp.pops, dt.exp.pops.chrX, dt.exp.miRNA.pops, dt.exp.miRNA.pops.chrX, dt.deg.chrX, dt.meth.exp, dt.meth.exp.hr, dt.meth.diff, dt.meth.diff.exp, dt.meth.exp.deg, dt.balaton.deg? loaded\n")
}else{
	if(RNASEQ_TYPE=="TOTAL"){
		# CG methylation RoadMap tissue and PT at chrX 
		# see ~/Pipelines/bin/R/RoadMap/local.R
		#load("~/results/RoadMap/X-inactivation/RData/10X/dt.meth.chrX.per.tissue.CG.Genes_20.20.WGoxBS.RData") # WGoxBS only (PT)
		load("~/results/RoadMap/X-inactivation/RData/10X/dt.meth.chrX.per.tissue.CG.Genes_20.20.RData") # WGoxBS and SureSelect (PT, PT.SS.Tech, PT.SS.Bio)
		## dt.meth.all.chrX: Note that ensembl_gene_id can be null <NA> because of left-join (all.x=TRUE)
		## ensembl_gene_id from RnBeads (Ensembl 73 or gencode 18)
		dt.meth.chrX=rbindlist(dt.meth.chrX.per.tissue)[grepl("PT",tissue) & !is.na(ensembl_gene_id),
															.(f.met=sum(female.c)/sum(female.cov)*100, m.met=sum(male.c)/sum(male.cov)*100, diff.met=sum(female.c)/sum(female.cov)*100-sum(male.c)/sum(male.cov)*100, num.sites=sum(num.sites)),
															by="tissue,ensembl_gene_id,region"]
		############################
		## get PT expression      ##
		############################
		# 1. FG.JD.AGA (GRCh37)
		library(DESeq2) # to load DESeq2
		deseq.RData <- "/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/DESeq2/AGA.withoutChrY/deseq.AGA.RData"
		load(deseq.RData) 

		my.res<-res
		my.res$ensembl_gene_id=rownames(my.res)
		dt.exp.pops=as.data.table(as.data.frame(my.res))
		dt.chr<-as.data.table(getBM(attributes = c("ensembl_gene_id", "chromosome_name"), filters = "ensembl_gene_id", values = dt.exp.pops[,unique(ensembl_gene_id)], mart = myMart)) # is a data.frame
		dt.exp.pops=merge(dt.exp.pops, dt.chr, by="ensembl_gene_id")
		dt.exp.pops.chrX<-dt.exp.pops[chromosome_name=="X"][,`new.padj` := p.adjust(pvalue, method="BH")]
	}else if(RNASEQ_TYPE=="SMALL"){
		# CG methylation RoadMap tissue and PT at chrX of miR
		# see ~/Pipelines/bin/R/RoadMap/local.R
		load("~/results/RoadMap/X-inactivation/RData/dt.meth.chrX.mir.per.tissue.CG.miR_Genes_10.10.RData")
		dt.meth.chrX=rbindlist(dt.meth.chrX.mir.per.tissue)[grepl("PT",tissue) & !is.na(mirbase_id),
															.(f.met=sum(female.c)/sum(female.cov)*100, m.met=sum(male.c)/sum(male.cov)*100, diff.met=sum(female.c)/sum(female.cov)*100-sum(male.c)/sum(male.cov)*100, num.sites=sum(num.sites)),
															by="tissue,mirbase_id,region"]
		############################
		## get PT expression      ##
		############################
		# 2. FG.JD.pre.miRNA (GRCh37)
		library(DESeq2) # to load DESeq2
		deseq.RData <- "/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.pre.miRNA.GRCh38/DESeq2/AGA/deseq.AGA.RData"
		load(deseq.RData)

		my.res<-res
		my.res$mirbase_id=rownames(my.res)
		dt.exp.pops=as.data.table(as.data.frame(my.res))

		# get chr infor of mirbase_id from the biomart
		grch38 = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset = "hsapiens_gene_ensembl")
		dt.chr<-as.data.table(getBM(attributes =c("mirbase_id", "ensembl_gene_id","chromosome_name"), filters = "mirbase_id", values = dt.exp.pops[,unique(mirbase_id)], mart = grch38)) # is a data.frame
		dt.exp.pops=merge(dt.exp.pops, dt.chr, by="mirbase_id", all.x=TRUE, allow.cartesian=TRUE)
		dt.exp.pops.chrX<-dt.exp.pops[chromosome_name=="X"]
		dt.exp.pops.chrX[,`new.padj` := p.adjust(pvalue, method="BH")]
	}else{
		stop(paste(RNASEQ_TYPE, "Not supported\nEither TOTAL or SMALL"))
	}

	#####################
	## Merge Met & Exp ##
	#####################
	## ensembl_gene_id from RnBeads (Ensembl 73 or gencode 18)
	dt.meth.exp=merge(dt.meth.chrX, dt.exp.pops.chrX, by=my.id)

	dt.meth.diff<-as.data.table(reshape2::dcast(dt.meth.chrX[num.sites>=min.num.sites], my.formula, value.var="diff.met")) # dcast: long(row-based) to wide(column-based)
	dt.meth.female<-as.data.table(reshape2::dcast(dt.meth.chrX[num.sites>=min.num.sites], my.formula,  value.var="f.met")) # dcast: long(row-based) to wide(column-based)
	dt.meth.male<-as.data.table(reshape2::dcast(dt.meth.chrX[num.sites>=min.num.sites], my.formula,  value.var="m.met")) # dcast: long(row-based) to wide(column-based)
	dt.meth.exp.hr<-merge(merge(dt.meth.female, dt.meth.male, by=c("tissue",my.id)), dt.exp.pops.chrX, by=my.id)

	## dt.meth.diff[,ensembl_gene_id] from RnBeads (Ensembl 73 or gencode 18)
	## dt.exp.pops[,ensembl_gene_id] from DESeq2 (Ensembl 75 or gencode 19?)
	dt.meth.diff.exp<-merge(dt.meth.diff, dt.exp.pops.chrX, by=my.id)
	dt.meth.diff.exp[, `Expression?` := ifelse(log2FoldChange>0, 'female>male', 'female<male')]

	if(RNASEQ_TYPE=="TOTAL"){
		# one-to-many: (1-ensg-to-2-hgnc)
		fields <- c("ensembl_gene_id", "hgnc_symbol","gene_biotype", "description")
		dt.gene<-as.data.table(getBM(attributes = fields, filters = my.id, values = dt.meth.exp.hr[new.padj<0.01, unique(ensembl_gene_id)], mart = myMart)) # is a data.table
		dt.meth.exp.deg<-merge(dt.meth.exp.hr[new.padj<0.01], dt.gene, by=my.id)
		dt.meth.exp.deg[, `Expression?` := ifelse(log2FoldChange>0, 'female>male', 'female<male')]

		dt.gene<-as.data.table(getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.exp.pops.chrX[new.padj<0.01, unique(ensembl_gene_id)], mart = myMart)) # is a data.table
		dt.deg.chrX<-merge(dt.exp.pops.chrX[new.padj<0.01], dt.gene, by="ensembl_gene_id")

		dt.balaton.deg<-merge(dt.deg.chrX, dt.balaton, by="hgnc_symbol", all.x=TRUE)[,list(hgnc_symbol,ensembl_gene_id,gene_biotype,baseMean,log2FoldChange,new.padj,Balaton,Carrel,other_literature,Y_homology,description)]
		#dt.balaton.deg[,XCI:=ifelse(grepl("PAR",Y_homology),"PAR", ifelse(is.na(Balaton)|Balaton=="No call",'Unknown',ifelse(grepl("VE$",Balaton),"Variable", ifelse(grepl("E$",Balaton),"Escaped",ifelse(grepl("S$",Balaton),"Inactivated",Balaton)))))]
		#dt.balaton.deg[hgnc_symbol=="XIST",XCI:="Escaped"]
		#> dt.balaton.deg[XCI=="Discordant",list(hgnc_symbol,ensembl_gene_id,Balaton,Carrel)]
		#   hgnc_symbol ensembl_gene_id    Balaton Carrel
		#1:         CHM ENSG00000188419 Discordant   0.78
		#2:        PIN4 ENSG00000102309 Discordant   0.33
		#3:     ZC3H12B ENSG00000102053 Discordant   0.11
		#dt.balaton.deg[XCI=="Discordant", XCI:=ifelse(Carrel==1,"Escaped",ifelse(Carrel==0,"Inactivated","Variable"))] # Classify Discondordant case
		dt.balaton.deg[, `Expression?` := ifelse(log2FoldChange>0, 'female>male', 'female<male')]

		write.csv(dt.balaton.deg[order(new.padj)], file=paste0("~/results/chrX/CSV/dt.deg.chrX.balaton.",RNASEQ_TYPE,".csv"))
		# save it
		save(dt.meth.chrX, dt.exp.pops, dt.exp.pops.chrX, dt.exp.miRNA.pops, dt.exp.miRNA.pops.chrX, dt.deg.chrX, dt.meth.exp, dt.meth.exp.hr, dt.meth.diff, dt.meth.diff.exp, dt.meth.exp.deg, dt.balaton.deg, file=my.chrX.RData)
	}else if(RNASEQ_TYPE=="SMALL"){
		dt.meth.exp.deg<-dt.meth.exp.hr[new.padj<0.1]
		# save it
		save(dt.meth.chrX, dt.exp.pops, dt.exp.pops.chrX, dt.exp.miRNA.pops, dt.exp.miRNA.pops.chrX, dt.deg.chrX, dt.meth.exp, dt.meth.exp.hr, dt.meth.diff, dt.meth.diff.exp, dt.meth.exp.deg, file=my.chrX.RData)
	}else{
		stop(paste(RNASEQ_TYPE, "Not supported\nEither TOTAL or SMALL"))
	}
	###########################
	## Save Meth & Exp Data  ##
	###########################
	# all chrX genes
	#write.csv(dt.exp.pops.chrX, file=paste0("~/results/chrX/CSV/dt.exp.pops.chrX.",RNASEQ_TYPE,".csv"))
	#write.csv(dt.meth.exp.deg, file=paste0("~/results/chrX/CSV/dt.meth.exp.deg.",RNASEQ_TYPE,".csv"))
	#write.csv(dt.meth.exp.hr, file=paste0("~/results/chrX/CSV/dt.meth.exp.deg.",RNASEQ_TYPE,".csv"))
	#write.csv(dt.meth.diff.exp, file=file.path("~/results/chrX/CSV/",paste("PT.meth.exp",RNASEQ_TYPE,"csv",sep=".")))
}

blankPlot <- ggplot()+geom_blank(aes(1,1))+
	theme(
		plot.background = element_blank(), 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		panel.border = element_blank(),
		panel.background = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.text.x = element_blank(), 
		axis.text.y = element_blank(),
		axis.ticks = element_blank(),
		axis.line = element_blank()
		)

################################################
# used by bin/R/chrX/plot.met.exp.by.sex.R   ###
# "dt.meth.diff.exp" must be pre-defined     ###
# dt.balaton must be pre-defined             ###
# min.num.site already applied >4            ###
################################################
plotMethExp<-function(METH_TYPE, min.depth=10, my.padj, my.target=c("CPGi","Gene-body")){
	## At least 10% difference and at least 10x coverage
	my.dt.meth.exp=dt.meth.diff.exp[tissue==meth.type[[METH_TYPE]] & !is.na(`Expression?`) & !is.na(eval(parse(text=my.region[[my.target[1]]]))) & !is.na(eval(parse(text=my.region[[my.target[2]]]))) ]

	if(my.padj<1){
		my.dt.meth.exp=my.dt.meth.exp[new.padj<=my.padj]
	}

	if(my.padj==0.01){
		# get hgnc_symbol - it could be one-to-many: (1-ensg-to-2-hgnc)
		my.dt.meth.exp=my.dt.meth.exp[abs(log2FoldChange)>log2(1.1) & baseMean>min.depth]
		my.dt.meth.exp<-merge(my.dt.meth.exp, getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = my.dt.meth.exp[,unique(ensembl_gene_id)], mart = myMart), by="ensembl_gene_id")
		# attach Balaton annotation
		my.dt.meth.exp<-merge(my.dt.meth.exp, dt.balaton, by="hgnc_symbol", all.x=TRUE)
		# classify XCI status based on Balaton et al.
		my.dt.meth.exp[,XCI:=ifelse(grepl("PAR",Y_homology),"PAR", ifelse(is.na(Balaton)|Balaton=="No call",'Unknown',ifelse(grepl("VE$",Balaton),"Variable", ifelse(grepl("E$",Balaton),"Escaped",ifelse(grepl("S$",Balaton),"Inactivated",Balaton)))))]
		my.dt.meth.exp[hgnc_symbol=="XIST",XCI:="Escaped"]
		my.dt.meth.exp[XCI=="Discordant", XCI:=ifelse(Carrel==1,"Escaped",ifelse(Carrel==0,"Inactivated","Variable"))] # Classify Discondordant case

		p <-ggplot(my.dt.meth.exp, aes_string(x=my.region[[my.target[1]]], y=my.region[[my.target[2]]], label="hgnc_symbol")) +
			geom_text(check_overlap = TRUE, hjust =1.2, nudge_x=0.05, size=4.5) +
			geom_point(aes(col=`Expression?`, shape=XCI),size=7,alpha=0.7)
		if(packageVersion("ggplot2")<="2.1.0"){
			p <-p + facet_wrap(~XCI,nrow=2)
		}
	}else{
		p <-ggplot(my.dt.meth.exp, aes_string(x=my.region[[my.target[1]]], y=my.region[[my.target[2]]])) + 
			geom_point(aes(col=`Expression?`, size=abs(log2FoldChange)),alpha=.8)
	}

	p <-p +
		labs(x=paste("% Methylation Difference (Female - Male)\n",my.target[1]),y=paste("% Methylation Difference (Female - Male)\n",my.target[2]) ) +
		geom_hline(yintercept=0) + geom_vline(xintercept=0) +
		scale_colour_manual(values=my.col[["Expression?"]]) +
		theme_Publication()

	if(packageVersion("ggplot2")<="2.1.0"){
		print(p + ggtitle(paste("Methylation Difference (padj<=",my.padj,")")) )
	}else{
		q<-ggplot(my.dt.meth.exp, aes_string("`Expression?`", my.region[[my.target[1]]])) + 
			geom_boxplot(aes(col=`Expression?`)) + 
			labs(y='') + 
			scale_x_discrete(labels=c('F<M','F>M')) + 
			scale_colour_manual(values=my.col[["Expression?"]]) +
			coord_flip() + 
			theme_Publication() + 
			theme(axis.ticks.y = element_blank(), legend.position="none" )
		r<-ggplot(my.dt.meth.exp, aes_string("`Expression?`", my.region[[my.target[2]]])) + 
			geom_boxplot(aes(col=`Expression?`)) + 
			labs(y='') + 
			scale_x_discrete(labels=c('F<M','F>M')) + 
			scale_colour_manual(values=my.col[["Expression?"]]) +
			theme_Publication() + 
			theme(legend.position="none")
		gridExtra::grid.arrange(r, p + theme(legend.position="none"), blankPlot, q, ncol=2, nrow=2, widths=c(1, 4), heights=c(4,1))
	}
}


##############################
## X-axis: % Meth of female ##
## Y-axis: % Meth of male   ##
## DEG (FDR<0.01) only      ##
## dt.balaton.deg.meth      ## 
## needs to be pre-defined  ##
##############################
plotMethXY<-function(METH_TYPE, my.target, my.padj, min.depth=10, min.num.sites=2){
	if(!my.target %in% c("Promoter","miR Promoter", "CPGi", "miR CPGi", "Gene-body", "miR Gene-body", "CPGi-shore","Enhancer")){
		cat(paste(my.target, "Not supported\n"))
		cat("One of following supported: Promoter, CPGi, Gene-body\n")
		return(1)
	}
	# DEG only
	if(my.padj==0.01){
		p<-ggplot(dt.balaton.deg.meth[tissue==meth.type[[METH_TYPE]] & region==my.region[[my.target]] & abs(log2FoldChange)>log2(1.1) & baseMean>min.depth & num.sites>=min.num.sites], aes(f.met, m.met, label=hgnc_symbol)) + 
			#geom_point(aes(col=`Expression?`, size=abs(log2FoldChange), shape=XCI)) + 
			geom_point(aes(col=`Expression?`, shape=XCI),size=7,alpha=0.7) + 
			scale_colour_manual(values=my.col[["Expression?"]]) +
			geom_text(check_overlap = TRUE, hjust =1.2, nudge_x=0.05, size=4.5) +
			ggtitle(paste("% Methylation at ",my.target,"(No. of CpGs>=",min.num.sites,")")) + 
			labs(x="% Methylation (female)", y="% Methylation (male)") +
			geom_abline(slope=1, lty=3) +
			theme_Publication()
		#print(p+theme(legend.position=c(0.1,0.8)))
		print(p+facet_wrap(~XCI))
	}else{
		dt.dummy<-dt.meth.exp
		if(my.padj<1){
			dt.dummy<-dt.meth.exp[new.padj<my.padj]
		}
		dt.dummy[, `Expression?` := ifelse(log2FoldChange>0, 'female>male', 'female<male')]
		p<-ggplot(dt.dummy[tissue==meth.type[[METH_TYPE]] & region==my.region[[my.target]] & baseMean>=min.depth & num.sites>=min.num.sites], aes(f.met, m.met)) +
			geom_point(aes(col=`Expression?`, shape=`Expression?`),size=5, alpha=0.7) + 
			scale_colour_manual(values=my.col[["Expression?"]]) +
			ggtitle(paste("% Methylation at ",my.target,"(No. of CpGs>=",min.num.sites,"padj<",my.padj,")")) + 
			labs(x="% Methylation (female)", y="% Methylation (male)") +
			geom_abline(slope=1, lty=3) +
			theme_Publication()
		#print(p+theme(legend.position=c(0.1,0.8)))
		print(p)
	}
}

##################################
## X-axis: XCI                  ##
## Y-axis: Boxplot of meth.diff ##
## DEG (FDR<0.01) only          ##
## dt.balaton.deg.meth         ## 
## needs to be pre-defined      ##
##################################
plotBoxBalatonMethDiff<-function(METH_TYPE, my.target, min.depth=10, min.num.sites=2){
	if(!my.target %in% c("Promoter","CPGi","Gene-body","CPGi-shore","Enhancer")){
		cat(paste(my.target, "Not supported\n"))
		cat("One of following supported: Promoter, CPGi, Gene-body\n")
		return(1)
	}   
	p<-ggplot(dt.balaton.deg.meth[tissue==meth.type[[METH_TYPE]] & region==my.region[[my.target]] & abs(log2FoldChange)>log2(1.1) & baseMean>min.depth & num.sites>=min.num.sites], aes(XCI, diff.met)) + 
		geom_boxplot(aes(col=`Expression?`),alpha = .7) + 
		geom_hline(yintercept=0) +
		scale_colour_manual(values=my.col[["Expression?"]]) +
		ggtitle(paste("Methylation difference at",my.target,"(No. of CpGs>=",min.num.sites,")")) + 
		labs(x="XCI status", y="% Methylation Difference\n(Female - Male)") +
		theme_Publication()
	print(p)
}

##################################
## X-axis: XCI                  ##
## Y-axis: % methylation by sex ##
## DEG (FDR<0.01) only          ##
## dt.balaton.deg.meth         ## 
## needs to be pre-defined      ##
##################################
plotBoxBalatonMethbySex<-function(METH_TYPE, my.target, min.depth=10, min.num.sites=2){
	dt.balaton.deg.meth.gender=rbind(
				dt.balaton.deg.meth[,.(tissue,ensembl_gene_id, baseMean, log2FoldChange, XCI, region, `Gender`="Female",`meth`=f.met, num.sites)],
				dt.balaton.deg.meth[,.(tissue,ensembl_gene_id, baseMean, log2FoldChange, XCI, region, `Gender`="Male",`meth`=m.met, num.sites)]
			)
	if(!my.target %in% c("Promoter","CPGi","Gene-body","CPGi-shore","Enhancer")){
		cat(paste(my.target, "Not supported\n"))
		cat("One of following supported: Promoter, CPGi, Gene-body\n")
		return(1)
	}   
	p<-ggplot(dt.balaton.deg.meth.gender[tissue==meth.type[[METH_TYPE]] & region==my.region[[my.target]] & abs(log2FoldChange)>log2(1.1) & baseMean>min.depth & num.sites>=min.num.sites], aes(XCI, meth)) + 
		geom_boxplot(aes(col=Gender),alpha = .7) + 
		scale_colour_manual(values=my.col[["Gender"]]) +
		ggtitle(paste("% Methylation at",my.target,"(No. of CpGs>=",min.num.sites,")")) + 
		labs(x="XCI status", y="% Methylation") +
		theme_Publication()
	print(p)
}

#################
## Carrel 2005 ##
## Use Balaton ##
#################
if(FALSE){
	gr.xci<-with(read.table("~/data/X-inactivation/Carrel.2005/Carrel.science.2005.hg19.bed"), GRanges(V1, IRanges(V2+1, V3), `strand`=V6, `score`=V5 ,`gene`=V4, `strata`=V7))
	correct.gene.name<-getBM(attributes = c("chromosome_name", "ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = gr.xci$gene, mart = myMart) # is a data.frame
	# gr.xci having the correct.gene.name or not
	gr.xci$is.hgnc<-gr.xci$gene %in% unique(correct.gene.name[correct.gene.name$chromosome_name=="X","hgnc_symbol"])
	# from bin/R/RoadMap/plot.roadmap.gender.meth.diff.by.XCI.Carrel.science.2005.R
	# See ~/data/X-inactivation/README for detail
	dt.xci.ensg<-list()
	#########################
	# method1: by gene name #
	#########################
	dt.xci.ensg[['ensg']]<-as.data.table(
							with(
								merge(as.data.frame(gr.xci[gr.xci$is.hgnc,]), unique(correct.gene.name[correct.gene.name$chromosome_name=="X",c("ensembl_gene_id","hgnc_symbol")]), by.x="gene",by.y="hgnc_symbol"),
								data.frame(gene, overlap=width, seqnames, start, end, width, strand, score, ensembl_gene_id, strata, is.hgnc)
							))
	################################################
	# method2: by coordinate (lifted over to hg19) #
	# gene not having the correct.gene.name only   #
	################################################
	# find the overlap between the gene (grch37) and XCI
	gr.no.hgnc<-gr.xci[!gr.xci$is.hgnc,]
	gr.overlap<-mergeByOverlaps(gr.ensg, gr.no.hgnc, minoverlap=5L) # gr.ensg # from Annotation.R

	# find the overlap size between the gene and XCI
	overlap.size=c()
	for(i in 1:nrow(gr.overlap)){overlap.size[i]<-width(intersect(gr.overlap[i,]$gr.ensg, gr.overlap[i,]$gr.no.hgnc, ignore.strand=T))}

	bar<-as.data.table(
			as.data.frame(
				with(gr.overlap, GRanges(seqnames(gr.overlap$gr.no.hgnc), ranges(gr.overlap$gr.no.hgnc), `strand`=strand(gr.overlap$gr.no.hgnc), `score`=score, `ensembl_gene_id`=gene_id, `gene`=gene, `strata`=strata, overlap=overlap.size, is.hgnc))
		))
	# remove duplicated ensg_id (pick the largest overlap)
	dt.xci.ensg[['no.ensg']]<-unique(merge(bar, bar[,list(overlap=max(overlap)),gene], by=c("gene","overlap")))
	##############################
	## merge method 1 & method 2 #
	##############################
	#rbindlist(dt.xci.ensg)

	#######################################
	# chromosome location to identify PAR #
	#######################################
	#dt.query=as.data.table(getBM(attributes =c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"), filters = "ensembl_gene_id", values = dt.deg.chrX[,ensembl_gene_id], mart = myMart))
	#dt.par<-as.data.table(with(as.data.frame(gr.par), data.frame(par=names(gr.par),seqnames,start,end,strand)))
	#setkeyv(dt.par,c("start","end"))
	#foverlaps(dt.query, dt.par, by.x=c("start_position","end_position"), type="any", nomatch=0L)
}

plotBoxMethDiffAllTissue<-function(min.num.sites=4){
	if(RNASEQ_TYPE=="TOTAL"){
		load("~/results/RoadMap/X-inactivation/RData/dt.meth.chrX.per.tissue.CG.Genes_20.20.RData") # WGoxBS and SureSelect (PT, PT.SS.Tech, PT.SS.Bio)
		dt.foo<-dt.meth.chrX.per.tissue
		rm(dt.meth.chrX.per.tissue)
	}else if(RNASEQ_TYPE=="SMALL"){
		load("~/results/RoadMap/X-inactivation/RData/dt.meth.chrX.mir.per.tissue.CG.miR_Genes_10.10.RData")
		dt.foo<-dt.meth.chrX.mir.per.tissue
		rm(dt.meth.chrX.mir.per.tissue)
	}else{
		stop(paste(RNASEQ_TYPE, "Not supported\nEither TOTAL or SMALL"))
	}
	## dt.meth.all.chrX: Note that ensembl_gene_id can be null <NA> because of left-join (all.x=TRUE)
	## ensembl_gene_id from RnBeads (Ensembl 73 or gencode 18)
	dt.meth.chrX.all.tissue=rbindlist(dt.foo)[!grepl(".SS.",tissue) & !is.na(ensembl_gene_id),
														.(f.met=sum(female.c)/sum(female.cov)*100, m.met=sum(male.c)/sum(male.cov)*100, diff.met=sum(female.c)/sum(female.cov)*100-sum(male.c)/sum(male.cov)*100, num.sites=sum(num.sites)),
														by="tissue,ensembl_gene_id,region"]
	p1<-ggplot(dt.meth.chrX.all.tissue[region %in% cpg.contexts & num.sites>=min.num.sites], aes(tissue, diff.met)) + 
		geom_boxplot(aes(fill=region),outlier.shape=NA,alpha=0.9) +
		geom_hline(yintercept=0) +
		scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Genes"]],my.region.rev[["Promo_10.10"]])) +	
		coord_cartesian(ylim=c(-50,70)) +
		scale_y_continuous(breaks=seq(from=-50,to=70,by=10)) + 
		scale_x_discrete(limits=c("AD","AO","EG","FT","GA","PO","SB","SX","PT")) + 
		ggtitle("Methylation Difference by Sex") +
		labs(x="Tissues", y="% Methylation Difference\n(Female - Male)") +
		theme_Publication() 
	print(p1)

	if(FALSE){
	p2<-ggplot(dt.meth.chrX.all.tissue[region %in% cpg.contexts & num.sites>=min.num.sites], aes(tissue, round((f.met+m.met)/2,2))) + 
		geom_boxplot(aes(fill=region),outlier.shape=NA,alpha=0.9) +
		scale_fill_manual(values=cbPalette, name="Region", labels=c(my.region.rev[["CPGi"]],my.region.rev[["Genes"]],my.region.rev[["Promo_10.10"]])) +	
		scale_x_discrete(limits=c("AD","AO","EG","FT","GA","PO","SB","SX","PT")) + 
		ggtitle("Methylation Difference by Sex") +
		labs(x="Tissues", y="% Methylation") +
		theme_Publication() 
	print(p2)
	}

}# end of plotBoxMethDiffAllTissue

