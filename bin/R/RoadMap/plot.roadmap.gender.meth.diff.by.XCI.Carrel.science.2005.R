#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First crated 14/Mar/2016
# Last modified 23/Aug/2016

TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/RoadMap/local.R")

my.cpg.type="CG" # CG, CHH, CHG
#############################
## Get XCI Genes           ##
## from Carrel Nature 2005 ##
#############################
fields <- c("chromosome_name", "ensembl_gene_id", "hgnc_symbol", "gene_biotype")
# method1: by gene name
if(TRUE){
	my.RData="~/results/RoadMap/X-inactivation/RData/carrel.xci.genes.grch37.RData"
	if(file.exists(my.RData)){
		load(my.RData)
		cat("xci.anno and dt.xci.anno loaded\n")
	}else{
		fields <- c("chromosome_name", "ensembl_gene_id", "hgnc_symbol", "gene_biotype")
		# method1: by gene name
		carrel.xci=read.csv("~/data/X-inactivation/Carrel.science.2005.csv")
		xci<-getBM(attributes = fields, filters = "hgnc_symbol", values = carrel.xci$gene, mart = grch37) # is a data.frame
		xci.anno<-merge(carrel.xci, xci, by.x="gene", by.y="hgnc_symbol") # isa 'data.frame'
		xci.anno<-xci.anno[xci.anno$chromosome_name=="X",]
		dt.xci.anno<-as.data.table(xci.anno)
		save(xci.anno, dt.xci.anno, file="~/results/RoadMap/X-inactivation/RData/carrel.xci.genes.grch37.RData")
	}
}

# method2: by coordinate (lifted over from hg18) 
# incorrect lifting TKTL1 for example
if(FALSE){
	gr.genes<-get.region("Genes")
	mcols(gr.genes)<-DataFrame(mcols(gr.genes), ensembl_gene_id=names(gr.genes)) # add ensg id

	# See ~/data/X-inactivation/README for detail
	gr.xci<-rtracklayer::import.bed("~/data/X-inactivation/carrel.science.2005.coordinate.hg19.bed", genome="hg19")
	# find the overlap between the gene and XCI
	dt.overlap<-as.data.table(with(mergeByOverlaps(gr.genes, gr.xci, minoverlap=50L), data.frame(ensembl_gene_id=ensembl_gene_id, score=score)))
	carrel.xci<-dt.overlap[,round(mean(score)),by="ensembl_gene_id"] # isa 'data.table' (group by gene name and merge score)
	carrel.xci[,type:=ifelse(carrel.xci$V1==0, "XI", ifelse(carrel.xci$V1==9, "XIE", "XIC"))] # isa 'data.table'
	setnames(carrel.xci, c("ensembl_gene_id","xci","type"))

	xci<-getBM(attributes = fields, filters = "ensembl_gene_id", values = carrel.xci$ensembl_gene_id, mart = grch37) # is a data.frame
	xci.anno<-merge(carrel.xci, xci, by="ensembl_gene_id") # isa 'data.table'
	rm(gr.genes,gr.xci)
	write.csv(xci.anno, file="~/data/X-inactivation/Carrel.science.2005.grch37.coord.based.anno.csv", row.names=F)

	#TKTL1 (ENSG00000007350, grch37 ensembl)
	#chrX 153524024 153558700 TKTL1 1 +

	#TKTL1 (ENSG00000268013, grch37 ensembl)
	#chromosome HG1497_PATCH: 153464139 153498818
	
	#TKTL1 (lifted over from hg18 to hg19 locally)
	#chrX 153391803 153426352 TKTL1 1 +

	#TKTL1 (lifted over from hg18 to hg19 from the web)
	#chrX 153391803 153426352 TKTL1 1 +

	#TKTL1 (hg18 based on Carrel Nature 2005)
	#chrX 153044997 153079546 TKTL1 1 +
}

# XI : X-chr inactivated genes
# XIC: between XI and XIE
# XIE: X-chr inactivation escaped genes
######################################################################
## JOB1: re-produce 'extended data figure 9' of Schultz 2015 Nautre ##
######################################################################
if(TRUE){
	my.cpg.type="CG" # CG, CHH, CHG
	# see the code above
	my.RData="/home/ssg29/results/RoadMap/X-inactivation/RData/XCI.CG.2016-04-20_08_37_AM.RData" #: Genes, promoter, CPGi etc
	#my.RData="/home/ssg29/results/RoadMap/X-inactivation/RData/XCI.CG.2016-08-01_11_59_AM.RData" #: various promoters only
	if(file.exists(my.RData)){
		cat("loading...\n")
		load(my.RData)
		cat("cpg.tissue.list and dt.cpg.xci loaded\n")
	}else{
		xci.list<-split(dt.xci.anno, dt.xci.anno[,type]) # isa 'list'
		cpg.tissue.list=list() # to save methylation level by XCI 
		cpg.xci.list=list()
		#cpg.contexts=c("Promo_15.15", "Promo_15.10", "Promo_15.05", "Promo_15.00", "Promo_10.15","Promo_10.10","Promo_10.05","Promo_10.00") # ~/results/RoadMap/X-inactivation/RData/XCI.CG.2016-08-01_11_59_AM.RData
		cpg.contexts=c("Tiling5K","Genes","Promo_15.05", "Promo_10.00", "CPGi","CPGi-no-promo","CPGi-with-promo","Enhancer")
		#cpg.contexts=c("Tiling5K","Genes","Promo_15.05", "Promo_10.00", "CPGi","CPGi-no-promo","CPGi-with-promo","CPGi-no-gene","CPGi-with-gene","CPGi-shores","Enhancer")
		# 1. for this tissue
		for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){
			dummy=list()
			## Load this RoadMap tissue
			dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)
			cat("Filtering out autosomes...\n")
			dt.query=dt.query[V1=="X"] #chrX only
			# 2. for this xci type (e.g. "XI", "XIC", "XIE") and "chrX"
			for(my.xci in c(names(xci.list),"chrX")){
				cat(paste0("Processing ", my.xci, "...\n"))
				# 3. for this genomic region
				for(my.context in cpg.contexts){
					# call a user function defined within local.R
					cpg.tissue.list[[my.tissue]][[my.xci]][[my.context]]=get.meth.for.xci(my.tissue, my.xci, my.context)
				}
				dummy[[my.tissue]][[my.xci]]=rbindlist(cpg.tissue.list[[my.tissue]][[my.xci]]) # concat all cpg context types 
			}
			cat(paste0("Done for ", my.tissue, "...\n"))
			cpg.xci.list[[my.tissue]]=rbindlist(dummy[[my.tissue]]) # concat all xci types 
		}# end of my.tissue 
		dt.cpg.xci=rbindlist(cpg.xci.list) # concat all xci types 
		save(cpg.tissue.list,dt.cpg.xci, file=file.path("~/results/RoadMap/X-inactivation/RData",paste("XCI",my.cpg.type,time.stamp,"RData",sep=".")))
		rm(cpg.xci.list,cpg.tissue.list)
	}
	# make publication figures
	file.name<-file.path("~/results/RoadMap/X-inactivation",paste("XCI",my.cpg.type,time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=11.7, height=8.3, title="methylation difference by XCI")

	for(i in dt.cpg.xci[,unique(cpg.region)]){
		if(my.cpg.type=="CG"){
			if(grepl('CPGi',i) || grepl('Promo',i)){
				my.lim=c(-25,43)
			}else{
				my.lim=c(-25,25)
			}
			my.y.tick=seq(from=my.lim[1],to=my.lim[2], by=10)  # Ticks every .1
		}else{
			if(grepl('CPGi',i)){
				my.lim=c(-1.5,1)
			}else{
				my.lim=c(-1,1)
			}
			my.y.tick=seq(from=my.lim[1],to=my.lim[2], by=0.5)  # Ticks every .005
		}
		# Boxplot by CpG context (chrX only) 
		r<-ggplot(dt.cpg.xci[cpg.region==i], aes(XCI, diff.met*100)) + 
			geom_boxplot(aes(fill=tissue),  outlier.shape=NA) + 
			labs(x="X-Chromosome Inactivation Class", y="% methylation difference (female - male)") + 
			geom_hline(yintercept=0) + 
			ggtitle(i) +
			scale_fill_manual(values=cbPalette, name="Tissue types") +
			coord_cartesian(ylim=my.lim) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
			theme_Publication() +
			scale_y_continuous(breaks=my.y.tick) 
		print(r)
	}
	dev.off()
}

#print for a test
if(FALSE){
	sapply(cpg.tissue.list[["SX"]], function(i) sapply(i, nrow))
	sum(rowSums(sapply(cpg.tissue.list[["SX"]], function(i) sapply(i, nrow))))
	dt.cpg.xci[cpg.region=="Genes",list(total.cpg=sum(num.sites),cnt=.N),by="tissue,XCI,cpg.region"]

	# Hitogram-line by CpG context (chrX only) 
	p<-ggplot(dt.cpg.xci, aes(diff.met)) + 
		geom_density(aes(colour=type)) + 
		labs(x="Methylation Difference (female - male)") + 
		geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
		ggtitle(paste0(tissue.list[[my.tissue]], ": chrX")) +
		scale_colour_manual(values=cbPalette, name="CpG\nContext")
	print(p)

	# Hitogram-fill by CpG context (chrX only) 
	q<-ggplot(dt.cpg.xci, aes(diff.met)) + 
		geom_density(aes(fill=type),alpha=.2) + 
		labs(x="Methylation Difference (female - male)") + 
		geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
		ggtitle(paste0(tissue.list[[my.tissue]], ": chrX")) +
		scale_fill_manual(values=cbPalette, name="CpG\nContext")
	print(q)
}

if(FALSE){
	## initialise dt.meth.chrX.per.tissue
	## Note that ensembl_gene_id can be null <NA> because of left-join (all.x=TRUE)
	## dt.meth.all.chrX[is.na(ensembl_gene_id),.N,"tissue,region"]
	prep.dt.meth.chrX.per.tissue(my.cpg.type="CG")

}

####################################################################
## JOB3: count the number of diff meth (>=10%) regions by tissues ##
####################################################################
if(FALSE){
	# ensembl_gene_id from RnBeads (Ensembl 73 or gencode 18)
	load("~/results/RoadMap/X-inactivation/RData/dt.meth.chrX.per.tissue.CG.Genes_20.20.RData")# cpg.contexts=c("Genes","Promo_15.05", "Promo_10.03","CPGi","CPGi-shores")
	cat("dt.meth.chrX.per.tissue loaded\n")
	dt.meth.all.chrX<-rbindlist(dt.meth.chrX.per.tissue)

	min.diff=.1
	dt.meth.summary<-rbind(
		dt.meth.all.chrX[male.c/male.cov - female.c/female.cov >=min.diff,.(cnt.region=.N, num.cpg=sum(num.sites),class="female<male"), by="tissue,region"],  # male 10% less methylated than female
		dt.meth.all.chrX[female.c/female.cov - male.c/male.cov >=min.diff,.(cnt.region=.N, num.cpg=sum(num.sites),class="female>male"), by="tissue,region"] # female 10% more methylated than male
		)

	file.name<-file.path("~/results/RoadMap/X-inactivation",paste("number.diff.meth.region.by.tissue",my.cpg.type,time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=11.7, height=8.3, title="Number of Diff Meth by Tissue")
	cpg.contexts=c("Genes","Promo_15.05", "Promo_10.03","CPGi","CPGi-shores")
	for(i in cpg.contexts){
		p<-ggplot(dt.meth.summary[region==i], aes(tissue, cnt.region, fill=class)) + 
			geom_bar(stat="identity") + 
			scale_fill_manual(values=cbPalette) + 
			ggtitle(paste0("No. of Differentially Methylated (>=",min.diff*100,"%) ",i, " (chrX only)"))
		print(p)
	}
	dev.off()
}# end of JOB3

##############################################################
## JOB4: 
## XIE where female-male are negative at CPGi
##############################################################
if(FALSE){
	# ensembl_gene_id from RnBeads (Ensembl 73 or gencode 18)
	load("~/results/RoadMap/X-inactivation/RData/dt.meth.chrX.per.tissue.CG.Genes_20.20.RData")# cpg.contexts=c("Genes","Promo_15.05", "Promo_10.03","CPGi","CPGi-shores")
	cat("dt.meth.chrX.per.tissue loaded\n")
	dt.meth.all.chrX<-rbindlist(dt.meth.chrX.per.tissue)
	dt.meth.xci<-merge(dt.meth.all.chrX, xci.anno, by="ensembl_gene_id")
	# example
	# CPGi of PT where male>female, get gene and XCI
	dt.meth.xci[(type=="XI" | type=="XIE") & tissue=="PT" & region=="CPGi" & female.c/female.cov-male.c/male.cov<0, .(num.region=.N), "gene,type"]
	# CPGi of PT where male 10% more than female, get gene and XCI  
	dt.meth.xci[(type=="XI" | type=="XIE") & tissue=="PT" & region=="CPGi" & male.c/male.cov-female.c/female.cov>=.1, .(ensembl_gene_id,gene,type,female.met=female.c/female.cov, male.met=male.c/male.cov)]
	#    ensembl <- gene <- id   gene type female.met  male.met
	# 1: ENSG00000006756   ARSD  XIE 0.56190476 0.8086957
	# 2: ENSG00000006756   ARSD  XIE 0.01501976 0.3454545
	# 3: ENSG00000063587 ZNF275   XI 0.82646593 0.9320221
	# 4: ENSG00000063601  MTMR1   XI 0.81391304 0.9370370
	# CPGi of PT where female 10% more than male, get gene and XCI  
	dt.meth.xci[(type=="XI" | type=="XIE") & tissue=="PT" & region=="CPGi" & female.c/female.cov-male.c/male.cov>=.1, .(ensembl_gene_id,gene,type,female.met=female.c/female.cov, male.met=male.c/male.cov)]

	# per region,tissue,type,gene
	dt.dummy1<-dt.meth.xci[(type=="XI" | type=="XIE") & (female.c/female.cov - male.c/male.cov >=min.diff),.(cnt.region=.N, num.cpg=sum(num.sites), class="female>male"), by="region,tissue,type,gene"]   # female 10% more methylated than male
	dt.dummy2<-dt.meth.xci[(type=="XI" | type=="XIE") & (male.c/male.cov - female.c/female.cov >=min.diff),.(cnt.region=.N, num.cpg=sum(num.sites), class="female<male"), by="region,tissue,type,gene"]  # male 10% less methylated than female

	# per region,tissue,type,class (aggregate gene)

	dt.meth.xci.summary<-rbind(
		dt.dummy1[,.(num.gene=.N,cnt.region=sum(cnt.region),num.cpg=sum(num.cpg)),"region,tissue,type,class"], # female 10% more methylated than male
		dt.dummy2[,.(num.gene=.N,cnt.region=sum(cnt.region),num.cpg=sum(num.cpg)),"region,tissue,type,class"]  # male 10% less methylated than female
			)

	cpg.contexts=c("Genes","Promo_15.05", "Promo_10.03","CPGi","CPGi-shores")

	file.name<-file.path("~/results/RoadMap/X-inactivation",paste("number.xci.gene.by.tissue",my.cpg.type,time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=11.7, height=8.3, title="Number of Diff Meth by Tissue")
	for(i in cpg.contexts){
		for(j in c("female<male", "female>male")){
			p<-ggplot(dt.meth.xci.summary[region==i & class==j], aes(tissue,num.gene)) + 
				geom_bar(aes(fill=type),stat="identity",position="dodge") +
				scale_fill_manual(values=cbPalette) + 
				ggtitle(paste0("No. of genes where methylation of ",i," ",j," (at least ",min.diff*100,"%, chrX only)"))
			print(p)
		}
	}
	dev.off()
}# end of JOB4

cat("All done\n")
