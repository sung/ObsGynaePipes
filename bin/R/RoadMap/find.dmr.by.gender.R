#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 3/Mar/2016
# Last modified 9/May/2016

TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/RoadMap/local.R")
#source("~/Pipelines/lib/methylkit.R") # defines fast.fisher

############
# Make DMR #
############
# prep.dmr defined within local.R
#for(my.tissue in c(avail.tissues[!avail.tissues %in% "PA"])){print(system.time(prep.dmr(my.tissue,"CHG")))}
#for(my.tissue in c(avail.tissues[!avail.tissues %in% "PA"])){print(system.time(prep.dmr(my.tissue,"CHH")))}

##############################
## Update DMR Combined Rank ##
##############################
# update.dmr.rank defined within local.R
#for(my.tissue in c(avail.tissues[!avail.tissues %in% "PA"])){print(system.time(update.dmr.rank(my.tissue,"CHG")))}
#for(my.tissue in c(avail.tissues[!avail.tissues %in% "PA"])){print(system.time(update.dmr.rank(my.tissue,"CHH")))}

if(FALSE){
	my.tissue="PT"
	my.cpg.type="CHG" # CG, CHH, CHG

	#my.RData=file.path("~/results/RoadMap/CH/RData",paste(my.tissue,my.cpg.type,"dt.dmr.RData",sep="."))
	my.RData=file.path("~/results/RoadMap/CH/RData",paste(my.tissue,my.cpg.type,"dt.dmr.rank.RData",sep="."))
	if(!file.exists(my.RData)){
		stop(paste0(my.RData," not found\n"))
	}else{
		print(system.time(load(my.RData)))
		cat("dt.dmr.per.region loaded\n")
	}

	#########################
	## Filter by num.sites ##
	#########################
	#new.dt.dmr.per.region=lapply(dt.dmr.per.region, function(i) i[num.sites>=summary(i[,num.sites])[3]]) # at least median num.sites
	new.dt.dmr.per.region=lapply(dt.dmr.per.region, function(i) i[num.sites>=i[,median(num.sites)]]) # at least median num.sites

	sapply(dt.dmr.per.region, nrow) # number of region before filtering 
	sapply(dt.dmr.per.region, function(i) i[,median(num.sites)]) # number of median num.sites

	sapply(new.dt.dmr.per.region, nrow) # number of region after filtering
	sapply(new.dt.dmr.per.region, function(i) quantile(i[,met],.999)*100) # top 0.1% methylation level
	sapply(new.dt.dmr.per.region, function(i) nrow(i[met>=quantile(i[,met],.999)])) # number of top 0.1% region


	# boxplot of meth.diff by chr
	ggplot(rbindlist(new.dt.dmr.per.region)[met.f>0.0001 & met.m>0.0001], aes(chr, met.f-met.m, col=region)) + geom_boxplot(outlier.shape=NA) + coord_cartesian(ylim=c(-0.01,0.01))
	ggplot(rbindlist(dt.dmr.per.region), aes(chr, met.f-met.m, fill=region)) + geom_boxplot(outlier.shape=NA) + coord_cartesian(ylim=c(-0.01,0.01))

	# Boxplot of % meth by region
	ggplot(rbindlist(new.dt.dmr.per.region)[met.f>0.0001 & met.m>0.0001], aes(region, met.f-met.m)) + geom_boxplot(aes(fill=region),outlier.shape=NA) + coord_cartesian(ylim=c(-0.01,0.01))

	#################################
	## Selected Methylated regions ##
	#################################
	file.name<-file.path("~/results/RoadMap/CH/Figures",paste("top.meth.boxplot",my.cpg.type,my.tissue,time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=11.7, height=8.3, title="Percentage of Occurrence by Methylation Level")

	my.region=sapply(new.dt.dmr.per.region, nrow)[sort(names(sapply(new.dt.dmr.per.region, nrow)))] # number of region (sorted by region name)
	my.up.to=ceiling(max(sapply(new.dt.dmr.per.region, function(i) i[,max(met)*100])))

	ggplot(rbindlist(new.dt.dmr.per.region), aes(region, met*100)) + 
	geom_boxplot(aes(fill=region)) + 
	labs(y="% met") + 
	scale_x_discrete(labels=paste0(names(my.region), "(n=", my.region, ")")) +
	scale_y_continuous(breaks=seq(from=0,to=my.up.to,by=1)) +
	ggtitle(paste(my.cpg.type,"Methylated Regions from",my.tissue))

	################################
	## Top 100 methylated regions ##
	################################
	top.dt.met.per.region=lapply(new.dt.dmr.per.region, function(i) i[order(-met)][1:100]) 

	my.top.region=sapply(top.dt.met.per.region, nrow)[sort(names(sapply(top.dt.met.per.region, nrow)))] # number of region (sorted by region name)
	my.up.to=ceiling(max(sapply(top.dt.met.per.region, function(i) i[,max(met)*100])))

	ggplot(rbindlist(top.dt.met.per.region), aes(region, met*100)) + 
	geom_boxplot(aes(fill=region)) + 
	labs(y="% met") + 
	scale_x_discrete(labels=paste0(names(my.top.region), "(n=", my.top.region, ")")) +
	scale_y_continuous(breaks=seq(from=0,to=my.up.to,by=1)) +
	ggtitle(paste("Top 100 Highly Methylated",my.cpg.type,"Regions from",my.tissue))

	################################
	## Top 0.1% methylated region ##
	################################
	top.dt.met.per.region=lapply(new.dt.dmr.per.region, function(i) i[met>=quantile(i[,met],.999)]) 

	my.top.region=sapply(top.dt.met.per.region, nrow)[sort(names(sapply(top.dt.met.per.region, nrow)))] # number of region (sorted by region name)
	my.up.to=ceiling(max(sapply(top.dt.met.per.region, function(i) i[,max(met)*100])))

	sapply(top.dt.met.per.region, function(i) i[,min(met)*100]) # min meth 
	sapply(top.dt.met.per.region, function(i) i[,max(met)*100]) # max meth

	ggplot(rbindlist(top.dt.met.per.region), aes(region, met*100)) + 
	geom_boxplot(aes(fill=region)) + 
	labs(y="% met") + 
	scale_x_discrete(labels=paste0(names(my.top.region), "(n=", my.top.region, ")")) +
	scale_y_continuous(breaks=seq(from=0,to=my.up.to,by=1)) +
	ggtitle(paste("Top 0.1% Highly Methylated",my.cpg.type,"Regions from",my.tissue))

	dev.off()
} # Boxplot of % meth by region

##########################################
## Where the top methylated regions are ##
## nearest gene annotations             ##
##########################################
if(FALSE){
	gr.genes<-get.region("Genes") # RnBeads gene definitions (use 'gr.ensg' for UCSC ensGene)
	mcols(gr.genes)<-DataFrame(mcols(gr.genes), ensembl_gene_id=names(gr.genes)) # add ensg id

	dt.top.meth.list=list()
	for(i in names(top.dt.met.per.region)){
		#top.dt.met.per.region[["Tiling500"]][order(-met)]
		gr.target=makeGRangesFromDataFrame(top.dt.met.per.region[[i]],keep.extra.columns=TRUE) 
		seqlevels(gr.target)=mapSeqlevels(seqlevels(gr.target), style="UCSC") # X=>chrX

		# find.nearest.genes is defined in bin/R/RoadMap/local.R
		# be sure to pre-defined 'gr.genes'
		dt.target.with.gene=find.nearest.genes(gr.target)

		dt.top.meth.list[[i]]=dt.target.with.gene

		# save top meth file
		my.file=file.path("~/results/RoadMap/CH/TopMeth",paste("top.meth",my.cpg.type,my.tissue,i,"csv",sep="."))
		write.csv(dt.target.with.gene[order(-met)], file=my.file, row.names=F)
	}

	############################
	## Gene Ontology Analysis ##
	############################
	cat("doing GO analysis...\n")
	library(goseq) # for goseq
	library(GO.db) # for ,GOBPPARENTS
	library(GOstats) # for GOGraph

	p.cuff.off=.01
	i="Tiling500"
	# target: highly methylated genes
	# background 1: all genes
	myGenes<-as.integer(names(gr.genes) %in% dt.top.meth.list[[i]][,unique(ensembl_gene_id)])
	#myGenes<-as.integer(names(gr.genes) %in% unique(gr.target$ensembl_gene_id))
	names(myGenes)<-names(gr.genes)

	myGenome="hg19"
	pwf=goseq::nullp(myGenes,myGenome,"ensGene") # isa data.frame  (check available ID via supportedGeneIDs())
	GO.wall=goseq::goseq(pwf,myGenome,"ensGene") # isa data.frame
	GO.wall$over_represented_qvalue=p.adjust(GO.wall$over_represented_pvalue, method="BH")
	GO.wall$under_represented_qvalue=p.adjust(GO.wall$under_represented_pvalue, method="BH")
	enriched.GO.wall=list('over'=GO.wall[GO.wall$over_represented_qvalue<p.cuff.off,], # isa data.frame
						'under'=GO.wall[GO.wall$under_represented_qvalue<p.cuff.off,] # isa data.frame
					)
	# Find the level (depth) from the root (depth-first or shortest) given a GO term
	for(j in names(enriched.GO.wall)){
		enriched.GO.wall[[j]]$depth<-as.numeric(sapply(enriched.GO.wall[[j]]$category, function(i) {
					if(Ontology(i)=="BP"){
						graph::DFS(GOGraph(i,GOBPPARENTS), "all")[i]
					}else if(Ontology(i)=="MF"){
						graph::DFS(GOGraph(i,GOMFPARENTS), "all")[i]
					}else{
						graph::DFS(GOGraph(i,GOCCPARENTS), "all")[i]
					}
				} # end of function(i)
			)# end of sapply
		) # end of as.numeric
	}

	cat("writing GO analysis results...\n")
	write.csv(enriched.GO.wall[["over"]], file=file.path("~/results/RoadMap/CH/TopMeth", paste("go.over.top.meth",my.cpg.type,my.tissue,i,"csv",sep=".")), row.names=F)
	#write.csv(enriched.GO.wall[["over"]], file=file.path("~/results/RoadMap/CH/TopMeth", paste("go.under.top.meth",my.cpg.type,my.tissue,i,"csv",sep=".")), row.names=F)

	################################
	## GO terms to Word Cloud     ##
	## for top methylated regions ##
	################################
	file.name<-file.path("~/results/RoadMap/CH/Figures",paste("go.word.cloud.",my.cpg.type,my.tissue,time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=11.7, height=8.3, title="Enriched GO terminology for highly mehtyled non-CG")

	cat("making WordCloud...\n")
	library(tm)
	library(wordcloud)
	library(SnowballC)

	terms<-as.character(unlist(sapply(enriched.GO.wall[["over"]][enriched.GO.wall[["over"]]$depth-1>=3,"term"], function(i) strsplit(i, " ")))) # isa character vector
	#terms<-as.character(unlist(sapply(enriched.GO.wall[["over"]][,"term"], function(i) strsplit(i, " ")))) # isa character vector
	terms<-tm::Corpus(VectorSource(terms)) # isa 'VCorpus'

	terms<- tm_map(terms, content_transformer(tolower))
	terms<- tm_map(terms, removePunctuation)
	terms<- tm_map(terms, removeNumbers)
	terms<- tm_map(terms, removeWords, stopwords("english"))
	terms<- tm_map(terms, stripWhitespace)
	#terms<- tm_map(terms, stemDocument) #install.packages("SnowballC")
	wordcloud(terms, scale=c(5,0.5), max.words=100, random.order=FALSE, rot.per=0.35, use.r.layout=FALSE, colors=brewer.pal(8, "Dark2")) # brewer.pal from RColorBrewer

	dev.off()
}

################################
## DMR met.f>met.m frequency  ##
################################
if(FALSE){
	my.cpg.type="CHG" # CG, CHH, CHG

	for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){
		my.RData=file.path("~/results/RoadMap/CH/RData",paste(my.tissue,my.cpg.type,"dt.dmr.rank.RData",sep="."))
		if(!file.exists(my.RData)){
			stop(paste0(my.RData," not found\n"))
		}else{
			print(system.time(load(my.RData)))
			cat("dt.dmr.per.region loaded\n")
		}
		#########################
		## Filter by num.sites ##
		#########################
		new.dt.dmr.per.region=lapply(dt.dmr.per.region, function(i) i[num.sites>=i[,median(num.sites)]]) # at least median num.sites
		# topX
		for(j in c(100, 500, 1000, 2000)){
			top.dt.dmr.per.region=lapply(new.dt.dmr.per.region, function(i) i[order(rank)][1:j]) # top X DMR by combined rank
			lapply(top.dt.dmr.per.region, function(i) i[,list(proportion=round(.N/nrow(i)*100,2)),met.f-met.m>=0])
		}
	}
}


###############################
## Where the DMR regions are ##
## nearest gene annotations  ##
###############################
if(FALSE){
	dt.top.dmr.list=list()
	top.dt.dmr.per.region=lapply(new.dt.dmr.per.region, function(i) i[order(rank)][1:100]) # top 100 DMR by combined rank
	for(i in names(top.dt.dmr.per.region)){
		gr.target=makeGRangesFromDataFrame(top.dt.dmr.per.region[[i]],keep.extra.columns=TRUE)
		seqlevels(gr.target)=mapSeqlevels(seqlevels(gr.target), style="UCSC") # X=>chrX (style="NCBI", chrX=>X)
		# find.nearest.genes is defined in bin/R/RoadMap/local.R
		# be sure to pre-defined 'gr.genes'
		dt.target.with.gene=find.nearest.genes(gr.target)
		dt.top.dmr.list[[i]]=dt.target.with.gene

		# save top DMR file
		my.file=file.path("~/results/RoadMap/CH/DMR/gender/",paste("top.dmr",my.cpg.type,my.tissue,i,"csv",sep="."))
		write.csv(dt.target.with.gene[order(rank)], file=my.file, row.names=F)
	}
	############################
	## GO terms to Word Cloud ##
	## for DMR                ##
	############################
}


###########################
## Manhattan plot of DMR ##
###########################
if(FALSE){
	library(qqman)
	for(i in names(new.dt.dmr.per.region)){
		file.name<-file.path("~/results/RoadMap/CH/Figures",paste("dmr.manhattan",i,my.cpg.type,my.tissue,time.stamp,"tiff",sep="."))
		tiff(filename=file.name, width=11.7, height=8.27, units="in", res=300, compression='lzw')

		my.dmr<-new.dt.dmr.per.region[[i]]
		upper<-log(my.dmr[,max(rank)])

		top.diff.met<-my.dmr[1:10,paste(chr,strand,start,end,sep=".")]
		dummy<-with(my.dmr, 
				data.frame(
					SNP=paste(chr,strand,start,end,sep="."),
					CHR=sapply(chr, function(i) if(i=="X"){i<-23}else if(i=="Y"){i<-24}else if(i=="M"){i<-25}else{i<-as.numeric(i)}),
					BP=start,
					P=p.value,
					RANK=upper-log(rank),
					NUM.SITE=num.sites
				))
		main.title = paste("Manhattan plot of ",i," (",nrow(my.dmr),"regions of >=",my.dmr[,min(num.sites)],my.cpg.type, "at", my.tissue,")")
		manhattan(dummy, p="RANK", logp=FALSE, ylab=paste0(round(upper,2),"-log(rank)"), genomewideline=FALSE, suggestiveline=FALSE, main=main.title, highlight=top.diff.met, chrlabs=c(1:22, "X")) 
		dev.off()
	}
}

cat("All done\n")
