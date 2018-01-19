#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# Last modified 3/Mar/2016

source ("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/bin/R/RoadMap/local.R")

####################
# Load RnBeads     #
# AGA Boy vs. Girl #
####################
if(FALSE){
	library(RnBeads)
	logger.start("Logger started", fname = NA )

	options(fftempdir=temp.dir) # enable 'ff' package to write on disk
	unixtools::set.tempdir(temp.dir) # unixtools downloaded from https://rforge.net/doc/packages/unixtools
	cat("Temp dir:\n")
	print(tempdir())

	rnb.set<-RnBeads::load.rnb.set(path=preprocessed.rnb.set.dir) #isa 'RnBSet'
	rnb.gl<-rnb.RnBSet.to.GRangesList(rnb.set) # isa 'GRangeList'

	destroy(rnb.set) #The ff files behind an RnBeads object can be deleted completely from the hard disk by executing the destructor method:
	warnings() # print out the warnings
	logger.completed()
	logger.close()
}

################################
# Load processed diff.met data # 
################################
#rnb.dm.rdata<-"~/results/RnBeads/AGA.Boy.Girl/CpG/RData/rnb.meth.diff.all.AGA.Boy.Girl.reports-2015-03-30_04_47_PM.RData"
cat(paste0("Loading rnb.dm.rdata: ",rnb.dm.rdata,"...\n"))
load(rnb.dm.rdata) # loading "rnb.diff.met"
cat("rnb.diff.met loaded\n")

###################################
## Methylation Difference by Sex ##
## Chr-wise                      ##
## WGBS                          ##
## density plot and BoxPlot
###################################
# pdf output filename 
file.name<-file.path(run.dir,paste("rnb.meth.diff",my.project,my.run.date,"pdf",sep=".")) # this will store 'rnb.diff.met' which is a list of data.frame for sites, promoters...
pdf(file=file.name, width=11.7, height=8.3, title=my.project)

# Note that 'rnb.diff.met' coordinate for the "sites" for Boy.Girl analysis, which is not converged to BED format yet
# + strand: Start: 1-based C-position, End:   1-based G-position, 
# - strand: Start: 1-based C-position, End:   1-based G-position
# therefore 'Start' for +strand, 'End' for -strand
# GRange is supposed to have 1-based 
# for each cpg context
for(i in names(rnb.diff.met)){
	select=rnb.diff.met[[i]]$Chromosome!="chrY"

	my.case=ifelse(i=="sites","mean.g1","mean.mean.g1") # mean.g1: meth(girl), mean.g2: meth(boy)
	my.ctrl=ifelse(i=="sites","mean.g2","mean.mean.g2")

	my.x=ifelse(i=="sites","mean.diff","mean.mean.diff") # meth(girl)-meth(boy)
	my.y=ifelse(i=="sites","log10(1/diffmeth.p.val)","log10(1/comb.p.val)")

	# meth diff colour/fill by chr
	g<-ggplot(rnb.diff.met[[i]][select,], aes_string(x=my.x,y="..density..")) + geom_freqpoly(aes(color=Chromosome),alpha=.7) + labs(x="Methylation Difference (Girl - Boy)") + ggtitle(i)
	g1<-ggplot(rnb.diff.met[[i]][select,], aes_string(x=my.x)) + geom_density(aes(fill=Chromosome),alpha=.7) + labs(x="Methylation Difference (Girl - Boy)") + ggtitle(i)

	# meth diff by chr facet
	#g<-ggplot(rnb.diff.met[[i]][select,], aes_string(x=my.x,y="..density..")) + geom_freqpoly() + labs(x="Methylation Difference (Girl - Boy)") + facet_grid(Chromosome ~ .) + ggtitle(i)

	# meth profile over the start position
	#g<-ggplot(rnb.diff.met[[i]], aes(x=Start, y=mean.mean.diff)) + geom_point(aes(colour=mean.mean.diff),alpha=0.8) + scale_colour_gradient2() + facet_grid(Chromosome ~ .)

	# volcano plot by chromosome facet
	#g<-ggplot(rnb.diff.met[[i]][select,], aes_string(x=my.x,y=my.y)) + geom_point(alpha=0.3) + labs(x="Methylation Difference (Girl - Boy)") + facet_grid(Chromosome ~ .) + ggtitle(i)

	# beta value distribution by sex and chr facet
	if(FALSE){
		dummy<-rbind(
			setNames(cbind(rnb.diff.met[[i]][select,c("Chromosome","Start","End",my.ctrl)], "Boy"),  c("Chromosome","Start","End","Meth","Sex")),
			setNames(cbind(rnb.diff.met[[i]][select,c("Chromosome","Start","End",my.case)], "Girl"), c("Chromosome","Start","End","Meth","Sex"))
			)
		g2<-ggplot(dummy, aes(Meth)) + geom_density(aes(colour=Sex), alpha=.7) + facet_grid(Chromosome ~ .)
	}
	print(g)
	print(g1)

	# save as PT.dt.merged
	# to make the format compatible with the RoadMap data
	# mean.covg.g1: female
	# mean.covg.g2: male
	if(i=="sites"){
		query<-rnb.diff.met[[i]] #isa 'data.frame'
		query$Chromosome<-simplify2array(strsplit(as.character(query$Chromosome), "chr"))[2,] # chrX = > X
		PT.CG.dt.merged<-with(query, data.table(V1=Chromosome,V2=ifelse(Strand=='+',Start,End), V3=Strand, V5.x=round(mean.g1*mean.covg.g1,1), V6.x=mean.covg.g1, V5.y=round(mean.g2*mean.covg.g2,1), V6.y=mean.covg.g2))
		my.new.RData<-file.path("~/data/RoadMap/BS-Seq/Schultz2015/PT.CG.dt.merged.RData")
		save(PT.CG.dt.merged, file=my.new.RData)
	}
}
dev.off()

##################################################
# proportion of sites/regions                   ##
# differentially methylated in case and control ##
##################################################
file.name<-file.path(run.dir,paste("diff.ratio.by.chr",my.project,my.run.date,"pdf",sep=".")) # this will store 'rnb.diff.met' which is a list of data.frame for sites, promoters...
pdf(file=file.name, width=11.7, height=8.3, title=my.project)

bin.stat<-list() # a list of list 
my.diffs=c(0, .05, .10, .15, .20, .25) # meth diff more than or equal to 10%
for(my.diff in my.diffs){
	for(i in names(rnb.diff.met)){
		select=rnb.diff.met[[i]]$Chromosome!="chrY"
		my.x=ifelse(i=="sites","mean.diff","mean.mean.diff") # meth(girl)-meth(boy)

		dummy<-rnb.diff.met[[i]][select,]
		dummy$Chromosome<-droplevels(dummy$Chromosome) # adjust level (drop unused factor)

		bin.stat[[paste0("diff>",my.diff)]][[i]]<-t(
						 sapply(split(dummy, dummy$Chromosome), 
								function(j){
									filter=abs(j[,my.x])>=my.diff
									table(factor(j[filter,my.x]>0, levels=c(F,T))) # TRUE if girls>boys, FALSE OW
									#table(j[filter,my.x]>0)
								} )
						)

		if(TRUE){
			g<-ggplot(melt(bin.stat[[paste0("diff>",my.diff)]][[i]]/rowSums(bin.stat[[paste0("diff>",my.diff)]][[i]])), aes(x=Var1,y=value,fill=Var2)) + 
				geom_bar(stat="identity") + geom_hline(aes(yintercept=.5)) + labs(x="Chromosome",y="Ratio") + 
				#coord_polar(theta = "y") +
				scale_fill_manual(values=cbPalette[1:2], name="Methylation\nDifference", breaks=c("FALSE", "TRUE"),labels=c("Boy>=Girl", "Boy<Girl")) +
				ggtitle(paste0(i,": diff>",my.diff))
			print(g)
		}
	}
}
dev.off()

#####################################################################
# count sites/regions differentially methylated in case and control #
# below for chrX only                                               #
#####################################################################
file.name<-file.path(run.dir,paste("cnt.cpg.by.meth.diff",my.project,my.run.date,"pdf",sep=".")) # this will store 'rnb.diff.met' which is a list of data.frame for sites, promoters...
pdf(file=file.name, width=11.7, height=8.3, title=my.project)
for(j in names(rnb.diff.met)){
	q<-ggplot(melt(sapply(bin.stat, function(i) i[[j]]["chrX",])), aes(x=Var2,y=value,fill=Var1)) + 
		geom_bar(stat="identity") +
		labs(x="Minimum Methylation Difference",y="Number of sites/regions in chrX") + 
		scale_fill_manual(values=cbPalette[1:2], name="Methylation\nDifference", breaks=c("FALSE", "TRUE"),labels=c("Boy>=Girl", "Boy<Girl")) +
		ggtitle(j)
	print(q)
}
dev.off()

#sapply(bin.stat, function(i) i[["cpgislands"]]["chrX",])

# Median methylation difference between case and control
#lapply(rnb.diff.met[1], function(j) summary(j$mean.g2)["Median"] - summary(j$mean.g1)["Median"])
#lapply(rnb.diff.met[2:6], function(j) summary(j$mean.mean.g2)["Median"] - summary(j$mean.mean.g1)["Median"])

################################################################
## Enrichment of CPGi where boys|girls are more methylated    ##
## in Chromosome X                                            ##
## Note that rnb.diff.met stores 0-based coordiate (bed file) ##
################################################################
fields <- c("ensembl_gene_id", "hgnc_symbol","chromosome_name","strand", "description","gene_biotype") # get the gene names

i="cpgislands"
my.upflank=500 #10^4 # 10kb up from the start 
my.downflank=500 #10^4 # 10kb down from the end
my.diff=.1 # meth diff more than or equal to 10%
my.x=ifelse(i=="sites","mean.diff","mean.mean.diff") # meth(girl)-meth(boy)

if(FALSE){
	# extends a bit (isa GRanges)
	# gr.ensg from config/Annotation.R
	gr.ext.ensg<-unlist(GRangesList(lapply(gr.ensg, function(i) promoters(i, upstream=2000, downstream=end(i)-start(i)+1+100)))) # ncludes 2000 +TSS and 100 +TES
}
#
# 1. genes of interests: genes overlapped with cpgi and chrX where boys are more methylated
#
select=rnb.diff.met[[i]]$Chromosome=="chrX" & abs(rnb.diff.met[[i]][,my.x]) >= my.diff# chrX only & at least this difference
dummy<-rnb.diff.met[[i]][select,]
dummy$Chromosome<-droplevels(dummy$Chromosome) # adjust level (drop unused factor)
filter=dummy[,my.x] < 0 # girl-boy < 0 In other words, boys are more methylated  
#my.gr<-data.frame2GRanges(dummy[filter,], chrom.column=1,start.column=2,end.column=3) # data.frame2GRanges is from 'RnBeads'
my.gr=makeGRangesFromDataFrame(dummy[filter,],keep.extra.columns=TRUE) # from 'GenomicRanges' package
my.gr.ext<-unlist(
				  GRangesList(lapply(my.gr, function(i) promoters(i, upstream=my.upflank, downstream=end(i)-start(i)+1+my.downflank))) #includes upflank +TSS and downflank from +TES 
				 )
my.overlap=as.matrix( findOverlaps(gr.ensg, my.gr.ext, ignore.strand=TRUE, minoverlap=10L) ) #By default, any overlap is accepted

my.new.gr<-my.gr.ext[my.overlap[,"subjectHits"]]
mcols(my.new.gr)<-DataFrame(mcols(my.new.gr), `ensembl_gene_id`=names(gr.ensg[my.overlap[,"queryHits"]]))

my.target.genes<-unique(mcols(my.new.gr)$ensembl_gene_id )
my.gene.info=getBM(attributes = fields, filters = "ensembl_gene_id", values=my.target.genes, mart = grch37) # is a data.frame (grch37 from config/Annotation.R)

file.name<-file.path(run.dir,paste(i,"girls.more.meth",my.project,my.run.date,"csv",sep=".")) # this will store 'rnb.diff.met' which is a list of data.frame for sites, promoters...
write.csv(dummy[filter,], file=file.name, row.names=F)
file.name<-file.path(run.dir,paste(i,"near.gene.girls.more.meth",my.project,my.run.date,"csv",sep=".")) # this will store 'rnb.diff.met' which is a list of data.frame for sites, promoters...
write.csv(my.gene.info, file=file.name, row.names=F)

# merge methylation and expression
# res from 'bin/R/DEG/DESeq2.complex.design.R'
my.exp<-res[rownames(res) %in% mcols(my.new.gr)$ensembl_gene_id,]
my.exp<-data.frame(my.exp, `ensembl_gene_id`=rownames(my.exp))
dummy<-merge(mcols(my.new.gr), my.exp) 
select<-!is.na(dummy$log2FoldChange) & !is.na(dummy$padj) & dummy$baseMean>=minRead
ggplot(as.data.frame(dummy[select,]), aes(-mean.mean.diff,log2FoldChange)) + geom_point(aes(color=baseMean),size=5, type=1) + labs(x="met(M)-met(F)",y="exp(M)-exp(F)")
#getTopDeseq(dummy[select,], "met(boy) > met(girl)")

#
# 2. background: genes overlapped with cpgi
#
select=rnb.diff.met[[i]]$Chromosome!="chrY" # not chrY
dummy<-rnb.diff.met[[i]][select,]
dummy$Chromosome<-droplevels(dummy$Chromosome) # adjust level (drop unused factor)
filter=TRUE # all cpgi region
#my.gr<-data.frame2GRanges(dummy[filter,], chrom.column=1,start.column=2,end.column=3) # data.frame2GRanges is from 'RnBeads'
my.gr=makeGRangesFromDataFrame(dummy[filter,],keep.extra.columns=TRUE)
my.gr.ext<-unlist(
				  GRangesList(lapply(my.gr, function(i) promoters(i, upstream=my.upflank, downstream=end(i)-start(i)+1+my.downflank))) #includes upflank +TSS and downflank from +TES 
				 )
my.overlap=as.matrix( findOverlaps(gr.ensg, my.gr.ext, ignore.strand=TRUE, minoverlap=10L) ) #By default, any overlap is accepted
my.bg.ensg<-gr.ensg[my.overlap[,"queryHits"]] # isa 'GRanges'

# parpare genes for GO analysis
myGenes<-as.integer(unique(mcols(my.bg.ensg)$gene_id) %in% my.target.genes)
names(myGenes)<-unique(mcols(my.bg.ensg)$gene_id)

myGenome="hg19"
pwf=goseq::nullp(myGenes,myGenome,"ensGene") # isa data.frame  (check available ID via supportedGeneIDs())
GO.wall=goseq::goseq(pwf,myGenome,"ensGene") # isa data.frame
enriched.GO.wall=list('over'=GO.wall[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.001,], # isa data.frame
					  'under'=GO.wall[p.adjust(GO.wall$under_represented_pvalue, method="BH")<.001,] # isa data.frame
				)
# Find the level from the root (depth-first or shortest) given a GO term
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
# list to data.frame
cat("writing GO analysis results...\n")
write.table(enriched.GO.wall[["over"]][enriched.GO.wall[["over"]]$depth-1>=5,], file=file.path(top.dir, paste0("Cuffcompare/",myType,"/",myType, ".goseq.coding.95.vs.5.over.csv")))
write.table(enriched.GO.wall[["under"]][enriched.GO.wall[["under"]]$depth-1>=5,], file=file.path(top.dir, paste0("Cuffcompare/",myType,"/",myType, ".goseq.coding.95.vs.5.under.csv")))

############################
## GO terms to Word Cloud ##
## for those genes above  ##
############################
cat("making WordCloud...\n")
library(tm)
library(wordcloud)
library(SnowballC)

terms<-as.character(unlist(sapply(enriched.GO.wall[["over"]][enriched.GO.wall[["over"]]$depth-1>=5,"term"], function(i) strsplit(i, " ")))) # isa character vector
terms<-tm::Corpus(VectorSource(terms)) # isa 'VCorpus'

terms<- tm_map(terms, content_transformer(tolower))
terms<- tm_map(terms, removePunctuation)
terms<- tm_map(terms, removeNumbers)
terms<- tm_map(terms, removeWords, stopwords("english"))
terms<- tm_map(terms, stripWhitespace)
#terms<- tm_map(terms, stemDocument) #install.packages("SnowballC")
wordcloud(terms, scale=c(5,0.5), max.words=100, random.order=FALSE, rot.per=0.35, use.r.layout=FALSE, colors=brewer.pal(8, "Dark2")) # brewer.pal from RColorBrewer

#############################
## Promoters & Cpgislands ###
#############################
#gr.promoters=data.frame2GRanges(rnb.diff.met[["promoters"]], chrom.column=1,start.column=2,end.column=3) # from RnBeads
gr.promoters=makeGRangesFromDataFrame(rnb.diff.met[["promoters"]],keep.extra.columns=TRUE) # from GenomicRanges
#gr.cpgi=data.frame2GRanges(rnb.diff.met[["cpgislands"]], chrom.column=1,start.column=2,end.column=3)
gr.cpgi=makeGRangesFromDataFrame(rnb.diff.met[["cpgislands"]],keep.extra.columns=TRUE)
my.overlap=as.matrix( findOverlaps(gr.promoters, gr.cpgi, ignore.strand=TRUE, minoverlap=10L) ) #By default, any overlap is accepted
my.promoters<-gr.promoters[my.overlap[,"queryHits"]] # isa 'GRanges' (overlap between promoters and cpgi)

bin.stat<-list() # a list of list 
my.diffs=c(0, .05, .10, .15, .20, .25) # meth diff more than or equal to 10%
for(my.diff in my.diffs){
	select=my.promoters$Chromosome!="chrY"
	my.x="mean.mean.diff"

	dummy<-my.promoters[select,]
	dummy$Chromosome<-droplevels(dummy$Chromosome) # adjust level (drop unused factor)

	bin.stat[[paste0("diff>",my.diff)]][[i]]<-t(
						sapply(split(dummy, dummy$Chromosome), 
							function(j){
								filter=abs(j[,my.x])>=my.diff
								table(factor(j[filter,my.x]>0, levels=c(F,T))) # TRUE if girls>boys, FALSE OW
								#table(j[filter,my.x]>0)
							} )
					)

	if(TRUE){
		g<-ggplot(melt(bin.stat[[paste0("diff>",my.diff)]][[i]]/rowSums(bin.stat[[paste0("diff>",my.diff)]][[i]])), aes(x=Var1,y=value,fill=Var2)) + 
			geom_bar(stat="identity") + geom_hline(aes(yintercept=.5)) + labs(x="Chromosome",y="Ratio") + 
			#coord_polar(theta = "y") +
			scale_fill_manual(values=cbPalette[1:2], name="Methylation\nDifference", breaks=c("FALSE", "TRUE"),labels=c("Boy>=Girl", "Boy<Girl")) +
			ggtitle(paste0(i,": diff>",my.diff))
		print(g)
	}
}

my.target.genes<-unique(mcols(my.target.ensg)$gene_id) 
my.gene.info=getBM(attributes = fields, filters = "ensembl_gene_id", values=my.target.genes, mart = grch37) # is a data.frame

