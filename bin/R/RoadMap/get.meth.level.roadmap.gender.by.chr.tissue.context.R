#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 3/Mar/2016
# Last modified 17/Aug/2016
# Modified from bin/R/RoadMap/plot.roadmap.gender.meth.diff.by.chr.tissue.context.R
# Modified from bin/R/RoadMap/plot.roadmap.gender.meth.distribution.by.chr.tissue.context.R
# use bin/R/RoadMap/count.roadmap.cpg.R for CX and region counts

TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/RoadMap/local.R")

################################
## RoadMap                    ##
## Schultz et al. Nature 2015 ##
## Tissue comparison          ##
## Female-Male %Met Difference##
################################
my.cpg.type="CG" # CG, CHH, CHG
# cpg.context is fixed for the RData
#cpg.contexts=c("Sites","Genes","Promo_15.05","Promo_15.00","Promo_10.10","Promo_10.03", "Promo_10.00", "CPGi","CPGi-with-promo","CPGi-no-promo","CPGi-shores","Enhancer") # 18/Aug/2016 (do not change)
cpg.contexts=c("Sites","Genes","Promo_15.05","CPGi")
# 9 tissues ("PT", "AD","AO","EG","FT","GA","PO","SB","SX"): (7/March/2016)
# PA: not available yet
# PT: for our Placenta
# ~/results/RoadMap/BS-Seq/RData/10X/CG.dt.meth.per.tissue.per.region.RData stores following:
# a summary of methylation-level and count by chromosome, region and tissue (18/Oct/2016)
#my.RData=file.path("~/results/RoadMap/BS-Seq/RData/10X",paste(my.cpg.type,"dt.meth.per.tissue.per.region.RData",sep="."))
my.RData=file.path("~/results/RoadMap/BS-Seq/RData/10X",paste(my.cpg.type,"PT.SS.dt.meth.per.tissue.per.region.RData",sep="."))
if(file.exists(my.RData)){
	load(my.RData)
	cat("dt.meth.per.tissue loaded\n")
}else{
	dt.meth.per.tissue=list() # to save %met by cpg contexts
	#for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){
	for(my.tissue in c("PT.SS.Tech","PT.SS.Bio")){
		## Load this RoadMap tissue
		dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)
		dt.query=dt.query[V6.x>=10 & V6.y>=10] # apply CpG depth for PT.SS.Tech or PT.SS.Bio

		dt.meth.per.region=list() # to save %met by cpg contexts
		for(i in cpg.contexts){
			cat(paste0("Processing ", i, "...\n"))
			if(i=='Sites'){
				dt.meth<-rbind(
							dt.query[,list(tissue=my.tissue, context=i, num.sites=.N, f.met=round(sum(V5.x)/sum(V6.x)*100,2), m.met=round(sum(V5.y)/sum(V6.y)*100,2), met=round(sum(V5.x+V5.y)/sum(V6.x+V6.y)*100,2), f.c=sum(V5.x), f.t=sum(V6.x)-sum(V5.x), m.c=sum(V5.y), m.t=sum(V6.y)-sum(V5.y)), by="V1"],
							dt.query[,list(V1="genome",tissue=my.tissue, context=i, num.sites=.N, f.met=round(sum(V5.x)/sum(V6.x)*100,2), m.met=round(sum(V5.y)/sum(V6.y)*100,2), met=round(sum(V5.x+V5.y)/sum(V6.x+V6.y)*100,2), f.c=sum(V5.x), f.t=sum(V6.x)-sum(V5.x), m.c=sum(V5.y), m.t=sum(V6.y)-sum(V5.y))]
							)
			}else{
				subject<-with(as.data.frame(get.region(i)), data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand)) #isa 'data.frame'

				cat("\tRemoving leading 'chr'...\n")
				subject$Chromosome<-simplify2array(strsplit(as.character(subject$Chromosome), "chr"))[2,] # chrX = > X
				cat("\tconvert to data.table...\n")
				dt.subject<-as.data.table(subject) # isa data.table
				rm(subject)

				if(grepl('Gene',i) || grepl('Promo',i)){
					subject.keys=c("Chromosome","Strand","Start","End") # regions of interests 
					query.key=c("V1","V3","V2","End")                   # columns from dt.query (roadmap merged data)
				# cpgislands, Tiling, genomeTiling500
				}else{
					subject.keys=c("Chromosome","Start","End")
					query.key=c("V1","V2","End")
				}
				setkeyv(dt.subject,subject.keys)

				cat("\tFinding CpGs overlapping with ",i,"...\n")
				system.time(dt.overlap<-foverlaps(dt.query, dt.subject, by.x=query.key, type="any", nomatch=0L))

				cat("\tAggregating %met...\n")
				dt.meth<-rbind(
							dt.overlap[,list(tissue=my.tissue, context=i, num.sites=.N, f.met=round(sum(V5.x)/sum(V6.x)*100,2), m.met=round(sum(V5.y)/sum(V6.y)*100,2), met=round(sum(V5.x+V5.y)/sum(V6.x+V6.y)*100,2), f.c=sum(V5.x), f.t=sum(V6.x)-sum(V5.x), m.c=sum(V5.y), m.t=sum(V6.y)-sum(V5.y)), by="V1"],
							dt.overlap[,list(V1="genome",tissue=my.tissue, context=i, num.sites=.N, f.met=round(sum(V5.x)/sum(V6.x)*100,2), m.met=round(sum(V5.y)/sum(V6.y)*100,2), met=round(sum(V5.x+V5.y)/sum(V6.x+V6.y)*100,2), f.c=sum(V5.x), f.t=sum(V6.x)-sum(V5.x), m.c=sum(V5.y), m.t=sum(V6.y)-sum(V5.y))]
							)
				rm(dt.overlap)
			}# end if i=='Sites'
			#       V1 tissue context     num f.met m.met   met       f.c        f.t       m.c        m.t
			# 1:      1     PT    CPGi  170182 11.91 12.84 12.38  419366.0  3100743.5  466377.8  3166981.7
			# 2:      7     PT    CPGi  103374 19.32 20.26 19.80  407558.3  1701437.7  441174.4  1736020.1
			# 3:      4     PT    CPGi   81875 14.55 15.62 15.11  248553.7  1459244.3  285036.5  1539211.5
			#24: genome     PT    CPGi 1918675 15.71 16.62 16.17 6205138.9 33299524.6 6741399.3 33832692.7
			dt.meth.per.region[[i]]=dt.meth
		}# end of cpg.contexts
		dt.meth.per.tissue[[my.tissue]]<-rbindlist(dt.meth.per.region) # rbind by tissue
		cat(paste0("Done for ", my.tissue, "...\n"))
	}# end of my.tissue 
	cat("saving dt.meth.per.tissue...\n")
	save(dt.meth.per.tissue, file=my.RData)
}

if(FALSE){
	cat("Merging all samples...\n")
	dt.meth.all<-data.table::rbindlist(dt.meth.per.tissue)

	file.name<-file.path("~/results/RoadMap/BS-Seq",paste(my.cpg.type,"meth.level",time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=11.7, height=8.3, title=paste0(my.cpg.type, "Methylation Level of Various Tissues by Gender"))

	# Boxplot of meth level by tissue
	p<-ggplot(dt.meth.all[V1=="genome" & context %in% c("Sites","Genes","Promo_15.05","CPGi","CPGi-shores","Enhancer")], aes(tissue, met)) + 
		geom_bar(aes(fill=context), stat="identity",position="dodge") + 
		ggtitle("genome-wide methylation in various context") +
		theme_bw() +
		scale_fill_manual(values=cbPalette, name=paste0(my.cpg.type,"\nregional contexts"))
	print(p)

	# Boxplot of meth level @ promoter by tissue
	p<-ggplot(dt.meth.all[V1=="genome" & context %in% c("CPGi", "Promo_15.05","Promo_15.00","Promo_10.10","Promo_10.03","Promo_10.00")], aes(tissue, met)) + 
		geom_bar(aes(fill=context), stat="identity",position="dodge") + 
		ggtitle("genome-wide methylation @ various promoters") +
		theme_bw() +
		scale_fill_manual(values=cbPalette, name=paste0(my.cpg.type,"\nregional contexts"))
	print(p)

	# Boxplot of genome-wide meth level @ CpGi by tissue
	p<-ggplot(dt.meth.all[V1=="genome" & context %in% c("CPGi", "CPGi-with-promo", "CPGi-no-promo", "Promo_10.03")], aes(tissue, met)) + 
		geom_bar(aes(fill=context), stat="identity",position="dodge") + 
		ggtitle("genome-wide methylation @ CPGi and nearby") +
		theme_bw() +
		scale_fill_manual(values=cbPalette, name=paste0(my.cpg.type,"\nregional contexts"))
	print(p)

	# Boxplot of genome-wide meth level diff by gender
	p<-ggplot(dt.meth.all[V1=="genome" & context %in% c("Sites","Genes","Promo_15.05","CPGi"),list(tissue,context,meth.diff=f.met-m.met)], aes(tissue, meth.diff)) +
		geom_bar(aes(fill=context), stat="identity",position="dodge") + 
		ggtitle("genome-wide methylation difference by gender at various regions") +
		theme_bw() +
		scale_fill_manual(values=cbPalette, name=paste0(my.cpg.type,"\nregional contexts"))
	print(p)

	# Boxplot of chrX meth level diff by gender
	p<-ggplot(dt.meth.all[V1=="X" & context %in% c("Sites","Genes","Promo_15.05","CPGi","CPGi-shores","Enhancer"),list(tissue,context,meth.diff=f.met-m.met)], aes(tissue, meth.diff)) +
		geom_bar(aes(fill=context), stat="identity",position="dodge") + 
		ggtitle("methylation difference (ChrX only) by gender at various regions") +
		theme_bw() +
		scale_fill_manual(values=cbPalette, name=paste0(my.cpg.type,"\nregional contexts"))
	print(p)

	# Boxplot of meth level for PT 
	p<-ggplot(dt.meth.all[V1=="genome" & tissue=="PT" & context %in% c("Sites","Genes","Promo_15.05","CPGi","CPGi-shores","Enhancer")], aes(context, met)) + 
		geom_bar(aes(fill=context), stat="identity",position="dodge") + 
		ggtitle("genome-wide methylation of PT at various regions") +
		theme_bw() +
		scale_fill_manual(values=cbPalette, name=paste0(my.cpg.type,"\nregional contexts"))
	print(p)

	# Boxplot of meth level for PT by gender
	foo<-rbind(dt.meth.all[,list(tissue, Chr=V1, context, Gender="Female",`% methylation`=f.met)], dt.meth.all[,list(tissue, Chr=V1, context, Gender="Male",`% methylation`=m.met)])
	my.tissue="PT"; my.chr="genome"
	p<-ggplot(foo[tissue==my.tissue & Chr==my.chr & context %in% c("Sites","Genes","Promo_10.03","CPGi","CPGi-shores","Enhancer","CPGi-with-promo","CPGi-no-promo")], aes(context, `% methylation`)) +
		geom_bar(aes(fill=Gender), stat="identity",position="dodge") + 
		ggtitle(paste(my.tissue, my.chr, "methylation difference by gender at various regions")) +
		theme_bw() +
		scale_fill_manual(values=my.col[["Gender"]])
	print(p)

	# Boxplot of genome-wide meth level by tissue and gender
	for(my.context in cpg.contexts){
		p<-ggplot(foo[Chr=="genome" & context==my.context], aes(tissue, `% methylation`)) + 
			geom_bar(aes(fill=Gender),stat="identity",position="dodge") +
			ggtitle(paste("Genome-wide methylation at",my.context)) +
			theme_bw() +
			scale_fill_manual(values=my.col[["Gender"]])
		print(p)
	}

	# Boxplot of chrX meth level by tissue and gender
	for(my.context in cpg.contexts){
		p<-ggplot(foo[Chr=="X" & context==my.context], aes(tissue, `% methylation`)) + 
			geom_bar(aes(fill=Gender),stat="identity",position="dodge") +
			ggtitle(paste("ChrX methylation at",my.context)) +
			theme_bw() +
			scale_fill_manual(values=my.col[["Gender"]])
		print(p)
	}

	# Boxplot of meth level by gender by chr
	#for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){
	#	for(my.context in cpg.contexts){
	#		#my.tissue="PT"; my.context="Sites"
	#		p<-ggplot(foo[Chr!="genome" & tissue==my.tissue & context==my.context], aes(Chr, `% methylation`)) +
	#			geom_bar(aes(fill=Gender), stat="identity",position="dodge") + 
	#			ggtitle(paste(my.tissue, "at", my.context)) +
	#			theme_bw() +
	#			scale_fill_manual(values=my.col[["Gender"]])
	#		print(p)
	#	}
	#}
	dev.off()
}

#write.csv(rbindlist(dt.meth.per.tissue)[V1=="genome"], file=file.path("~/results/RoadMap/BS-Seq",paste(my.cpg.type,"PT.SS.meth.level.csv",sep=".")), row.names=FALSE)
#write.csv(rbindlist(dt.meth.per.tissue)[V1=="X"], file=file.path("~/results/RoadMap/BS-Seq",paste(my.cpg.type,"PT.SS.meth.level.chrX.csv",sep=".")), row.names=FALSE)

cat("All done\n")
