#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 3/Mar/2016
# Last modified 17/Aug/2016
# Modified from bin/R/RoadMap/plot.roadmap.gender.meth.diff.by.chr.tissue.context.R

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
plot.per.tissue=TRUE
# 9 tissues ("PT", "AD","AO","EG","FT","GA","PO","SB","SX"): (7/March/2016)
# PA: not available yet
# PT: for our Placenta
# ~/results/RoadMap/BS-Seq/RData/CG.dt.meth.per.tissue.RData have following regions by gender-specific way:
# "Sites","Tiling5K","Genes","Promo_10.05","CPGi","CPGi-shores" (16/Oct/2016)
my.RData=file.path("~/results/RoadMap/BS-Seq/RData",paste(my.cpg.type,"dt.meth.per.tissue.RData",sep="."))
if(file.exists(my.RData)){
	cat("loading dt.meth.per.tissue.RData...\n")
	load(my.RData)
	cat("dt.meth.per.tissue loaded\n")
}else{
	dt.meth.per.tissue=list() # to save %met by cpg contexts
	for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){
		## Load this RoadMap tissue
		dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)

		#cpg.contexts=c("Sites","Tiling5K","Genes","Promo_10.10","CPGi","Enhancer", "CPGi-shores")
		cpg.contexts=c("Sites","Tiling5K","Genes","Promo_15.05","CPGi","Enhancer", "CPGi-shores")
		dt.meth.per.region=list() # to save %met by cpg contexts
		for(i in cpg.contexts){
			cat(paste0("Processing ", i, "...\n"))
			if(i=='Sites'){
				dt.meth<-rbind(dt.query[,.(V1,V3,V2,End,1,V5.x/V6.x,'cnt.C'=V5.x,'cnt.T'=V6.x-V5.x,"Female",i,my.tissue)], 
								dt.query[,.(V1,V3,V2,End,1,V5.y/V6.y,'cnt.C'=V5.y,'cnt.T'=V6.y-V5.y,"Male",i,my.tissue)])
			}else{
				subject<-with(as.data.frame(get.region(i)), data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand)) #isa 'data.frame'

				cat("\tRemoving leading 'chr'...\n")
				#subject$Chromosome<-simplify2array(strsplit(as.character(subject$Chromosome), "chr"))[2,] # chrX = > X
				subject$Chromosome<-substr(subject$Chromosome, 4, 5) #chrX=>X, chr13=>13
				cat("\tconvert to data.table...\n")
				dt.subject<-as.data.table(subject) # isa data.table
				rm(subject)

				if(grepl('Gene',i) || grepl('Promo',i)){
					subject.keys=c("Chromosome","Strand","Start","End") # regions of interests 
					query.key=c("V1","V3","V2","End")                   # columns from dt.query (roadmap merged data)
					group.keys=c("V1","V3","Start","End")               # V3: strand
					group.keys.list="V1,V3,Start,End"                   # same as above 
				# cpgislands, Tiling, genomeTiling500
				}else{
					subject.keys=c("Chromosome","Start","End")
					query.key=c("V1","V2","End")
					group.keys=c("V1","Strand","Start","End")
					group.keys.list="V1,Strand,Start,End"
				}
				setkeyv(dt.subject,subject.keys)

				cat("\tFinding CpGs overlapping with ",i,"...\n")
				print(system.time(dt.overlap<-foverlaps(dt.query, dt.subject, by.x=query.key, type="any", nomatch=0L)))

				cat("\tAggregating %met...\n")
				dt.meth<-rbind(
							dt.overlap[,list(.N,sum(V5.x)/sum(V6.x),'cnt.C'=sum(V5.x),'cnt.T'=sum(V6.x)-sum(V5.x),"Female",i,my.tissue),by=group.keys.list],
							dt.overlap[,list(.N,sum(V5.y)/sum(V6.y),'cnt.C'=sum(V5.y),'cnt.T'=sum(V6.y)-sum(V6.y),"Male",i,my.tissue),by=group.keys.list]
						)
				#> dt.meth 
				#V1 V3        V2       End        V5     V6     i my.tissue
				#1: 10  +     66269     66269 1.0000000 Female Sites        AD
				#2: 10  +     66273     66273 1.0000000 Female Sites        AD
				#---
				#58485806:  X  - 154929107 154929107 1.0000000   Male Sites        AD
				#58485807:  X  + 154929126 154929126 1.0000000   Male Sites        AD
				rm(dt.overlap)
			}# end if i=='sites'
			setnames(dt.meth,c("chr","strand","start","end","num.sites","meth","cnt.C","cnt.T","Gender","Region","Tissue"))

			# avoid storing single CH sites
			# this is to avoid 'too large for hashing' error of rbindlist
			if(grepl("CH",my.cpg.type) & i=="Sites"){cat("pass storing dt.meth\n")}else{dt.meth.per.region[[i]]=dt.meth}
		}# end of cpg.contexts

		# Error in data.table::rbindlist(dt.meth.per.region) :
		# length 1220193705 is too large for hashing
		# see https://github.com/Rdatatable/data.table/blob/master/src/rbindlist.c
		# there seems to be a limit up to 1073741824
		dt.meth.merged<-data.table::rbindlist(dt.meth.per.region)
		dt.meth.per.tissue[[my.tissue]]<-dt.meth.merged
		cat(paste0("Done for ", my.tissue, "...\n"))
	}# end of my.tissue 
	cat("saving dt.meth.per.tissue...\n")
	save(dt.meth.per.tissue,file=my.RData) 
}

my.tissue="PT"
if(FALSE){
	file.name<-file.path("~/results/RoadMap/BS-Seq/Beta-distribution/Tissue/",paste(my.tissue,my.cpg.type,time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=11.7, height=8.3, title=paste("methylation value distrubution by sex in ",my.tissue))
	for(i in dt.meth.per.tissue[[my.tissue]][,unique(Region)]){
		dt.meth=dt.meth.per.tissue[[my.tissue]][Region==i]
		# beta value (%met) distribution of this cpg context
		p1<-ggplot(dt.meth, aes(meth, col=Gender)) + 
			geom_density(size=1.2) + 
			ggtitle(paste(my.tissue, i,sep=":")) + 
			theme_bw() +
			labs(x="Methylation level", y="Density") +
			scale_colour_manual(values=my.col[["Gender"]])
		print(p1+theme_Publication())

		# beta value (%met) distribution of this cpg context in autosome only (no chrX)
		p2<-ggplot(dt.meth[chr!="X"], aes(meth, col=Gender)) + 
			geom_density(size=1.2) + 
			ggtitle(paste(my.tissue, "autosomes",i, sep=":")) + 
			theme_bw() +
			labs(x="Methylation level", y="Density") +
			scale_colour_manual(values=my.col[["Gender"]])
		print(p2+theme_Publication())

		# beta value (%met) distribution of this cpg context in chrX only
		p3<-ggplot(dt.meth[chr=="X"], aes(meth, col=Gender)) + 
			geom_density(size=1.2) + 
			ggtitle(paste(my.tissue, "chrX",i, sep=":")) + 
			theme_bw() +
			labs(x="Methylation level", y="Density") +
			scale_colour_manual(values=my.col[["Gender"]])
		#print(p3+theme_Publication())
	}
	# beta value (%met) distribution of this cpg context by chromosome and cpg context
	p4<-ggplot(dt.meth.per.tissue[[my.tissue]], aes(meth, col=Gender)) + 
		geom_density() + 
		ggtitle(my.tissue) +
		scale_colour_manual(values=my.col[["Gender"]]) + 
		theme_bw() + 
		theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
		facet_grid(Region ~ chr)
	#print(p4)
	# a matrix of beta-distribution (tissue by region)
	dev.off()
	if(!plot.per.tissue){
		cat("Mering all samples...\n")
		dt.meth.all<-data.table::rbindlist(dt.meth.per.tissue) # too big (not memory efficient)?

		file.name<-file.path("~/results/RoadMap/BS-Seq/Beta-distribution",paste("all.tissue.beta.dist",my.cpg.type,time.stamp,"pdf",sep="."))
		pdf(file=file.name, width=11.7, height=8.3, title=paste("methylation value distrubution by sex in ",my.tissue))

		p5<-ggplot(dt.meth.all, aes(meth, col=Gender)) + 
			geom_density() + 
			ggtitle("Methylation Level Distribution by Tissue") +
			scale_colour_manual(values=my.col[["Gender"]]) + 
			theme_bw() + 
			theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
			facet_grid(Region ~ Tissue) # Row: Region, Column:Tissue 
		print(p5)

		p6<-ggplot(dt.meth.all[chr=="X"], aes(meth, col=Gender)) + 
			geom_density() + 
			ggtitle("Methylation Level Distribution of ChrX by Tissue") +
			scale_colour_manual(values=my.col[["Gender"]]) + 
			theme_bw() + 
			theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
			facet_grid(Region ~ Tissue) # Row: Region, Column:Tissue 
		print(p6)
		dev.off()
	}
} # FALSE

	########################
	## Beta Distribution  ##
	## via plot_grid      ##
	########################
	library("cowplot") # only for R >= 3.1.2
	dt.meth=dt.meth.per.tissue[[my.tissue]][Region=="Sites"]
	p.site.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_density(size=1.2) + 
			labs(x="",y="Density") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			#ylim(c(0,0.062)) +
			theme_Publication() +
			theme(legend.position="none")

	dt.meth=dt.meth.per.tissue[[my.tissue]][Region=="Genes"]
	p.gene.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_density(size=1.2) + 
			labs(x="",y="") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			#ylim(c(0,0.062)) +
			theme_Publication() +
			#theme(legend.position="top") 
			theme(legend.position=c(0.91,0.91)) 
	
	dt.meth=dt.meth.per.tissue[[my.tissue]][Region=="Promo_15.05"]
	p.promoter.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_density(size=1.2) + 
			labs(x="% Methylation Level", y="Density") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			#ylim(c(0,0.062)) +
			theme_Publication() +
			theme(legend.position="none")

	dt.meth=dt.meth.per.tissue[[my.tissue]][Region=="CPGi"]
	p.cpgi.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_density(size=1.2) + 
			labs(x="% Methylation Level", y="") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			#ylim(c(0,0.062)) +
			theme_Publication() + 
			theme(legend.position="none")

	file.name<-file.path("~/results/RoadMap/BS-Seq/Beta-distribution/Tissue/",paste("Fig1",my.tissue,my.cpg.type,time.stamp,"tiff",sep="."))
	tiff(filename=file.name,width=10, height=8,units="in",res=300, compression = 'lzw') #A4 size 
	cowplot::plot_grid(p.site.beta, p.gene.beta, p.promoter.beta, p.cpgi.beta, labels=c("A","B","C","D"), ncol = 2, nrow = 2)
	dev.off()

	#######################
	## Meth Diff Boxplot ##
	## via plot_grid     ##
	#######################
	library("cowplot") # only for R >= 3.1.2
	dt.meth=dt.meth.per.tissue[[my.tissue]][Region=="Sites"]
	p.site.box<-ggplot(dcast(dt.meth, chr+strand+start+end~Gender, value.var="meth"), aes(chr, (Female-Male)*100)) + 
			geom_boxplot() + 
			labs(x="",y="Density") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			#ylim(c(0,0.062)) +
			theme_Publication() +
			theme(legend.position="none")

	dt.meth=dt.meth.per.tissue[[my.tissue]][Region=="Genes"]
	p.gene.box<-ggplot(dcast(dt.meth, chr+strand+start+end~Gender, value.var="meth"), aes(chr, (Female-Male)*100)) + 
			geom_boxplot() + 
			labs(x="% methylation difference (female - male)") + 
			geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
			ylim(c(-100,100)) +
			theme_Publication()

	
	dt.meth=dt.meth.per.tissue[[my.tissue]][Region=="Promo_15.05"]
	p.promoter.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_density(size=1.2) + 
			labs(x="% Methylation Level", y="Density") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			#ylim(c(0,0.062)) +
			theme_Publication() +
			theme(legend.position="none")

	dt.meth=dt.meth.per.tissue[[my.tissue]][Region=="CPGi"]
	p.cpgi.beta<-ggplot(dt.meth, aes(meth*100, col=Gender)) + 
			geom_density(size=1.2) + 
			labs(x="% Methylation Level", y="") +
			scale_colour_manual(name="",values=my.col[["Gender"]]) +
			#ylim(c(0,0.062)) +
			theme_Publication() + 
			theme(legend.position="none")

# median meth level by gender and tissue
#lapply(dt.meth.per.tissue, function(i) i[Region=="Sites",median(meth),"Gender"])

#write.csv(data.table::rbindlist(dt.meth.per.tissue.csv), file="~/results/RoadMap/BS-Seq/Beta-distribution/non.weighted.mean.meth.level.by.cpg.context.all.tissues.csv")
cat("All done\n")
