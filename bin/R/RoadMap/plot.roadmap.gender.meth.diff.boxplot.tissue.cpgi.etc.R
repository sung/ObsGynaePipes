#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 24/Mar/2016
# Last modified 23/Aug/2016
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
# 8 tissues ("AD","AO","EG","FT","GA","PO","SB","SX"): (7/March/2016)
# PA: not available yet
# PT: for our Placenta
dt.meth.chrX.per.tissue<-list()
#for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){
for(my.tissue in c("PT","PT.SS.Tech","PT.SS.Bio")){
	## Load this RoadMap tissue
	dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)
	cat("Filtering out autosomes...\n")
	dt.query=dt.query[V1=="X"] #chrX only
	################################
	# 1. aggregate by CPG context ##
	################################
	cpg.contexts=c("Sites","Tiling5K","Genes","Promo_10.10","CPGi","CPGi-shores","Enhancer")
	#cpg.contexts=c("miR_Promo", "miR_Promo_01.01", "miR_Genes", "miR_Genes_01.01")
	dt.meth.chrX.per.region=list() # to save chrX only (key: cpg context)
	for(i in cpg.contexts){
		cat(paste0("Processing ", i, "...\n"))
		if(i=='Sites'){
			group.keys.list="V1,V3,V2,End"                   # same as above 
			dt.meth<-dt.query[,list(diff.met=V5.x/V6.x-V5.y/V6.y,i,my.tissue),by=group.keys.list]
		}else{
			subject<-with(as.data.frame(get.region(i)), data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand)) #isa 'data.frame'
			subject<-subject[subject$Chromosome=="chrX",] #chrX only
			cat("\tRemoving leading 'chr'...\n")
			#subject$Chromosome<-simplify2array(strsplit(as.character(subject$Chromosome), "chr"))[2,] # chrX = > X
			subject$Chromosome<-substr(subject$Chromosome, 4, 5) #chrX=>X, chr13=>13
			cat("\tconvert to data.table...\n")
			dt.subject<-as.data.table(subject) # isa data.table
			rm(subject)

			if(grepl('Gene',i) || grepl('Promo',i) || grepl('enst',i)){
				subject.keys=c("Chromosome","Strand","Start","End") # regions of interests 
				query.key=c("V1","V3","V2","End")                   # columns from dt.query (roadmap merged data)
				group.keys=c("V1","V3","Start","End")               # V3: strand
				group.keys.list="V1,V3,Start,End"                   # same as above 
			}else{
				subject.keys=c("Chromosome","Start","End")
				query.key=c("V1","V2","End")
				group.keys=c("V1","Strand","Start","End")
				group.keys.list="V1,Strand,Start,End"
			}
			setkeyv(dt.subject,subject.keys)

			cat("\tFinding CpGs overlapping with ",i,"...\n")
			system.time(dt.overlap<-foverlaps(dt.query, dt.subject, by.x=query.key, type="any", nomatch=0L))

			cat("\tAggregating %met...\n")
			dt.meth<-dt.overlap[,list(diff.met=sum(V5.x)/sum(V6.x)-sum(V5.y)/sum(V6.y), i, my.tissue),by=group.keys.list]
		}# end if i=='Sites'
		setnames(dt.meth, c("chr","strand","start","end","diff.met","Region","Tissue"))
		dt.meth.chrX.per.region[[i]]=dt.meth
	}# end of cpg.contexts
	dt.meth.chrX.per.tissue[[my.tissue]]<-rbindlist(dt.meth.chrX.per.region)
	cat(paste0("Done for ", my.tissue, "...\n"))
}# end of my.tissue 

cat("Merging all samples...\n")
dt.meth.all<-data.table::rbindlist(dt.meth.chrX.per.tissue)

#file.name<-file.path("~/results/RoadMap/BS-Seq/MethDiff",paste("all.tissue.chrX",my.cpg.type,time.stamp,"pdf",sep="."))
#pdf(file=file.name, width=11.7, height=8.3, title=paste0("methylation difference by gender in various ",my.cpg.type," context"))

#tiff(filename="~/results/chrX/meth.diff.by.tissue.PROmiR.tiff", width=11.7, height=8.27, units="in", res=300, compression='lzw')
#tiff(filename="~/results/chrX/meth.diff.PT.PROmiR.tiff", width=11.7, height=8.27, units="in", res=300, compression='lzw')
tiff(filename="~/results/chrX/meth.diff.PT.tiff", width=11.7, height=8.27, units="in", res=300, compression='lzw')

if(my.cpg.type=="CG"){
	my.lim=c(-50,80)
	my.y.tick=seq(from=my.lim[1],to=my.lim[2], by=20)  # Ticks every .2
}else{
	my.lim=c(-1.5,1)
	my.y.tick=seq(from=my.lim[1],to=my.lim[2], by=0.5)  # Ticks every .005
}
# Boxplot by tissue (chrX only) 
a3<-ggplot(dt.meth.all, aes(Tissue,diff.met*100)) + 
	geom_boxplot(aes(fill=Region), outlier.shape=NA,alpha=0.7) + 
	labs(y="% methylation difference\n(Female - Male)", x="Tissues") + 
	geom_hline(yintercept=0) + 
	ggtitle(paste0(my.cpg.type," Methylation Difference in ChrX")) +
	scale_fill_manual(values=cbPalette, name=paste0(my.cpg.type,"\nContexts")) +
	#scale_fill_Publication(name="CpG\nContext") +
	coord_cartesian(ylim=my.lim) + # zoom-in this range, rather than limiting (http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/)
	scale_y_continuous(breaks=my.y.tick) 
print(a3+theme_Publication())

dev.off()
cat("All done\n")
