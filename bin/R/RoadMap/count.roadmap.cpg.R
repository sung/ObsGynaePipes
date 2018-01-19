#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# Last modified 3/Mar/2016

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
dt.count.region=list() # to store the number of 1) cpg-regions and 2) cpg within the region
dt.per.tissue.count=list() # aggregation of dt.count.region by tissue 
for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){
	## Load this RoadMap tissue
	dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)
	################################
	# 1. aggregate by CPG context ##
	################################
	#cpg.contexts=names(rnb.diff.met)
	cpg.contexts=c("Sites","Tiling5K","Genes","Promo_15.05","Promo_10.03","CPGi","CPGi-shores")
	for(i in cpg.contexts){
		cat(paste0("Processing ", i, "...\n"))
		if(i=='Sites'){
			dt.count.region[[my.tissue]][[i]]=dt.query[,list(tissue=my.tissue, region=i, cpg.sum=.N, region.cnt=.N, density=1),V1]
		}else{
			subject<-with(as.data.frame(get.region(i)), data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand)) #isa 'data.frame'

			cat("\tRemoving leading 'chr'...\n")
			#subject$Chromosome<-simplify2array(strsplit(as.character(subject$Chromosome), "chr"))[2,] # chrX = > X
			subject$Chromosome<-substr(subject$Chromosome, 4, 5) #chrX=>X
			dt.subject<-as.data.table(subject) # isa data.table
			cat("\tdone...\n")

			if(grepl('Gene',i) || grepl('Promo',i) || grepl('enst',i)){
				subject.keys=c("Chromosome","Strand","Start","End") # columns from rnb.diff.met 
				object.key=c("V1","V3","V2","End")                  # columsn from dt.query (roadmap merged data)
				group.keys="V1,V3,Start,End"                        # V3: strand
			# cpgislands, tiling, genomeTiling500
			}else{
				subject.keys=c("Chromosome","Start","End")
				object.key=c("V1","V2","End")
				group.keys="V1,Strand,Start,End"
			}
			setkeyv(dt.subject,subject.keys)
			cat("\tFinding CpGs overlapping with ",i,"...\n")
			system.time(dt.overlap<-foverlaps(dt.query, dt.subject, by.x=object.key, type="any", nomatch=0L))

			cat("\tAggregating %met...\n")
			dt.cpg.region<-dt.overlap[,list(sum(V5.x)/sum(V6.x)-sum(V5.y)/sum(V6.y),num.sites=.N),by=group.keys]
			dt.count.region[[my.tissue]][[i]]=dt.cpg.region[,list(tissue=my.tissue, region=i, cpg.sum=sum(num.sites), region.cnt=.N, density=sum(num.sites)/.N),V1]

			####################
			## GRanges way    ##
			## via findOvelap ##
			## should be same ##
			####################
			if(FALSE){
				df.query<-as.data.frame(dt.query)	
				colnames(df.query)<-c("chr","start","strand","context","flag","num.C.f","depth.f","num.C.m","depth.m","end")

				gr.query=makeGRangesFromDataFrame(df.query,keep.extra.columns=TRUE)
				gr.subject=makeGRangesFromDataFrame(subject,keep.extra.columns=TRUE)

				gr.overlap<-findOverlaps(gr.query, gr.subject, ignore.strand=TRUE)
			}
		}
	}# end of cpg.contexts
	dt.per.tissue.count[[my.tissue]]<-rbindlist(dt.count.region[[my.tissue]]) # merge all cpg context of this tissue
	cat(paste0("Done for ", my.tissue, "...\n"))
}# end of my.tissue 

my.RData=paste("~/results/RoadMap/BS-Seq/RData/cpg.cnt.by.region",my.cpg.type,"RData",sep=".")
save(dt.per.tissue.count, file=my.RData)
write.csv(rbindlist(dt.per.tissue.count), file=paste("~/results/RoadMap/BS-Seq/RData/cpg.cnt.by.region",my.cpg.type,"csv",sep="."))

cat("All done\n")
