#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 18/Oct/2016
# Last modified 18/Oct/2016
# Modified from bin/R/RoadMap/plot.roadmap.gender.meth.distribution.by.chr.tissue.context.R 
# Modified from bin/R/RoadMap/plot.roadmap.gender.meth.diff.by.chr.tissue.context.R

TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/RoadMap/local.R")

my.cpg.type="CG" # CG, CHH, CHG
#for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){
for(my.tissue in c("PT.SS.Tech","PT.SS.Bio")){
	my.RData=file.path("~/results/RoadMap/BS-Seq/RData/10X",paste(my.tissue,my.cpg.type,"dt.meth.region.RData",sep="."))
	if(file.exists(my.RData)){
		cat(paste(my.RData, " already there. Skipping...\n"))
	}else{
		## Load this RoadMap tissue
		dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)
		dt.query=dt.query[V6.x>=10 & V6.y>=10] # apply CpG depth for PT.SS.Tech or PT.SS.Bio

		#cpg.contexts=c("Sites","Tiling5K","Genes","Promo_10.10","CPGi","Enhancer", "CPGi-shores")
		cpg.contexts=c("Sites","Genes","Promo_15.05","CPGi")
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
							dt.overlap[,list(.N,sum(V5.y)/sum(V6.y),'cnt.C'=sum(V5.y),'cnt.T'=sum(V6.y)-sum(V5.y),"Male",i,my.tissue),by=group.keys.list]
						)
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
		dt.meth.region<-data.table::rbindlist(dt.meth.per.region)

		save(dt.meth.region, file=my.RData)
		cat(paste0("Done for ", my.tissue, "...\n"))

	}# end of my.tissue 
}

cat("All done\n")
