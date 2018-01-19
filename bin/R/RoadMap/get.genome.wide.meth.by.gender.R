#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 14/Jul/2016
# Last modified 14/Jul/2016

#stop("use bin/R/RoadMap/get.meth.level.roadmap.gender.by.chr.tissue.context.R")

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
# 9 tissues ("PT", "AD","AO","EG","FT","GA","PO","SB","SX"): (7/March/2016)
# PA: not available yet
# PT: for our Placenta
# use following RData file instead:
# ~/results/RoadMap/BS-Seq/RData/CG.dt.meth.per.tissue.per.region.RData
my.RData=file.path("~/results/RoadMap/BS-Seq/RData",paste(my.cpg.type,"dt.genome.wide.meth.RData",sep="."))
if(file.exists(my.RData)){
	load(my.RData)
	cat("dt.genome.wide.meth loaded\n")
}else{
	dt.genome.wide.meth=list() # to save %met by cpg contexts
	for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){
		## Load this RoadMap tissue
		dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)

		# genome-wide meth level of this tissue by gender
		dt.genome.wide.meth[[my.tissue]]=rbind(
			dt.query[,list(tissue=my.tissue,num=.N,f.met=round(sum(V5.x)/sum(V6.x)*100,2), m.met=round(sum(V5.y)/sum(V6.y)*100,2), met=round(sum(V5.x+V5.y)/sum(V6.x+V6.y)*100,2), f.c=sum(V5.x), f.t=sum(V6.x)-sum(V5.x), m.c=sum(V5.y), m.t=sum(V6.y)-sum(V5.y)),"V1"],
			dt.query[,list(V1="genome",tissue=my.tissue,num=.N,f.met=round(sum(V5.x)/sum(V6.x)*100,2), m.met=round(sum(V5.y)/sum(V6.y)*100,2), met=round(sum(V5.x+V5.y)/sum(V6.x+V6.y)*100,2), f.c=sum(V5.x), f.t=sum(V6.x)-sum(V5.x), m.c=sum(V5.y), m.t=sum(V6.y)-sum(V5.y))]
			)

		cat(paste0("Done for ", my.tissue, "...\n"))
	}
	cat("saving dt.genome.wide.meth...\n")
	save(dt.genome.wide.meth, file=my.RData)
}

dt.genome.meth=rbindlist(dt.genome.wide.meth)

cat("All done\n")
