################################
## RoadMap                    ##
## Schultz et al. Nature 2015 ##
## Tissue comparision         ##
################################
my.cpg.context <- "CG" # either CG, CHH, CHG or all
source("~/Pipelines/bin/R/RoadMap/local.R")

#all.tissues<-c("AD","AO","EG","FT","GA","PO","SB","SX") # where both STL2 and STL3 are available
#all.tissues<-c("PT",avail.tissues[!avail.tissues %in% "PA"])
all.tissues<-avail.tissues[!avail.tissues %in% "PA"]
# PA: not available yet

#####################
## Process RoadMap ##
#####################
#my.tissue='EG'
#dt.merge.roadmap(my.tissue, my.cpg.context)
for(my.tissue in all.tissues){dt.merge.roadmap(my.tissue,my.cpg.context)}

################
## Process PT ##
################
# process.PT.bed("CHH")

cat("All is done...\n")
