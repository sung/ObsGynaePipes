################################
## RoadMap                    ##
## Schultz et al. Nature 2015 ##
## Tissue comparision         ##
################################
my.cpg.context <- "CH" # CpG context
source("~/Pipelines/bin/R/RoadMap/local.R")

# 9 tissues ("AD" "AO" "EG" "FT" "GA" "PA" "PO" "SB" "SX") where both STL2 and STL3 are available

is.filter<-FALSE
prep.roadmap("AO",is.filter, my.cpg.context)
#for(my.tissue in all.tissues){prep.roadmap(my.tissue,is.filter)}
cat("All is done\n")
