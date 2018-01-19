#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# CpG>=10x
# First created 9/Nov/2016
# Last modified 9/Nov/2016
# This scripts is to make BedGraph for % methylation by sex

source("~/Pipelines/bin/R/RoadMap/local.R")
options(scipen=999)

################## 
## PT: Placenta ##
################## 
dt.query=load.my.tissue.dt.merged(my.tissue="PT",my.cpg.type="CG")

#######################################
## Ensembl Format (e.g. X, not chrX) ##
#######################################
# Female
my.file="~/results/RoadMap/POPS/BedGraph/WGBS/PT.WGBS.Female.bedgraph"
write.table(dt.query[order(V1,V2)][,.(V1,V2-1,V2,round(V5.x/V6.x*100,1))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=my.file)
# Female-fwd
my.file="~/results/RoadMap/POPS/BedGraph/WGBS/PT.WGBS.Female.fwd.bedgraph"
write.table(dt.query[order(V1,V2)][V3=="+",.(V1,V2-1,V2,round(V5.x/V6.x*100,1))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=my.file)
# Female-rev
my.file="~/results/RoadMap/POPS/BedGraph/WGBS/PT.WGBS.Female.rev.bedgraph"
write.table(dt.query[order(V1,V2)][V3=="-",.(V1,V2-1,V2,round(V5.x/V6.x*100,1))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=my.file)

# Male
my.file="~/results/RoadMap/POPS/BedGraph/WGBS/PT.WGBS.Male.bedgraph"
write.table(dt.query[order(V1,V2)][,.(V1,V2-1,V2,round(V5.y/V6.y*100,1))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=my.file)
# Male-fwd
my.file="~/results/RoadMap/POPS/BedGraph/WGBS/PT.WGBS.Male.fwd.bedgraph"
write.table(dt.query[order(V1,V2)][V3=="+",.(V1,V2-1,V2,round(V5.y/V6.y*100,1))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=my.file)
# Male-rev
my.file="~/results/RoadMap/POPS/BedGraph/WGBS/PT.WGBS.Male.rev.bedgraph"
write.table(dt.query[order(V1,V2)][V3=="-",.(V1,V2-1,V2,round(V5.y/V6.y*100,1))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=my.file)

####################################
## UCSC Format (e.g. chrX, not X) ##
####################################
# Female
my.file=gzfile("~/results/RoadMap/POPS/BedGraph/WGBS/PT.WGoxBS.Female.hg19.bedgraph.gz")
write.table(dt.query[order(V1,V2)][,.(paste0("chr",V1),V2-1,V2,round(V5.x/V6.x*100,1))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=my.file)
# Female-fwd
my.file=gzfile("~/results/RoadMap/POPS/BedGraph/WGBS/PT.WGoxBS.Female.fwd.hg19.bedgraph.gz")
write.table(dt.query[order(V1,V2)][V3=="+",.(paste0("chr",V1),V2-1,V2,round(V5.x/V6.x*100,1))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=my.file)
# Female-rev
my.file=gzfile("~/results/RoadMap/POPS/BedGraph/WGBS/PT.WGoxBS.Female.rev.hg19.bedgraph.gz")
write.table(dt.query[order(V1,V2)][V3=="-",.(paste0("chr",V1),V2-1,V2,round(V5.x/V6.x*100,1))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=my.file)

# Male
my.file=gzfile("~/results/RoadMap/POPS/BedGraph/WGBS/PT.WGoxBS.Male.hg19.bedgraph.gz")
write.table(dt.query[order(V1,V2)][,.(paste0("chr",V1),V2-1,V2,round(V5.y/V6.y*100,1))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=my.file)
# Male-fwd
my.file=gzfile("~/results/RoadMap/POPS/BedGraph/WGBS/PT.WGoxBS.Male.fwd.hg19.bedgraph.gz")
write.table(dt.query[order(V1,V2)][V3=="+",.(paste0("chr",V1),V2-1,V2,round(V5.y/V6.y*100,1))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=my.file)
# Male-rev
my.file=gzfile("~/results/RoadMap/POPS/BedGraph/WGBS/PT.WGoxBS.Male.rev.hg19.bedgraph.gz")
write.table(dt.query[order(V1,V2)][V3=="-",.(paste0("chr",V1),V2-1,V2,round(V5.y/V6.y*100,1))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=my.file)
