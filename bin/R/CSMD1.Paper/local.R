#!/usr/bin/Rscript --vanilla
# used by bin/R/CSMD1.Paper/CSMD1.paper.figures.R
library(data.table)
time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p') #2015-01-07_01_30PM
TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/RoadMap/local.R") # time.stamp, load.my.tissue.dt.merged

my.cpg.type="CG" # CG, CHH, CHG
my.tissue="PT"
min.doc <- 10 # depth of coverage
my.flank=1000
my.ensg='ENSG00000183117' # CSMD1
my.first.exon_id="ENSE00002127078" # the first exon at the putative TSS

gr.csmd<-import("~/data/genome/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/Homo_sapiens.GRCh37.75.CSMD1.exon.union.sorted.ucsc.gtf") # ucsc style (e.g. chr8)
gr.csmd.ext<-promoters(range(reduce(gr.csmd)), upstream=my.flank, downstream=end(range(reduce(gr.csmd)))-start(range(reduce(gr.csmd)))+1+my.flank) # includes 1kb +TSS and 1kb +TES 
