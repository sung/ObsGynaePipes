
source("~/Pipelines/bin/R/CSMD1.Paper/local.R")
my.RData=file.path("~/results/RoadMap/BS-Seq/RData",paste(my.tissue,my.cpg.type,"dt.meth.region.RData",sep="."))
cat("loading dt.meth.region.RData...\n")
load(my.RData)
cat("dt.meth.region loaded\n")

#SAT1 (spermine acetyltransferase 1): ENSG00000130066
#SMS (spermine synthase): ENSG0000010217
my.hgnc=c("SAT1","SMS")
my.hgnc=c("XIST")
fields <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description")
df.mart<-getBM(attributes = fields, filters = "hgnc_symbol", values = my.hgnc, mart=myMart)

# GRanges
gr.target<-gr.ensg[gr.ensg$gene_id %in% df.mart$ensembl_gene_id,]

# GRangeList
gr.target.list<-split(gr.target, gr.target$gene_id)

# % Methylation within 
# Genes,Promo_15.05,CPGi
lapply(gr.target.list, function(i) dt.meth.region[Region %in% c("Genes","Promo_15.05","CPGi") & chr==substr(as.character(seqnames(i)),4,5) & start(i)>=start & start(i)<=end] )
lapply(gr.target.list[1], function(i) {my.gene=df.mart[df.mart$ensembl_gene_id==names(i),]$hgnc_symbol; dt.meth.region[Region %in% c("Genes","Promo_15.05","CPGi") & chr==substr(as.character(seqnames(i)),4,5) & start(i)>=start & start(i)<=end,.(`gene_id`=my.gene, Region, Gender, `meth`=round(meth*100,2), `No_CpG`=num.sites),]} )
