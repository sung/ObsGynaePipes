#!/usr/bin/Rscript --vanilla
# to parse '.tracking' file from the cuffcompare
# Sung Gong <sung@bio.cc>

minFpkm=0.1
TR_PREFIX='GRCh38' # GRCh37|GRCh38
mySource="Placentome" # "FG"
top.dir=file.path("~/results/RNA-Seq",mySource) #top.dir="~/results/RNA-Seq/Placentome"
#myCohort="SGA" # "CTLM", "CTLF", "POPS", "CTL", "PET", "SGA"

source("~/Pipelines/config/Annotation.R")
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/config/cufflink.R")
source("~/Pipelines/bin/R/Placentome/local.R")

dl.frag=list()
for(myCohort in c("CTL","PET","SGA")){
	myProject=paste(mySource,myCohort,sep=".") # e.g. Placentome.CTL
	cuffRData=file.path(top.dir, paste0("RData/",myProject,".",TR_PREFIX,".",ENS_VER,".cuffcompare.RData"))

	cat(paste("loading",cuffRData,"...\n"))
	load (cuffRData) # 
	cat("A list named ", myProject, "loaded\n")

	myCuff<-get(myProject) # myCuff<-FG.AGA
	eval(parse(text=paste('rm(', myProject, ')'))) # remove the original object from memory

	dl.frag[[myCohort]]=myCuff$Transfrag[,`cohort`:=myCohort]   # isa 'data.table' (or data.frame?) (per each transcript (TCONS_))
}

myRData=file.path(top.dir, paste0("RData/",mySource,".ALL.",TR_PREFIX,".",ENS_VER,"dt.frag.RData"))
save(dl.frag, file=myRData)

venn.list<-list()
venn.list[["FPKM_0.1"]]<-lapply(dl.frag, function(i) i[fpkm.iqr>0.1 & evi.ratio>=.95 & class_code=="=",ensembl_transcript_id])
venn.list[["FPKM_1"]]<-lapply(dl.frag, function(i) i[fpkm.iqr>1 & evi.ratio>=.95 & class_code=="=",ensembl_transcript_id])
venn.list[["FPKM_10"]]<-lapply(dl.frag, function(i) i[fpkm.iqr>10 & evi.ratio>=.95 & class_code=="=",ensembl_transcript_id])
venn.list[["FPKM_100"]]<-lapply(dl.frag, function(i) i[fpkm.iqr>100 & evi.ratio>=.95 & class_code=="=",ensembl_transcript_id])
venn.list[["NO_filter"]]<-lapply(dl.frag, function(i) i[class_code=="=",ensembl_transcript_id])

library(VennDiagram)
for(my.type in names(venn.list)){
	file.name<-file.path(paste0(top.dir,"/Figures/Venn/venn.complete.match.",my.type,".tiff"))
	venn.diagram(
		x=venn.list[[my.type]],
		filename = file.name,
		height=1500,
		width=1500,
		main=my.type,
		#main.cex = 1.5,
		main.fontfamily = "sans",
		main.fontface = "bold",
		col = "black",
		#fill = c("white","white","white"),
		alpha = 0.50,
		#cat.col = c("dodgerblue", "goldenrod1", "darkorange1"),
		cex = 1.4,
		fontfamily = "sans",
		#cat.cex = 1.3,
		cat.fontfamily = "sans",
		cat.fontface = "bold",
		margin = 0.05
	)
}

dt.frag.cohort<-rbindlist(dl.frag)
dt.uniq<-dt.frag.cohort[fpkm.iqr>100 & evi.ratio>=.95 & class_code=="=",list(.N,paste(cohort,collapse=".")),ensembl_transcript_id][N==1]
#    ensembl_transcript_id N  V2
# 1:       ENST00000278175 1 CTL
# 2:       ENST00000360467 1 PET

fields <- c("ensembl_transcript_id", "ensembl_gene_id","hgnc_symbol","gene_biotype")
dt.uniq.enst<-data.table(getBM(attributes =fields, filters = fields[1], values = dt.uniq$ensembl_transcript_id, mart = myMart))

merge(dt.frag.cohort, dt.uniq)
