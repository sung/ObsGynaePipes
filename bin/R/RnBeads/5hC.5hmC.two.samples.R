#!/usr/bin/Rscript --vanilla
cat("Loading configuation...\n")
source("~/Pipelines/config/RnBeads.sbs.R")
source("~/Pipelines/lib/methylkit.R") # defines faster.fisher
cat("Loading configuration done...\n")
rnb.dm.rdata<-file.path(rdata.dir,paste0("rnb.meth.diff.",my.region,".RData")) # this will store 'rnb.diff.met' which is a list of data.frame for sites, promoters...
    cat(paste0("Loading rnb.dm.rdata: ",rnb.dm.rdata,"...\n"))
    load(rnb.dm.rdata) # loading "rnb.diff.met"
    cat("rnb.diff.met loaded\n")
my.project
i<-"sites"
5hmC.cpg<-rnb.diff.met[[i]]
5hmC.cpg=rnb.diff.met[[i]]
head(rnb.diff.met[[i]])
5hmC.sites<-rnb.diff.met[[i]]
5hmC.cpg<-rnb.diff.met[["sites"]]
5hmC.cpg
5hmc<-rnb.diff.met[[i]]
dummy<-rnb.diff.met[[i]]
5hmC.cpg<-dummy
5hmC.cpg <- dummy
dim(dummy)
5hmC.cpg=dummy
dummy->5hmC.cpg
5hmC.cpg<-data.frame()
5hmCcpg<-data.frame()
5hm<-data.frame()
cpg.5hmc<-dummy
cat("Loading configuation...\n")
source("~/Pipelines/config/RnBeads.sbs.R")
source("~/Pipelines/lib/methylkit.R") # defines faster.fisher
cat("Loading configuration done...\n")
rnb.dm.rdata<-file.path(rdata.dir,paste0("rnb.meth.diff.",my.region,".RData")) # this will store 'rnb.diff.met' which is a list of data.frame for sites, promoters...
rnb.dm.rdata
    cat(paste0("Loading rnb.dm.rdata: ",rnb.dm.rdata,"...\n"))
    load(rnb.dm.rdata) # loading "rnb.diff.met"
    cat("rnb.diff.met loaded\n")
i
cpg.5mc<-rnb.diff.met[[i]]
dim(cpg.5hmc)
dim(cpg.5mc)
head(cpg.5hmc)
select<-cpg.5hmc$Chromosome=="chr10" & cpg.5hmc$Start>=126135592 & cpg.5hmc$End<=126140753
select.5hmc<-select
select.5mc<-cpg.5mc$Chromosome=="chr10" & cpg.5mc$Start>=126135592 & cpg.5mc$End<=126140753
dim(cpg.5hmc[select.5hmc,])
dim(cpg.5mc[select.5mc,])
nkx.5hmc<-cpg.5hmc[select.5hmc,]
nkx.5mc<-cpg.5mc[select.5mc,]
dim(nkx.5hmc)
dim(nkx.5mc)
head(nkx.5hmc)
head(nkx.5mc)
nkx<-list()
nkx[["5mC"]]<-nkx.5mc
nkx[["5hmC"]]<-nkx.5hmc
lapply(nkx, length)
lapply(nkx, nrow)
names(nkx)
join(nkx[["5hmC"]],nkx[["5hmC"]])
merge(nkx[["5hmC"]],nkx[["5hmC"]])
head(merge(nkx[["5hmC"]],nkx[["5hmC"]]))
head( merge(nkx[["5hmC"]],nkx[["5hmC"]], by=c("Chromosome","Start","End")) )
head( merge(nkx[["5hmC"]],nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand")) )
head( merge(nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand")) )
merge(nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"))
merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"))
nrow(merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand")))
merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"))
merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"), all=T)
nrow(merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"), all=T))
nrow(merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"), all.x=T))
merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"), all.x=T)
merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"), all.x=T)
merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"), all.y=T)
merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"), all=F)
merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"), all=F)
merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"), all=FT
merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"), all=T)
library(VennDiagram)
venn.diagram(nkx, filename="nkx.5hmC.5mC")
nkx
names(nkx)
head(nkx[["5hmC"]])
lapply(nkx, function(i) subset(i, select(c("Chromosome","Start","End","Strand","combinedRank"))))
lapply(nkx, function(i) subset(i, select=c("Chromosome","Start","End","Strand","combinedRank")))
nkx<-lapply(nkx, function(i) subset(i, select=c("Chromosome","Start","End","Strand","combinedRank")))
nkx
lapply(nkx,nrow)
venn.diagram(nkx, filename="nkx.5mc.5hmc.count.tiff")
nrow(merge(nkx[["5mC"]][,c("Chromosome","Start","End","Strand","combinedRank")],nkx[["5hmC"]][,c("Chromosome","Start","End","Strand","combinedRank")], by=c("Chromosome","Start","End","Strand"), all=T))
nrow(merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand"), all=T))
nrow(merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand")))
nrow(merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand")))
lapply(nkx,head)
lapply(nkx, function(i) i<-i[,c("Start")] )
venn.diagram(lapply(nkx, function(i) i<-i[,c("Start")] ), filename="nkx.5mc.5hmc.count.tiff")
merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand"), all=T)
lapply(nkx, function(i) i<-i[,paste(i$Chromosome,i$Start,i$End, i$Start, sep=".")] )
lapply(nkx, function(i) i<-paste(i$Chromosome,i$Start,i$End, i$Start, sep=".") )
lapply(nkx, function(i) i<-paste(i$Chromosome,i$Start,i$End, i$Str, sep=".") )
venn.diagram(lapply(nkx, function(i) i<-paste(i$Chromosome,i$Start,i$End, i$Str, sep=".") ), filename="nkx.5mc.5hmc.count.tiff")
merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand"), all=T)
merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand"), all=T, all.x=T)
merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand"), all=T, all.x=F)
merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand"), all.x=F)
merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand"), all.x=T)
merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand"), all.x=T)
merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand"), all.y=T)
nrow(merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand"), all.y=T))
nrow(merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand"), all.y=T))-nrow(merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand")))
nrow(nkx[["5hmC"]]) - nrow(merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand")))
nrow(nkx[["5mC"]]) - nrow(merge(nkx[["5mC"]], nkx[["5hmC"]], by=c("Chromosome","Start","End","Strand")))
my.project
preprocessed.rnb.set.dir
 source("~/Pipelines/bin/R/RnBeads/set.user.annotation.sbs.R") 
rnb.set<-load.rnb.set(path="~/results/RnBeads/SGA.AGA/CpG/all/reports-2015-01-21_04_30_PM/rnbSet_unnormalized")
i<-"sites"
rnb.meth[[i]] <-meth(rnb.set, type=i, row.names=TRUE)
rnb.set<-load.rnb.set(path="~/results/RnBeads/SGA.AGA/CpG/all/reports-2015-01-21_04_30_PM/rnbSet_unnormalized")
i
meth.5mC.cpg<-rnb.meth[[i]]
str(rnb.set)
meth.5mC.cpg<-meth(rnb.set, type=i, row.names=TRUE) 
head(meth.5mC.cpg)
anno.5hmC.cpg <-annotation(rnb.set, type=i, add.names=TRUE)
raw.cpg<-data.frame(anno.5hmC.cpg, meth.5mC.cpg)
dim(raw.cpg)
head(raw.cpg)
select<-raw.cpg$Chromosome=="chr10" & raw.cpg$Start>=126135592 & raw.cpg$End<=126140753
nkx.raw.cpg<-raw.cpg[select,]
dim(nkx.raw.cpg)
head(nkx.raw.cpg)
nkx[["CpG"]]<-nkx.raw.cpg
names(nkx)
nkx[["5hmC"]]
nkx[["CpG"]]
head(nkx[["CpG"]])
merge(nkx[["5hmC"]], nkx[["CpG"]], by=c("Chromosome","Start","End","Strand"))
merge(nkx[["5hmC"]], nkx[["CpG"]], by=c("Chromosome","Start","End","Strand"), all.x=T)
dummy<-merge(nkx[["5hmC"]], nkx[["CpG"]], by=c("Chromosome","Start","End","Strand"), all.x=T)
head(dummy)
is.na(dummy$AGA2F)
is.na(dummy$SGA2F)
table(is.na(dummy$SGA2F))
rownames(dummy)
colnames(dummy)
colnames(dummy)[,ncol(dummy)-8:ncol(dummy)]
colnames(dummy)[ncol(dummy)-8:ncol(dummy)]
colnames(dummy)[1:ncol(dummy)]
colnames(dummy)[12:ncol(dummy)]
lapply(colnames(dummy)[12:ncol(dummy)], function(i) table(is.na(dummy$i))
)
lapply(colnames(dummy)[12:ncol(dummy)], function(i) table(is.na(dummy[,c(i)]))
)
table(is.na(dummy$SGA2F))
table(is.na(dummy$AGA2F))
table(is.na(dummy$AGA3F))
table(is.na(dummy$SGA3F))
colnames(dummy)
table(is.na(dummy$AGA8M))
table(is.na(dummy$SGA8M))
dummy
head(dummy)
head(raw.cpg)
history
history()
ls()
meth.5mC.cpg
samples(rnb.set)
rnb.sample.groups(rnb.set)
samples(rnb.set)
rnb.set.2F<-remove.samples(rnb.set, c("AGA3F","SGA3F","AGA5M","SGA5M","AGA8M","SGA8M"))
samples(rnb.set.2F)
savehistory(file="bin/R/RnBeads/5hC.5hmC.two.samples.R")
