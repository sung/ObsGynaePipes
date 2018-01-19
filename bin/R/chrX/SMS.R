TR_PREFIX="GRCh38"
library(data.table)
source("~/Pipelines/config/Annotation.R")
dt.samples.all<-fread("/home/ssg29/results/RNA-Seq/Placentome/Meta/reconstruct.stringtie.csv", na.strings="")

dt.sms.ase<-fread("~/Pipelines/data/ASE/SMS.HC.ASE.GRCh38.txt")
dt.sms.ase.samples<-merge(dt.sms.ase, unique(dt.samples.all[,.(Library,BarCode,Condition,Sex,CRN)]))
dt.sms.ase.samples[,.(No_affected=length(unique(CRN)), Ratio=length(unique(CRN))/209),"Pos,Ref,Alt"][order(-No_affected)]

# exonic hetSNPs
dt.dummy<-lapply(dt.sms.ase.samples[,.(No_affected=length(unique(CRN)), Ratio=length(unique(CRN))/209),"Pos,Ref,Alt"][order(Pos)]$Pos, function(i) dt.exon[hgnc_symbol=="SMS" & i>=start_position & i<=end_position])
names(dt.dummy)<-dt.sms.ase.samples[,.(No_affected=length(unique(CRN)), Ratio=length(unique(CRN))/209),"Pos,Ref,Alt"][order(Pos)]$Pos
dt.sms.ase.samples[Pos %in% names(dt.dummy[sapply(dt.dummy, nrow)>0])]
dt.sms.ase.samples[,.(No_affected=length(unique(CRN)), Ratio=length(unique(CRN))/209),"Pos,Ref,Alt"][order(-No_affected)][Pos %in% names(dt.dummy[sapply(dt.dummy, nrow)>0])]
