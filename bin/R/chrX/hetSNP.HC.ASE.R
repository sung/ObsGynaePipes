#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 3/Nov/2017
# Last modified 3/Nov/2017
library(data.table)
TR_PREFIX="GRCh38"
ENS_VER=82 # 82 (for manual downlaod) 90 (for Cholestatis project)
load(paste0("~/data/Annotation/Ensembl/",TR_PREFIX,".",ENS_VER,"/dt.ensg.RData")) # load 'gr.ensg, dt.ensg'
source("~/Pipelines/bin/R/GTEx/local.R")# defines my.tissue dt.tissue, dt.exp.all.chrX
dt.samples.all<-fread("/home/ssg29/results/RNA-Seq/Placentome/Meta/reconstruct.stringtie.csv", na.strings="")

######################
## Inactivated gene ##
## GRCh38           ##
######################
dt.target<-merge(dt.exp.all.chrX[baseMean>100 & new.padj>0.4  & abs(log2FoldChange) <= log2(1.1),.N,"ensembl_gene_id"][N==20], 
				 dt.ensg[,.(`ensembl_gene_id`=as.character(ensembl_gene_id),`hgnc_symbol`=as.character(hgnc_symbol),chromosome_name,strand,start_position,end_position)])[order(hgnc_symbol)]

################################
# for Irving (+ve controls)
my.target.gene=merge(dl.exp.gtex$Placenta, dt.ensg)[order(-baseMean)][gene_biotype=="protein_coding"][1:500,hgnc_symbol]
dt.target<-dt.ensg[hgnc_symbol %in% my.target.gene]; dt.target$hgnc_symbol<-as.character(dt.target$hgnc_symbol)
################################

# SLX*tophat.HC.ASE.txt: from JD samples (GRCh38)
print(system.time(dt.ase.list<-lapply(split(dt.target, dt.target$hgnc_symbol), 
									function(i){
										my.cmd<-paste0("for i in `ls /scratch/obsgynae/POPS/results/SLX*/GATK/D*/SLX*tophat.HC.ASE.txt`; do a=`echo $i | cut -d/ -f9`; SLX=`echo $a| cut -d. -f1`; Barcode=`echo $a| cut -d. -f2`; awk -v S=$SLX -v B=$Barcode '$1==",i$chromosome_name," && $2>=",i$start_position," && $2<=",i$end_position,"{print S,B,$1,$2,$4,$5,$6,$7}' $i; done\n")
										cat(my.cmd)
										dummy<-system(my.cmd,intern=T)
										if(length(dummy)>0){
											data.table(cbind(i$hgnc_symbol, as.data.table(t(simplify2array(strsplit(dummy, " "))))))
										}
									}
								)
))

dt.ase<-rbindlist(dt.ase.list); setnames(dt.ase, c("Gene","Library","BarCode","Chr","Pos","Ref","Alt","RefCount","AltCount"))
print(system.time(dt.ase[, (c("Pos","RefCount","AltCount")) := lapply(.SD,as.numeric), .SDcols=c("Pos","RefCount","AltCount")]))
dt.ase.samples<-merge(dt.ase, unique(dt.samples.all[,.(Library,BarCode,Condition,Sex,CRN)]))

##
## dbSNP
##
system.time(dt.dbsnp<-fread("zcat ~/data/dbSNP/Homo_sapiens/Ensembl/GRCh38/all.vcf.gz",skip="#CHROM",select=c(1:5),col.names=c("CHROM","POS","ID","REF","ALT"),key=c("CHROM","POS")))

##################
# exonic hetSNPs #
# GRCh38         #
##################
No.sample<-dt.samples.all[Source!="FG",length(unique(CRN))]
dt.query=dt.ase.samples[,.(No_affected=length(unique(CRN)), Ratio=length(unique(CRN))/No.sample),"Gene,Chr,Pos,Ref,Alt"][,.(Chr,`start`=Pos,`end`=Pos)]
setkeyv(dt.exon,c("chromosome_name","start_position","end_position"))
# exonic variant only
dt.exonic.snp=foverlaps(dt.query,dt.exon,by.x=c("Chr","start","end"),type="within",nomatch=0L)[,.N,"Chr,start,end"] # exonic 
# exonic varaint + reasonble frequency (e.g. >.5)
dt.top.exonic.snp<-merge(dt.ase.samples,dt.exonic.snp,by.x=c("Chr","Pos"), by.y=c("Chr","start"))[,.(No_affected=length(unique(CRN)), Ratio=length(unique(CRN))/No.sample),"Chr,Pos,Ref,Alt"][Ratio>0.5][order(-Ratio)] # most frequent exonic hetSNPs

# this top exonic snp should be in these three CRN samples
dt.foo<-merge(dt.ase.samples[CRN %in% c(1302386,1321127,919968)],dt.top.exonic.snp,by.x=c("Chr","Pos","Ref","Alt"),by.y=c("Chr","Pos","Ref","Alt"))[,.N,"Chr,Pos,Ref,Alt"][N==3] 

# varaints of known rsID only
dt.bar<-merge(dt.foo[,.(Chr,Pos,Ref,Alt)], dt.dbsnp, by.x=c("Chr","Pos","Ref","Alt"),by.y=c("CHROM","POS","REF","ALT"))
# write in VCF
write.table(dt.bar[,.(Chr,Pos,ID,Ref,Alt,".",".",".")], sep="\t",file="~/results/chrX/CSV/positive.control.candidates.HC.ASE.GRCh38.vcf", col.names=F, row.names=F, quote=F)
# all samples with these varaints of known rsID only
dt.final<-merge(dt.ase.samples, dt.bar, by.x=c("Chr","Pos","Ref","Alt"), by.y=c("Chr","Pos","Ref","Alt"))
# write all samples and varaints in CSV
write.csv(dt.final[,.(No_affected=length(unique(CRN)), Ratio=length(unique(CRN))/No.sample),"Gene,Chr,Pos,Ref,Alt,ID"][order(-No_affected)], file="~/results/chrX/CSV/positive.control.candidates.HC.ASE.GRCh38.csv", row.names=F, quote=F)


##########
## SMS  ##
##########
dt.sms.ase<-fread("~/Pipelines/data/ASE/SMS.HC.ASE.GRCh38.txt")
dt.sms.ase.samples<-merge(dt.sms.ase, unique(dt.samples.all[,.(Library,BarCode,Condition,Sex,CRN)]))
dt.sms.ase.samples[,.(No_affected=length(unique(CRN)), Ratio=length(unique(CRN))/No.sample),"Gene,Chr,Pos,Ref,Alt"][order(-No_affected)]

# exonic hetSNPs
dt.dummy<-lapply(dt.sms.ase.samples[,.(No_affected=length(unique(CRN)), Ratio=length(unique(CRN))/No.sample),"Gene,Chr,Pos,Ref,Alt"][order(Pos)]$Pos, function(i) dt.exon[hgnc_symbol=="SMS" & i>=start_position & i<=end_position])
names(dt.dummy)<-dt.sms.ase.samples[,.(No_affected=length(unique(CRN)), Ratio=length(unique(CRN))/No.sample),"Gene,Chr,Pos,Ref,Alt"][order(Pos)]$Pos
dt.sms.ase.samples[Pos %in% names(dt.dummy[sapply(dt.dummy, nrow)>0])]
dt.sms.ase.samples[,.(No_affected=length(unique(CRN)), Ratio=length(unique(CRN))/No.sample),"Gene,Chr,Pos,Ref,Alt"][order(-No_affected)][Pos %in% names(dt.dummy[sapply(dt.dummy, nrow)>0])]

# inactive genes from the same individual found in SMS
dt.ase.samples[CRN %in% dt.sms.ase.samples[Pos %in% names(dt.dummy[sapply(dt.dummy, nrow)>0]),CRN]]
write.csv(dt.ase.samples[CRN %in% dt.sms.ase.samples[Pos %in% names(dt.dummy[sapply(dt.dummy, nrow)>0]),CRN]][order(Gene,Pos,CRN)], file="~/results/chrX/CSV/Inactive.genes.Placenta.GTEx.HC.ASE.GRCh38.csv", row.names=F, quote=F)
