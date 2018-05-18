#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# circRNA prediced by CIRI2
# https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbx014/3058729

#########################
# FG CIRI2 result files #
#########################
my.cmd<-"for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-916*.Homo_sapiens.SE125.v2/CIRI2/D*/SLX-916*.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done"
fg.meta<-system(my.cmd,intern=T)

########################
# JD CIRI2 result file #
########################
my.cmd<-c(
          "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-9792.Homo_sapiens.v1/CIRI2/D*/SLX-9792.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done",
        "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10281.Homo_sapiens.v1/CIRI2/D*/SLX-10281.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done",
        "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10283.Homo_sapiens.v1/CIRI2/D*/SLX-10283.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done",
        "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10284.Homo_sapiens.v1/CIRI2/D*/SLX-10284.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done",
        "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10285.Homo_sapiens.v1/CIRI2/D*/SLX-10285.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done",
        "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10287.Homo_sapiens.v1/CIRI2/D*/SLX-10287.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done",
        "for i in `ls ~/rcs/rcs-ssg29-obsgynae/POPS/results/SLX-10402.Homo_sapiens.v1/CIRI2/D*/SLX-10402.D*.merged.bwa.sam.CIRI2.txt`; do SLX=`echo $i | cut -d'/' -f11 | cut -d'.' -f1`; BARCODE=`echo $i | cut -d'/' -f11 | cut -d'.' -f2`; printf \"$SLX $BARCODE $i\n\"; done"
        )
jd.meta<-unlist(
            sapply(my.cmd, system, intern=T, USE.NAMES=F)
            )

dt.pops<-data.table(t(simplify2array(strsplit(c(fg.meta,jd.meta), " ")))) # n=324 (114+210)

# read all the CIRI2 result files
dl.pops<-apply(dt.pops, 1, function(i){fread(i[3])[,-"junction_reads_ID"][,`:=`(SLX=i[1],Barcode=i[2])]})
dt.pops.ciri<-rbindlist(dl.pops)
#write.csv(dt.pops.ciri, file=gzfile("~/results/RNA-Seq/Placentome/CSV/circRNA.CIRI2.txt.gz"), row.names=F)
dt.pops.ciri[,circRNA_end-circRNA_start+1,circRNA_ID]

#############
# No Filter #
#############
dt.pops.ciri[,.(`sum`=sum(`#junction_reads`),`avg`=sum(`#junction_reads`)/nrow(dt.pops),`freq`=length(unique(paste(SLX,Barcode,sep=".")))/nrow(dt.pops),`cnt`=length(unique(paste(SLX,Barcode,sep=".")))),circRNA_ID]

############################
# Filter (at least MIN>FREQ of samples should have at least MIN.JNC.CNT)
############################
i<-0; dl.dummy<-list()
for(MIN.JNC.CNT in seq(1,5)){
    for(MIN.FREQ in seq(.1,1,by=0.1)){
        i<-i+1
        cnt<-dt.pops.ciri[`#junction_reads`>=MIN.JNC.CNT,.(`freq`=length(unique(paste(SLX,Barcode,sep=".")))/nrow(dt.pops)),circRNA_ID][freq>=MIN.FREQ,.N]
        dl.dummy[[i]]<-data.table(`min junction count`=MIN.JNC.CNT,`min freq`=MIN.FREQ,`Number of circRNA`=cnt)
    }
}
rbindlist(dl.dummy)
library(ggplot2)
ggplot(rbindlist(dl.dummy), aes(`min freq`,`Number of circRNA`,group=`min junction count`)) + geom_point(aes(color=factor(`min junction count`)),size=5) + geom_line()

MIN.JNC.CNT=1
MIN.FREQ=1/10
my.circRNA<-dt.pops.ciri[`#junction_reads`>=MIN.JNC.CNT,.(`sum`=sum(`#junction_reads`),`avg`=sum(`#junction_reads`)/nrow(dt.pops),`freq`=length(unique(paste(SLX,Barcode,sep=".")))/nrow(dt.pops),`cnt`=length(unique(paste(SLX,Barcode,sep=".")))),circRNA_ID][freq>=MIN.FREQ,circRNA_ID]

dt.foo<-dcast.data.table(dt.pops.ciri[circRNA_ID %in% my.circRNA],
                 circRNA_ID+strand~SLX+Barcode, 
                 value.var="#junction_reads",fill=0,sep=":") # long (row-based) to wide (column-based)

dt.foo[,1:5]
dim(dt.foo)
