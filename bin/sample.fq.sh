#!/bin/bash 
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 4/May/2017
# Last modified: 21/Jul/2017
source $HOME/lib/sung.sh #defines 'mkdir_unless'    
alias seqtk='~/Install/seqtk/seqtk'

SLX="SLX-11369"

# sample 10K reads
# -s100: initial seeds to preserve the same order between r_1 and r_2

#time seqtk sample -s100 ~/data/fastq/SLX-11369/SLX-11369.NoIndex.HJ75CBBXX.s_2.r_1.fq.gz 10000 | awk 'BEGIN{OFS="\n"}NR%4==1{n=$0}NR%4==2{s=substr($0,1,50)}NR%4==0{q=substr($0,1,50); print n,s,"+",q}' | gzip -9 > ~/data/fastq/SLX-11369-PE50/SLX-11369-PE50.NoIndex.HJ75CBBXX.s_2.r_1.fq.gz
#time seqtk sample -s100 ~/data/fastq/SLX-11369/SLX-11369.NoIndex.HJ75CBBXX.s_2.r_2.fq.gz 10000 | awk 'BEGIN{OFS="\n"}NR%4==1{n=$0}NR%4==2{s=substr($0,1,50)}NR%4==0{q=substr($0,1,50); print n,s,"+",q}' | gzip -9 > ~/data/fastq/SLX-11369-PE50/SLX-11369-PE50.NoIndex.HJ75CBBXX.s_2.r_2.fq.gz

#echo -e "20K"
#mkdir_unless ~/data/fastq/$SLX-20K
#time seqtk sample -s100 ~/data/fastq/$SLX/$SLX.NoIndex.HJ75CBBXX.s_2.r_1.fq.gz 20000 |gzip -9 > ~/data/fastq/$SLX-20K/$SLX-20K.NoIndex.HJ75CBBXX.s_2.r_1.fq.gz
#time seqtk sample -s100 ~/data/fastq/$SLX/$SLX.NoIndex.HJ75CBBXX.s_2.r_2.fq.gz 20000 |gzip -9 > ~/data/fastq/$SLX-20K/$SLX-20K.NoIndex.HJ75CBBXX.s_2.r_2.fq.gz

#echo -e "40K"
#mkdir_unless ~/data/fastq/$SLX-40K
#time seqtk sample -s100 ~/data/fastq/$SLX/$SLX.NoIndex.HJ75CBBXX.s_2.r_1.fq.gz 40000 |gzip -9 > ~/data/fastq/$SLX-40K/$SLX-40K.NoIndex.HJ75CBBXX.s_2.r_1.fq.gz
#time seqtk sample -s100 ~/data/fastq/$SLX/$SLX.NoIndex.HJ75CBBXX.s_2.r_2.fq.gz 40000 |gzip -9 > ~/data/fastq/$SLX-40K/$SLX-40K.NoIndex.HJ75CBBXX.s_2.r_2.fq.gz
#
#echo -e "60K"
#mkdir_unless ~/data/fastq/$SLX-60K
#time seqtk sample -s100 ~/data/fastq/$SLX/$SLX.NoIndex.HJ75CBBXX.s_2.r_1.fq.gz 60000 |gzip -9 > ~/data/fastq/$SLX-60K/$SLX-60K.NoIndex.HJ75CBBXX.s_2.r_1.fq.gz
#time seqtk sample -s100 ~/data/fastq/$SLX/$SLX.NoIndex.HJ75CBBXX.s_2.r_2.fq.gz 60000 |gzip -9 > ~/data/fastq/$SLX-60K/$SLX-60K.NoIndex.HJ75CBBXX.s_2.r_2.fq.gz
#
#echo -e "80K"
#mkdir_unless ~/data/fastq/$SLX-80K
#time seqtk sample -s100 ~/data/fastq/$SLX/$SLX.NoIndex.HJ75CBBXX.s_2.r_1.fq.gz 80000 |gzip -9 > ~/data/fastq/$SLX-80K/$SLX-80K.NoIndex.HJ75CBBXX.s_2.r_1.fq.gz
#time seqtk sample -s100 ~/data/fastq/$SLX/$SLX.NoIndex.HJ75CBBXX.s_2.r_2.fq.gz 80000 |gzip -9 > ~/data/fastq/$SLX-80K/$SLX-80K.NoIndex.HJ75CBBXX.s_2.r_2.fq.gz
#
echo -e "100K"
mkdir_unless ~/data/fastq/$SLX-100K
time seqtk sample -s100 ~/data/fastq/$SLX/$SLX.NoIndex.HJ75CBBXX.s_2.r_1.fq.gz 100000 |gzip -9 > ~/data/fastq/$SLX-100K/$SLX-100K.NoIndex.HJ75CBBXX.s_2.r_1.fq.gz
time seqtk sample -s100 ~/data/fastq/$SLX/$SLX.NoIndex.HJ75CBBXX.s_2.r_2.fq.gz 100000 |gzip -9 > ~/data/fastq/$SLX-100K/$SLX-100K.NoIndex.HJ75CBBXX.s_2.r_2.fq.gz
