#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 26/Jan/2015
# Last modified: 26/Jan/2015
# Optimised and customised to run at the Darwin HPC

#source $HOME/RNA-Seq/config/rna_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.bash #defines PATH 

export SLX="MY_SLX" # e.g. SLX-8080 
export PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1

mkdir_unless $PROJECT_DIR

echo "making $PROJECT_DIR/$SLX.fq.read.base.cnt.txt"
#for i in `ls ~/data/fastq/$SLX/$SLX.*.fq.gz`; do printf "$i "; zcat $i | awk 'BEGIN{sum=0;cnt=0}{if(NR%4==2){cnt++;sum+=length($0);}}END{print cnt,sum;}' ; done >  $PROJECT_DIR/$SLX.fq.read.base.cnt.txt

for i in `ls $HOME/data/fastq/$SLX/$SLX.*.fq.gz`; do UNIQ=`zcat $i | awk 'NR%4==1{print $0}' | sort -S 60G -T $HOME/scratch/temp -u | wc -l`; printf "$i $UNIQ "; zcat $i | awk 'BEGIN{sum=0;cnt=0}NR%4==2{cnt++;sum+=length($0);}END{print cnt,sum;}'; done > $PROJECT_DIR/$SLX.fq.read.base.cnt.txt 

#for i in `ls /home/ssg29/data/fastq/SLX-9345.PE75/SLX-9345.D*.fq.gz`; do FQ=`echo $i | cut -d '/' -f7`; TOTAL=`zcat $i | wc -l`; UNIQ=`zcat $i | awk 'NR%4==1{print $0}' | sort | uniq -c | awk '$1>1{print 0}' | wc -l`; echo "$FQ $((TOTAL/4)) $UNIQ"; done > ~/results/SLX-9345.Homo_sapiens.PE75.v1/SLX-9345.fq.read.uniq.cnt.txt
