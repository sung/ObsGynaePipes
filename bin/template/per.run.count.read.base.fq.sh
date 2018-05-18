#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 26/Jan/2015
# Last modified: 26/Jan/2015
# Optimised and customised to run at the Darwin HPC

source $HOME/Pipelines/config/rna_seq.config # export envrionment variables
source $HOME/Pipelines/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.hpc.bash #defines PATH 

export SLX="MY_SLX" # e.g. SLX-8080 
export PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1

mkdir_unless $PROJECT_DIR

echo "making $PROJECT_DIR/$SLX.fq.read.base.cnt.txt"

# FORMAT
# FASTQ_FILE CNT_READ
#for i in `ls $HOME/data/fastq/$SLX/$SLX.*.fq.gz`; do printf "$i "; zcat $i | echo $((`wc -l`/4)); done > $PROJECT_DIR/$SLX.fq.read.base.cnt.txt 

# FORMAT
# FASTQ_FILE CNT_READ CNT_BASE
for i in `ls $HOME/data/fastq/$SLX/$SLX.*.fq.gz`; do printf "$i "; zcat $i | awk 'BEGIN{sum=0;cnt=0}NR%4==2{cnt++;sum+=length($0);}END{print cnt,sum}'; done > $PROJECT_DIR/$SLX.fq.read.base.cnt.txt 

# FORMAT
# FASTQ_FILE UNIQ_CNT_READ CNT_READ CNT_BASE
#for i in `ls $HOME/data/fastq/$SLX/$SLX.*.fq.gz`; do UNIQ=`zcat $i | awk 'NR%4==1{print $0}' | sort -S $SORT_MEM -T $HOME/scratch/temp -u | wc -l`; printf "$i $UNIQ "; zcat $i | awk 'BEGIN{sum=0;cnt=0}NR%4==2{cnt++;sum+=length($0);}END{print cnt,sum;}'; done > $PROJECT_DIR/$SLX.fq.read.base.cnt.txt 

