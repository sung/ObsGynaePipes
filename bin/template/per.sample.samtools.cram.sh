#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 5/Mar/2018
# Last modified: 5/Mar/2018
# Optimised and customised to run at the Darwin HPC

source $HOME/Pipelines/config/rna_seq.config # export envrionment variables
source $HOME/Pipelines/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.hpc.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

mkdir_unless $PROJECT_DIR	
mkdir_unless $PROJECT_DIR/scratch
mkdir_unless $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Cell.$Barcode\e[0m\n"

###################################
## FASTQ FROM CRAM USING SAMTOOLS #
###################################
My_cram=/home/ssg29/rds/rds-ssg29-obsgynae/POPS/results/WTSI/Xten/$Barcode.cram

#############
## samtools #
#############

FASTQ_R1=$HOME/data/fastq/$SLX/$SLX.$Barcode.XXX.s_x.r_1.fq.gz
FASTQ_R2=$HOME/data/fastq/$SLX/$SLX.$Barcode.XXX.s_x.r_2.fq.gz
if [ -s $My_cram ];then
    echo -e "samtools fastq -F 0x200 $My_cram -1 $FASTQ_R1 -2 $FASTQ_R2"
    time samtools fastq -F 0x200 $My_cram -1 $FASTQ_R1 -2 $FASTQ_R2
else
    echo -e "\e[031m$My_cram not found\e[0m\n"
    exit
fi

echo -e "\e[32m$Barcode done for per.sample.samtools.cram\e[0m"
