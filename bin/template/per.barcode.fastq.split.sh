#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 21/Jul/2014
# Last modified: 21/Jul/2014 
# Optimised and customised to run at the Darwin HPC

source $HOME/Methyl-Seq/config/methyl_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

export SLX="MY_SLX" # e.g. SLX-8080 
export PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1
export Barcode="MY_BARCODE" # e.g. A001

FASTQ_DIR=$HOME/scratch/data/fastq/$SLX
mkdir_unless $FASTQ_DIR/split
cd $FASTQ_DIR/split


FASTQ_FILES=`ls $FASTQ_DIR/$SLX.$Barcode*.fq.gz | grep -v lost`
for FastQ_file in $FASTQ_FILES	# per fastq file
do
	#NUM_CHUNK within config/methyl_seq.config
	#/home/ssg29/scratch/data/fastq/SLX-8075/SLX-8075.A008.C39MDACXX.s_1.r_1.fq.gz
	MY_PREFIX=`echo $FastQ_file | cut -d/ -f8 | cut -d. -f1,2,3,4,5`
	echo -e "fastqutils split -gz $FastQ_file $MY_PREFIX $NUM_CHUNK &" 
	time fastqutils split -gz $FastQ_file $MY_PREFIX $NUM_CHUNK &
done
wait 
