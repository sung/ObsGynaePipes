#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 15/Jul/2014
# Last modified: 15/Jul/2014 
# Optimised and customised to run at the Darwin HPC

source $HOME/Methyl-Seq/config/methyl_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1
Barcode="MY_BARCODE" # e.g. A001

	# This BAM file should be indexed to split by reference (chromosome)
	# SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.bam
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.bam 

	if [ ! -s ${My_bam%.bam}.sorted.bam ];then
		if [ -s $My_bam ];then
			echo -e "samtools sort $My_bam ${My_bam%.bam}.sorted"
			time samtools sort -m 5000000000 $My_bam ${My_bam%.bam}.sorted
		fi
	fi
	if [ ! -s ${My_bam%.bam}.sorted.bam.bai ];then
		if [ -s ${My_bam%.bam}.sorted.bam ];then
			echo -e "samtools index ${My_bam%.bam}.sorted.bam "
			time samtools index ${My_bam%.bam}.sorted.bam 
		else
			echo -e "\e[031m${My_bam%.bam}.sorted.bam not found\e[0m\n"
			exit
		fi
	fi

	#####################################
	# 6. split bam file by chromosome  ##
	# this is to feed 'bismark_methylation_extractor' by chromosome to speed up
	# input: position-based sorted bam
	# output: per chr bam file
	#####################################
	if [ $RUN_SPLIT_BAM_BY_CHR -eq 1 ]; then
		for Chr in ${UCSC_CHR[*]}
		do
			echo -e "samtools view -b ${My_bam%.bam}.sorted.bam $Chr > ${My_bam%.bam}.$Chr.bam &"
			time samtools view -b ${My_bam%.bam}.sorted.bam $Chr > ${My_bam%.bam}.$Chr.bam &
		done
		wait
	fi
	echo -e "\e[32m$SLX.$Barcode done for splitting bam by chromosome\e[0m" 
