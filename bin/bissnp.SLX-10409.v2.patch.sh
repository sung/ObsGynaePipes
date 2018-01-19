#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# Called by bin/template/bissnp.sh
# First created: 14/Jun/2016
# Last modified: 14/Jun/2016
#######################
# [PATCH] SLX-10409.v2
if [ $IS_PE -eq 1 ]; then
	if [ $IS_FLD -eq 1 ]; then
		My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.bam # Fludigm
	else
		My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.bam # WGBS/SureSelect
	fi
else
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.bam 
fi

if [ ! -s $My_bam ];then
	echo -e "\e[031m$My_bam not found\e[0m\n"
	exit
fi
# sort it
if [ ! -s ${My_bam%.bam}.sorted.bam ];then
	echo -e "samtools sort -m 10000000000 $My_bam ${My_bam%.bam}.sorted"
	time samtools sort -m 10000000000 $My_bam ${My_bam%.bam}.sorted # 10GB mem
fi
# index it
if [ ! -s ${My_bam%.bam}.sorted.bam.bai ];then
	if [ -s ${My_bam%.bam}.sorted.bam ];then
		echo -e "samtools index ${My_bam%.bam}.sorted.bam "
		time samtools index ${My_bam%.bam}.sorted.bam 
	else
		echo -e "\e[031m${My_bam%.bam}.sorted.bam not found\e[0m\n"
		exit
	fi
fi
# check the result
if [ ! -s ${My_bam%.bam}.sorted.bam ]; then
	echo -e "\e[031m${My_bam%.bam}.sorted.bam not found. samtools sorting failed?\e[0m\n"
	exit
fi
if [ ! -s ${My_bam%.bam}.sorted.bam.bai ]; then
	echo -e "\e[031m${My_bam%.bam}.sorted.bam.bai not found. samtools indexing failed?\e[0m\n"
	exit
fi
echo -e "\e[32m$SLX.$Barcode done for sorting bam by position\e[0m" 
# END OF PATCH
#######################

