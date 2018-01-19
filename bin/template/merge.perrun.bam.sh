#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 29/Jul/2014
# Last modified: 29/Jul/2014 
# Optimised and customised to run at the Darwin HPC

source $HOME/Methyl-Seq/config/methyl_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

SLX="MY_SLX" # e.g. "SLX-8074.SLX-8080"  "SLX-8075.SLX-8077.SLX-8081" 
SLX_ARRAY=(${SLX//./ }) # replace . to a space (e.g. "SLX-8075.SLX-8077.SLX-8081" => "SLX-8075 SLX-8077 SLX-8081")
PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. ~/scratch/results/SLX-8074.SLX-8080.v1
Barcode="MY_BARCODE" # e.g. A001

# define output after merging
if [ $IS_PE -eq 1 ]; then
	# SLX-8074.SLX-8080.A010.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.bam
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.bam
else
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.bam 
fi
#######################################
# 1. MERGE BAM FILE FROM PREVIOUS RUN #
# IF YOU HAVE MORE THAN TWO RUNS
#######################################
if [ $RUN_MERGE_OTHER_BAM -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/Bismark
	mkdir_unless $PROJECT_DIR/Bismark/Alignment
	mkdir_unless $PROJECT_DIR/Bismark/Alignment/$Barcode

	My_bam_array=() # initialise it
	for THIS_SLX in ${SLX_ARRAY[*]} 
	do
		echo $THIS_SLX
		# input
		THIS_PROJECT_DIR=$RESULT_DIR/$THIS_SLX.MY_VERSION # e.g. SLX-8074.v1
		if [ $IS_PE -eq 1 ]; then
			My_dummy_bam=$THIS_PROJECT_DIR/Bismark/Alignment/$Barcode/$THIS_SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.sorted.bam 
		else
			My_dummy_bam=$THIS_PROJECT_DIR/Bismark/Alignment/$Barcode/$THIS_SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.sorted.bam 
		fi
		if [ -s $My_dummy_bam ]; then
			My_bam_array+=($My_dummy_bam) #push it
		else
			echo -e "\e[031m $My_dummy_bam not found\e[0m\n"
			exit
		fi
	done
	Merged_perrun_bam=$(printf " %s" "${My_bam_array[@]}")

	# if there are more than 1 bam file
	if [ ${#My_bam_array[@]} -gt 1 ]; then
		echo -e "samtools merge $My_bam $Merged_perrun_bam &"
		time samtools merge $My_bam $Merged_perrun_bam &
	else
		echo -e "\e[031m There is only one bam file: $Merged_perrun_bam\e[0m\n"
		exit
	fi
	echo -e "\e[32m$SLX.$Barcode done for merge.perrun.bam \e[0m" 
fi #end of RUN_MERGE_OTHER_BAM

##########################################################
# 2. MERGE BAM FILE FROM PREVIOUS RUN BY CHROMOSOME WISE #
# IF YOU HAVE MORE THAN TWO RUNS
##########################################################
if [ $RUN_MERGE_OTHER_BAM_BY_CHR -eq 1 ]; then
	for Chr in ${UCSC_CHR[*]}
	do
		My_bam_array=() # initialise it
		for THIS_SLX in ${SLX_ARRAY[*]} 
		do
			# input
			THIS_PROJECT_DIR=$RESULT_DIR/$THIS_SLX.MY_VERSION # e.g. SLX-8074.v1
			if [ $IS_PE -eq 1 ]; then
				My_dummy_chr_bam=$THIS_PROJECT_DIR/Bismark/Alignment/$Barcode/$THIS_SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.$Chr.bam 
			else
				My_dummy_chr_bam=$THIS_PROJECT_DIR/Bismark/Alignment/$Barcode/$THIS_SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.$Chr.bam 
			fi
			if [ -s $My_dummy_chr_bam ]; then
				My_bam_array+=($My_dummy_chr_bam) #push it
			else
				echo -e "\e[031m $My_dummy_chr_bam not found\e[0m\n"
				exit
			fi
		done # end of for SLX
		Merged_perchr_bam=$(printf " %s" "${My_bam_array[@]}")

		# output
		if [ $IS_PE -eq 1 ]; then
			# SLX-8074.SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.chr1.bam
			My_Chr_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.$Chr.bam 
		else
			My_Chr_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.$Chr.bam 
		fi
		# if there are more than 1 bam file
		if [ ${#My_bam_array[@]} -gt 1 ]; then
			echo -e "samtools merge $My_Chr_bam $Merged_perchr_bam &"
			time samtools merge $My_Chr_bam $Merged_perchr_bam &
		else
			echo -e "\e[031m There is only one bam file: $Merged_perchr_bam\e[0m\n"
			exit
		fi
	done # end of for Chr
	wait
	echo -e "\e[32m$SLX.$Barcode done for merge.perrun.bam chromosome-wise\e[0m" 
fi #end of RUN_MERGE_OTHER_BAM_BY_CHR

##############################################
## 3. SORT THE FINAL BAM FILE THEN INDEX IT ##
## Input: unsorted bam file 
## Output: sorted.bam
##############################################
wait # wait for the previous process
if [ $RUN_SAM_SORT_BY_COORD2 -eq 1 ]; then
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam not found\e[0m\n"
		exit
	fi
	# output
	# SLX-8074.SLX-8080.A010.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.sorted.bam
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
fi

##################################### 
## DE-DUPLICATE MERGED BAM         ##
## INPUT: q$MAPQ.bam               ##
## OUTPUT: q$MAPQ.deduplicated.bam ##
#####################################
## This step will be omitted
## as it only de-dup ~1.5% of the total reads
if [ $RUN_BM_DEDUP2 -eq 1 ]; then
	# input
	if [ $IS_PE -eq 1 ]; then
		# SLX-8074.SLX-8080.A010.r_1_val_1.fq.gz_bismark_bt2_pe.q1.bam
		My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.bam
	else
		My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.bam 
	fi
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam not found\e[0m\n"
		exit
	fi
	# output
	# SLX-8074.SLX-8080.A010.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.deduplicated.bam
	if [ ! -s ${My_bam%.bam}.deduplicated.bam ]; then
		echo -e "deduplicate_bismark $My_bam &"
		time deduplicate_bismark --paired --bam $My_bam &
	fi

	for Chr in ${UCSC_CHR[*]}
	do
		# input 
		if [ $IS_PE -eq 1 ]; then
			# SLX-8074.SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.chr1.bam
			My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.$Chr.bam 
		else
			My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.$Chr.bam 
		fi
		if [ ! -s $My_bam ];then
			echo -e "\e[031m$My_bam not found\e[0m\n"
			exit
		fi
		# output
		# SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.Chr1.deduplicated.bam
		if [ ! -s ${My_bam%.bam}.deduplicated.bam ]; then
			echo -e "deduplicate_bismark $My_bam &"
			time deduplicate_bismark --paired --bam $My_bam &
		fi
	done
	wait
	echo -e "\e[32m$SLX.$Barcode done for deduplicate_bismark\e[0m" 
fi # end of RUN_BM_DEDUP2
