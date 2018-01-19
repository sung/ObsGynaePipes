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
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

	####################################################
	## 1. Get per lane MAPQ filtered bam files
	## e.g. SLX-8081.A002.C40D1ANXX.s_5.r_1_val_1.fq.gz_bismark_bt2_pe.q1.bam
	####################################################
	BAM_FILES=`ls $PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode*q$MAPQ.bam`
	for My_bam in $BAM_FILES # per bam file (which is per-lane)
	do
		My_bam_array+=(${My_bam}) #push it
	done
	Merged_perlane_bam=$(printf " %s" "${My_bam_array[@]}")

	# [sorted]
	# per-lane bam files need to be sorted before merging?
	# samtools merge is to 'Merge multiple sorted alignments' (http://samtools.sourceforge.net/samtools.shtml)
	# even without per-lane bam sorting, NO of reads are same from the merged BAM

	################################
	## 2. Merge per-lane bam file  #
	################################
	if [ $IS_PE -eq 1 ]; then
		My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.bam 
	else
		# e.g. PNAS-2013.SRR530647.CELL.s_1.r_1_trimmed.fq.gz_bismark_bt2.bam
		My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.bam 
	fi

	if [ $RUN_MERGE_PER_LANE_BAM -eq 1 ]; then
		# if there is only one lane
		if [ ${#My_bam_array[@]} -eq 1 ]; then
			echo -e "ln -s $Merged_perlane_bam ${My_bam%.bam}.q$MAPQ.bam"
			time ln -s $Merged_perlane_bam ${My_bam%.bam}.q$MAPQ.bam 
		else
		# only merge if there are more than two bam files
			echo -e "samtools merge ${My_bam%.bam}.q$MAPQ.bam  $Merged_perlane_bam"
			time samtools merge ${My_bam%.bam}.q$MAPQ.bam  $Merged_perlane_bam
			echo -e "\e[32m$SLX.$Barcode done for per-lane bam merging\e[0m" 
		fi
	fi
	# output
	# This will be made after merging
	# e.g. SLX-8080.v1/Bismark/Alignment/A012/SLX-8080.A012.r_1_val_1.fq.gz_bismark_bt2_pe.q1.bam

	##################################################### 
	# 3. Bismark De-duplicate reads
	## input: q$MAPQ.bam
	## output: q$MAPQ.deduplicated.bam
	##################################################### 
	if [ $RUN_BM_DEDUP -eq 1 ]; then
		if [ -s ${My_bam%.bam}.q$MAPQ.bam ]; then
			echo deduplicate_bismark ${My_bam%.bam}.q$MAPQ.bam 
			if [ $IS_PE -eq 1 ]; then
				time deduplicate_bismark --paired --bam ${My_bam%.bam}.q$MAPQ.bam
			else
				time deduplicate_bismark --single --bam ${My_bam%.bam}.q$MAPQ.bam
			fi
		else
			echo -e "\e[031m ${My_bam%.bam}.q$MAPQ.bam not found\e[0m\n"
			exit
		fi

		if [ -s ${My_bam%.bam}.q$MAPQ.deduplicated.bam ]; then
			# expecting  A011.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.bam
			# flagstat for q$MAPQ.deduplicated.bam
			echo -e "samtools view -u -F 0x100 ${My_bam%.bam}.q$MAPQ.deduplicated.bam | samtools flagstat - > ${My_bam%.bam}.q$MAPQ.deduplicated.bam.flagstat"
			time samtools view -u -F 0x100 ${My_bam%.bam}.q$MAPQ.deduplicated.bam | samtools flagstat - > ${My_bam%.bam}.q$MAPQ.deduplicated.bam.flagstat
		else
			echo -e "\e[031m${My_bam%.bam}.q$MAPQ.deduplicated.bam not found\e[0m\n"
			exit
		fi
		echo -e "\e[32m$SLX.$Barcode done for deduplicate_bismark\e[0m" 
		My_bam=${My_bam%.bam}.q$MAPQ.deduplicated.bam # e.g. A011.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.bam
	else
		My_bam=${My_bam%.bam}.q$MAPQ.bam # e.g. A011.r_1_val_1.fq.gz_bismark_bt2_pe.q1.bam
	fi # end of RUN_BM_DEDUP

	##################################
	# 4. MERGE Previous BAM FILE #
	# IF YOU HAVE MORE THAN TWO
	# moved to bin/template/merge.perrun.bam.sh
	##################################
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam not found\e[0m\n"
		exit
	fi
	##########################################################################
	# 5. sort the bam file then index it
	# Input:  A011.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.bam 
	# Output: A011.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.sorted.bam
	# This file will be used for BisSNP which expects a sorted BAM file
	# [NB] Make sure you've got enough memory allocated for this job! (check $NT)
	##########################################################################
	if [ $RUN_SAM_SORT_BY_COORD -eq 1 ]; then
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
	fi
	#####################################
	# 6. split bam file by chromosome  ##
	# this is to feed 'bismark_methylation_extractor' by chromosome to speed up
	# input: read-name sorted bam (this BAM file should be indexed)
	# output: per chr bam file
	#####################################
	if [ $RUN_SPLIT_BAM_BY_CHR -eq 1 ]; then
		for Chr in ${UCSC_CHR[*]}
		do
			echo -e "samtools view -b ${My_bam%.bam}.sorted.bam $Chr > ${My_bam%.bam}.$Chr.bam &"
			time samtools view -b ${My_bam%.bam}.sorted.bam $Chr > ${My_bam%.bam}.$Chr.bam &
		done
		wait
		echo -e "\e[32m$SLX.$Barcode done for splitting bam by chromosome\e[0m" 
	fi

echo -e "\e[32m$Barcode done for merging/deduplicating\e[0m" 
