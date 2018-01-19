#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 15/Jul/2014
# Last modified: 25/Jul/2014 
# Optimised and customised to run at the Darwin HPC

source $HOME/Methyl-Seq/config/methyl_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1
Barcode="MY_BARCODE" # e.g. A001

#########################################################################
# 10. Bismark Methyl Calls 
# --comprehensive: merge different strand (OT, OB, CTOT, CTOB)
# --merge_non_CpG: CpG vs non-CpG (rather than CpG, CHG, CGG)
# --bedGraph calls 'bismark2bedGraph'
# --cytosine_report calls 'bedGraph2cytosine', which is now renamed as 'coverage2cytosine' 
# The option '--cytosine_report' above run this in one go, but if you prefer run separately:
# alterantively, time coverage2cytosine --genome_folder $GENOME_DIR --dir $PROJECT_DIR/Bismark/MethylCall/$Barcode --output $Barcode.cov2cytosine.txt $Barcode.cov
# this produces genome-wide methylation report for all cytosines in the genome 
# --gzip: This option does not work on bedGraph and genome-wide cytosine reports as they are 'tiny' anyway
#######################
if [ $RUN_BM_METHYL -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/Bismark/MethylCall
	mkdir_unless $PROJECT_DIR/Bismark/MethylCall/$Barcode
	for Chr in ${UCSC_CHR[*]}
	do 
		mkdir_unless $PROJECT_DIR/Bismark/MethylCall/$Barcode/$Chr
		#####################################################################################################################
		# IF paired-end, BAM file should be sorted by the read name (samtools sort -n) before running bismark_methylextract #
		# if bismark was run, the BAM is already sorted by the read name
		# to run by chromosoime, do: run_pipeline "sort.per.chr.bam" 1  
		# this will sort chr.bam by the read name
		#####################################################################################################################
		# input
		# SLX-8080.A012.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.chr2.sorted.bam
		if [ $IS_PE -eq 1 ]; then
			My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.$Chr.sorted.bam 
			BISMARK_OPT="--paired-end --no_overlap --ignore_r2 $BM_IGNORE2"
		else
			My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.$Chr.sorted.bam 
			BISMARK_OPT="--single-end"
		fi
		if [ ! -s $My_bam ];then
			echo -e "\e[031m$My_bam not found\e[0m\n"
			exit
		fi
		echo -e "bismark_methylation_extractor $SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.$Chr.sorted.bam"
			#--comprehensive \ # merge strand information
		time bismark_methylation_extractor \
			$BISMARK_OPT \
			--merge_non_CpG \
			--report --gzip \
			--buffer_size $((${BM_BUFFER_SIZE%G}/10))G \
			--bedGraph \
			--CX_context \
			--cytosine_report \
			--output $PROJECT_DIR/Bismark/MethylCall/$Barcode/$Chr \
			--genome_folder $GENOME_DIR \
			$My_bam &
	done
	wait
	echo -e "\e[32m$SLX.$Barcode done for Bismark calling\e[0m" 
fi
