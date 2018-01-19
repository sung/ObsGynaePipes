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

#####################################################################################################################
# IF paired-end, BAM file should be sorted by the read name (samtools sort -n) before running bismark_methylextract #
# if bismark was used, the BAM must have been sorted by the read name
#####################################################################################################################
if [ $IS_PE -eq 1 ]; then
	if [ $IS_FLD -eq 1 ]; then
		My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.bam  # for Fluidigm
	else
		My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.bam 
	fi
	BISMARK_OPT="--paired-end --no_overlap --ignore $BM_IGNORE1 --ignore_r2 $BM_IGNORE2"
else
	# PNAS-2013.SRR530647.r_1_trimmed.fq.gz_bismark_bt2.q1.deduplicated.bam
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.bam 
	BISMARK_OPT="--single-end"
fi

if [ ! -s $My_bam ];then
	echo -e "\e[031m$My_bam not found\e[0m\n"
	exit
fi

	#########################################################################
	# 10. Bismark Methyl Calls 
	# time:
	# input: A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.q20.bam 
	# output: results/SGA_v1/Bismark/MethylCall/A001/ 
	# 1) CpG_context_A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.txt.gz
	# 2) Non_CpG_context_A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.txt.gz
	# 3) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.M-bias_R2.png
	# 4) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.bam_splitting_report.txt
	# 5) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.M-bias.txt
	# 6) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.M-bias_R1.png
	# 7) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.bismark.cov
	# 8) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.bedGraph
	# 9) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.CpG_report.txt
	#
	# --comprehensive: merge different strand (OT, OB, CTOT, CTOB)
	# --merge_non_CpG: CpG vs non-CpG (rather than CpG, CHG, CGG)
	# --bedGraph calls 'bismark2bedGraph'
	# --cytosine_report calls 'bedGraph2cytosine', which is now renamed as 'coverage2cytosine' 
	# --CX_context: both CG (default) and CH context
	# The option '--cytosine_report' above run this in one go, but if you prefer run separately:
	# alterantively, time coverage2cytosine --genome_folder $GENOME_DIR --dir $PROJECT_DIR/Bismark/MethylCall/$Barcode --output $Barcode.cov2cytosine.txt $Barcode.cov
	# this produces genome-wide methylation report for all cytosines in the genome 
	# --gzip: This option does not work on bedGraph and genome-wide cytosine reports as they are 'tiny' anyway
	# [todo] other parameters? (e.g. --ignore after M-bias plot)
	#######################
	if [ $RUN_BM_METHYL -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/Bismark/MethylCall
		mkdir_unless $PROJECT_DIR/Bismark/MethylCall/$Barcode
		echo bismark_methylation_extractor $Barcode
		time bismark_methylation_extractor \
			$BISMARK_OPT \
			--comprehensive \
			--merge_non_CpG \
			--report --gzip \
			--buffer_size $BM_BUFFER_SIZE \
			--bedGraph \
			--cytosine_report \
			--output $PROJECT_DIR/Bismark/MethylCall/$Barcode \
			--genome_folder $GENOME_DIR \
			$My_bam
		echo -e "\e[32m$SLX.$Barcode done for Bismark calling\e[0m" 
	fi
