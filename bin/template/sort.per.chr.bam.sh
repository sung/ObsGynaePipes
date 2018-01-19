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

###############################
# sort bam file by read name ##
# to run bismark methyl extract
###############################
for Chr in ${UCSC_CHR[*]}
do
	###############################
	# input
	# SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.Chr1.bam
	# SLX-8074.SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.Chr1.bam
	###############################
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.$Chr.bam 
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam not found\e[0m\n"
		exit
	fi

	###############################
	# output: name-based sorted chr bam file
	# e.g. SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.Chr1.sorted.bam
	# e.g. SLX-8074.SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.Chr1.sorted.bam
	###############################
	if [ $RUN_SAM_SORT_BY_NAME -eq 1 ]; then
		if [ ! -s ${My_bam%.bam}.sorted.bam ]; then
			echo -e "time samtools sort -n -m 3000000000 $My_bam ${My_bam%.bam}.sorted &"
			time samtools sort -n -m 3000000000 $My_bam ${My_bam%.bam}.sorted & # 3GB mem
		fi
		# NB. you cannot index name-based sorted bam
	fi

	wait
	if [ $RUN_SAM_SORT_BY_NAME -eq 1 ]; then
		echo -e "\e[32m$SLX.$Barcode done for sort by name\e[0m" 
	fi

	############################################################
	# 1. Add RG (Read Group) to bam header
	# 'CREATE_INDEX=true' will make unknown error - disabled
	# input
	# SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.Chr1.bam
	# SLX-8074.SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.Chr1.bam
	# output: my.RG.bam
	# this step is to run BisSNP per chromosome later
	############################################################
	if [ $RUN_PICARD_RG_BY_CHR -eq 1 ]; then
		if [ ! -s ${My_bam%.bam}.RG.bam ]; then
			echo java -Xmx$RG_MEM -jar $PICARD/AddOrReplaceReadGroups.jar $My_bam
			time java -Xmx$RG_MEM -jar $PICARD/AddOrReplaceReadGroups.jar \
				INPUT=$My_bam \
				OUTPUT=${My_bam%.bam}.RG.bam \
				RGID=$SLX \
				RGLB=$SLX \
				RGPL=illumina \
				RGPU=run \
				RGSM=$Barcode \
				RGCN=$RG_CN \
				VALIDATION_STRINGENCY=$RG_STRINGENCY \
				SORT_ORDER=$RG_ORDER
		else
			echo -e "\e[31m${My_bam%.bam}.RG.bam not found\e[0m" 
			exit
		fi
		echo -e "\e[32m$SLX.$Barcode done for Picard RG\e[0m" 
	fi

	if [ $RUN_PICARD_RG_INDEX_BY_CHR -eq 1 ]; then
		if [ -s ${My_bam%.bam}.RG.bam ]; then
			echo -e "samtools index ${My_bam%.bam}.RG.bam "
			time samtools index ${My_bam%.bam}.RG.bam 
		else
			echo -e "\e[031m${My_bam%.bam}.RG.bam not found\e[0m\n"
			exit
		fi
		echo -e "\e[32m$SLX.$Barcode done for Picard RG indexing\e[0m" 
	fi

	echo -e "\e[32m$SLX.$Barcode done for sort.per.chr.bam\e[0m" 
done




