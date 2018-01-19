#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 30/Jan/2015
# Last modified: 30/Jan/2015 
# Optimised and customised to run at the SBS machine

source $HOME/Pipelines/config/target.seq.sbs.config # calls global.config, slurm.config, samples.config 
source $HOME/lib/sung.sh

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

mkdir_unless $PROJECT_DIR	

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample.Lane=$SLX.$Barcode.$Lane\e[0m\n"

#/disk1/ssg29/data/fastq/SLX-9338/SLX-9338.HAL21.000000000-ACU88.s_1.r_1.fq.gz
#/disk1/ssg29/data/fastq/SLX-9338/SLX-9338.HAL21.000000000-ACU88.s_1.r_2.fq.gz
FASTQ_FILES=`ls $FASTQ_DIR/$SLX.$Barcode.$Cell.$Lane.*.fq.gz | grep -v lost`

##################################################
# 2. run initial fastqc 
# output: results/SLX123/FastQC/A001/
##################################################
if [ $RUN_FASTQC -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/FastQC
	mkdir_unless $PROJECT_DIR/FastQC/$Barcode
	for FastQ_file in $FASTQ_FILES
	do
		echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode"
		time fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode 
	done
fi

############################################################
# 3. Trim FastQ file per lane 	
# input: .fastq.gz
# output: results/SGA_v6/Trim/A001/SLX-8074.A001.C2N01ACXX.s_1.r_1.1_val_1.fq.gz
# output: results/SGA_v6/Trim/A001/SLX-8074.A001.C2N01ACXX.s_1.r_1.1_val_2.fq.gz
# output: results/SGA_v6/Trim/A001/SLX-8074.A001.C2N01ACXX.s_1.r_1.1.fastq.gz_trimming_report.txt
############################################################
if [ $RUN_TRIM -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/Trim
	mkdir_unless $PROJECT_DIR/Trim/$Barcode
	for FastQ_file in $FASTQ_FILES
	do
		FastQ_array+=($FastQ_file) #push it
	done
	dummy_fastq=$(printf " %s" "${FastQ_array[@]}")

	if [ $IS_PE -eq 1 ]; then
		echo -e "time trim_galore --paired --adapter $TR_ADAPTOR1 --adapter2 $TR_ADAPTOR2 --output_dir $PROJECT_DIR/Trim/$Barcode --fastqc $dummy_fastq\n"
		time trim_galore \
			--paired \
			--length $TR_LEN \
			--fastqc \
			--adapter $TR_ADAPTOR1 \
			--adapter2 $TR_ADAPTOR2 \
			--output_dir $PROJECT_DIR/Trim/$Barcode \
			--clip_R1 $TR_CLIP_LEN \
			--clip_R2 $TR_CLIP_LEN \
			$dummy_fastq
	else
		echo -e "time trim_galore --adapter $TR_ADAPTOR1 --adapter2 $TR_ADAPTOR2 --output_dir $PROJECT_DIR/Trim/$Barcode --fastqc $dummy_fastq\n"
		time trim_galore \
			--fastqc \
			--adapter $TR_ADAPTOR1 \
			--adapter2 $TR_ADAPTOR2 \
			--output_dir $PROJECT_DIR/Trim/$Barcode \
			--clip_R1 $TR_CLIP_LEN \
			$dummy_fastq
	fi
fi # end of run trim

if [ $IS_PE -eq 1 ]; then
	Trimmed_Fw_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_1_val_1.fq.gz 
	Trimmed_Rv_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_2_val_2.fq.gz 
	if [ ! -s $Trimmed_Fw_FastQ ];then
		echo -e "\e[031m$Trimmed_Fw_FastQ not found\e[0m\n"
		exit
	fi
	if [ ! -s $Trimmed_Rv_FastQ ];then
		echo -e "\e[031m$Trimmed_Rv_FastQ not found\e[0m\n"
		exit
	fi
	TRIMMED_FASTQ_FILES="$Trimmed_Fw_FastQ $Trimmed_Rv_FastQ"
else
	#SLX-8764.A001.C4EYHANXX.s_7.r_1_trimmed.fq.gz
	TRIMMED_FASTQ_FILES="$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_1"_trimmed.fq.gz
	if [ ! -s $TRIMMED_FASTQ_FILES ];then
		echo -e "\e[031m$TRIMMED_FASTQ_FILES not found\e[0m\n"
		exit
	fi
fi

MY_BWA_OUT=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.$Lane.bwa.sam
if [ $RUN_BWA -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/BWA
	mkdir_unless $PROJECT_DIR/BWA/$Barcode
	#-R  STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
	echo -e "bwa $TRIMMED_FASTQ_FILES \n"
	time bwa mem \
		-M \
		-t $NT \
		-R "@RG\tID:$SLX\tPL:illumina\tPU:run\tLB:$SLX\tSM:$Barcode\tCN:CamObsGynae" \
		$BWA_INDEX_BASE".fa" \
		$TRIMMED_FASTQ_FILES > $MY_BWA_OUT

	if [ -s $MY_BWA_OUT ];then
		echo -e "samtools view -bS $MY_BWA_OUT  > ${MY_BWA_OUT%.sam}.bam"
		time samtools view -bS $MY_BWA_OUT  > ${MY_BWA_OUT%.sam}.bam 

		if [ -s ${MY_BWA_OUT%.sam}.bam ];then
			echo -e "rm -f $MY_BWA_OUT \n"
			time rm -f $MY_BWA_OUT 
		else
			echo -e "\e[031m${MY_BWA_OUT%.sam}.bam not found\e[0m\n"
			exit
		fi
	else
		echo -e "\e[031m$MY_BWA_OUT not found\e[0m\n"
		exit
	fi
fi
My_bam=${MY_BWA_OUT%.sam}.bam

echo -e "\e[32m$Barcode done for mapping\e[0m" 
