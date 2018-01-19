#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 11/Apr/2014
# Last modified: 18/Jul/2014 
# Optimised and customised to run at the Darwin HPC

source $HOME/Methyl-Seq/config/methyl_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

export SLX="MY_SLX" # e.g. SLX-8080 
export PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1
export Barcode="MY_BARCODE" # e.g. A001
export Cell="MY_CELL" # e.g. C48CWACXX  
export Lane="MY_LANE" # e.g. s_1 
export Chunk="MY_CHUNK" # e.g. 1 

FASTQ_DIR=$HOME/scratch/data/fastq/$SLX

# '--maxins' from bismark (Default: 500), Michelle identifies this from BioAnaliser
# see data/README_INFO to get the insertion size
BM_MAXINS=`grep $Barcode $LIB_SIZE_FILE | grep ${SLX/SLX-/} | cut -f7` 

if [[ ! $BM_MAXINS =~ ^-?[0-9]+$ ]]; then
	echo "No max lib insertion size found"
	exit
fi

mkdir_unless $PROJECT_DIR	
mkdir_unless $PROJECT_DIR/scratch
mkdir_unless $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Barcode\e[0m\n"
echo -e  "Max insertion size=$BM_MAXINS"


<<fastq
SLX-8075.A011.000000000-A7DHM.s_1.r_1.1.fastq.gz
SLX-8075.A011.000000000-A7DHM.s_1.r_2.1.fastq.gz
SLX-8075.A011.000000000-A7DHM.s_1.r_1.2.fastq.gz
SLX-8075.A011.000000000-A7DHM.s_1.r_2.2.fastq.gz
SLX-8075.A011.000000000-A7DHM.s_1.r_1.3.fastq.gz
SLX-8075.A011.000000000-A7DHM.s_1.r_2.3.fastq.gz
fastq

# assuming paired-end
# [todo] what if single end?
if [ $IS_PE -eq 1 ]; then
	##################################################
	# 1. Get Flowcell & Lane info
	##################################################
	Fw_FastQ=$FASTQ_DIR/$SLX.$Barcode.$Cell.$Lane.r_1.$Chunk.fastq.gz # SLX-8074.A001.C2N01ACXX.s_1.r_1.1.fastq.gz
	Rv_FastQ=$FASTQ_DIR/$SLX.$Barcode.$Cell.$Lane.r_2.$Chunk.fastq.gz # SLX-8074.A001.C2N01ACXX.s_1.r_2.1.fastq.gz

	if [ ! -s $Fw_FastQ ];then
		echo -e "\e[031m$Fw_FastQ not found\e[0m\n"
		exit

	fi
	if [ ! -s $Rv_FastQ ];then
		echo -e "\e[031m$Rv_FastQ not found\e[0m\n"
		exit
	fi

	##################################################
	# 2. run initial fastqc 
	# output: results/SGA_v1/FastQC/A001/
	##################################################
	if [ $RUN_FASTQC -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/FastQC
		mkdir_unless $PROJECT_DIR/FastQC/$Barcode
		echo -e "time fastqc $Fw_FastQ -o $PROJECT_DIR/FastQC/$Barcode"
		time fastqc $Fw_FastQ -o $PROJECT_DIR/FastQC/$Barcode 

		echo -e "time fastqc $Rv_FastQ -o $PROJECT_DIR/FastQC/$Barcode"
		time fastqc $Rv_FastQ -o $PROJECT_DIR/FastQC/$Barcode 
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
		echo -e "\ntrim_galore $Fw_FastQ $Rv_FastQ"
		time trim_galore \
			--paired --trim1 \
			--output_dir $PROJECT_DIR/Trim/$Barcode \
			--stringency $TR_STRNCY \
			$Fw_FastQ $Rv_FastQ
			#--fastqc \
	fi # end of run trim

	# e.g. SLX-8075.A011.000000000-A7DHM.s_1.r_1.1_val_1.fq.gz (SLX-8075.A011.000000000-A7DHM.s_1.r_2.1_val_2.fq.gz)
	Trimmed_Fw_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_1.$Chunk"_val_1.fq.gz"
	Trimmed_Rv_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_2.$Chunk"_val_2.fq.gz"
	if [ ! -s $Trimmed_Fw_FastQ ];then
		echo -e "\e[031m$Trimmed_Fw_FastQ not found\e[0m\n"
		exit

	fi
	if [ ! -s $Trimmed_Rv_FastQ ];then
		echo -e "\e[031m$Trimmed_Rv_FastQ not found\e[0m\n"
		exit
	fi

	###################
	## 4. BsExpress
	###################
	if [ $RUN_BS_EXPRESS -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/BsExpress
		source $BS_EXP_WRAPPER
	fi # end of BsExpress

	###########################
	# 5. Mapping with bismark #
	# input: .fq.gz (.fastq)
	# output: results/SGA_v1/Bismark/Alignment/A001/A001.r_1_val_1.fq.gz_bismark_bt2_PE_report.txt
	# output: results/SGA_v1/Bismark/Alignment/A001/A001.r_1_val_1.fq.gz_bismark_bt2_pe.bam
	# [todo] other parameters? (e.g. -L for seed length)
	###########################
	if [ $RUN_BM_ALN -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/Bismark
		mkdir_unless $PROJECT_DIR/Bismark/Alignment
		mkdir_unless $PROJECT_DIR/Bismark/Alignment/$Barcode
		echo -e "\nbismark -1 $Trimmed_Fw_FastQ -2 $Trimmed_Rv"
		time bismark \
		--bam --bowtie2 --maxins $BM_MAXINS -p $((NT/2)) \
		--temp_dir $PROJECT_DIR/scratch/$Barcode \
		--output_dir $PROJECT_DIR/Bismark/Alignment/$Barcode \
		$GENOME_DIR \
		-1 $Trimmed_Fw_FastQ \
		-2 $Trimmed_Rv_FastQ
	fi #end of bismark alignment

	####################################################
	# 7. Remove reads of poor mapping quality (MAPQ) ##
	# QC: by mapping quality (MAPQ) - remove poorly mapped reads based on MAPQ
	# do not sort bam at this stage as this breaks methyl_extractor later
	# input: BAM from bismark alignment
	# output: A011.r_1_val_1.fq.gz_bismark_bt2_pe.q20.bam
	####################################################
	# e.g. SLX-8077.A002.C48CWACXX.s_1.r_1_val_1.fq.gz_bismark_bt2_pe.bam
	# e.g. SLX-8075.A011.000000000-A7DHM.s_1.r_1.1_val_1.fq.gz_bismark_bt2_pe.bam
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_1.$Chunk"_val_1.fq.gz_bismark_bt2_pe.bam" # this will be generated after mapping 
	if [ $RUN_FLT_BY_MAPQ -eq 1 ]; then
		if [ -s $My_bam ]; then
			# filter
			echo -e "samtools view -bh -q $MAPQ $My_bam -o ${My_bam%.bam}.q$MAPQ.bam "
			time samtools view -bh -q $MAPQ $My_bam -o ${My_bam%.bam}.q$MAPQ.bam 
		else
			echo -e "\e[031m$My_bam not found\e[0m\n"
			exit
		fi

	fi

	##################################################
	# run preseq 
	##################################################
	if [ $RUN_PRESEQ -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/Preseq
		mkdir_unless $PROJECT_DIR/Preseq/$Barcode
		# sort inital bam file
		if [ -s $My_bam ]; then
			echo -e "samtools sort $My_bam ${My_bam%.bam}.sorted"
			time samtools sort $My_bam ${My_bam%.bam}.sorted
		else
			echo -e "\e[031m$My_bam not found\e[0m\n"
			exit
		fi
		if [ -s ${My_bam%.bam}.sorted.bam ]; then
			echo -e "preseq lc_extrap -v -B -P -o $PROJECT_DIR/Preseq/$Barcode/$Barcode.preseq.txt ${My_bam%.bam}.sorted.bam"
			time preseq lc_extrap -v -B -P -o $PROJECT_DIR/Preseq/$Barcode/$Barcode.preseq.txt ${My_bam%.bam}.sorted.bam 
		else
			echo -e "\e[031m${My_bam%.bam}.sorted.bam not found\e[0m\n"
			exit
		fi
	fi

fi # end of IS_PE

echo -e "\e[32m$Barcode done for mapping\e[0m" 
