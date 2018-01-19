#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 21/Jul/2014
# Last modified: 21/Jul/2014 
# assumes 1 case and 1 control separated by ',' (e.g. D709_D501,D709_D502)
# Optimised and customised to run at the Darwin HPC

source $HOME/RNA-Seq/config/rna_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

export SLX="MY_SLX" # e.g. SLX-8080 
export PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1
export Barcode="MY_BARCODE" # e.g. A001
export Cell="MY_CELL" # e.g. C48CWACXX  
export Lane="MY_LANE" # e.g. s_1 
export Chunk="MY_CHUNK" # e.g. 1 

if [ $NUM_CHUNK -eq 1 ];then
	FASTQ_DIR=$HOME/scratch/data/fastq/$SLX
else
	FASTQ_DIR=$HOME/scratch/data/fastq/$SLX/split
fi

mkdir_unless $PROJECT_DIR	
mkdir_unless $PROJECT_DIR/scratch
mkdir_unless $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Barcode\e[0m\n"
echo -e  "Max insertion size=$BM_MAXINS"

<<fastq
SLX-8546.D709_D501.C4FN0ANXX.s_4.r_1.fq.gz 
SLX-8546.D709_D501.C4FN0ANXX.s_5.r_1.fq.gz 
SLX-8546.D709_D501.C4FN0ANXX.s_6.r_1.fq.gz 
fastq

# assuming single-end
# [todo] what if paired end?
if [ $IS_SE -ne 1 ]; then
	echo -e "\e[031This script is for signle-end only\e[0m\n"
	exit
fi

# per each barcode (sample) 
# do 1).FastQC 2).trim 3.tophat 4.cufflink 

	if [ $NUM_CHUNK -eq 1 ];then
		FastQ_file=$FASTQ_DIR/$SLX.$Barcode.$Cell.$Lane.r_1.fq.gz # SLX-8546.D709_D501.C4FN0ANXX.s_4.r_1.fq.gz
	else
		FastQ_file=$FASTQ_DIR/$SLX.$Barcode.$Cell.$Lane.r_1.$Chunk.fastq.gz 
	fi

	if [ ! -s $FastQ_file ]; then
		echo -e "\e[031m$FastQ_file not found\e[0m\n"
		exit
	fi
	##################################################
	# 1. Run the inital fastqc 
	##################################################
	if [ $RUN_FASTQC -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/FastQC
		mkdir_unless $PROJECT_DIR/FastQC/$Barcode
		echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode"
		time fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode 
	fi
	############################################################
	# 2. Trim 	
	############################################################
	if [ $RUN_TRIM -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/Trim
		mkdir_unless $PROJECT_DIR/Trim/$Barcode
		echo -e "\ntrim_galore -a $TR_ADAPTOR1 -o $PROJECT_DIR/Trim/$Barcode --fastqc --quality $TR_QUAL --stringency $TR_STRNCY $FastQ_file"
		time trim_galore \
			--output_dir $PROJECT_DIR/Trim/$Barcode \
			--adapter $TR_ADAPTOR1 \
			--quality $TR_QUAL \
			--stringency $TR_STRNCY \
			$FastQ_file
			#--fastqc \
	fi # end of run trim

	if [ $NUM_CHUNK -eq 1 ];then
		Trimmed_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_1_trimmed.fq.gz 
	else
		Trimmed_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_1.$Chunk"_trimmed.fq.gz"
	fi

	if [ ! -s $Trimmed_FastQ ]; then
		echo -e "\e[031m$Trimmed_FastQ not found\e[0m\n"
		exit
	fi

	####################################################################
	# 2. tophat: Map the reads for each sample to the reference genome #
	# real    61m11.381s (for R3 with 2 threads) 82m22.024s (for R4)
	# output: $TOPHAT_OUT/accepted_hits.bam
	####################################################################
	TOPHAT_OUT=$PROJECT_DIR/TopHat/$Barcode/$Lane
	if [ $RUN_TOPHAT -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/TopHat
		mkdir_unless $PROJECT_DIR/TopHat/$Barcode
		mkdir_unless $TOPHAT_OUT
		echo -e "\ntophat2 -p $NT --library-type $LIB_TYPE --output-dir $TOPHAT_OUT --max-multihits $TH_MH $TH_PREFILTER -transcriptome-index=$TR_INDEX_DIR/$TR_PREFIX $BOWTIE2_INDEX_BASE $Trimmed__FastQ"

		time tophat2 \
			--num-threads $NT \
			--library-type $LIB_TYPE \
			--output-dir $TOPHAT_OUT \
			--max-multihits $TH_MH \
			$TH_PREFILTER \
			--transcriptome-index=$TR_INDEX_DIR/$TR_PREFIX \
			$BOWTIE2_INDEX_BASE $Trimmed_FastQ
	fi
