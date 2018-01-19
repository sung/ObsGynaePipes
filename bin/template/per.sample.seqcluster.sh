#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 21/Jul/2014
# Last modified: 21/Jul/2014 
# Optimised and customised to run at the Darwin HPC
# Based on http://seqcluster.readthedocs.org/mirna_annotation.html

source $HOME/RNA-Seq/config/rna_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001
FASTQ_DIR=$HOME/data/fastq/$SLX

mkdir_unless $PROJECT_DIR	
mkdir_unless $PROJECT_DIR/scratch
mkdir_unless $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Barcode\e[0m\n"

<<fastq
SLX-9169.D704_D501.C6A20ANXX.s_1.r_1.fq.gz
SLX-9169.D704_D501.C6GNPANXX.s_5.r_1.fq.gz
SLX-9169.D704_D501.C6GNPANXX.s_6.r_1.fq.gz
SLX-9169.D704_D501.C6GW5ANXX.s_1.r_1.fq.gz
SLX-9169.D704_D501.C6GW5ANXX.s_2.r_1.fq.gz
SLX-9169.D704_D501.C6GW5ANXX.s_3.r_1.fq.gz
fastq

# assuming single-end
# [todo] what if paired end?
if [ $IS_SE -eq 1 ]; then
# per each barcode (sample) 
# do 1)FastQC 2)trim 3)tophat 4)bedtool, 5)htseq 
	Trimmed_FastQ_array=() #initialise
	FASTQ_FILES=`ls $FASTQ_DIR/$SLX.$Barcode*.fq.gz | grep -v lost`
	for FastQ_file in $FASTQ_FILES	# per fastq file (per lane)
	do
		##################################################
		# 0. Run the inital fastqc 
		##################################################
		if [ $RUN_FASTQC -eq 1 ]; then
			mkdir_unless $PROJECT_DIR/FastQC
			mkdir_unless $PROJECT_DIR/FastQC/$Barcode
			echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode"
			time fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode 
		fi
		#/home/ssg29/data/fastq/SLX-9176/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1.fq.gz
		#/home/ssg29/results/SLX-9176.Homo_sapiens.v1/Trim/NEBsmRNA01.v1/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1_trimmed.fq.gz
		#/home/ssg29/results/SLX-9176.Homo_sapiens.v1/Trim/NEBsmRNA01.v1/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1.fq.gz_trimming_report.txt
		Trimmed_FastQ=${FastQ_file/data\/fastq\/$SLX/results\/$PROJECT\/Trim\/$Barcode}
		Trimmed_FastQ=${Trimmed_FastQ%.fq.gz}_trimmed.fq.gz #SLX-8546.D709_D501.C4FN0ANXX.s_6.r_1_trimmed.fq.gz 
		############################################################
		# 1. Trimming via cutadapt 
		############################################################
		if [ $RUN_TRIM_NEB -eq 1 ]; then
			mkdir_unless $PROJECT_DIR/Trim
			mkdir_unless $PROJECT_DIR/Trim/$Barcode
			echo -e "\ncutadapt $FastQ_file"
			time cutadapt \
				--trimmed-only \
				-q $TR_QUAL \
				-O $TR_STRNCY \
				-m $TR_LEN \
				-a $TR_ADAPTOR1 \
				-g $TR_ADAPTOR2 \
				-o $Trimmed_FastQ \
				$FastQ_file &> ${Trimmed_FastQ%_trimmed.fq.gz}.fq.gz_trimming_report.txt
			if [ ! -s $Trimmed_FastQ ];then
				echo -e "\e[031m$Trimmed_FastQ not found\e[0m\n"
				exit
			else
				echo -e "\nfastqc $Trimmed_FastQ -o $PROJECT_DIR/Trim/$Barcode"
				time fastqc $Trimmed_FastQ -o $PROJECT_DIR/Trim/$Barcode 
			fi
		fi # end of run trim

		if [ -s $Trimmed_FastQ ]; then
			Trimmed_FastQ_array+=($Trimmed_FastQ) #push it multi-laned fastq
		else
			echo -e "\e[031m$Trimmed_FastQ not found\e[0m\n"
			exit
		fi
	done # end of per FastQ_file 

	##
	## Merged FastQ
	## 
	Merged_FastQ=$(printf ",%s" "${Trimmed_FastQ_array[@]}")
	Merged_FastQ=${Merged_FastQ:1} # remove first comma 
	echo -e "\nzcat $Merged_FastQ | zip -9 > $PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.trim.merged.fq.gz"
	zcat $Merged_FastQ | zip -9 > $PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.trim.merged.fq.gz 

	##
	## 2. Collapse Reads
	##
	mkdir_unless $PROJECT_DIR/Seqcluster
	mkdir_unless $PROJECT_DIR/Seqcluster/collapsed
	mkdir_unless $PROJECT_DIR/Seqcluster/collapsed/$Barcode
	PYTHONPATH=$HOME/Install/seqcluster/anaconda/lib/python2.7/site-packages:$PYTHONPATH
	time seqcluster collapse -f $PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.trim.merged.fq.gz -o $PROJECT_DIR/Seqcluster/collapsed/$Barcode/ 

	##
	## 3. Aligner 
	##
	mkdir_unless $PROJECT_DIR/Seqcluster
	mkdir_unless $PROJECT_DIR/Seqcluster/mirAligner
	mkdir_unless $PROJECT_DIR/Seqcluster/mirAligner/$Barcode
	Collapsed_FastQ=$PROJECT_DIR/Seqcluster/collapsed/$Barcode/$SLX.$Barcode.trim.merged_trimmed.fastq
	#SLX-9176.NEBsmRNA01.trim.merged_trimmed.fastq
	if [ -s $Collapsed_FastQ ]; then
		time java -jar ~/Install/miraligner.jar \
			-freq \
			-pre \
			-sub 1 \
			-trim 3 \
			-add 3 \
			-s hsa \
			-i $Collapsed_FastQ \
			-db $HOME/data/Annotation/smallRNA \
			-o $PROJECT_DIR/Seqcluster/mirAligner/$Barcode/$SLX.$Barcode \
			2> $PROJECT_DIR/Seqcluster/mirAligner/$Barcode/$SLX.$Barcode.miraligner.log
	else
		echo -e "\e[031m$Collapsed_FastQ not found\e[0m\n"
		exit
	fi
else
	echo -e "\e[031This script is for signle-end only\e[0m\n"
	exit
fi

echo -e "\e[32m$Barcode done for mapping\e[0m" 
