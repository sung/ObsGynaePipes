#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 17/Aug/2015
# Last modified: 17/Aug/2015
# Optimised and customised to run at the Darwin HPC

source $HOME/RNA-Seq/config/rna_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001

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
# do 1)FastQC 2)trim 3)star (2-pass mode) 4)bedtool, 5)htseq 
	Trimmed_FastQ_array=() #initialise
	FASTQ_FILES=`ls $HOME/data/fastq/$SLX/$SLX.$Barcode*.fq.gz | grep -v lost`
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
		############################################################
		# 1. Trimming via trim_galore 	
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

		#/home/ssg29/data/fastq/SLX-9176/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1.fq.gz
		#/home/ssg29/results/SLX-9176.Homo_sapiens.v1/Trim/NEBsmRNA01.v1/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1_trimmed.fq.gz
		#/home/ssg29/results/SLX-9176.Homo_sapiens.v1/Trim/NEBsmRNA01.v1/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1.fq.gz_trimming_report.txt
		Trimmed_FastQ=${FastQ_file/data\/fastq\/$SLX/results\/$PROJECT\/Trim\/$Barcode}
		Trimmed_FastQ=${Trimmed_FastQ%.fq.gz}_trimmed.fq.gz #SLX-8546.D709_D501.C4FN0ANXX.s_6.r_1_trimmed.fq.gz 

		if [ -s $Trimmed_FastQ ]; then
			Trimmed_FastQ_array+=($Trimmed_FastQ) #push it
		else
			echo -e "\e[031m$Trimmed_FastQ not found\e[0m\n"
			exit
		fi
	done # end of per FastQ_file 
	###################################################################
	# 2. STAR 
	# Map the reads for each sample to the reference genome           #
	###################################################################
	STAR_OUT=$PROJECT_DIR/STAR/$Barcode 
	STAR_TMP=$PROJECT_DIR/scratch/$Barcode/STAR # will be created by STAR
	if [ $RUN_STAR -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/STAR
		mkdir_unless $STAR_OUT
		Merged_FastQ=$(printf ",%s" "${Trimmed_FastQ_array[@]}")
		Merged_FastQ=${Merged_FastQ:1} # remove the first comma 

		time STAR \
			--runThreadN $NT \
			--genomeDir $STAR_INDEX_DIR \
			--readFilesIn $Merged_FastQ \
			--readFilesCommand zcat \
			--outFileNamePrefix $STAR_OUT \
			--outTmpDir $STAR_TMP \
			--outReadsUnmapped Fastx \
			--twopassMode Basic \
			--sjdbGTFfile $STAR_GTF \
			--outSAMtype BAM SortedByCoordinate \
			--outSAMattrRGline ID:$SLX PL:illumina PU:run LB:$SLX SM:$Barcode CN:CamObsGynae 
			#--outFilterIntronMotifs RemoveNoncanonical \ # for cufflinks
	fi # end of RUN_STAR

	######################
	# 3. CoverageBed
	# Input: bam from tophat
	# (NB, to use 'bedtools genomecov', BAM _must_ be sorted by position, 
	# which is already done by tophat)
	# Output: $sample.genomecov.txt
	######################

	######################
	# 4. HTSeq
	# For paired-end data, the alignment have to be sorted either by read name or by alignment position.
	# 'sort -T $SCRATCH_OUT -k 1,1' will do a sort by read name, or
	# samtools sort -n $TOPHAT_OUT/accepted_hits.bam TOPHAT_OUT/accepted_hits.sorted.bam
	# output:
	######################

else
	echo -e "\e[031This script is for signle-end only\e[0m\n"
	exit
fi

echo -e "\e[32m$Barcode done for mapping\e[0m" 
