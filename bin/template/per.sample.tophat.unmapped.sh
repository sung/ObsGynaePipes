#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 7/Jul/2015
# Last modified: 29/Mar/2017
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

RAW_JUNCS=$PROJECT_DIR/$PROJECT.merged.junctions.junc # also defined within the config file for project-wise

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$PROJECT.$Barcode\e[0m\n"

TOPHAT_OUT=$PROJECT_DIR/TopHat/$Barcode
TOPHAT_OUT2=$PROJECT_DIR/TopHat/$Barcode/unmapped
if [ $IS_SE -eq 1 ]; then
	###################################################################
	# 1. TOPHAT 
	# Map the unmapped reads 
	# output: $TOPHAT_OUT2/accepted_hits.bam
	###################################################################
	Unmapped_BAM=$TOPHAT_OUT/unmapped.bam
	Unmapped_FastQ=$TOPHAT_OUT/unmapped.fastq.gz
	if [ ! -s $Unmapped_FastQ ]; then
		if [ -s $Unmapped_BAM ]; then
			echo -e "\nbamutils tofastq $Unmapped_BAM | gzip > $Unmapped_FastQ"
			time bamutils tofastq $Unmapped_BAM | gzip -9 -f > $Unmapped_FastQ
		else
			echo -e "\e[031m$Unmapped_BAM not found\e[0m\n"
			exit
		fi
	fi

	if [ $RUN_TOPHAT_UNMAPPED -eq 1 ]; then
		if [ -s $Unmapped_FastQ ] && [ -s $RAW_JUNCS ]; then
			mkdir_unless $TOPHAT_OUT2
			echo -e "\ntophat2 -p $NT --library-type $LIB_TYPE --output-dir $TOPHAT_OUT2 --raw-juncs $RAW_JUNCS $Unmapped_FastQ"
			time tophat2 \
				--num-threads $NT \
				--library-type $LIB_TYPE \
				--output-dir $TOPHAT_OUT2 \
				--max-multihits $TH_MH \
				--raw-juncs $RAW_JUNCS \
				--no-novel-junc \
				$BOWTIE2_INDEX_BASE $Unmapped_FastQ
				#$TH_PREFILTER \
				#--transcriptome-index $TR_INDEX_DIR/$TR_PREFIX \
		else
			echo -e "\e[031m$Unmapped_FastQ or $RAW_JUNCS not found\e[0m\n"
			exit
		fi

		## Merge 1) inital mapped bam + 2) rescued from un-mapped bam files
		if [ -s $TOPHAT_OUT2/accepted_hits.bam ];then 
			echo -e "\nsamtools merge $TOPHAT_OUT/merged_accepted_hits.bam $TOPHAT_OUT2/accepted_hits.bam $TOPHAT_OUT/accepted_hits.bam"
			time samtools merge $TOPHAT_OUT/merged_accepted_hits.bam $TOPHAT_OUT2/accepted_hits.bam $TOPHAT_OUT/accepted_hits.bam
		else
			echo -e "\e[031m$TOPHAT_OUT2/accepted_hits.bam not found. tophat failed? \e[0m\n"
			exit
		fi

		if [ -s $TOPHAT_OUT/merged_accepted_hits.bam ];then
			echo -e "samtools index $TOPHAT_OUT/merged_accepted_hits.bam\n"
			time samtools index $TOPHAT_OUT/merged_accepted_hits.bam
		fi

	fi # end of RUN_TOPHAT_UNMAPPED

	if [ -s $TOPHAT_OUT/merged_accepted_hits.bam ];then 
		######################
		# 4. HTSeq
		# For paired-end data, the alignment have to be sorted either by read name or by alignment position.
		# 'sort -T $SCRATCH_OUT -k 1,1' will do a sort by read name, or
		# samtools sort -n $TOPHAT_OUT/merged_accepted_hits.bam TOPHAT_OUT/merged_accepted_hits.sorted.bam
		# output:
		######################
		if [ $RUN_HTSEQ_UNMAPPED -eq 1 ]; then
			HTSEQ_OUT=$PROJECT_DIR/HTSeq/$Barcode
			mkdir_unless $PROJECT_DIR/HTSeq
			mkdir_unless $HTSEQ_OUT

			# count based on illumina iGenome GTF
			echo -e "\nsamtools view $TOPHAT_OUT/merged_accepted_hits.bam | htseq-count --stranded=$HTSEQ_STRAND --type=$HTSEQ_TYPE --quiet --idattr=$HTSEQ_IDATTR - $GTF > $HTSEQ_OUT/$SLX.$Barcode.igenome.HTSeq.count.unmap.rescued.txt"
			time samtools view $TOPHAT_OUT/merged_accepted_hits.bam | htseq-count \
				--stranded=$HTSEQ_STRAND \
				--type=$HTSEQ_TYPE \
				--quiet \
				--idattr=$HTSEQ_IDATTR \
				- $GTF > $HTSEQ_OUT/$SLX.$Barcode.igenome.HTSeq.count.unmap.rescued.txt
		fi

		#####################
		## 5. DEXSeq Count ##
		#####################
		DEXSEQ_OUT=$PROJECT_DIR/DEXSeq/$Barcode
		if [ $RUN_DEXSEQ_CNT -eq 1 ]; then
			mkdir_unless $PROJECT_DIR/DEXSeq
			mkdir_unless $DEXSEQ_OUT
			echo -e "\npython ~/R/x86_64-unknown-linux-gnu-library/3.1/DEXSeq/python_scripts/dexseq_count.py $TOPHAT_OUT/merged_accepted_hits.bam"
			time python ~/R/x86_64-unknown-linux-gnu-library/3.1/DEXSeq/python_scripts/dexseq_count.py \
				--paired=no \
				--stranded=reverse \
				--format=bam \
				--order=pos \
				$DEX_GFF \
				$TOPHAT_OUT/merged_accepted_hits.bam \
				$DEXSEQ_OUT/$SLX.$Barcode.merged_accepted_hits.dexseq.cnt.txt 
		fi
	else
		echo -e "\e[031m$TOPHAT_OUT/merged_accepted_hits.bam not found. tophat failed? \e[0m\n"
		exit
	fi
else
	echo -e "\e[031m Paired-end mode\e[0m\n"
	######################
	## 0. Prepare Reads ##
	######################
	Unmapped_Fw_FastQ=$TOPHAT_OUT/unmapped.sorted_1.fastq.gz # paired read1
	Unmapped_Rv_FastQ=$TOPHAT_OUT/unmapped.sorted_2.fastq.gz # paired read2
	Unmapped_Un_FastQ=$TOPHAT_OUT/unmapped.sorted.fastq.gz   # un-paired 
	if [ ! -s $Unmapped_Fw_FastQ ] && [ ! -s $Unmapped_Rv_FastQ ]; then
		# sort unmapped.bam by readname
		if [ ! -s $TOPHAT_OUT/unmapped.sorted.bam ];then  
			# 2. assuming readnames of unmapped.bam are uniq 
			echo -e "samtools sort -n -@ $NT -m 60G -o $TOPHAT_OUT/unmapped.sorted.bam -O bam -T $SLX.$Barcode $TOPHAT_OUT/unmapped.bam"
			time samtools sort -n -@ $NT -m 60G -o $TOPHAT_OUT/unmapped.sorted.bam -O bam -T $SLX.$Barcode $TOPHAT_OUT/unmapped.bam
		fi
		echo -e "bam bam2FastQ --readname --gzip --in $TOPHAT_OUT/unmapped.sorted.bam &> /dev/null"
		time bam bam2FastQ --readname --gzip --in $TOPHAT_OUT/unmapped.sorted.bam &> /dev/null
	fi
	###################################################################
	# 1. TOPHAT 
	# Map the unmapped reads 
	# output: $TOPHAT_OUT2/accepted_hits.bam
	###################################################################
	if [ $RUN_TOPHAT_UNMAPPED -eq 1 ]; then
		if [ -s $Unmapped_Fw_FastQ ] && [ -s $Unmapped_Rv_FastQ ] && [ -s $RAW_JUNCS ]; then
			mkdir_unless $TOPHAT_OUT2
			echo -e "\ntophat2 -p $NT --library-type $LIB_TYPE --output-dir $TOPHAT_OUT2 --max-multihits $TH_MH --raw-juncs $RAW_JUNCS --no-novel-junc $BOWTIE2_INDEX_BASE $Unmapped_Fw_FastQ $Unmapped_Rv_FastQ"
			time tophat2 \
				--num-threads $NT \
				--library-type $LIB_TYPE \
				--output-dir $TOPHAT_OUT2 \
				--max-multihits $TH_MH \
				--raw-juncs $RAW_JUNCS \
				--no-novel-junc \
				$BOWTIE2_INDEX_BASE \
				$Unmapped_Fw_FastQ $Unmapped_Rv_FastQ
		else
			echo -e "\e[031m$Unmapped_Fw_FastQ or $Unmapped_Rv_FastQ or $RAW_JUNCS not found\e[0m\n"
			exit
		fi

		## Merge 1) inital mapped bam + 2) rescued from un-mapped bam files
		if [ -s $TOPHAT_OUT2/accepted_hits.bam ];then 
			echo -e "\nsamtools merge -@ $NT $TOPHAT_OUT/merged_accepted_hits.bam $TOPHAT_OUT2/accepted_hits.bam $TOPHAT_OUT/accepted_hits.bam"
			time samtools merge -@ $NT $TOPHAT_OUT/merged_accepted_hits.bam $TOPHAT_OUT2/accepted_hits.bam $TOPHAT_OUT/accepted_hits.bam
		else
			echo -e "\e[031m$TOPHAT_OUT2/accepted_hits.bam not found. tophat failed? \e[0m\n"
			exit
		fi

		if [ -s $TOPHAT_OUT/merged_accepted_hits.bam ];then
			echo -e "samtools index $TOPHAT_OUT/merged_accepted_hits.bam\n"
			time samtools index $TOPHAT_OUT/merged_accepted_hits.bam
		fi

	fi # end of RUN_TOPHAT_UNMAPPED

	if [ -s $TOPHAT_OUT/merged_accepted_hits.bam ];then 
		######################
		# 4. HTSeq
		# For paired-end data, the alignment have to be sorted either by read name or by alignment position.
		# 'sort -T $SCRATCH_OUT -k 1,1' will do a sort by read name, or
		# samtools sort -n $TOPHAT_OUT/merged_accepted_hits.bam TOPHAT_OUT/merged_accepted_hits.sorted.bam
		# TopHat retruns BAM sorted by position, hence "--order=pos" for htseq-count in case of pair-read 
		# output:
		######################
		if [ $RUN_HTSEQ_UNMAPPED -eq 1 ]; then
			HTSEQ_OUT=$PROJECT_DIR/HTSeq/$Barcode
			mkdir_unless $PROJECT_DIR/HTSeq
			mkdir_unless $HTSEQ_OUT
			# count based on illumina iGenome GTF
			# -f 1: pair-end only (there could be single-end not paired from the previous rescuing process of unmapped.bam)
			echo -e "\nsamtools view -f 1 $TOPHAT_OUT/merged_accepted_hits.bam | htseq-count --order=pos --stranded=$HTSEQ_STRAND --type=$HTSEQ_TYPE --quiet --idattr=$HTSEQ_IDATTR - $GTF > $HTSEQ_OUT/$SLX.$Barcode.HTSeq.count.unmap.rescued.$TR_PREFIX.$ENS_VER.txt"
			time samtools view -f 1 $TOPHAT_OUT/merged_accepted_hits.bam | htseq-count \
				--order=pos \
				--stranded=$HTSEQ_STRAND \
				--type=$HTSEQ_TYPE \
				--quiet \
				--idattr=$HTSEQ_IDATTR \
				- $GTF > $HTSEQ_OUT/$SLX.$Barcode.HTSeq.count.unmap.rescued.$TR_PREFIX.$ENS_VER.txt
		fi
		######################
		# 5. featureCount
		# featureCounts automatically sorts reads by name 
		# if paired reads are not in consecutive positions in the SAM or BAM file
		# -T: No. of core
		# -a: GTF file
		# -Q: min MQ
		# --primary: Count primary alignments only
		# -s: 0 (unstranded), 1 (stranded) and 2 (reversely stranded).
		# -p: If specified, fragments (or templates) will be counted
		# -B: Count read pairs that have both ends successfully aligned only.
		# -C: Do not count read pairs that have their two ends mapping to different chromosomes
		######################	
		if [ $RUN_FTCNT_UNMAPPED -eq 1 ]; then
			FTCNT_OUT=$PROJECT_DIR/featureCount/$Barcode
			mkdir_unless $PROJECT_DIR/featureCount
			mkdir_unless $FTCNT_OUT
			FTCNT_FILE=$FTCNT_OUT/$SLX.$Barcode.featureCount.unmap.rescued.$TR_PREFIX.$ENS_VER.txt
			#FTCNT_FILE=$FTCNT_OUT/$SLX.$Barcode.featureCount.unmap.rescued.$TR_PREFIX.$ENS_VER.known.novel.10.tr.reconstruction.txt
			echo -e "featureCounts -T $NT -a $GTF -Q $FTCNT_MQ -s $FTCNT_STRAND -p -C -o $FTCNT_FILE $TOPHAT_OUT/merged_accepted_hits.bam\n"
			time featureCounts \
				-T $NT \
				-a $GTF \
				-Q $FTCNT_MQ \
				-s $FTCNT_STRAND \
				-p -C \
				-o $FTCNT_FILE \
				$TOPHAT_OUT/merged_accepted_hits.bam
			time awk 'BEGIN{OFS="\t"}NR>2{print $1,$7}' $FTCNT_FILE | sort > ${FTCNT_FILE%.txt}.clean.txt
		fi
	else
		echo -e "\e[031m$TOPHAT_OUT/merged_accepted_hits.bam not found. tophat failed? \e[0m\n"
		exit
	fi
fi

echo -e "\e[32m$Barcode done for tophat.persample.unmapped \e[0m" 
