#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 11/Apr/2014
# Last modified: 21/Jul/2014 
# Optimised and customised to run at the Darwin HPC

source $HOME/RNA-Seq/config/rna_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Barcode\e[0m\n"

	#####################################################
	# 5. cufflinks: Assemble transcripts for each sample #
	# The SAM file supplied to Cufflinks must be sorted by reference position. 
	# If you aligned your reads with TopHat, your alignments will be properly sorted already
	# o/w, sort -k 3,3 -k 4,4n hits.sam > hits.sam.sorted
	# computation time? real    352m42.053s (for R4)
	# output: $PROJECT_DIR/Cufflinks/$Barcode/transcripts.gtf 
	#####################################################
	#My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.bam # SE125
	My_bam=$PROJECT_DIR/TopHat/$Barcode/final_merged_accepted_hits.bam # SE125 + SE50
	if [ $RUN_CUFFLINK -eq 1 ]; then
		if [ -s $My_bam ]; then
			mkdir_unless $PROJECT_DIR/Cufflinks
			mkdir_unless $PROJECT_DIR/Cufflinks/$Barcode
			# --GTF: reference transcript only (http://cufflinks.cbcb.umd.edu/manual.html)
			# --GTF-guide: reference transcript + novel (http://cufflinks.cbcb.umd.edu/howitworks#hrga)
			# http://cole-trapnell-lab.github.io/cufflinks/cufflinks/#advanced-reference-annotation-based-transcript-rabt-assembly-options
			echo -e "cufflinks $My_bam\n"
			time cufflinks \
				--num-threads $NT \
				--library-type $LIB_TYPE \
				--frag-len-mean $CUFF_FRAG_LEN_MEAN \
				--frag-len-std-dev $CUFF_FRAG_LEN_SD \
				--GTF-guide $GTF \
				--frag-bias-correct $GENOME \
				--label $SLX.$Barcode.CUFF \
				--multi-read-correct \
				--output-dir $PROJECT_DIR/Cufflinks/$Barcode \
				$My_bam
				#--library-norm-method classic-fpkm # default

			if [ ! -s $PROJECT_DIR/Cufflinks/$Barcode/transcripts.gtf ]; then
				echo -e "\e[031m$CUFFLINK_OUT/transcripts.gtf not found. cufflinks failed? \e[0m\n"
				exit
			fi
		else
			echo -e "\e[031m$My_bam not found\e[0m\n"
			exit
		fi
	fi # end of RUN_CUFFLINK
