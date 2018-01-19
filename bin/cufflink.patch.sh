#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 28/Oct/2015
# Last modified: 28/Oct/2015
# Optimised and customised to run at the Darwin HPC

source $HOME/lib/sung.sh #defines 'mkdir_unless'

SPECIES='Homo_sapiens' # 'Homo_sapiens' or 'Mus_musculus' 
PROJECT="$SLX.$SPECIES.SE125.v2"
PROJECT_DIR="/scratch/ssg29/results/$PROJECT"
TR_PREFIX=GRCh38 # GRCH37.gff from /scratch/genome/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf ($GTF)
ENS_VER=82 # 75 (for igenome), 82 (for manual downlaod)
GENOME=$HOME/data/genome/$SPECIES/Ensembl/$TR_PREFIX/Sequence/WholeGenomeFasta/genome.fa
GTF=$HOME/data/genome/$SPECIES/Ensembl/$TR_PREFIX/Annotation/Genes/genes.gtf

mkdir_unless $PROJECT_DIR/Cufflinks
mkdir_unless $PROJECT_DIR/Cufflinks/$Barcode

My_bam=$PROJECT_DIR/TopHat/$Barcode/final_merged_accepted_hits.bam # SE125 + SE50
if [ -s $My_bam ];then 
	echo -e "cufflinks $My_bam\n"
	time cufflinks \
		--num-threads 6 \
		--library-type fr-firststrand \
		--GTF-guide $GTF \
		--frag-bias-correct $GENOME \
		--label $SLX.$Barcode.CUFF \
		--multi-read-correct \
		--output-dir $PROJECT_DIR/Cufflinks/$Barcode \
		$My_bam
		#--library-norm-method classic-fpkm # default
else
	echo -e "\e[031m$My_bam not found. tophat failed? \e[0m\n"
	exit
fi
echo -e "\e[32m$Barcode done for Cufflinks\e[0m" 
