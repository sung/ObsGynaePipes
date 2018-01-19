#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 21/Oct/2015
# Last modified: 21/Oct/2015
# Optimised and customised to run at the Darwin HPC
# based on StringTie <https://github.com/gpertea/stringtie>

source $HOME/lib/sung.sh #defines 'mkdir_unless'

SPECIES='Homo_sapiens' # 'Homo_sapiens' or 'Mus_musculus' 
PROJECT="$SLX.$SPECIES.SE125.v2"
PROJECT_DIR="/scratch/ssg29/results/$PROJECT"
TR_PREFIX=GRCh38 # GRCH37.gff from /scratch/genome/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf ($GTF)
ENS_VER=82 # 75 (for igenome), 82 (for manual downlaod)
GTF=$HOME/data/genome/$SPECIES/Ensembl/$TR_PREFIX/Annotation/Genes/genes.gtf

mkdir_unless $PROJECT_DIR/StringTie
mkdir_unless $PROJECT_DIR/StringTie/$Barcode

My_bam=$PROJECT_DIR/TopHat/$Barcode/final_merged_accepted_hits.bam # SE125 + SE50
if [ -s $My_bam ];then 
	echo -e "stringtie $My_bam\n"
	time stringtie \
		$My_bam \
		-p 6 \
		-G $GTF \
		-l $SLX.$Barcode.STRG \
		-c 10 \
		-j 5 \
		-B \
		-o $PROJECT_DIR/StringTie/$Barcode/$PROJECT.$Barcode.stringtie.$TR_PREFIX.$ENS_VER.assembly.gtf \
		-C $PROJECT_DIR/StringTie/$Barcode/$PROJECT.$Barcode.stringtie.$TR_PREFIX.$ENS_VER.coverage.txt \
		-A $PROJECT_DIR/StringTie/$Barcode/$PROJECT.$Barcode.stringtie.$TR_PREFIX.$ENS_VER.gene.cnt.txt
else
	echo -e "\e[031m$My_bam not found. tophat failed? \e[0m\n"
	exit
fi
echo -e "\e[32m$Barcode done for StringTie\e[0m" 
