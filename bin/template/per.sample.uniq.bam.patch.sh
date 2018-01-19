#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 4/May/2017
# Last modified: 4/May/2017
# Optimised and customised to run at the Darwin HPC

source $HOME/RNA-Seq/config/rna_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

mkdir_unless $PROJECT_DIR	
mkdir_unless $PROJECT_DIR/scratch
mkdir_unless $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Cell.$Barcode\e[0m\n"

#My_bam=$PROJECT_DIR/TopHat/$Barcode/accepted_hits.bam                # TopHat 1-pass
My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.bam          # TopHat 2-pass

if [ ! -s ${My_bam%.bam}.uniq.bam ]; then
	echo -e "samtools view $My_bam | sort | uniq"
	time samtools view $My_bam \
		| awk '{print $0}' | sort -T $PROJECT_DIR/scratch/$Barcode \
		| uniq | sort -T $PROJECT_DIR/scratch/$Barcode -k3,3 -k4,4n | samtools view -bT $GENOME - > ${My_bam%.bam}.uniq.bam
fi

My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.uniq.bam          # TopHat 2-pass (uniq entry only)
if [ $RUN_BU_MD -eq 1 ]; then
	echo -e "bam dedup --in $My_bam --out ${My_bam%.bam}.dedup.bam"
	time bam dedup \
		--in $My_bam \
		--out ${My_bam%.bam}.dedup.bam \
		--force --verbose --log
fi

echo -e "\e[32m$Barcode done for per.sample.bamUtil\e[0m"
