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

########################################
## DE-DUPLICATE SAME ENTRIES FROM BAM ## 
########################################
My_bam=$PROJECT_DIR/TopHat/$Barcode/accepted_hits.9.bam          # TopHat 1-pass x% downsampled
#My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.bam          # TopHat 2-pass
# sort by position
#if [ ! -s ${My_bam%.bam}.uniq.sorted.bam ]; then
#	echo -e "samtools view $My_bam | awk | sort -u | smtools view -bT $GENOME - > ${My_bam%.bam}.uniq.sorted.bam"
#	time samtools view $My_bam \
#		| awk '{print $0}' | sort -k 3,3 -k 4,4n -S 60G -T $PROJECT_DIR/scratch/$Barcode -u \
#		| samtools view -bT $GENOME - \
#		| samtools sort -@ $NT -m 60G -o ${My_bam%.bam}.uniq.sorted.bam -O bam -T $SLX.$Barcode -
#		#| samtools view -bT $GENOME - > ${My_bam%.bam}.uniq.sorted.bam
#fi

####################
## MARK DUPLICATE ##
####################
#Note: The input file must be sorted by coordinate.
#My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.sorted.bam
#merged_accepted_hits.uniq.sorted.dedup.bam
if [ ! -s ${My_bam%.bam}.dedup.bam ]; then
	echo -e "bam dedup --in $My_bam --out ${My_bam%.bam}.dedup.bam --force --verbose --log"
	time bam dedup \
		--in $My_bam \
		--out ${My_bam%.bam}.dedup.bam \
		--force --verbose --log
fi

#############
## dupRadar #
#############
#mkdir -p $PROJECT_DIR/dupRadar/$Barcode
#echo -e "Rscript /home/ssg29/Pipelines/bin/template/per.sample.dupRadar.R $SLX $Barcode ${My_bam%.bam}.dedup.bam $NT\n"
#time Rscript /home/ssg29/Pipelines/bin/template/per.sample.dupRadar.R $SLX $Barcode ${My_bam%.bam}.dedup.bam $NT

echo -e "\e[32m$Barcode done for per.sample.bamUtil\e[0m"
