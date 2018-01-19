#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>

TR_PREFIX=GRCh37
SPECIES='Homo_sapiens' # 'Homo_sapiens' or 'Mus_musculus' 

#ALL_SLX="SLX-9168 SLX-9169"
#VERSION=$SPECIES."SE125.v1"
ALL_SLX="SLX-10281 SLX-10402 SLX-9792 SLX-10284 SLX-10285 SLX-10283 SLX-10287"
VERSION=$SPECIES."v2"

FILE="/home/ssg29/results/RNA-Seq/Placentome/Meta/FG.JD.SE125."$TR_PREFIX".mapping.txt"
#echo "SLX Barcode Category Read" > $FILE 
#touch $FILE 

for SLX in $ALL_SLX; do
	PROJECT=$SLX.$VERSION
	PROJECT_DIR=~/results/$RESULT_DIR/$PROJECT

	BARCODE=$(for i in `ls /home/ssg29/scratch/data/fastq/$SLX/$SLX*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f8 | cut -d'.' -f2 ; done | uniq)

	for j in $BARCODE; do printf "$SLX $j Mapped "; awk '/Mapped/{print $3}' $PROJECT_DIR/TopHat/$j/align_summary.txt; done >> $FILE 
	for j in $BARCODE; do printf "$SLX $j Rescued "; awk '/Mapped/{print $3}' $PROJECT_DIR/TopHat/$j/unmapped/align_summary.txt; done >> $FILE 
	for j in $BARCODE; do printf "$SLX $j Unmapped "; awk '/Input/{input=$3};/Mapped/{mapped=$3}END{print input-mapped}' $PROJECT_DIR/TopHat/$j/unmapped/align_summary.txt; done >> $FILE 
done

echo -e "\e[32mALL done for $PROJECT\e[0m" 
