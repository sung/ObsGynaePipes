#!/bin/bash

#printf "Library,BarCode,HtseqFile\n"
#SLXS="SLX-9176 SLX-9791 SLX-9794" # FG
SLXS="SLX-10288 SLX-10290 SLX-11365 SLX-11366 SLX-11367" # JD
for SLX in $SLXS;
do 
	BARCODE=$(for i in `ls /home/ssg29/scratch/data/fastq/$SLX/$SLX*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f8 | cut -d'.' -f2 ; done | uniq)
	for i in $BARCODE;
	do
		MIRNA=/home/ssg29/results/$SLX.Homo_sapiens.v1/miRDeep2/core/$i/$SLX.Homo_sapiens.v1.$i.miRNA.cnt.txt
		#MIRNA=/home/ssg29/results/$SLX.Homo_sapiens.v3/miRDeep2/core/$i/$SLX.Homo_sapiens.v3.$i.miRNA.cnt.txt
		printf "$SLX,$i,$MIRNA\n"
	done
done
