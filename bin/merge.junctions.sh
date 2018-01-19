#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 15/Sep/2015
# Last modified: 15/Sep/2015
# Optimised and customised to run at the Darwin HPC

source $HOME/RNA-Seq/config/rna_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.bash #defines PATH 

PROJECT=SGA.AGA.total.SE125
PROJECT_DIR=$HOME/results/RNA-Seq/$PROJECT
BARCODE_FILE=$PROJECT_DIR/Meta/FG.barcode.txt
COHORT=ALL

RAW_BEDS=$PROJECT_DIR/Junctions/FG.SE125.SE50.$COHORT.merged.junctions.bed
RAW_JUNCS=$PROJECT_DIR/Junctions/FG.SE125.SE50.$COHORT.merged.junctions.junc

echo -e "merging juction file $RAW_JUNCS\n"

#time cat $PROJECT_DIR/TopHat/*/junctions.bed > $RAW_BEDS
if [ $COHORT == 'AGA' ];then
	time cat `awk 'BEGIN{FS=",";OFS=" "}$3==0{a="/home/ssg29/results/"$1".Homo_sapiens.SE125.v2/TopHat/"$2"/junctions.bed"; b="/home/ssg29/results/"$1".Homo_sapiens.v5/TopHat/"$2"/junctions.bed"; print a,b}' $BARCODE_FILE` > $RAW_BEDS
elif [ $COHORT == 'SGA' ];then
	time cat `awk 'BEGIN{FS=",";OFS=" "}$3==1{a="/home/ssg29/results/"$1".Homo_sapiens.SE125.v2/TopHat/"$2"/junctions.bed"; b="/home/ssg29/results/"$1".Homo_sapiens.v5/TopHat/"$2"/junctions.bed"; print a,b}' $BARCODE_FILE` > $RAW_BEDS
elif [ $COHORT == 'ALL' ];then
	time cat `awk 'BEGIN{FS=",";OFS=" "}{a="/home/ssg29/results/"$1".Homo_sapiens.SE125.v2/TopHat/"$2"/junctions.bed"; b="/home/ssg29/results/"$1".Homo_sapiens.v5/TopHat/"$2"/junctions.bed"; print a,b}' $BARCODE_FILE` > $RAW_BEDS
else
	echo -e "\e[031m$COHORT not supported\e[0m\n"
	exit
fi

time bed_to_juncs < $RAW_BEDS | sort -k 1,4 -u | sort -k 1,1 > $RAW_JUNCS

echo -e "ALL done\n"
