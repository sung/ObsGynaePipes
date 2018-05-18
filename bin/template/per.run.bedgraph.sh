#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 10/Jun/2016
# Last modified: 21/Jun/2016
# Optimised and customised to run at the Darwin HPC

source $HOME/Pipelines/config/rna_seq.config # export envrionment variables
source $HOME/Pipelines/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.hpc.bash #defines PATH 

PROJECT=Placentome #SGA.AGA.total.SE125
PROJECT_DIR=$HOME/results/RNA-Seq/$PROJECT

COHORT="MY_COHORT" # CTLF, CTLM, POPS, CTL, SGA, PET
echo -e "COHORT=$COHORT"

mkdir_unless $PROJECT_DIR/BedGraph
mkdir_unless $PROJECT_DIR/BedGraph/$COHORT

####################
## Heahtly Female ##
####################
#awk 'BEGIN{FS=","}NR>1 && $6==0 && $1!="84C" && $1!="88" && $1!="80C" && $1!="93C" && $1!="94C" && $7=="F"{split($2,a,"/"); myFile="/home/ssg29/results/"a[5]"/Coverage/"a[7]"/"$3"."a[7]".genomecov.split.unmap.rescued.bedgraph"; print myFile}' ~/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/Meta/meta.ALL.csv > ~/results/RNA-Seq/Placentome/Meta/CTLF.GRCh37.75.bedgraph.list
##################
## Heahtly Male ##
##################
#awk 'BEGIN{FS=","}NR>1 && $6==0 && $1!="84C" && $1!="88" && $1!="80C" && $1!="93C" && $1!="94C" && $7=="M"{split($2,a,"/"); myFile="/home/ssg29/results/"a[5]"/Coverage/"a[7]"/"$3"."a[7]".genomecov.split.unmap.rescued.bedgraph"; print myFile}' ~/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/Meta/meta.ALL.csv > ~/results/RNA-Seq/Placentome/Meta/CTLM.GRCh37.75.bedgraph.list

if [ -s $PROJECT_DIR/Meta/$COHORT.$TR_PREFIX.$ENS_VER.bedgraph.list ];then
	LIST=`awk 'BEGIN{ORS=" "}{print;}' $PROJECT_DIR/Meta/$COHORT.$TR_PREFIX.$ENS_VER.bedgraph.list`
else
	echo -e "\e[031$PROJECT_DIR/Meta/$COHORT.$TR_PREFIX.$ENS_VER.bedgraph.list not found\e[0m\n"
	exit
fi

echo -e "bedtools unionbedg -i $LIST"
time bedtools unionbedg -i $LIST | awk '{avg=0;sum=0;for (i=4;i<=NF;i++) sum+=$i; avg=sum/(NF-3); printf "%s\t%s\t%s\t%.1f\n", $1,$2,$3,avg}' > $PROJECT_DIR/BedGraph/$COHORT/$COHORT.$TR_PREFIX.$ENS_VER.bedgraph

echo -e "bedGraphToBigWig $PROJECT_DIR/BedGraph/$COHORT/$COHORT.$TR_PREFIX.$ENS_VER.bedgraph"
time bedGraphToBigWig $PROJECT_DIR/BedGraph/$COHORT/$COHORT.$TR_PREFIX.$ENS_VER.bedgraph $PROJECT_DIR/BedGraph/$TR_PREFIX.samtool.idxstats.chr.size.txt $PROJECT_DIR/BedGraph/$COHORT/$COHORT.$TR_PREFIX.$ENS_VER.bw 

echo -e "gzip -f -9 --rsyncable $PROJECT_DIR/BedGraph/$COHORT/$COHORT.$TR_PREFIX.$ENS_VER.bedgraph"
time gzip -f -9 --rsyncable $PROJECT_DIR/BedGraph/$COHORT/$COHORT.$TR_PREFIX.$ENS_VER.bedgraph

echo -e "bedtools unionbedg done for $COHORT"
