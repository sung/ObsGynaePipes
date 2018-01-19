#!/bin/bash
COHORT="CTL"

echo -e "COHORT=$COHORT"

mkdir -p  ~/results/RNA-Seq/Placentome/BedGraph/$COHORT
LIST=`awk 'BEGIN{ORS=" "}{print;}' /scratch/ssg29/results/RNA-Seq/Placentome/Meta/$COHORT.GRCh38.82.bedgraph.list`

echo -e "bedtools unionbedg -i $LIST"
time bedtools unionbedg -i $LIST | awk '{avg=0;sum=0;for (i=4;i<=NF;i++) sum+=$i; avg=sum/(NF-3); printf "%s\t%s\t%s\t%.1f\n", $1,$2,$3,avg}' > ~/results/RNA-Seq/Placentome/BedGraph/$COHORT/$COHORT.GRCh38.82.bedgraph

echo -e "bedGraphToBigWig ~/results/RNA-Seq/Placentome/BedGraph/$COHORT/$COHORT.GRCh38.82.bedgraph"
time bedGraphToBigWig ~/results/RNA-Seq/Placentome/BedGraph/$COHORT/$COHORT.GRCh38.82.bedgraph ~/results/RNA-Seq/Placentome/BedGraph/GRCh38.samtool.idxstats.chr.size.txt ~/results/RNA-Seq/Placentome/BedGraph/$COHORT/$COHORT.GRCh38.82.bw 

echo -e "gzip -9 ~/results/RNA-Seq/Placentome/BedGraph/$COHORT/$COHORT.GRCh38.82.bedgraph"
time gzip -9 ~/results/RNA-Seq/Placentome/BedGraph/$COHORT/$COHORT.GRCh38.82.bedgraph
