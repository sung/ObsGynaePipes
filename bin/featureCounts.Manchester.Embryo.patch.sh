#!/bin/bash

# Read Count at Exon Level of CSMD1 from RoadMap RNA-Seq data (BED format)

REF="GRCh37" # GRCh38
ENS=75 #82
ENS_VER=$REF.$ENS
NT=4
MY_GENE="CSMD1"
# see bin/fetch.exon.union.sh
#/home/ssg29/data/genome/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/Homo_sapiens.GRCh37.75.CSMD1.exon.union.sorted.ucsc.gtf
GFF=$HOME/data/genome/Homo_sapiens/Ensembl/$REF/Annotation/Genes/Homo_sapiens.$ENS_VER.$MY_GENE.exon.union.sorted.ucsc.gtf 

#for i in `ls ~/data/RoadMap/Consortium/RNA-Seq/Placenta/GSM*.bed.gz`; do
#for i in `ls ~/data/RoadMap/Consortium/RNA-Seq/Fetal/GSM*.bed.gz`; do
for i in `ls ~/data/ManchesterEmryro/RNA-Seq/*hg19_.unique.bed12.gz`; do
	echo $i; 
	time bedtools coverage -a $i -b $GFF -s -counts | \
		sort -k1,1 -k4,4n -k5,5n | \
		awk 'BEGIN{FS="\t";OFS="\t"}{split($9,a,";"); for(j in a) if(a[j]~/exon_id/){split(a[j],c," ")}; print c[2],$10}' | \
		sed 's/"//g' > ${i%.bed12.gz}.$MY_GENE.exon_id.cnt.txt; 
done
echo -e "ALL done\n"
