#!/bin/bash
# Author: Sung Gong <sung@bio.cc>
# Last modified in 7/Mar/2018

SLX=SLX-ILM-Plasma2017
REF="GRCh37" # illumina mapped on hg19 version 
ENS=82 #82
ENS_VER=$REF.$ENS
NT=16
GTF=$HOME/data/genome/Homo_sapiens/Ensembl/$REF/Annotation/Genes/Homo_sapiens.$ENS_VER.UCSC.gtf # Homo_sapiens.GRCh37.82.UCSC.gtf 

FTCNT_DIR=$HOME/results/$SLX.Homo_sapiens.v1/featureCount
BARCODE=$(for i in `ls $HOME/data/fastq/$SLX/$SLX*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f7 | cut -d'.' -f2 ; done | uniq)
for i in $BARCODE
do 
    echo $SLX"."$i
    mkdir -p $FTCNT_DIR/$i
    MY_CNT_FILE=$FTCNT_DIR/$i/$SLX.$i.featureCount.$ENS_VER.UCSC.txt
    MY_BAM=$HOME/rcs/rcs-ssg29-obsgynae/POPS/results/Illumina/Plasma-RNA-2017/Bam/$i"_accepted_hits.bam"

    # -T: No. of core
    # -a: GTF file
    # -s: 2 (reversely stranded)
    # -f: Perform read counting at feature level (eg. counting reads for exons rather than genes).
    # -O: Assign reads to all their overlapping meta-features (or features if -f is specified).
    if [ -s $MY_BAM ]; then
        echo -e "featureCounts -o $MY_CNT_FILE $MY_BAM\n"
        time featureCounts \
            -T $NT \
            -a $GTF \
            -Q 10 \
            -s 2 \
            -p \
            -C \
            -o $MY_CNT_FILE \
            $MY_BAM
        time awk 'BEGIN{OFS="\t"}NR>2{print $1,$7}' $MY_CNT_FILE | sort > ${MY_CNT_FILE%.txt}.clean.txt
    else
        echo -e "\e[031m$MY_BAM not found\e[0m\n"
        exit
    fi
done
