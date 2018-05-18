#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 2/May/2018
# Last modified: 2/May/2018 
# Optimised and customised to run at the UIS CSD3 (eg. peta4 skylake) machines

#ALL_SLX="SLX-9168 SLX-9169" # FG 
#ALL_SLX="SLX-9792 SLX-10281 SLX-10284 SLX-10283 SLX-10285 SLX-10287 SLX-10402" # JD

#ALL_SLX="SLX-10283" # e-5 D712_D502 D712_D503

NT=20
## ASSUME SINGLE-END ##

for SLX in $ALL_SLX;do
    #PROJECT=$SLX."Homo_sapiens.SE125.v2" # FG Placenta total RNA-Seq
    PROJECT=$SLX."Homo_sapiens.v1" # JD Placenta total RNA-Seq
    RESULT_DIR=$HOME/rcs/rcs-ssg29-obsgynae/POPS/results
    PROJECT_DIR=$RESULT_DIR/$PROJECT
    ALL_BARCODES=($(for i in `ls $HOME/rcs/rcs-ssg29-obsgynae/POPS/data/fastq/$SLX/$SLX.*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f10 | cut -d'.' -f2 ; done | uniq))
    MY_BARCODES=($(for i in ${ALL_BARCODES[*]}; do CNT=`grep -c 'CIRI finished' $PROJECT_DIR/CIRI2/$i/$SLX.$i.merged.bwa.sam.CIRI2.txt.log`; if [ $CNT -eq 0  ];then echo -e "$i"; fi; done))
    for Barcode in ${MY_BARCODES[*]};do
        ########################
        ## Merge per-lane SAM ## 
        ########################
        MERGED_SAM=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.merged.bwa.sam
        #dummy=`ls $PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode*s_*bwa.sam` # isa character
        if [ ! -s $MERGED_SAM ];then
            # avoid globbing (expansion of *).
            #set -f                      
            #SAM_FILES=(${dummy// / }) # isa array
            SAM_FILES=($(for i in `ls $PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode*s_*bwa.sam`; do echo $i; done))
            HEADER=${SAM_FILES[0]}
            echo -e "merging sam for $SLX:$Barcode"
            time (grep ^@SQ $HEADER; grep ^@RG $HEADER; for f in ${SAM_FILES[@]}; do grep ^@PG $f; done; for f in ${SAM_FILES[@]}; do grep -v @ $f; done) > $MERGED_SAM 
        fi
        ###############
        ## RUN CIRI2  #
        ###############
        mkdir -p $PROJECT_DIR/CIRI2/$Barcode
        CIRI_OUT=$PROJECT_DIR/CIRI2/$Barcode/$SLX.$Barcode.merged.bwa.sam.CIRI2.txt
        CIRI_LOG=$PROJECT_DIR/CIRI2/$Barcode/$SLX.$Barcode.merged.bwa.sam.CIRI2.log
        #if [ ! -s $CIRI_OUT ];then
            echo -e "perl ~/Install/CIRI_v2.0.6/CIRI2.pl -T $NT -I $MERGED_SAM -F ~/data/genome/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa -A ~/data/genome/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.90.gtf -M MT -O $CIRI_OUT"
            time perl ~/Install/CIRI_v2.0.6/CIRI2.pl \
                -T $NT \
                -I $MERGED_SAM \
                -F ~/data/genome/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
                -A ~/data/genome/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.90.gtf \
                -M MT \
                -O $CIRI_OUT > $CIRI_LOG
            echo -e "\e[32m$SLX:$Barcode done for CIRI2\e[0m"
        #else
        #    echo -e "\e[32mskip $SLX:$Barcode for CIRI2 (done previously)\e[0m"
        #fi
    done
done
