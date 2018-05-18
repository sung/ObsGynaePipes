#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 9/May/2018
# Last modified: 9/May/2018 
# Optimised and customised to run at the UIS CSD3 (eg. peta4 skylake) machines

source $HOME/Pipelines/config/rna_seq.config # export envrionment variables
source $HOME/Pipelines/lib/sung.sh # defines user defined functions (e.g. make_run_script)
source $HOME/config/sung.hpc.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

mkdir -p $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Barcode\e[0m\n"

########################
## Merge per-lane SAM ##
########################
MERGED_SAM=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.merged.bwa.sam
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
if [ $RUN_CIRI -eq 1 ]; then
    mkdir -p $PROJECT_DIR/CIRI2/$Barcode
    CIRI_OUT=$PROJECT_DIR/CIRI2/$Barcode/$SLX.$Barcode.merged.bwa.sam.CIRI2.txt
    CIRI_LOG=$PROJECT_DIR/CIRI2/$Barcode/$SLX.$Barcode.merged.bwa.sam.CIRI2.log
    if [ ! -s $CIRI_OUT ];then
        echo -e "perl ~/Install/CIRI_v2.0.6/CIRI2.pl -I $MERGED_SAM -O $CIRI_OUT"
        time perl ~/Install/CIRI_v2.0.6/CIRI2.pl \
            -T $NT \
            -I $MERGED_SAM \
            -F $GENOME \
            -A $GTF \
            -M MT \
            -O $CIRI_OUT > $CIRI_LOG
    else
        echo -e "\e[32mskip $SLX:$Barcode for CIRI2 (done previously)\e[0m"
    fi
fi
echo -e "\e[32m$Barcode done for CIRI\e[0m" 
