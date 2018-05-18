#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 2/May/2018
# Last modified: 2/May/2018 
# Optimised and customised to run at the UIS CSD3 (eg. peta4 skylake) machines

source $HOME/Pipelines/config/rna_seq.config # export envrionment variables

#ALL_SLX="SLX-9168 SLX-9169" # 9168 e-10, 9169: e-13
#ALL_SLX="SLX-9169" # 9168 e-10, 9169: e-13
#ALL_SLX="SLX-10281 SLX-10402 SLX-9792" # login-e-11
#ALL_SLX="SLX-10284 SLX-10285 SLX-10283 SLX-10287" # login-e-9
NT=16
## ASSUME SINGLE-END ##

for SLX in $ALL_SLX;do
    #PROJECT=$SLX."Homo_sapiens.SE125.v2" # FG Placenta total RNA-Seq
    PROJECT=$SLX."Homo_sapiens.v1" # JD Placenta total RNA-Seq
    RESULT_DIR=$HOME/rcs/rcs-ssg29-obsgynae/POPS/results
    PROJECT_DIR=$RESULT_DIR/$PROJECT
    ALL_BARCODES=($(for i in `ls $HOME/rcs/rcs-ssg29-obsgynae/POPS/data/fastq/$SLX/$SLX.*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f10 | cut -d'.' -f2 ; done | uniq))
    for Barcode in ${ALL_BARCODES[*]};do
        mkdir -p $PROJECT_DIR/BWA/$Barcode
        FASTQ_FILES=`ls $PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode*trimmed.fq.gz`
        #################################################
        ## Mapping via BWA                             ##
        ## For each fastq files (distinct cell-lane)   ##
        #################################################
        for FastQ_file in $FASTQ_FILES;do	# per fastq file (per lane)
            Cell=`echo $FastQ_file | cut -d/ -f11 | cut -d. -f3`
            Lane=`echo $FastQ_file | cut -d/ -f11 | cut -d. -f4`
            MY_BWA_PREFIX=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.$Lane.bwa
            #FLAG=`grep -c -v ^@ $MY_BWA_PREFIX.sam`
            FLAG=$(stat -c%s "$MY_BWA_PREFIX.sam") # file size
            if [ $FLAG -lt 1200 ];then
                # -t INT        number of threads [1]
                # -T INT        minimum score to output [30]
                # -U INT        penalty for an unpaired read pair [17]
                # -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
                echo -e "bwa mem -t $NT $BWA_INDEX_BASE.fa $FastQ_file > $MY_BWA_PREFIX.sam\n"
                time bwa mem \
                    -t $NT \
                    -T 19 \
                    -U 15 \
                    -R "@RG\tID:$SLX\tPL:illumina\tPU:run\tLB:$SLX\tSM:$Barcode\tCN:CamObsGynae" \
                    $BWA_INDEX_BASE.fa \
                    $FastQ_file \
                    1> $MY_BWA_PREFIX.sam \
                    2> $MY_BWA_PREFIX.log
            else
                echo -e "\e[32mskip $SLX:$Barcode:$Cell:$Lane (done previously)\e[0m" 
            fi
        done
        echo -e "\e[32m$SLX:$Barcode done for bwa-mem\e[0m" 
    done
done
