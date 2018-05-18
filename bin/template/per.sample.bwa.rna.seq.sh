#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 2/May/2018
# Last modified: 2/May/2018 
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

##################################################
# 2. run initial fastqc 
# output: results/SLX123/FastQC/A001/
##################################################
if [ $RUN_FASTQC -eq 1 ]; then
    mkdir -p $PROJECT_DIR/FastQC/$Barcode
    FASTQ_FILES=`ls $HOME/data/fastq/$SLX/$SLX.$Barcode*.fq.gz | grep -v lost`
    for FastQ_file in $FASTQ_FILES	# per fastq file (per lane, both forward and reverse)
    do
        echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode"
        time fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode 
    done
fi

############################################################
# 3. Trim FastQ file per lane 	
# input: .fastq.gz
# output: results/SGA_v6/Trim/A001/SLX-8074.A001.C2N01ACXX.s_1.r_1.1_val_1.fq.gz
# output: results/SGA_v6/Trim/A001/SLX-8074.A001.C2N01ACXX.s_1.r_1.1_val_2.fq.gz
# output: results/SGA_v6/Trim/A001/SLX-8074.A001.C2N01ACXX.s_1.r_1.1.fastq.gz_trimming_report.txt
############################################################
if [ $IS_SE -eq 1 ]; then
    echo -e "\e[031m Single-end mode\e[0m\n"
    #############################
    # 1-b. Triming by cutadapt ##
    #############################
	if [ $RUN_TRIM -eq 1 ]; then
        mkdir -p $PROJECT_DIR/Trim/$Barcode
        FASTQ_FILES=`ls $HOME/data/fastq/$SLX/$SLX.$Barcode*.fq.gz | grep -v lost`
        for FastQ_file in $FASTQ_FILES	# per fastq file (per lane)
        do
            Trimmed_FastQ=${FastQ_file/data\/fastq\/$SLX/results\/$PROJECT\/Trim\/$Barcode}
            Trimmed_FastQ=${Trimmed_FastQ%.fq.gz}_trimmed.fq.gz #SLX-8546.D709_D501.C4FN0ANXX.s_6.r_1_trimmed.fq.gz 
            echo -e "\ncutadapt $FastQ_file"
            time python3.6 $HOME/.local/bin/cutadapt \
                -j $NT \
                -a $TR_ADAPTOR1 \
                -q $TR_QUAL \
                -O $TR_STRNCY \
                -m $TR_LEN \
                -u $TR_CUT_LEN \
                -o $Trimmed_FastQ \
                $FastQ_file &> ${Trimmed_FastQ%_trimmed.fq.gz}.fq.gz_trimming_report.txt
        done
    fi
    #################################################
    ## Mapping via BWA                             ##
    ## result BAM files sorted by read-name for PE #
    #################################################
    if [ $RUN_BWA -eq 1 ]; then
        mkdir -p $PROJECT_DIR/BWA/$Barcode
        FASTQ_FILES=`ls $PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode*trimmed.fq.gz`
        ## For each fastq files (distinct cell-lane)
        for FastQ_file in $FASTQ_FILES;do	# per fastq file (per lane)
            Cell=`echo $FastQ_file | cut -d/ -f8 | cut -d. -f3`
            Lane=`echo $FastQ_file | cut -d/ -f8 | cut -d. -f4`

            MY_BWA_PREFIX=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.$Lane.bwa
            # -t INT        number of threads [1]
            # -T INT        minimum score to output [30]
            # -U INT        penalty for an unpaired read pair [17]
            # -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
            echo -e "bwa mem -t $NT $BWA_INDEX_BASE.fa $FastQ_file\n"
            time bwa mem \
                -t $NT \
                -T 19 \
                -U 15 \
                -R "@RG\tID:$SLX\tPL:illumina\tPU:run\tLB:$SLX\tSM:$Barcode\tCN:CamObsGynae" \
                $BWA_INDEX_BASE.fa \
                $FastQ_file \
                1> $MY_BWA_PREFIX.sam \
                2> $MY_BWA_PREFIX.log
        done
    fi
else
    echo -e "\e[031m Pair-end mode\e[0m\n"
    #############################
    # 1-b. Triming by cutadapt ##
    #############################
	if [ $RUN_TRIM -eq 1 ]; then
        mkdir -p $PROJECT_DIR/Trim/$Barcode
        FASTQ_FILES=`ls $HOME/data/fastq/$SLX/$SLX.$Barcode*r_1.fq.gz | grep -v lost`
        ## For each fastq files (distinct cell-lane)
        for FastQ_file in $FASTQ_FILES;do	# per fastq file (per lane)
            Cell=`echo $FastQ_file | cut -d/ -f7 | cut -d. -f3`
            Lane=`echo $FastQ_file | cut -d/ -f7 | cut -d. -f4`

            Fw_FastQ=$HOME/data/fastq/$SLX/$SLX.$Barcode.$Cell.$Lane.r_1.fq.gz # SLX-8074.A001.C2N01ACXX.s_1.r_1.fq.gz
            Rv_FastQ=$HOME/data/fastq/$SLX/$SLX.$Barcode.$Cell.$Lane.r_2.fq.gz # SLX-8074.A001.C2N01ACXX.s_1.r_2.fq.gz
            if [ ! -s $Fw_FastQ ];then
                echo -e "\e[031m$Fw_FastQ not found\e[0m\n"
                exit
            fi
            if [ ! -s $Rv_FastQ ];then
                echo -e "\e[031m$Rv_FastQ not found\e[0m\n"
                exit
            fi
            Trimmed_Fw_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_1_val_1.fq.gz 
            Trimmed_Rv_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_2_val_2.fq.gz 

            echo -e "\ncutadapt $Fw_FastQ $Rv_FastQ"
            time python3.6 $HOME/.local/bin/cutadapt \
                -j $NT \
                -a $TR_ADAPTOR1 \
                -A $TR_ADAPTOR2 \
                -q $TR_QUAL \
                -O $TR_STRNCY \
                -m $TR_LEN \
                -o $Trimmed_Fw_FastQ \
                -p $Trimmed_Rv_FastQ \
                $Fw_FastQ $Rv_FastQ &> ${Trimmed_Fw_FastQ%.fq.gz}.trimming_report.txt
        done
    fi
    #################################################
    ## Mapping via BWA                             ##
    ## result BAM files sorted by read-name for PE #
    #################################################
    if [ $RUN_BWA -eq 1 ]; then
        mkdir -p $PROJECT_DIR/BWA/$Barcode
        FASTQ_FILES=`ls $PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode*r_1*.fq.gz`
        ## For each fastq files (distinct cell-lane)
        for Fwd_FastQ_file in $FASTQ_FILES;do	# per fastq file (per lane)
            Cell=`echo $Fwd_FastQ_file | cut -d/ -f8 | cut -d. -f3`
            Lane=`echo $Fwd_FastQ_file | cut -d/ -f8 | cut -d. -f4`
            Rv_FastQ_file=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_2_val_2.fq.gz 
            if [ ! -s $Rv_FastQ_file ];then
                echo -e "\e[031m$Rv_FastQ_file not found\e[0m\n"
                exit
            fi

            MY_BWA_PREFIX=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.$Lane.bwa
            # -t INT        number of threads [1]
            # -T INT        minimum score to output [30]
            # -U INT        penalty for an unpaired read pair [17]
            # -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
            echo -e "bwa mem -t $NT $BWA_INDEX_BASE.fa $Fwd_FastQ_file\n"
            time bwa mem \
                -t $NT \
                -T 19 \
                -U 15 \
                -R "@RG\tID:$SLX\tPL:illumina\tPU:run\tLB:$SLX\tSM:$Barcode\tCN:CamObsGynae" \
                $BWA_INDEX_BASE.fa \
                $Fwd_FastQ_file $Rev_FastQ_file \
                1> $MY_BWA_PREFIX.sam \
                2> $MY_BWA_PREFIX.log
        done
    fi
fi
###############
## Merge SAM ##
###############
MERGED_SAM=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.merged.bwa.sam
if [ ! -s $MERGED_SAM ];then
    dummy=`ls $PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.*.bwa.sam` # isa character
    # avoid globbing (expansion of *).
    set -f                      
    SAM_FILES=(${dummy// / }) # isa array
    HEADER=${SAM_FILES[0]}
    time (grep ^@SQ $HEADER; grep ^@RG $HEADER; for f in ${SAM_FILES[@]}; do grep ^@PG $f; done; for f in ${SAM_FILES[@]}; do grep -v @ $f; done) > $MERGED_SAM 
fi

#perl ~/Install/CIRI_v2.0.6/CIRI2.pl -T 16 -I /home/ssg29/results/SLX-9168.Homo_sapiens.SE125.v2/BWA/D701_D501/SLX-9168.D701_D501.merged.bwa.sam -F ~/data/genome/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa -A ~/data/genome/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.90.gtf -O /home/ssg29/results/SLX-9168.Homo_sapiens.SE125.v2/BWA/D701_D501/SLX-9168.D701_D501.merged.bwa.CIRI2.txt -G /home/ssg29/results/SLX-9168.Homo_sapiens.SE125.v2/BWA/D701_D501/SLX-9168.D701_D501.merged.bwa.CIRI2.log -M MT

echo -e "\e[32m$Barcode done for mapping\e[0m" 
