#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 30/Jan/2015
# Last modified: 19/Mar/2018 
# Optimised and customised to run at the UIS CSD3 (eg. peta4 skylake) machines

source $HOME/Pipelines/config/target.seq.sbs.config # calls global.config, slurm.config, samples.config 
source $HOME/Pipelines/lib/sung.sh # defines user defined functions (e.g. make_run_script)
source $HOME/config/sung.hpc.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

mkdir -p $PROJECT_DIR	

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Barcode\e[0m\n"

##################################################
# 2. run initial fastqc 
# output: results/SLX123/FastQC/A001/
##################################################
if [ $RUN_FASTQC -eq 1 ]; then
    mkdir -p $PROJECT_DIR/FastQC
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
if [ $IS_PE -eq 1 ]; then
    Trimmed_Fw_FastQ_array=() #initialise
    Trimmed_Rv_FastQ_array=() #initialise

    FASTQ_FILES=`ls $HOME/data/fastq/$SLX/$SLX.$Barcode*r_1.fq.gz | grep -v lost`
    Cell_lane_array=`for i in $FASTQ_FILES; do echo $i | cut -d/ -f7 | cut -d. -f3,4 ; done`
    ## For each lane 
    for CELL_LANE in ${Cell_lane_array[*]}
    do
        Cell=`echo $CELL_LANE | cut -d. -f1`
        Lane=`echo $CELL_LANE | cut -d. -f2`
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
        #############################
        # 1-b. Triming by cutadapt ##
        #############################
        if [ $RUN_TRIM -eq 1 ]; then
            mkdir -p $PROJECT_DIR/Trim
            mkdir -p $PROJECT_DIR/Trim/$Barcode

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
        fi

        # scrap trimmed reads
        if [ -s $Trimmed_Fw_FastQ ];then
            Trimmed_Fw_FastQ_array+=($Trimmed_Fw_FastQ) #push it
        else
            echo -e "\e[031m$Trimmed_Fw_FastQ not found\e[0m\n"
            exit
        fi

        if [ -s $Trimmed_Rv_FastQ ];then
            Trimmed_Rv_FastQ_array+=($Trimmed_Rv_FastQ) #push it
        else
            echo -e "\e[031m$Trimmed_Rv_FastQ not found\e[0m\n"
            exit
        fi
    done # endof CELL_LANE

    Merged_Fw_FastQ=$(printf ",%s" "${Trimmed_Fw_FastQ_array[@]}")
    Merged_Fw_FastQ=${Merged_Fw_FastQ:1} # remove first comma 

    Merged_Rv_FastQ=$(printf ",%s" "${Trimmed_Rv_FastQ_array[@]}")
    Merged_Rv_FastQ=${Merged_Rv_FastQ:1} # remove first comma 

	TRIMMED_FASTQ_FILES="$Merged_Fw_FastQ $Merged_Rv_FastQ"
else
    echo -e "\e[031m Single-end mode\e[0m\n"
    Trimmed_FastQ_array=() #initialise
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
        if [ ! -s $Trimmed_FastQ ];then
            echo -e "\e[031m$Trimmed_FastQ not found\e[0m\n"
            exit
        else
            echo -e "\nfastqc $Trimmed_FastQ -o $PROJECT_DIR/Trim/$Barcode"
            #time fastqc $Trimmed_FastQ -o $PROJECT_DIR/Trim/$Barcode 
        fi
    done

    Merged_FastQ=$(printf ",%s" "${Trimmed_FastQ_array[@]}")
    Merged_FastQ=${Merged_FastQ:1} # remove first comma 
	TRIMMED_FASTQ_FILES=$Merged_FastQ
fi

#################################################
## Mapping via BWA                             ##
## result BAM files sorted by read-name for PE #
#################################################
MY_BWA_PREFIX=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa
if [ $RUN_BWA -eq 1 ]; then
    mkdir -p $PROJECT_DIR/BWA
    mkdir -p $PROJECT_DIR/BWA/$Barcode
    # -t INT        number of threads [1]
    # -M            mark shorter split hits as secondary
    # -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
    echo -e "bwa mem -M $TRIMMED_FASTQ_FILES | samblaster | samtools \n"
    time bwa mem \
        -M \
        -t $NT \
        -R "@RG\tID:$SLX\tPL:illumina\tPU:run\tLB:$SLX\tSM:$Barcode\tCN:CamObsGynae" \
        $BWA_INDEX_BASE.fa \
        $TRIMMED_FASTQ_FILES 2> $MY_BWA_PREFIX.log | \
        time $HOME/Install/samblaster/samblaster -M 2> $MY_BWA_PREFIX.samblaster.log | \
        time samtools view -Sb - > $MY_BWA_PREFIX.samblaster.bam
fi

echo -e "\e[32m$Barcode done for mapping\e[0m" 
