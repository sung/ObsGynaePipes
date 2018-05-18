#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# Last modified: 22/Mar/2018
# First created: 30/Jan/2015
# Based on GATK4

source $HOME/Pipelines/config/target.seq.sbs.config # calls global.config, slurm.config, samples.config 
source $HOME/lib/sung.sh
source $HOME/config/sung.hpc.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

mkdir -p $PROJECT_DIR	
mkdir -p $PROJECT_DIR/GATK4
mkdir -p $PROJECT_DIR/GATK4/$Barcode
mkdir -p $PROJECT_DIR/scratch
mkdir -p $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Cell.$Barcode\e[0m\n"

#[todo] is this step necessary?
# bwa makes n-sorted BAM files
#######################
## 0. Sort BAM files ##
## by coordinate     ##
#######################
# SLX-7634.HALO1.000000000-A5NK2.bwa.bam
# SLX-Xten2017.23454_6.XXX.bwa.samblaster.bam
if [ $RUN_SORT_BAM -eq 1 ]; then
    My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.samblaster.bam 
    if [ -s $My_bam ];then
        if [ ! -s ${My_bam%.bam}.sorted.bam ];then
            echo -e "samtools sort -@ $NT -m $SAMTOOLS_MEM -o ${My_bam%.bam}.sorted.bam -T $PROJECT_DIR/scratch/$Barcode $My_bam\n"
            time samtools sort -@ $NT -m $SAMTOOLS_MEM -o ${My_bam%.bam}.sorted.bam -T $PROJECT_DIR/scratch/$Barcode $My_bam
        fi

        if [ -s ${My_bam%.bam}.sorted.bam ];then
            # remove unsorted bam file
            echo -e "rm -f $My_bam\n"
            #time rm -f $My_bam
        else
            echo -e "\e[031m${My_bam%.bam}.sorted.bam is expected but not found\e[0m\n"
            exit
        fi
    else
        #SLX-7634.HALO1.000000000-A5NK2.bwa.sorted.bam
        # merge per-lane BAM files
        if [ ! -s ${My_bam%.bam}.sorted.bam ];then
            #SLX-7634.HALO1.000000000-A5NK2.s_1.bwa.bam
            #SLX-Xten2017.23454_6.XXX.bwa.samblaster.bam
            My_bam_array=`ls $PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.s_*.bwa.bam` # per-lane bam files
            Merged_bam_files=$(printf " %s" "${My_bam_array[@]}")
            #My_bam_array is not null
            if [ -n $My_bam_array ]; then
                # two or more bam files
                if [ ${#My_bam_array[@]} -gt 1 ]; then
                    # merge bam files
                    echo -e "samtools merge $My_bam $Merged_bam_files"
                    time samtools merge $My_bam $Merged_bam_files 
                else
                    echo -e "cp $Merged_bam_files $My_bam\n"
                    time cp $Merged_bam_files $My_bam
                fi

                if [ -s $My_bam ];then
                    echo -e "samtools sort $My_bam ${My_bam%.bam}.sorted\n"
                    time samtools sort -@ $NT -m 190G -o ${My_bam%.bam}.sorted.bam -T $SLX.$Barcode $My_bam

                    if [ -s ${My_bam%.bam}.sorted.bam ];then
                        # remove unsorted bam file
                        echo -e "rm -f $My_bam\n"
                        #time rm -f $My_bam
                    fi
                else
                    echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
                    exit
                fi
            fi # end of My_bam_array is not null
        else
            echo -e "\e[032m${My_bam%.bam}.sorted.bam found\e[0m\n"
        fi
    fi # end of if [ -s $My_bam ];then

    #SLX-7634.HALO1.000000000-A5NK2.bwa.sorted.bam
    #SLX-Xten2017.23454_6.XXX.bwa.samblaster.sorted.bam
    My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.samblaster.sorted.bam 
    if [ ! -s $My_bam ];then
        echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
        exit
    fi
    #SLX-7634.HALO1.000000000-A5NK2.bwa.sorted.bam.bai
    if [ ! -s $My_bam".bai" ];then
        echo -e "samtools index -@ $NT $My_bam\n"
        time samtools index -@ $NT $My_bam
    fi
fi

##################################################################################
## 1. Indel Realignment                                                         ##
## Note that indel realignment is no longer necessary for variant discovery     ##
## if you plan to use a variant caller that performs a haplotype assembly step, ##
## such as HaplotypeCaller or MuTect2. However it is still required when using  ##
## legacy callers such as UnifiedGenotyper or the original MuTect.              ##
##################################################################################
#https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php

###################################
## 2. Base-quality Recalibration ##
###################################
if [ $RUN_GATK_BASE_RECAL -eq 1 ]; then
    #SLX-7634.HALO1.000000000-A5NK2.bwa.sorted.bam
    #SLX-Xten2017.23454_6.XXX.bwa.samblaster.sorted.bam
    My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.samblaster.sorted.bam 
    if [ ! -s $My_bam ];then
        echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
        exit
    fi
    ###BaseRecalibrator
    echo -e "$GATK_HOME/gatk BQSRPipelineSpark -I $My_bam -O ${My_bam%.bam}.recal.bam 2> $PROJECT_DIR/GATK4/$Barcode/$SLX.$Barcode.$Cell.BQSRPipelineSpark.log\n"
    time $GATK_HOME/gatk --java-options "-Xmx$GATK_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode" \
        BQSRPipelineSpark \
        --spark-master local[$NT] \
        -R $GENOME.2bit \
        -I $My_bam \
        --known-sites $DBSNP \
        --known-sites $GATK_MILLS \
        -O ${My_bam%.bam}.recal.bam 2> $PROJECT_DIR/GATK4/$Barcode/$SLX.$Barcode.$Cell.BQSRPipelineSpark.log

	## Index bam 
	if [ -s ${My_bam%.bam}.recal.bam ]; then
        echo -e "samtools index -@ $NT ${My_bam%.bam}.recal.bam\n"
        time samtools index -@ $NT ${My_bam%.bam}.recal.bam 
	else
		echo -e "\e[031m${My_bam%.bam}.recal.bam not found\e[0m\n"
		exit
	fi
fi

#########################
## Variant Calling     ##
## via HaplotypeCaller ##
#########################
#https://software.broadinstitute.org/gatk/documentation/article?id=11068
#https://software.broadinstitute.org/gatk/documentation/article?id=3893
#https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php
# Input BAM files for HaplotypeCallerSpark do not necessarity need to be sorted by chr
if [ $RUN_GATK_HC_SPARK -eq 1 ];then
    My_HC_VCF=$PROJECT_DIR/GATK4/$Barcode/$SLX.$Barcode.$Cell.bwa.samblaster.recal.HC.spark.g.vcf
    My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.samblaster.recal.bam 
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
		exit
    else
	    if [ ! -s ${My_bam%.bam}.sorted.bam ];then
            echo -e "samtools sort -@ $NT -m 10G -o ${My_bam%.bam}.sorted.bam -T $PROJECT_DIR/scratch/$Barcode $My_bam\n"
            time samtools sort -@ $NT -m 10G -o ${My_bam%.bam}.sorted.bam -T $PROJECT_DIR/scratch/$Barcode $My_bam
        fi
	fi

    My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.samblaster.recal.sorted.bam 
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
		exit
	fi

    # HaplotypeCallerSpark
    # cannot make a compressed vcf file (i.e. .g.vcf.gz)
    # input BAM does not necessarily coordinate sorted 
    echo -e "$GATK_HOME/gatk HaplotypeCallerSpark -I $My_bam -O $My_HC_VCF 2> $PROJECT_DIR/GATK4/$Barcode/$SLX.$Barcode.$Cell.HaplotypeCallerSpark.log\n"
    time $GATK_HOME/gatk --java-options "-Xmx$GATK_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode" \
        HaplotypeCallerSpark \
        --spark-master local[$NT] \
        -ERC GVCF \
        -R $GENOME.2bit \
        -I $My_bam \
        -O $My_HC_VCF 2> $PROJECT_DIR/GATK4/$Barcode/$SLX.$Barcode.$Cell.HaplotypeCallerSpark.log
	echo -e "\e[032mDone: GATK HaplotypeCallerSpark \e[0m\n"

	if [ -s $My_HC_VCF ];then
        echo -e "bgzip -f $My_HC_VCF\n"
        time bgzip -f $My_HC_VCF
    fi
	if [ -s $My_HC_VCF.gz ];then
        echo -e "tabix -f -p vcf $My_HC_VCF.gz\n"
        time tabix -f -p vcf $My_HC_VCF.gz
    fi
fi

# Input BAM files for HaplotypeCaller should be sorted by chr
if [ $RUN_GATK_HC -eq 1 ];then
    My_HC_VCF=$PROJECT_DIR/GATK4/$Barcode/$SLX.$Barcode.$Cell.bwa.samblaster.recal.HC.g.vcf.gz 
    My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.samblaster.sorted.recal.bam 
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
		exit
	fi

    # HaplotypeCaller
    # input BAM must be indexed; therefore coordinate-sorted 
    echo -e "$GATK_HOME/gatk HaplotypeCaller -I $My_bam -O $My_HC_VCF 2> $PROJECT_DIR/GATK4/$Barcode/$SLX.$Barcode.$Cell.HaplotypeCaller.log\n"
    time $GATK_HOME/gatk --java-options "-Xmx$GATK_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode" \
        HaplotypeCaller \
        -ERC GVCF \
        -R $GENOME \
        -I $My_bam \
        -O $My_HC_VCF 2> $PROJECT_DIR/GATK4/$Barcode/$SLX.$Barcode.$Cell.HaplotypeCaller.log
	echo -e "\e[032mDone: GATK HaplotypeCaller \e[0m\n"

	if [ -s $My_HC_VCF ];then
        echo -e "tabix -f -p vcf $My_HC_VCF\n"
        time tabix -f -p vcf $My_HC_VCF
    fi
fi

echo -e "\e[32m$Barcode done for per.barcode.gatk4\e[0m" 
