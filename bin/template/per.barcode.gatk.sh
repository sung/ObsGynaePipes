#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# Last modified: 22/Mar/2018
# First created: 30/Jan/2015
# Based on GATK3

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
mkdir -p $PROJECT_DIR/GATK
mkdir -p $PROJECT_DIR/GATK/$Barcode
mkdir -p $PROJECT_DIR/scratch
mkdir -p $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Cell.$Barcode\e[0m\n"

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
#https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php (only version>=3.6-0 )
My_interval=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.gatk.bam.intervals
if [ $RUN_GATK_INTERVAL -eq 1 ]; then
	#Creating Intervals
	echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator\n"
	time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-I $My_bam \
		-R $GENOME \
		-nt $NT \
		-L $TARGET \
		--filter_mismatching_base_and_quals \
		-known $GATK_MILLS \
		-known $GATK_10KIndel \
		--out $My_interval
	echo -e "\e[032mDone Creating Interval\e[0m\n"
fi

if [ $RUN_GATK_REALIGN -eq 1 ]; then
	if [ -s $My_interval ]; then
		echo -e "time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T IndelRealigner\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T IndelRealigner \
			-I $My_bam \
			-R $GENOME \
			-targetIntervals $My_interval \
			-known $GATK_MILLS \
			-known $GATK_10KIndel \
			-model USE_READS \
			--filter_mismatching_base_and_quals \
			--filter_bases_not_stored \
			--out ${My_bam%.bam}.realigned.bam
		echo -e "\e[032mDone Indel Realigning\e[0m\n"
	
		if [ -s ${My_bam%.bam}.realigned.bam ]; then
			if [ -s ${My_bam%.bam}.realigned.bai ]; then
				echo -e "mv ${My_bam%.bam}.realigned.bai ${My_bam%.bam}.realigned.bam.bai\n"
				mv ${My_bam%.bam}.realigned.bai ${My_bam%.bam}.realigned.bam.bai
			else
				echo -e "samtools index -@ $NT ${My_bam%.bam}.realigned.bam\n"
				time samtools index -@ $NT ${My_bam%.bam}.realigned.bam 
			fi
		else
			echo -e "\e[031m${My_bam%.bam}.realigned.bam not found\e[0m\n"
			exit
		fi
	else
		echo -e "\e[031m$My_interval not found\e[0m\n"
		exit
	fi
fi

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
    My_recalTable=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.samblaster.sorted.bam.BaseRecalibrator.grp
    ###BaseRecalibrator
	echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar --T BaseRecalibrator -I $My_bam \n"
	time java -Xmx$GATK_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode -jar $GATK_HOME/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -nct $NT \
        -I $My_bam \
        -R $GENOME \
        --knownSites $DBSNP \
        --knownSites $GATK_MILLS \
        -cov ReadGroupCovariate \
        -cov QualityScoreCovariate \
        -cov CycleCovariate \
        -cov ContextCovariate \
        --out $My_recalTable 2> $PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.BaseRecalibrator.log
        #-L $TARGET \
        #--knownSites $GATK_10KIndel \
    echo -e "\e[032mDone: BaseRecalibrator\e[0m\n"
	if [ ! -s $My_recalTable ]; then
		echo -e "\e[031m$My_bam or $My_recalTable not found\e[0m\n"
		exit
    fi

    ###PrintReads (it takes long)
    echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T PrintReads -I $My_bam -BQSR $My_recalTable 2> $PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.PrintReads.log \n"
    time java -Xmx$GATK_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode -jar $GATK_HOME/GenomeAnalysisTK.jar \
        -T PrintReads \
        -R $GENOME \
        -nct $NT \
        -BQSR $My_recalTable \
        -I $My_bam \
        --out ${My_bam%.bam}.recal.bam 2> $PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.PrintReads.log
        #-S LENIENT \
	echo -e "\e[032mDone: GATK PrintReads \e[0m\n"

	## Rename indexed bam 
    if [ -s ${My_bam%.bam}.recal.bai ];then
        mv ${My_bam%.bam}.recal.bai ${My_bam%.bam}.recal.bam.bai
        echo -e "samtools index -@ $NT $My_bam\n"
        time samtools index -@ $NT $My_bam
    fi
fi

#########################
## Variant Calling     ##
## via HaplotypeCaller ##
#########################
# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
#https://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode
if [ $RUN_GATK_HC -eq 1 ];then
    My_HC_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.samblaster.sorted.recal.HC.g.vcf.gz 
    My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.samblaster.sorted.recal.bam 
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
		exit
	fi

	## Indexed bam 
    if [ ! -s $My_bam".bai" ];then
        echo -e "samtools index -@ $NT $My_bam\n"
        time samtools index -@ $NT $My_bam
    fi

    # HaplotypeCaller
    # input BAM must be indexed; therefore coordinate-sorted (GATK4 only?)
	echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T HaplotypeCaller -ERC GVCF -I $My_bam 2> $PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.HaplotypeCaller.log\n"
	time java -Xmx$GATK_MEM -jar -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode $GATK_HOME/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
        -ERC GVCF \
		-nct $NT \
		-R $GENOME \
		-I $My_bam \
        --out $My_HC_VCF 2> $PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.HaplotypeCaller.log
	echo -e "\e[032mDone: GATK HaplotypeCaller \e[0m\n"

    # Index vcf
	#if [ -s $My_HC_VCF ];then
    #    echo -e "tabix -f -p vcf $My_HC_VCF\n"
    #    time tabix -f -p vcf $My_HC_VCF
    #fi
fi

echo -e "\e[32m$Barcode done for per.barcode.gatk\e[0m" 
